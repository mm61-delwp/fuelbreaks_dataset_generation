import arcpy
import os
import sys
from typing import Tuple, Optional, List, Dict
import re
from numbers import Number as num


# Configuration parameters
OUTPATH     = r"C:\Data\fuelbreaks"                         # Path to output geodatabase
GDB_NAME    = r"fuelbreaks_evaluation.gdb"                        # Preferred name for output geodatabase
SPATIAL_REF = arcpy.SpatialReference(3111)                  # Vicgrid94

# Fishnet parameters
X1, Y1 = 2125000, 2259000                                   # Lower-left coordinates
X2, Y2 = 2940000, 2828000                                   # Upper-right coordinates
GRID_EXTENT = arcpy.Extent(X1, Y1, X2, Y2, 3111)            # (MinX, MinY, MaxX, MaxY, SpatialReference)
CELL_HA = 100                                               # Cell size in hectares
CELL_SHAPE = "HEXAGON"                                      # Coz it looks cool
CELL_M = CELL_HA * 10000                                    # Cell size in metres (calculated)


# Datasets required outside of data matrix processing
FIREFMZ = r"C:\Data\CSDL\FIRE.GDB\\FIREFMZ"


""" DATA MATRIX FOR DUMMIES

    "name here": {
        "fieldname":        str (required)  - name of field to create and populate in fishnet
        "source":           str (required)  - full path to features we want to add data to fishnet from
        "method":           str (required)  - determines how to report joins
        "source_field":     str (dependent) - used to classify grid cells or for numeric summary methods
        "match_option":     str (optional)  - determines spatial join type. Limited set from ArcGIS match_options
        "group_by":         str (optional)  - field name to group joins by (reduce multiple joins to 1 per group)
        "filter":           SQL (optional)  - apply a filter to the data source
        "distance":         int (optional)  - linear buffer distance; used with match option "WITHIN_A_DISTANCE"
        }

    Valid options for method:
        "presence"              - returns 1 for any grid cell with valid join, else 0
        "count"                 - returns count of valid joins found to the grid cell 
        "sum"                   - returns sum of source_field values for all valid joins
        "average"               - returns mean of source_field values for all valid joins
        "area_weighted_average" - returns sum(area/CELL_M * source_field) for all valid joins
        "maximum"               - returns largest source_field value for all valid joins
        "minimum"               - returns smallest source_field value for all valid joins
        "area"                  - returns area of intersection if intersected, else 0
        "percentage"            - returns proportion of area covered by intersection if intersected, else 0
        "classify"              - transfers source_field value from joined features
"""

data_matrix = {
    "fire_fmz": {
        "fieldname":    "fire_fmz",
        "source":       r"C:\Data\CSDL\FIRE.GDB\FIREFMZ",
        "method":       "classify",
        "source_field": "X_ZONETYPE",
        "match_option": "CLOSEST", # Probably should be "LARGEST_OVERLAP"? But it's very slow.
        "group_by":     None,
        "filter":       "ZONETYPE > 0"
    },
    "control_line_freq": {
        "fieldname":    "control_line_freq",
        "source":       r"c:\Projects\20250725_StatewideFuelBreaks\spatial_data\input_data.gdb\HistoricalControlLines_2013to2025",
        "method":       "count",
        "match_option": "INTERSECT",
        "method":       "count",
        "group_by":     "SEASON",
        "filter":       None
    },
    "bu_id": {
        "fieldname":    "bu_id",
        "source":       r"C:\Projects\20250725_StatewideFuelBreaks\spatial_data\statewide_burn_units_20220224.shp",
        "method":       "classify",
        "source_field": "BUID",
        "match_option": "CLOSEST",
        "group_by":     None,
        "filter":       None
    },
    "is_bu_boundary": {
        "fieldname":    "is_bu_boundary",
        "source":       r"C:\Projects\20250725_StatewideFuelBreaks\spatial_data\statewide_burn_units_20220224.shp",
        "method":       "presence",
        "match_option": "CROSSED_BY_THE_OUTLINE_OF",
        "group_by":     None,
        "filter":       None
    },
    "is_apzbmz_boundary": {
        "fieldname":    "is_apzbmz_boundary",
        "source":       r"C:\Data\CSDL\FIRE.GDB\FIREFMZ",
        "method":       "presence",
        "match_option": "CROSSED_BY_THE_OUTLINE_OF",
        "group_by":     None,
        "filter":       "ZONETYPE IN (1, 2)"
    },
    "is_bmzlmz": {
        "fieldname":    "is_bmzlmz",
        "source":       r"C:\Data\CSDL\FIRE.GDB\FIREFMZ",
        "method":       "presence",
        "match_option": "INTERSECT",
        "group_by":     None,
        "filter":       "ZONETYPE IN (2, 3)"
    },
    "is_apzbmz": {
        "fieldname":    "is_apzbmz",
        "source":       r"C:\Data\CSDL\FIRE.GDB\FIREFMZ",
        "method":       "presence",
        "match_option": "INTERSECT",
        "group_by":     None,
        "filter":       "ZONETYPE IN (1, 2)"
    },
    "is_1k_buff": {
        "fieldname":    "is_1k_buff",
        "source":       r"C:\Data\CSDL\FIRE.GDB\FIREFMZ",
        "method":       "presence",
        "match_option": "HAVE_THEIR_CENTER_IN",
        "distance":     1000,
        "group_by":     None,
        "filter":       "ZONETYPE = 0"
    },
    "is_heavy_strategic": {
        "fieldname":    "is_heavy_strategic",
        "source":       r"C:\Data\CSDL\FORESTSROADS.GDB\RDB_FIREACCESS",
        "method":       "presence",
        "match_option": "INTERSECT",
        "filter":       "EXISTING_FIREACCESS = 'Heavy'"
    },
    "is_medium_strategic": {
        "fieldname":    "is_medium_strategic",
        "source":       r"C:\Data\CSDL\FORESTSROADS.GDB\RDB_FIREACCESS",
        "method":       "presence",
        "match_option": "INTERSECT",
        "filter":       "EXISTING_FIREACCESS = 'Medium'"
    },
    "catchment_bdy_minor": {
        "fieldname":    "catchment_bdy_minor",
        "source":       r"c:\data\CSDL\WATER.gdb\SDL_CATCH",
        "method":       "presence",
        "match_option": "CROSSED_BY_THE_OUTLINE_OF",
    },
    "catchment_bdy_major": {
        "fieldname":    "catchment_bdy_major",
        "source":       r"C:\Projects\20250725_StatewideFuelBreaks\spatial_data\Catchment_Basin.shp",
        "method":       "presence",
        "match_option": "CROSSED_BY_THE_OUTLINE_OF",
    },
    "burnability": {
        "fieldname":    "burnability",
        "source":       r"C:\Projects\20250725_StatewideFuelBreaks\spatial_data\Burnability_2024_final.shp",
        "method":       "presence",
        "source_field": "gridcode",
        "match_option": "HAVE_THEIR_CENTER_IN",
        "filter":       "gridcode = 1"
    },
    "burn_freq": {
        "fieldname":    "burn_freq",
        "source":       r"C:\Projects\20250725_StatewideFuelBreaks\spatial_data\Risk2_gridcell_TimesBurnt_FireFMZ.shp",
        "method":       "average",
        "source_field": "times_burn",
        "match_option": "INTERSECT",
    },
    "locality_hl": {
        "fieldname":    "locality_hl",
        "source":       r"C:\Projects\20250725_StatewideFuelBreaks\spatial_data\Risk2_Locality_HouseLoss.shp",
        "method":       "average",
        "source_field": "hl",
        "match_option": "INTERSECT",
    },
    "bu_bdy_hl": {
        "fieldname":    "bu_bdy_hl",
        "source":       r"C:\Projects\20250725_StatewideFuelBreaks\spatial_data\BurnUnit_Risk2_HouseLoss_FireSize.shp",
        "method":       "sum", # total of boundaries for burn units on either/all side (what about intersections? should this be avg?)
        "source_field": "avg_hl",
        "match_option": "CROSSED_BY_THE_OUTLINE_OF"
    },
    "bu_bdy_firesize": {
        "fieldname":    "bu_bdy_firesize",
        "source":       r"C:\Projects\20250725_StatewideFuelBreaks\spatial_data\BurnUnit_Risk2_HouseLoss_FireSize.shp",
        "method":       "sum", # total of boundaries for burn units on either/all side (what about intersections? should this be avg?)
        "source_field": "avg_ha",
        "match_option": "CROSSED_BY_THE_OUTLINE_OF"
    },
    "bu_bdy_firesize_2h": {
        "fieldname":    "bu_bdy_firesize_2h",
        "source":       r"C:\Projects\20250725_StatewideFuelBreaks\spatial_data\BurnUnit_Risk2_HouseLoss_FireSize.shp",
        "method":       "sum", # total of boundaries for burn units on either/all side (what about intersections? should this be avg?)
        "source_field": "avg_ha_2h",
        "match_option": "CROSSED_BY_THE_OUTLINE_OF"
    },
    "bu_bdy_hl_red": {
        "fieldname":    "bu_bdy_hl_red",
        "source":       r"C:\Projects\20250725_StatewideFuelBreaks\spatial_data\BurnUnit_Risk2_HouseLoss_FireSize_SuppnEffects.shp",
        "method":       "sum", # total of boundaries for burn units on either/all side (what about intersections? should this be avg?)
        "source_field": "red_avg_hl",
        "match_option": "CROSSED_BY_THE_OUTLINE_OF"
    },
    "bu_bdy_firesize_2h_red": {
        "fieldname":    "bu_bdy_firesize_2h_red",
        "source":       r"C:\Projects\20250725_StatewideFuelBreaks\spatial_data\BurnUnit_Risk2_HouseLoss_FireSize_SuppnEffects.shp",
        "method":       "sum", # total of boundaries for burn units on either/all side (what about intersections? should this be avg?)
        "source_field": "red_avg__1",
        "match_option": "CROSSED_BY_THE_OUTLINE_OF"
    },
    "wc_travel_time": {
        "fieldname":    "wc_travel_time",
        "source":       r"C:\Projects\20250725_StatewideFuelBreaks\spatial_data\Workcentre_Travel_Times.shp",
        "method":       "average", # total of boundaries for burn units on either/all side (what about intersections? should this be avg?)
        "source_field": "FromBreak",
        "match_option": "INTERSECT"
    },
    
}

class Fishnet:
    """
    A class to handle fishnet generation and data processing operations.
    """
    
    def __init__(self, workspace: str, gdb_name: str, cell_size: int, shape: str = "HEXAGON",
                 extent: Optional[arcpy.Extent] = None,
                 spatial_reference: Optional[arcpy.SpatialReference] = None):
        """
        Initialise the Fishnet.
        
        Args:
            workspace (str): Path to the geodatabase or workspace
            spatial_reference (arcpy.SpatialReference, optional): Coordinate system to use
        """
        self.spatial_reference = spatial_reference or arcpy.SpatialReference(3111)  # Default to Vicgrid94
        
        # Create output directory and geodatabase
        self.output_path = OUTPATH
        self.gdb_name = GDB_NAME
        os.makedirs(self.output_path, exist_ok=True)
        self.output_gdb = os.path.join(self.output_path, gdb_name)
        if arcpy.Exists(self.output_gdb):
            # arcpy.Delete_management(self.output_gdb)
            arcpy.AddMessage(f"Workspace exists: {self.output_gdb}")
        else:
            arcpy.CreateFileGDB_management(self.output_path, gdb_name, "10.0")
            arcpy.AddMessage(f"Created workspace: {self.output_gdb}")
        self.workspace = self.output_gdb
        self.layer_cache = [] # so we don't need to re-add and re-project layers that are re-used

        # Set up ArcPy environment
        arcpy.env.workspace = self.workspace
        arcpy.env.overwriteOutput = True
        arcpy.env.cartographicCoordinateSystem = self.spatial_reference
        arcpy.env.outputCoordinateSystem = self.spatial_reference
        arcpy.env.parallelProcessingFactor = 80
        
        print(f"Initialised workspace: {workspace}")
    
        # Create analysis grid    
        try:
            # Output path
            grid_path = os.path.join(self.workspace, "grid_full")

            # generate grid
            arcpy.management.GenerateTessellation(grid_path, extent, shape, cell_size)

            # Add fields for area and centroid coordinates
            fields_to_add = [
                ("CELL_AREA", "DOUBLE", "", "CELL_AREA"),
                ("CENTROID_X", "DOUBLE", "", "CENTROID_X"),
                ("CENTROID_Y", "DOUBLE", "", "CENTROID_Y"),
            ]
            
            for field_name, field_type, field_length, field_alias in fields_to_add:
                if field_length:
                    arcpy.management.AddField(grid_path, field_name, field_type, field_length=field_length, field_alias=field_alias)
                else:
                    arcpy.management.AddField(grid_path, field_name, field_type, field_alias=field_alias)
            
            # Calculate field values
            with arcpy.da.UpdateCursor(grid_path, ["SHAPE@", "CELL_AREA", "CENTROID_X", "CENTROID_Y"]) as cursor:
                for row in cursor:
                    # alias field names for sanity
                    (shape, cell_area, x, y) = row

                    # Calculate new values
                    cell_area = shape.area/10000         # Cell area in hectares
                    x = shape.centroid.X                 # Centroid X
                    y = shape.centroid.Y                 # Centroid Y
                    
                    # Reassemble row with calculated values
                    updated_row = (shape, cell_area, x, y)
                    cursor.updateRow(updated_row)
            
            # Add spatial index
            try:
                arcpy.AddSpatialIndex_management(grid_path)
            except:
                pass # index probably already exists

            print(f"Created analysis grid: {grid_path}")
            self.fishnet = grid_path

            # create in_memory version of fishnet
            self.fishnet_mem = "in_memory\\fishnet"
            arcpy.management.CopyFeatures(self.fishnet, self.fishnet_mem)
            
        except Exception as e:
            print(f"Error creating analysis grid: {str(e)}")
            raise
        
    def clip_to(self, features, query = None):
        """ Remove grid cells that don't overlap supplied features, with optional query"""

        if query:
            # apply selection query to selection features
            features = arcpy.management.SelectLayerByAttribute(features, "NEW_SELECTION", query)

        # Output path
        output_path = os.path.join(self.workspace, "grid_clip")

        # select grid cells from fishnet
        grid_select = arcpy.management.SelectLayerByLocation(self.fishnet, "INTERSECT", features, selection_type="NEW_SELECTION")

        # export selected cells to grid_clip
        arcpy.conversion.ExportFeatures(grid_select, output_path)

        # Re-assign fishnet
        self.fishnet = output_path

        # update in_memory version of fishnet
        arcpy.management.CopyFeatures(self.fishnet, self.fishnet_mem)

        #remove selections
        arcpy.management.SelectLayerByAttribute(features, "CLEAR_SELECTION")
        arcpy.management.SelectLayerByAttribute(self.fishnet, "CLEAR_SELECTION") # probably not required?

        print(f"Clipped analysis grid: {output_path}")
    
    def join_layer(self, dataset_name: str, config: dict) -> None:
        """ Process a single spatial join based on configuration. 
        
        Dear future me... 
        The Fishnet class stores two references to fishnet data:
            - self.fishnet points to the GDB dataset, at time of writing 'grid_clip' but that might change
            - self.fishnet_mem points to an in_memory copy of the dataset
        This method uses the following flow:
            - load the join data to memory
            - process in_memory join data against fishnet_mem for efficiency
            - join results back to GDB fishet
            - refresh fishnet_mem with a copy of the GDB fishnet
        """
        
        # Extract configuration
        fieldname = config.get("fieldname", dataset_name)       # Field to add to fishnet
        source_og = config["source"]                            # Data location
        match_option = config.get("match_option", "INTERSECT")  # Type of join required
        method = config.get("method", "classify")               # Type of data to add to fishnet
        group_by = config.get("group_by", None)                 # How to summarise joined features
        source_field = config.get("source_field", None)         # Field to transfer/count/sum/average
        query = config.get("filter", None)                      # Query to restrict join to certain values
        distance = config.get("distance", None)                 # Buffer distance for selection

        fishnet_id = "GRID_ID"
        basename = os.path.basename(source_og)
        
        # Validate source exists
        if not arcpy.Exists(source_og):
            raise ValueError(f"Source does not exist: {source_og}")
        
        # Load source to memory if it doesn't already exist, reprojecting if required
        source = f"in_memory\\{os.path.splitext(basename)[0]}"

        if not arcpy.Exists(source):
            self.layer_cache.append(source)

            # Check if projection is needed
            source_sr = arcpy.Describe(source_og).spatialReference
            env_sr = arcpy.env.outputCoordinateSystem

            if source_sr.factoryCode != env_sr.factoryCode: # Different SR - need to project
                arcpy.management.Project(source_og, source, env_sr)
                print(f"Reprojected {basename} and loaded to memory: {source}")
            else: # Same SR - just copy for speed
                arcpy.management.CopyFeatures(source_og, source)
                print(f"Loaded {basename} to memory: {source}")
            
            # Add spatial indexing if not present (improves performance)
            try:
                arcpy.AddSpatialIndex_management(source)
                # print(f"Indexed {basename}")
            except:
                # print(f"Index exists for {basename}")
                pass # index probably already exists
        else:
            print(f"Using {basename} from memory: {source}")
        
        # Apply filter if specified
        if query:
            source = arcpy.management.SelectLayerByAttribute(source, 'NEW_SELECTION', query)
            print(f"Applied query to {basename}: {query}")

        # Check dataset size and warn if large (numbers are arbitrary - time depends on size & complexity & method)
        result = arcpy.management.GetCount(source)
        feature_count = int(result.getOutput(0))
        if feature_count > 50000:
            print(f"{basename} is large ({feature_count} features); Processing join & {method} may be slow.")

        # Process based on method
        if method == "count":
            
            # create join
            temp_join = "in_memory\\temp_join"
            
            arcpy.analysis.SpatialJoin(
                target_features=self.fishnet_mem,
                join_features=source,
                out_feature_class=temp_join,
                join_operation="JOIN_ONE_TO_MANY",
                join_type="KEEP_ALL",
                match_option=match_option,
                search_radius=distance
                )
            
            if group_by:
                # Count unique values in group_by field per fishnet cell
                temp_stats = "in_memory\\temp_stats"
                
                # get unique combinations
                case_fields = [fishnet_id, group_by]
                arcpy.analysis.Statistics(
                    in_table=temp_join,
                    out_table=temp_stats,
                    statistics_fields=[[group_by, "COUNT"]],  # don't actually care about this count, just uniques
                    case_field=case_fields
                )
                
                # count unique group_by values per fishnet cell
                final_stats = "in_memory\\final_stats"
                arcpy.analysis.Statistics(
                    in_table=temp_stats,
                    out_table=final_stats,
                    statistics_fields=[[group_by, "COUNT"]], # count of grouped values
                    case_field=[fishnet_id]
                )
                
                # join back to fishnet
                self.join_and_calc_field(
                    target=self.fishnet,
                    join_table=final_stats,
                    target_field=fishnet_id,
                    join_field=fishnet_id,
                    transfer_field=f"COUNT_{group_by}", # probably INT type?
                    new_fieldname=fieldname,
                    default_value=0
                )
                
                print(f"Populated fishnet[{fieldname}] with counts from {basename} grouped by {group_by}")

            else:
                # simple count of all intersecting features
                temp_stats = "in_memory\\temp_stats"
                arcpy.analysis.Statistics(
                    in_table=temp_join,
                    out_table=temp_stats,
                    statistics_fields=[["Join_Count", "SUM"]],  # Join_Count is auto-generated
                    case_field=[fishnet_id]
                )
                
                # Join back to fishnet
                self.join_and_calc_field(
                    target=self.fishnet,
                    join_table=temp_stats,
                    target_field=fishnet_id,
                    join_field=fishnet_id,
                    transfer_field="SUM_Join_Count", # probably INT type?
                    new_fieldname=fieldname,
                    default_value=0
                )

                print(f"Populated fishnet[{fieldname}] with counts from {basename}")                     
        
        elif method == "classify":

            # If no source field specified, error
            if not source_field:
                raise ValueError("No field identified for classify")
            
            # perform spatial join (temporary)
            temp_join = "in_memory\\temp_join"
            arcpy.analysis.SpatialJoin(
                target_features=self.fishnet_mem,
                join_features=source,
                out_feature_class=temp_join,
                join_operation="JOIN_ONE_TO_ONE",
                join_type="KEEP_ALL",
                match_option=match_option,
                search_radius=distance
            )

            # Join back to fishnet
            self.join_and_calc_field(
                target=self.fishnet,
                join_table=temp_join,
                target_field=fishnet_id,
                join_field=fishnet_id,
                transfer_field=source_field, # Type = source_field Type
                new_fieldname=fieldname,
                default_value=None
            )

            print(f"Populated fishnet[{fieldname}] with values from {basename}[{source_field}]")
      
        
        elif method == "sum":
            # create join
            temp_join = "in_memory\\temp_join"
            
            arcpy.analysis.SpatialJoin(
                target_features=self.fishnet_mem,
                join_features=source,
                out_feature_class=temp_join,
                join_operation="JOIN_ONE_TO_MANY",
                join_type="KEEP_ALL",
                match_option=match_option,
                search_radius=distance
                )
            
            # simple sum of values field of all intersecting features
            temp_stats = "in_memory\\temp_stats"
            arcpy.analysis.Statistics(
                in_table=temp_join,
                out_table=temp_stats,
                statistics_fields=[[source_field, "SUM"]],
                case_field=[fishnet_id]
            )
            
            # Join back to fishnet
            self.join_and_calc_field(
                target=self.fishnet,
                join_table=temp_stats,
                target_field=fishnet_id,
                join_field=fishnet_id,
                transfer_field=f"SUM_{source_field}", # Type = source_field Type??
                new_fieldname=fieldname,
                default_value=0
            )

            print(f"Populated fishnet[{fieldname}] with sum values from {basename}[{source_field}]")
        
        elif method == "area":
            # create join
            temp_join = "in_memory\\temp_join"
            arcpy.analysis.PairwiseIntersect([self.fishnet_mem, source], out_feature_class=temp_join)

            # Add hectare field and populate
            arcpy.management.AddField(
                in_table=temp_join,
                field_name="Area_Ha",
                field_type="FLOAT", # Type = FLOAT always
                field_scale=2)
            
            arcpy.management.CalculateField(temp_join, "Area_Ha", f"!SHAPE@AREA!/10000", "PYTHON3")

            # simple sum of values field of all intersecting features
            temp_stats = "in_memory\\temp_stats"

            area_field = "Area_Ha"
            arcpy.analysis.Statistics(
                in_table=temp_join,
                out_table=temp_stats,
                statistics_fields=[[area_field, "SUM"]], 
                case_field=[fishnet_id]
            )
            
            # Join back to fishnet
            self.join_and_calc_field(
                target=self.fishnet,
                join_table=temp_stats,
                target_field=fishnet_id,
                join_field=fishnet_id,
                transfer_field=f"SUM_{area_field}",
                new_fieldname=fieldname,
                default_value=0
            )

            print(f"Populated fishnet[{fieldname}] with area values from {basename}[{source_field}]")

        elif method == "average":
            # create join
            temp_join = "in_memory\\temp_join"
            arcpy.analysis.PairwiseIntersect([self.fishnet_mem, source], out_feature_class=temp_join)
            
            # simple sum of values field of all intersecting features
            temp_stats = "in_memory\\temp_stats"
            arcpy.analysis.Statistics(
                in_table=temp_join,
                out_table=temp_stats,
                statistics_fields=[[source_field, "MEAN"]], 
                case_field=[fishnet_id]
            )
            
            # Join back to fishnet
            self.join_and_calc_field(
                target=self.fishnet,
                join_table=temp_stats,
                target_field=fishnet_id,
                join_field=fishnet_id,
                transfer_field=f"MEAN_{source_field}", # Type = source_field Type??
                new_fieldname=fieldname,
                default_value=0
            )

            print(f"Populated fishnet[{fieldname}] with sum values from {basename}[{source_field}]")
        
        elif method == "area_weighted_average":
            # create join
            temp_join = "in_memory\\temp_join"
            arcpy.analysis.PairwiseIntersect([self.fishnet_mem, source], out_feature_class=temp_join)
            
            # calculate weighted value field (source_field * area)
            weighted_field = f"WEIGHTED_{source_field}"
            arcpy.management.AddField(temp_join, weighted_field, "DOUBLE")
            arcpy.management.CalculateField(temp_join, weighted_field, f"!{source_field}! * !shape.area@squaremeters!", "PYTHON3")
            
            # sum of weighted values for each fishnet cell
            temp_stats = "in_memory\\temp_stats"
            arcpy.analysis.Statistics(
                in_table=temp_join,
                out_table=temp_stats,
                statistics_fields=[[weighted_field, "SUM"]], 
                case_field=[fishnet_id]
            )
            
            # calculate area-weighted average - note that areas with no intersection effectively get 0 value
            weighted_avg_field = f"AREA_WEIGHTED_AVG"
            arcpy.management.AddField(temp_stats, weighted_avg_field, "DOUBLE")
            arcpy.management.CalculateField(temp_stats, weighted_avg_field, f"!SUM_{weighted_field}! / {CELL_M}", "PYTHON3" )
            
            # join back to fishnet
            self.join_and_calc_field(
                target=self.fishnet,
                join_table=temp_stats,
                target_field=fishnet_id,
                join_field=fishnet_id,
                transfer_field=weighted_avg_field,
                new_fieldname=fieldname,
                default_value=0
            )

            print(f"Populated fishnet[{fieldname}] with area-weighted average values from {basename}[{source_field}]")
            
        elif method == "presence":
            # Add field to fishnet
            if fieldname not in [f.name for f in arcpy.ListFields(self.fishnet)]:
                arcpy.management.AddField(self.fishnet, fieldname, "SHORT")
            
            # Set all to 0 initially
            arcpy.management.CalculateField(self.fishnet, fieldname, "0", "PYTHON3")

            # Spatial selection
            intersected = arcpy.management.SelectLayerByLocation(
                in_layer=self.fishnet_mem, 
                overlap_type=match_option, 
                select_features=source,
                search_distance=distance,
                selection_type="NEW_SELECTION")

            # Get unique grid IDs from selected features
            unique_grid_ids = set()
            with arcpy.da.SearchCursor(intersected, [fishnet_id]) as cursor:
                for row in cursor:
                    unique_grid_ids.add(row[0])

            # Clear selection
            arcpy.management.SelectLayerByAttribute(self.fishnet, "CLEAR_SELECTION")

            # Update field values using cursor for selected grid IDs only
            with arcpy.da.UpdateCursor(self.fishnet, [fishnet_id, fieldname]) as cursor:
                for row in cursor:
                    if row[0] in unique_grid_ids:
                        row[1] = 1
                        cursor.updateRow(row)

            print(f"Populated fishnet[{fieldname}] with presence from {basename}")
        
        else:
            raise ValueError(f"Unknown method: {method}")
        
        # clear any selection
        arcpy.management.SelectLayerByAttribute(source, 'CLEAR_SELECTION')

        # update in_memory version of fishnet
        arcpy.management.CopyFeatures(self.fishnet, self.fishnet_mem)


    ### HELPER METHODS
    def join_and_calc_field(self, target, join_table, target_field, join_field, transfer_field, new_fieldname, default_value=None):
        """Join table and calculate field with null handling."""
        
        # Create dictionary from join table or feature class
        value_dict = {}
        with arcpy.da.SearchCursor(join_table, [join_field, transfer_field]) as cursor:
            for row in cursor:
                value_dict[row[0]] = row[1] if row[1] is not None else default_value
        
        # Add field if it doesn't exist
        if new_fieldname not in [f.name for f in arcpy.ListFields(target)]:
            # Determine field type from transfer field
            fields = arcpy.ListFields(join_table, transfer_field)
            field_type = fields[0].type
            arcpy.AddField_management(target, new_fieldname, field_type)
        
        # Update values
        with arcpy.da.UpdateCursor(target, [target_field, new_fieldname]) as cursor:
            for row in cursor:
                row[1] = value_dict.get(row[0], default_value)
                cursor.updateRow(row)

    def get_field_type(self, fc, fieldname):
        """Get the field type for a field in a feature class."""
        fields = arcpy.ListFields(fc)
        for field in fields:
            if field.name == fieldname:
                if field.type in ["String", "Text"]:
                    return "TEXT"
                elif field.type in ["Integer", "SmallInteger"]:
                    return "LONG"
                elif field.type in ["Double", "Float", "Single"]:
                    return "DOUBLE"
                elif field.type == "Date":
                    return "DATE"
        return "TEXT"
    
    def check_config(self, config: dict) -> None:
        # Define valid values for each parameter
        VALID_METHODS = ["presence", "count", "sum", "average", "area_weighted_average", "maximum", "minimum", "area", "classify"]
        
        VALID_MATCH_OPTS = ["INTERSECT", "WITHIN_A_DISTANCE", "CONTAINS", "COMPLETELY_CONTAINS", "WITHIN", "CLOSEST",
                            "COMPLETELY_WITHIN","HAVE_THEIR_CENTER_IN", "LARGEST_OVERLAP", "CROSSED_BY_THE_OUTLINE_OF"]
        
        METHODS_REQUIRING_FIELD = ["classify", "sum", "average", "mean", "maximum", "minimum"]
        
        issues = []
        warnings = []
        
        # Check required fields
        if "source" not in config:
            issues.append(f"Missing required 'source' field")
        else:
            # Check if source exists
            if not arcpy.Exists(config["source"]):
                issues.append(f"Source does not exist: {config['source']}")

        # Check method validity
        method = config.get("method")
        if method not in VALID_METHODS:
            issues.append(f"Invalid method '{method}'. Valid methods: {', '.join(VALID_METHODS)}")
        
        # Check join_type validity
        match_option = config.get("match_option")
        if match_option and match_option.upper() not in VALID_MATCH_OPTS:
            issues.append(f"Invalid match_option '{match_option}'. Valid types: {', '.join(VALID_MATCH_OPTS)}")
        
        # Check if source_field is required but missing
        if method in METHODS_REQUIRING_FIELD:
            if "source_field" not in config or config["source_field"] is None:
                issues.append(f"Method '{method}' requires 'source_field' to be specified")
            elif arcpy.Exists(config["source"]):
                # Check if field exists in source
                field_names = [f.name for f in arcpy.ListFields(config["source"])]
                if config["source_field"] not in field_names:
                    issues.append(f"Field '{config['source_field']}' not found in source. Available fields: {', '.join(field_names[:10])}")
                else:
                    # Check field type compatibility
                    field_type = None
                    for f in arcpy.ListFields(config["source"]):
                        if f.name == config["source_field"]:
                            field_type = f.type
                            break
                    
                    if method in ["sum", "average", "maximum", "minimum"]:
                        if field_type not in ["Double", "Float", "Single", "Integer", "SmallInteger", "Long", "Short"]:
                            issues.append(f"Method '{method}' requires numeric field, but '{config['source_field']}' is type '{field_type}'")
                    
            # Check fieldname validity
            fieldname = config.get("fieldname")
            if len(fieldname) > 64:  # ArcGIS field name length limit
                issues.append(f"Field name '{fieldname}' exceeds 64 character limit")
            
            if not fieldname[0].isalpha():
                issues.append(f"Field name '{fieldname}' must start with a letter")
            
            if not all(c.isalnum() or c == '_' for c in fieldname):
                issues.append(f"Field name '{fieldname}' contains invalid characters (use only letters, numbers, underscore)")
            
            # Check group_by field if specified
            if "group_by" in config and config["group_by"] is not None:
                if method not in ["count"]:
                    warnings.append(f"'group_by' is only used with 'count' method, will be ignored for method '{method}'")
                elif arcpy.Exists(config["source"]):
                    field_names = [f.name for f in arcpy.ListFields(config["source"])]
                    if config["group_by"] not in field_names:
                        issues.append(f"{config} group_by field '{config['group_by']}' not found in source")
            
        # Check for duplicate field names
        fieldnames = [config.get("fieldname", name) for name, config in data_matrix.items()]
        duplicates = [name for name in fieldnames if fieldnames.count(name) > 1]
        if duplicates:
            issues.append(f"Duplicate field names found: {', '.join(set(duplicates))}")
        
        # Check for Field Calculator syntax (!field!)
        filter_string = config.get("filter")
        if filter_string:
            calc_fields = re.findall(r'![^!]+!', filter_string)
            if calc_fields:
                issues.append(f"Field Calculator syntax in filter string; replace with SQL: {', '.join(calc_fields)}")
            
            # Check for Python syntax ([field])
            python_fields = re.findall(r'\[[^\]]+\]', filter_string)
            if python_fields:
                issues.append(f"Python field syntax in filter string; replace with SQL: {', '.join(python_fields)}")
            
            # Check for unmatched quotes
            if filter_string.count("'") % 2 != 0:
                issues.append(f"Unmatched single quotes found, check filter string")
            
            if filter_string.count('"') % 2 != 0:
                issues.append(f"Unmatched double quotes found, check filter string")
        
        # Print report
        fieldname = config.get("fieldname")
        if issues or warnings:
            print(f"Configuration validation for {fieldname} found {len(issues)} error(s) and {len(warnings)} warning(s):")
            if issues:
                print(" ERRORS (must fix):")
                for issue in issues:
                    print(f"  • {issue}")
            if warnings:
                print(" WARNINGS (consider fixing):")
                for warning in warnings:
                    print(f"  • {warning}")
        else:
            print(f"Configuration for {fieldname} appears valid")
        
    def normalise(self, field_name: str, dataset: Optional[str], group_by: Optional[str] = None,
                method: Optional[str] = "ZERO", output_field: Optional[str] = None) -> None:
        """
        Normalise field values between 0 and 1
        Parameters:
            field_name      - name of field to be normalised
            dataset         - feature class (optional, defaults to self.fishnet)
            group_by        - grouping field name (optional - if provided, values will be normalised against group)
            method          - "ZERO" or "MIN" normalise relative to zero (default) or use minimum value as zero
            output_field    - field name to output normalised values (optional - default will overwrite). Will be added if it doesn't exist.
        """
        
        cursor_fields = [field_name]

        # Set dataset to self.fishnet if not provided
        if not dataset:
            dataset = self.fishnet
        
        # Set target to output_field if provided
        if output_field:
            target_field = output_field
            cursor_fields.append(output_field)
            
            if target_field not in [f.name for f in arcpy.ListFields(dataset)]: # Add field to fishnet
                arcpy.management.AddField(dataset, target_field, "FLOAT")
        else: # Overwrite input field
            target_field = field_name

        if group_by:
            cursor_fields.append(group_by) 
            
            # Use Statistics tool with case field
            stats_table = "in_memory\\norm_stats"
            arcpy.analysis.Statistics(
                in_table=dataset,
                out_table=stats_table,
                statistics_fields=[[field_name, "MIN"], [field_name, "MAX"]],
                case_field=[group_by]
            )
            
            # iterate cursor over stats table
            with arcpy.da.SearchCursor(stats_table, [group_by, f"MIN_{field_name}", f"MAX_{field_name}"]) as stats_cursor:
                for row in stats_cursor:
                    group = row[0]
                    min_val = 0 if method == "ZERO" else row[1]
                    max_val = row[2]

                    # format where_clause to handle string or int data types
                    if isinstance(group, str):
                        where_clause = f"{group_by} = '{group}'"
                    else:
                        where_clause = f"{group_by} = {group}"
            
                    # Update field values using cursor selection
                    with arcpy.da.UpdateCursor(dataset, cursor_fields, where_clause=where_clause) as cursor:
                        for row in cursor:

                            if max_val == min_val:  # Handle division by zero
                                norm_val = 0
                            else:
                                norm_val = (row[0] - min_val) / (max_val - min_val)
                            
                            if output_field:
                                row[1] = norm_val  # Update output field
                            else:
                                row[0] = norm_val  # Update original field
                            cursor.updateRow(row)
            
            # Clean up
            arcpy.management.Delete(stats_table)

        else:
            # Simple approach for overall min/max
            values = [row[0] for row in arcpy.da.SearchCursor(dataset, [field_name]) if row[0] is not None]
            if not values:
                print(f"No valid values found in {field_name}")
                return
                
            min_val = 0 if method == "ZERO" else min(values)
            max_val = max(values)
            
            if max_val == min_val:  # Handle division by zero
                print(f"All values in {field_name} are the same ({max_val}). Setting normalised values to 0.")
                arcpy.management.CalculateField(self.fishnet, target_field, "0", "PYTHON3")
                return

            # Update field values using cursor
            with arcpy.da.UpdateCursor(dataset, cursor_fields) as cursor:
                for row in cursor:
                    norm_val = (row[0] - min_val) / (max_val - min_val)
                    
                    if output_field:
                        row[1] = norm_val  # Update output field
                    else:
                        row[0] = norm_val  # Update original field
                    cursor.updateRow(row)

        print(f"Normalised {field_name} using {method} method" + (f" grouped by {group_by}" if group_by else ""))


        
def main():
    """
    Main program for creating fishnet dataset and assigning evaluation data to it.
    """

    try:

        # Initialise fishnet
        AnalysisDataset = Fishnet(OUTPATH, GDB_NAME, CELL_M, CELL_SHAPE, GRID_EXTENT, SPATIAL_REF)

        # Remove grid cells not in zoned public land
        AnalysisDataset.clip_to(FIREFMZ, "ZONETYPE > 0")
        
        # Spatial join with data layers
        for config in data_matrix:
            AnalysisDataset.check_config(data_matrix[config])
            AnalysisDataset.join_layer(config, data_matrix[config])
        
        print(f"\nDataset generation completed successfully!")
        print(f"Dataset output: {AnalysisDataset.fishnet}")



        #############################################################################
        #   Index generation code
        #   Customise or replace below to combine analysis fields as required
        #############################################################################

        # Copy grid

        new_grid = "in_memory\\fishnet"
        arcpy.management.CopyFeatures(AnalysisDataset.fishnet, new_grid)

        # add extra fields for index calculation
        extra_fields = [
            "cpi_c_score",
            "cpi_s_score", 
            "bu_bdy_travel_time",
            "bu_bdy_freq",
            "containment_priority_index",
            "strategic_location_index",
            "landscape_protection_index",
            "fuel_management_index",
            "strategic_protection_index",
            "frq",
            "bushfire_risk_index", 
            "asset_protection_index"]

        for field_name in extra_fields:
            arcpy.management.AddField(new_grid, field_name, "FLOAT")

        # normalise fields for Landscape Protection Index
        normalise_fields = [
        ("bu_bdy_firesize", "mfs"),
        ("bu_bdy_firesize_2h", "mfs2"),
        ("bu_bdy_hl", "mhl"),
        ("bu_bdy_firesize_2h_red", "mfs2r"),
        ("bu_bdy_hl_red", "mhlr"),
        ("control_line_freq", "clf"),
        ("burnability", None)]

        for source_field, output_field in normalise_fields:
            AnalysisDataset.normalise(source_field, dataset=new_grid, output_field=output_field)

        # calculate Containment Priority Index - containment score
        expr = "!mfs! + !mfs2! + !mhl!"
        arcpy.management.CalculateField(new_grid, "cpi_c_score", expr, "PYTHON")
        AnalysisDataset.normalise("cpi_c_score", dataset=new_grid)

        # calculate Containment Priority Index - suppression score
        expr = "!mfs2r! + !mhlr!"
        arcpy.management.CalculateField(new_grid, "cpi_s_score", expr, "PYTHON")
        AnalysisDataset.normalise("cpi_s_score", dataset=new_grid)

        # limit certain values to burn unit boundaries
        boundary_calculations = [
            ("bu_bdy_travel_time", "!wc_travel_time! * !is_bu_boundary!"),
            ("bu_bdy_freq", "!burn_freq! * !is_bu_boundary!")
        ]
        for field_name, expression in boundary_calculations:
            arcpy.management.CalculateField(new_grid, field_name, expression, "PYTHON")

        # calculate Containment Priority Index - travel time score
        AnalysisDataset.normalise("bu_bdy_travel_time", dataset=new_grid, output_field="cpi_t_score")

        # calculate Containment Priority Index - bushfire frequency score
        AnalysisDataset.normalise("bu_bdy_freq", dataset=new_grid, output_field="cpi_f_score")

        # final Containment Priority Index
        containment_priority_expr = "!cpi_c_score! + !cpi_s_score! + !cpi_t_score! + !cpi_f_score!"
        arcpy.management.CalculateField(new_grid, "containment_priority_index", containment_priority_expr, "PYTHON")
        AnalysisDataset.normalise("containment_priority_index", dataset=new_grid)

        # calculate Strategic Location Index 
        strategic_expr = "!clf! + !is_bmzlmz! + !is_bu_boundary! + !is_heavy_strategic! + !is_medium_strategic!/2 + !burnability!"
        arcpy.management.CalculateField(new_grid, "strategic_location_index", strategic_expr, "PYTHON")
        AnalysisDataset.normalise("strategic_location_index", dataset=new_grid)

        # calculate Bushfire Frequency Index
        AnalysisDataset.normalise("burn_freq", dataset=new_grid, output_field="bushfire_frequency_index")

        # LANDSCAPE PROTECTION INDEX
        landscape_expr = "!containment_priority_index! + !strategic_location_index! + !bushfire_frequency_index!"
        arcpy.management.CalculateField(new_grid, "landscape_protection_index", landscape_expr, "PYTHON")
        AnalysisDataset.normalise("landscape_protection_index", dataset=new_grid)

        #

        # restrict some scores to 1km buffer and renormalise 
        buffer_calculations = [
            ("mfs2r", "!mfs2r! * !is_1k_buff!"),
            ("mhlr", "!mhlr! * !is_1k_buff!"),
            ("clf", "!clf! * !is_1k_buff!")
        ]
        for field_name, expression in buffer_calculations:
            arcpy.management.CalculateField(new_grid, field_name, expression, "PYTHON")

        for field_name, _ in buffer_calculations:
            AnalysisDataset.normalise(field_name, dataset=new_grid)

        # calculate Fuel Management Index
        fuel_mgmt_expr = "(!is_apzbmz! + (!mfs2r! + !mhlr!)/2) * !is_1k_buff!"
        arcpy.management.CalculateField(new_grid, "fuel_management_index", fuel_mgmt_expr, "PYTHON")
        AnalysisDataset.normalise("fuel_management_index", dataset=new_grid)

        # calculate Strategic Protection Index
        expr = "(!clf! + !is_bu_boundary! + !is_heavy_strategic! + !is_medium_strategic!/2) * !is_1k_buff!"
        arcpy.management.CalculateField(new_grid, "strategic_protection_index", expr, "PYTHON")
        AnalysisDataset.normalise("strategic_protection_index", dataset=new_grid)

        # calculate Bushfire Risk Index
        AnalysisDataset.normalise("locality_hl", dataset=new_grid, output_field="loc_rr")
        arcpy.management.CalculateField(new_grid, "frq", "!burn_freq! * !is_1k_buff!", "PYTHON")
        AnalysisDataset.normalise("frq", dataset=new_grid)

        expr = "(!loc_rr! + !frq!) * !is_1k_buff!"
        arcpy.management.CalculateField(new_grid, "bushfire_risk_index", expr, "PYTHON")
        AnalysisDataset.normalise("bushfire_risk_index", dataset=new_grid)

        # calculate ASSET PROTECTION INDEX
        asset_protection_expr = "!fuel_management_index! + !strategic_protection_index! + !bushfire_risk_index!"
        arcpy.management.CalculateField(new_grid, "asset_protection_index", asset_protection_expr, "PYTHON")
        AnalysisDataset.normalise("asset_protection_index", dataset=new_grid)

        ## PRODUCE CLEAN OUTPUT FEATURES
        dissolve_fields = ["GRID_ID", "asset_protection_index", "landscape_protection_index", 
                           "containment_priority_index", "strategic_location_index", "bushfire_frequency_index",    # LPI Indices
                           "fuel_management_index", "strategic_protection_index", "bushfire_risk_index"]            # API Indices
        arcpy.management.Dissolve(
            in_features=new_grid,
            out_feature_class="evaluation_dataset",
            dissolve_field=dissolve_fields
        )

        print(f"\nFuel Break Evaluation Dataset generation completed successfully!")
        print(f"Dataset output: {os.path.join(OUTPATH, GDB_NAME)}\evaluation_dataset")
        
    except Exception as e:
        print(f"Error in main process: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()