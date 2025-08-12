"""
Fishnet.py

A Python class for ArcPy that generates hexagonal analysis grids and performs spatial data processing
for geographic analysis workflows.

Author:
Date: 
Version: 1.0

Dependencies:
- ArcGIS Pro with ArcPy
- Python 3+

Usage:
    from fishnet import Fishnet
    
    fishnet = Fishnet(
        workspace=r"C:\Data\analysis.gdb",
        gdb_name="analysis.gdb", 
        cell_size=1000000,  # 100 hectares in square metres
        shape="HEXAGON",
        extent=arcpy.Extent(2125000, 2259000, 2940000, 2828000, 3111),
        spatial_reference=arcpy.SpatialReference(3111)
    )
"""

import arcpy
import os
import sys
from typing import Tuple, Optional, List, Dict
import re
from numbers import Number as num


class Fishnet:
    """
    A class to handle fishnet generation and data processing operations.
    
    This class creates hexagonal or square analysis grids and performs complex
    spatial data processing operations by overlaying multiple datasets onto
    the standardised grid system.
    """
    
    def __init__(self, workspace: str, gdb_name: str, cell_size: int, shape: str = "HEXAGON",
                 extent: Optional[arcpy.Extent] = None,
                 spatial_reference: Optional[arcpy.SpatialReference] = None):
        """
        Initialise the Fishnet.
        
        Args:
            workspace (str): Path to the geodatabase directory
            gdb_name (str): Name for the output file geodatabase
            cell_size (int): Grid cell size in square metres
            shape (str, optional): Grid shape - "HEXAGON" or "SQUARE". Default: "HEXAGON"
            extent (arcpy.Extent, optional): Spatial extent for grid generation
            spatial_reference (arcpy.SpatialReference, optional): Coordinate system. Default: EPSG:3111
        """
        self.spatial_reference = spatial_reference or arcpy.SpatialReference(3111)  # Default to Vicgrid94
        
        # Create output directory and geodatabase
        self.output_path = workspace
        self.gdb_name = gdb_name
        os.makedirs(self.output_path, exist_ok=True)
        self.output_gdb = os.path.join(self.output_path, gdb_name)
        
        if arcpy.Exists(self.output_gdb):
            arcpy.AddMessage(f"Workspace exists: {self.output_gdb}")
        else:
            arcpy.CreateFileGDB_management(self.output_path, gdb_name, "10.0")
            arcpy.AddMessage(f"Created workspace: {self.output_gdb}")
            
        self.workspace = self.output_gdb
        self.layer_cache = []  # Cache for re-used layers to avoid re-projection

        # Set up ArcPy environment
        arcpy.env.workspace = self.workspace
        arcpy.env.overwriteOutput = True
        arcpy.env.cartographicCoordinateSystem = self.spatial_reference
        arcpy.env.outputCoordinateSystem = self.spatial_reference
        arcpy.env.parallelProcessingFactor = 80
        
        print(f"Initialised workspace: {self.workspace}")
    
        # Create analysis grid    
        try:
            # Output path
            grid_path = os.path.join(self.workspace, "grid_full")

            # Generate grid
            arcpy.management.GenerateTessellation(grid_path, extent, shape, cell_size)

            # Add fields for area and centroid coordinates
            fields_to_add = [
                ("CELL_AREA", "CELL_AREA"),
                ("CENTROID_X", "CENTROID_X"),
                ("CENTROID_Y", "CENTROID_Y"),
            ]
            
            for field_name, field_alias in fields_to_add:
                arcpy.management.AddField(grid_path, field_name, "DOUBLE", "", field_alias=field_alias)
            
            # Calculate field values
            with arcpy.da.UpdateCursor(grid_path, ["SHAPE@", "CELL_AREA", "CENTROID_X", "CENTROID_Y"]) as cursor:
                for row in cursor:
                    shape, cell_area, x, y = row
                    # Calculate new values
                    cell_area = shape.area / 10000  # Cell area in hectares
                    x = shape.centroid.X            # Centroid X
                    y = shape.centroid.Y            # Centroid Y
                    # Update row with calculated values
                    updated_row = (shape, cell_area, x, y)
                    cursor.updateRow(updated_row)
            
            # Add spatial index
            try:
                arcpy.AddSpatialIndex_management(grid_path)
            except:
                pass  # Index probably already exists

            print(f"Created analysis grid: {grid_path}")
            self.fishnet = grid_path
            self.cell_size = cell_size

            # Create in_memory version of fishnet
            self.fishnet_mem = "in_memory\\fishnet"
            arcpy.management.CopyFeatures(self.fishnet, self.fishnet_mem)
            
        except Exception as e:
            print(f"Error creating analysis grid: {str(e)}")
            raise
        
    def clip_to(self, features, query=None):
        """
        Remove grid cells that don't overlap supplied features, with optional query.
        
        Args:
            features (str): Path to clipping feature class
            query (str, optional): SQL query to filter clipping features
        """
        if query:
            # Apply selection query to selection features
            features = arcpy.management.SelectLayerByAttribute(features, "NEW_SELECTION", query)

        # Output path
        output_path = os.path.join(self.workspace, "grid_clip")

        # Select grid cells from fishnet
        grid_select = arcpy.management.SelectLayerByLocation(
            self.fishnet, "INTERSECT", features, selection_type="NEW_SELECTION")

        # Export selected cells to grid_clip
        arcpy.conversion.ExportFeatures(grid_select, output_path)

        # Re-assign fishnet
        self.fishnet = output_path

        # Update in_memory version of fishnet
        arcpy.management.CopyFeatures(self.fishnet, self.fishnet_mem)

        # Remove selections
        arcpy.management.SelectLayerByAttribute(features, "CLEAR_SELECTION")
        arcpy.management.SelectLayerByAttribute(self.fishnet, "CLEAR_SELECTION")

        print(f"Clipped analysis grid: {output_path}")
    
    def join_layer(self, dataset_name: str, config: dict) -> None:
        """
        Process a single spatial join based on configuration.
        
        The method uses the following flow:
        - Load the join data to memory
        - Process in_memory join data against fishnet_mem for efficiency
        - Join results back to GDB fishnet
        - Refresh fishnet_mem with a copy of the GDB fishnet
        
        Parameters:
            dataset_name (str): Identifier for the dataset being processed
            config (dict): Configuration dictionary with join parameters
        """
        
        # Extract configuration
        fieldname = config.get("fieldname", dataset_name)
        source_og = config["source"]
        match_option = config.get("match_option", "INTERSECT")
        method = config.get("method", "classify")
        group_by = config.get("group_by", None)
        source_field = config.get("source_field", None)
        query = config.get("filter", None)
        distance = config.get("distance", None)

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

            if source_sr.factoryCode != env_sr.factoryCode:  # Different SR - need to project
                arcpy.management.Project(source_og, source, env_sr)
                print(f"Reprojected {basename} and loaded to memory: {source}")
            else:  # Same SR - just copy for speed
                arcpy.management.CopyFeatures(source_og, source)
                print(f"Loaded {basename} to memory: {source}")
            
            # Add spatial indexing if not present (improves performance)
            try:
                arcpy.AddSpatialIndex_management(source)
            except:
                pass  # Index probably already exists
        else:
            print(f"Using {basename} from memory: {source}")
        
        # Apply filter if specified
        if query:
            source = arcpy.management.SelectLayerByAttribute(source, 'NEW_SELECTION', query)
            print(f"Applied query to {basename}: {query}")

        # Check dataset size and warn if large
        result = arcpy.management.GetCount(source)
        feature_count = int(result.getOutput(0))
        if feature_count > 50000:
            print(f"{basename} is large ({feature_count} features); Processing may be slow.")

        # Process based on method
        process_func = f"_process_{method}"
        self.process_func(source, fieldname, fishnet_id, match_option, distance, group_by, basename)
        
        # Clear any selection
        arcpy.management.SelectLayerByAttribute(source, 'CLEAR_SELECTION')

        # Update in_memory version of fishnet
        arcpy.management.CopyFeatures(self.fishnet, self.fishnet_mem)

    def _process_count(self, source, fieldname, fishnet_id, match_option, distance, source_field, basename):
        """Process count method with optional grouping."""
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
            case_fields = [fishnet_id, group_by]
            arcpy.analysis.Statistics(
                in_table=temp_join,
                out_table=temp_stats,
                statistics_fields=[[group_by, "COUNT"]],
                case_field=case_fields
            )
            
            # Count unique group_by values per fishnet cell
            final_stats = "in_memory\\final_stats"
            arcpy.analysis.Statistics(
                in_table=temp_stats,
                out_table=final_stats,
                statistics_fields=[[group_by, "COUNT"]],
                case_field=[fishnet_id]
            )
            
            # Join back to fishnet
            self.join_and_calc_field(
                target=self.fishnet,
                join_table=final_stats,
                target_field=fishnet_id,
                join_field=fishnet_id,
                transfer_field=f"COUNT_{group_by}",
                new_fieldname=fieldname,
                default_value=0
            )
            
            print(f"Populated fishnet[{fieldname}] with counts from {basename} grouped by {group_by}")
        else:
            # Simple count of all intersecting features
            temp_stats = "in_memory\\temp_stats"
            arcpy.analysis.Statistics(
                in_table=temp_join,
                out_table=temp_stats,
                statistics_fields=[["Join_Count", "SUM"]],
                case_field=[fishnet_id]
            )
            
            # Join back to fishnet
            self.join_and_calc_field(
                target=self.fishnet,
                join_table=temp_stats,
                target_field=fishnet_id,
                join_field=fishnet_id,
                transfer_field="SUM_Join_Count",
                new_fieldname=fieldname,
                default_value=0
            )

            print(f"Populated fishnet[{fieldname}] with counts from {basename}")

    def _process_classify(self, source, fieldname, fishnet_id, match_option, distance, source_field, basename):
        """Process classify method."""
        if not source_field:
            raise ValueError("No field identified for classify")
        
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
            transfer_field=source_field,
            new_fieldname=fieldname,
            default_value=None
        )

        print(f"Populated fishnet[{fieldname}] with values from {basename}[{source_field}]")

    def _process_sum(self, source, fieldname, fishnet_id, match_option, distance, source_field, basename):
        """Process sum method."""
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
            transfer_field=f"SUM_{source_field}",
            new_fieldname=fieldname,
            default_value=0
        )

        print(f"Populated fishnet[{fieldname}] with sum values from {basename}[{source_field}]")

    def _process_area(self, source, fieldname, fishnet_id, match_option, distance, source_field, basename):
        """Process area method."""
        temp_join = "in_memory\\temp_join"
        arcpy.analysis.PairwiseIntersect([self.fishnet_mem, source], out_feature_class=temp_join)

        # Add hectare field and populate
        arcpy.management.AddField(
            in_table=temp_join,
            field_name="Area_Ha",
            field_type="FLOAT",
            field_scale=2
        )
        
        arcpy.management.CalculateField(temp_join, "Area_Ha", "!SHAPE@AREA!/10000", "PYTHON3")

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

        print(f"Populated fishnet[{fieldname}] with area values from {basename}")

    def _process_average(self, source, fieldname, fishnet_id, match_option, distance, source_field, basename):
        """Process average method."""
        temp_join = "in_memory\\temp_join"
        arcpy.analysis.PairwiseIntersect([self.fishnet_mem, source], out_feature_class=temp_join)
        
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
            transfer_field=f"MEAN_{source_field}",
            new_fieldname=fieldname,
            default_value=0
        )

        print(f"Populated fishnet[{fieldname}] with average values from {basename}[{source_field}]")

    def _process_area_weighted_average(self, source, fieldname, fishnet_id, match_option, distance, source_field, basename):
        temp_join = "in_memory\\temp_join"
        arcpy.analysis.PairwiseIntersect([self.fishnet_mem, source], out_feature_class=temp_join)
        
        # Calculate weighted value field (source_field * area)
        weighted_field = f"WEIGHTED_{source_field}"
        arcpy.management.AddField(temp_join, weighted_field, "DOUBLE")
        arcpy.management.CalculateField(
            temp_join, weighted_field, 
            f"!{source_field}! * !shape.area@squaremeters!", "PYTHON3"
        )
        
        temp_stats = "in_memory\\temp_stats"
        arcpy.analysis.Statistics(
            in_table=temp_join,
            out_table=temp_stats,
            statistics_fields=[[weighted_field, "SUM"]],
            case_field=[fishnet_id]
        )
        
        # Calculate area-weighted average
        weighted_avg_field = "AREA_WEIGHTED_AVG"
        arcpy.management.AddField(temp_stats, weighted_avg_field, "DOUBLE")
        arcpy.management.CalculateField(
            temp_stats, weighted_avg_field, 
            f"!SUM_{weighted_field}! / {self.cell_size}", "PYTHON3"
        )
        
        # Join back to fishnet
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

    def _process_presence(self, source, fieldname, fishnet_id, match_option, distance, source_field, basename):
        """Process presence method."""
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
            selection_type="NEW_SELECTION"
        )

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

    def normalise(self, field_name: str, dataset: Optional[str] = None, group_by: Optional[str] = None,
                  method: Optional[str] = "ZERO", output_field: Optional[str] = None) -> None:
        """
        Normalise field values between 0 and 1.
        
        Args:
            field_name (str): Name of field to be normalised
            dataset (str, optional): Feature class. Defaults to self.fishnet
            group_by (str, optional): Grouping field name for grouped normalisation
            method (str, optional): "ZERO" (default) or "MIN" normalisation method
            output_field (str, optional): Output field name. Default overwrites input field
        """
        
        cursor_fields = [field_name]

        # Set dataset to self.fishnet if not provided
        if not dataset:
            dataset = self.fishnet
        
        # Set target to output_field if provided
        if output_field:
            target_field = output_field
            cursor_fields.append(output_field)
            
            if target_field not in [f.name for f in arcpy.ListFields(dataset)]:
                arcpy.management.AddField(dataset, target_field, "FLOAT")
        else:
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
            
            # Iterate cursor over stats table
            with arcpy.da.SearchCursor(stats_table, [group_by, f"MIN_{field_name}", f"MAX_{field_name}"]) as stats_cursor:
                for row in stats_cursor:
                    group = row[0]
                    min_val = 0 if method == "ZERO" else row[1]
                    max_val = row[2]

                    # Format where_clause to handle string or int data types
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
                arcpy.management.CalculateField(dataset, target_field, "0", "PYTHON3")
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

    # Helper Methods
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
        """
        Validate configuration dictionary for common issues.
        
        Args:
            config (dict): Configuration dictionary to validate
        """
        # Define valid values for each parameter
        VALID_METHODS = ["presence", "count", "sum", "average", "area_weighted_average", 
                        "maximum", "minimum", "area", "classify"]
        
        VALID_MATCH_OPTS = ["INTERSECT", "WITHIN_A_DISTANCE", "CONTAINS", "COMPLETELY_CONTAINS", 
                           "WITHIN", "CLOSEST", "COMPLETELY_WITHIN", "HAVE_THEIR_CENTER_IN", 
                           "LARGEST_OVERLAP", "CROSSED_BY_THE_OUTLINE_OF"]
        
        METHODS_REQUIRING_FIELD = ["classify", "sum", "average", "mean", "maximum", "minimum"]
        
        issues = []
        warnings = []
        
        # Check required fields
        if "source" not in config:
            issues.append("Missing required 'source' field")
        else:
            # Check if source exists
            if not arcpy.Exists(config["source"]):
                issues.append(f"Source does not exist: {config['source']}")

        # Check method validity
        method = config.get("method")
        if method not in VALID_METHODS:
            issues.append(f"Invalid method '{method}'. Valid methods: {', '.join(VALID_METHODS)}")
        
        # Check match_option validity
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
                    issues.append(f"Field '{config['source_field']}' not found in source. "
                                f"Available fields: {', '.join(field_names[:10])}")

        # Print report
        fieldname = config.get("fieldname", "Unknown")
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

    def cleanup_memory(self):
        """Clean up in-memory datasets to free memory."""
        for layer in self.layer_cache:
            try:
                if arcpy.Exists(layer):
                    arcpy.management.Delete(layer)
            except:
                pass
        self.layer_cache.clear()
        
        # Clean up fishnet memory copy
        try:
            if arcpy.Exists(self.fishnet_mem):
                arcpy.management.Delete(self.fishnet_mem)
        except:
            pass
        
        print("Cleaned up in-memory datasets")

    def get_summary_stats(self) -> dict:
        """
        Get summary statistics for the current fishnet grid.
        
        Returns:
            dict: Summary statistics including cell count, total area, etc.
        """
        try:
            # Get basic counts
            result = arcpy.management.GetCount(self.fishnet)
            cell_count = int(result.getOutput(0))
            
            # Get area statistics
            total_area = 0
            with arcpy.da.SearchCursor(self.fishnet, ["CELL_AREA"]) as cursor:
                for row in cursor:
                    if row[0] is not None:
                        total_area += row[0]
            
            # Get field list
            fields = [f.name for f in arcpy.ListFields(self.fishnet) 
                     if f.type in ['Double', 'Float', 'Integer', 'SmallInteger'] 
                     and f.name not in ['OBJECTID', 'GRID_ID', 'CELL_AREA', 'CENTROID_X', 'CENTROID_Y']]
            
            stats = {
                'cell_count': cell_count,
                'total_area_ha': round(total_area, 2),
                'cell_size_ha': round(self.cell_size / 10000, 2),
                'data_fields': len(fields),
                'field_names': fields,
                'workspace': self.workspace,
                'fishnet_path': self.fishnet
            }
            
            return stats
            
        except Exception as e:
            print(f"Error getting summary stats: {e}")
            return {}

    def export_to_shapefile(self, output_path: str, fields: Optional[List[str]] = None):
        """
        Export fishnet to shapefile format.
        
        Args:
            output_path (str): Full path for output shapefile
            fields (List[str], optional): List of fields to include. If None, includes all fields.
        """
        try:
            if fields:
                # Create field mappings to include only specified fields
                field_mappings = arcpy.FieldMappings()
                for field_name in fields:
                    if field_name in [f.name for f in arcpy.ListFields(self.fishnet)]:
                        field_map = arcpy.FieldMap()
                        field_map.addInputField(self.fishnet, field_name)
                        field_mappings.addFieldMap(field_map)
                
                # Export with field mapping
                arcpy.conversion.FeatureClassToFeatureClass(
                    in_features=self.fishnet,
                    out_path=os.path.dirname(output_path),
                    out_name=os.path.basename(output_path),
                    field_mapping=field_mappings
                )
            else:
                # Export all fields
                arcpy.conversion.FeatureClassToFeatureClass(
                    in_features=self.fishnet,
                    out_path=os.path.dirname(output_path),
                    out_name=os.path.basename(output_path)
                )
            
            print(f"Exported fishnet to: {output_path}")
            
        except Exception as e:
            print(f"Error exporting to shapefile: {e}")
            raise