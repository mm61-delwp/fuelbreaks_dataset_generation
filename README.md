# Fuel breaks dataset generation script
Arcpy script to generate analysis dataset for fuel breaks



# Fishnet Class Documentation
The provided script was prepared for a specific purpose - to produce an evaluation dataset for fuel breaks on forested land.
Its engine is flexible though and may have other uses, so documentation is provided below.

## Table of Contents

- [Overview](#overview)
- [Installation & Requirements](#installation--requirements)
- [Quick Start](#quick-start)
- [Class Reference](#class-reference)
- [Configuration System](#configuration-system)
- [Spatial Join Methods](#spatial-join-methods)
- [Example](#complete-workflow-example)
- [Best Practices](#best-practices)

## Overview

The `Fishnet` class automates the creation of the analysis grid and performs complex spatial data processing operations. It's designed for analysing spatial relationships between different geographic datasets by overlaying them onto a standardised grid system.

### Key Features

- **Grid Generation**: Creates tessellated grid (including hexagonal but flexible) for spatial analysis
- **Flexible Spatial Joins**: Multiple methods for joining spatial data (presence, count, sum, average, etc.)
- **Data Matrix Processing**: Configuration-driven approach for processing multiple datasets
- **Performance Optimised**: Uses in-memory processing and spatial indexing
- **Coordinate System Management**: Automatic reprojection and spatial reference handling

## Installation & Requirements

### Dependencies
- ArcGIS Pro with ArcPy
- Python 3+
- Write access to output directory

### Import
```python
import arcpy
from fishnet import Fishnet  # If the class is in fishnet.py, otherwise can just be pasted into script
```

## Quick Start

```python
# Basic setup
fishnet = Fishnet(
    workspace=r"C:\Data\analysis.gdb",
    gdb_name="analysis.gdb", 
    cell_size=1000000,  # Cell size in square meters (100 hectares)
    shape="HEXAGON",
    extent=arcpy.Extent(2125000, 2259000, 2940000, 2828000, 3111),
    spatial_reference=arcpy.SpatialReference(3111)
)

# Clip grid to analysis area
fishnet.clip_to(r"C:\Data\study_area.shp")

# Process data using configuration
config = {
    "fieldname": "fire_stations",
    "source": r"C:\Data\fire_stations.shp",
    "method": "count",
    "match_option": "INTERSECT"
}
fishnet.join_layer("fire_stations", config)
```

## Class Reference

### Constructor

```python
Fishnet(workspace, gdb_name, cell_size, shape="HEXAGON", extent=None, spatial_reference=None)
```

**Parameters:**
- `workspace` (str): Path to output geodatabase directory
- `gdb_name` (str): Name for the output file geodatabase
- `cell_size` (int): Grid cell size in square meters
- `shape` (str, optional): Grid shape - "HEXAGON" or "SQUARE". Default: "HEXAGON"
- `extent` (arcpy.Extent, optional): Spatial extent for grid generation
- `spatial_reference` (arcpy.SpatialReference, optional): Coordinate system. Default: EPSG:3111 (VicGrid94)

**Returns:** Fishnet instance with generated grid

### Methods

#### `clip_to(features, query=None)`

Clips the analysis grid to intersecting features, reducing processing area.

**Parameters:**
- `features` (str): Path to feature class for clipping boundary
- `query` (str, optional): SQL query to filter clipping features

**Example:**
```python
# Clip to administrative boundary
fishnet.clip_to(r"C:\Data\admin_boundaries.shp", "STATE_NAME = 'Victoria'")
```

#### `join_layer(dataset_name, config)`

Performs spatial join operations based on configuration dictionary.

**Parameters:**
- `dataset_name` (str): Identifier for the dataset being processed
- `config` (dict): Configuration dictionary (see [Configuration System](#configuration-system))

#### `normalise(field_name, dataset=None, group_by=None, method="ZERO", output_field=None)`

Normalises field values to 0-1 range for analysis.

**Parameters:**
- `field_name` (str): Name of field to normalise
- `dataset` (str, optional): Target dataset. Default: fishnet grid
- `group_by` (str, optional): Field for grouped normalisation
- `method` (str, optional): "ZERO" (default) or "MIN" normalisation method
- `output_field` (str, optional): Output field name. Default: overwrites input field

**Example:**
```python
# Normalise fire frequency by fire management zone
fishnet.normalise("fire_frequency", group_by="fire_zone", output_field="fire_freq_norm")
```

#### `check_config(config)`

Provides check for common configuration errors and prints a report of issues to console if found.

**Parameters:**
- `config` (dict): Name configuration dictionary

**Example:**
```python
config = {
    "fieldname": "has_roads",
    "source": r"C:\Data\roads.shp", 
    "method": "presence",
    "match_option": "INTERSECTING"
}
# Check for issues in configuration
fishnet.check_config(config)
```
The function will print a short report stating that INTERSECTING isn't a valid match option (should be INTERSECT) etc.

## Configuration System

The data processing system uses configuration dictionaries to define spatial join operations. Each dataset is processed according to its configuration parameters.

### Required Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `fieldname` | str | Name of field to create in fishnet grid |
| `source` | str | Full path to source dataset |
| `method` | str | Processing method (see [Spatial Join Methods](#spatial-join-methods)) |

### Optional Parameters

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `source_field` | str | Field name for numeric operations | None |
| `match_option` | str | Spatial relationship type | "INTERSECT" |
| `group_by` | str | Field for grouping operations | None |
| `filter` | str | SQL filter expression | None |
| `distance` | int | Buffer distance in meters | None |

### Match Options

Valid spatial relationship types:

- `INTERSECT` - Features intersect
- `WITHIN_A_DISTANCE` - Features within specified distance
- `CONTAINS` - Target contains join features
- `COMPLETELY_CONTAINS` - Target completely contains join features  
- `WITHIN` - Target within join features
- `CLOSEST` - Nearest feature
- `COMPLETELY_WITHIN` - Target completely within join features
- `HAVE_THEIR_CENTER_IN` - Target centroid within join features
- `LARGEST_OVERLAP` - Feature with largest overlap area
- `CROSSED_BY_THE_OUTLINE_OF` - Target boundary intersects join features

## Spatial Join Methods

### `presence`
Returns 1 if any intersecting features exist, 0 otherwise.

```python
config = {
    "fieldname": "has_roads",
    "source": r"C:\Data\roads.shp", 
    "method": "presence",
    "match_option": "INTERSECT"
}
```

### `count`
Returns count of intersecting features.

```python
config = {
    "fieldname": "num_buildings",
    "source": r"C:\Data\buildings.shp",
    "method": "count",
    "match_option": "INTERSECT"
}
```

### `count` with grouping
Returns count of unique values in specified field.

```python
config = {
    "fieldname": "fire_seasons", 
    "source": r"C:\Data\fire_history.shp",
    "method": "count",
    "group_by": "SEASON",
    "match_option": "INTERSECT"
}
```

### `sum`
Returns sum of numeric field values from intersecting features.

```python
config = {
    "fieldname": "total_population",
    "source": r"C:\Data\census_blocks.shp",
    "method": "sum", 
    "source_field": "POPULATION",
    "match_option": "INTERSECT"
}
```

### `average`
Returns mean of numeric field values from intersecting features.

```python
config = {
    "fieldname": "avg_elevation",
    "source": r"C:\Data\elevation_points.shp",
    "method": "average",
    "source_field": "ELEVATION", 
    "match_option": "INTERSECT"
}
```

### `area_weighted_average`
Returns area-weighted average of field values, useful for continuous surfaces.

```python
config = {
    "fieldname": "weighted_risk",
    "source": r"C:\Data\risk_zones.shp", 
    "method": "area_weighted_average",
    "source_field": "RISK_VALUE",
    "match_option": "INTERSECT"
}
```

### `area`
Returns intersection area in hectares.

```python
config = {
    "fieldname": "forest_area_ha",
    "source": r"C:\Data\forest_cover.shp",
    "method": "area",
    "match_option": "INTERSECT" 
}
```

### `classify`
Transfers attribute value from intersecting feature (one-to-one relationship).

```python
config = {
    "fieldname": "land_use_type",
    "source": r"C:\Data\land_use.shp",
    "method": "classify",
    "source_field": "USE_CODE",
    "match_option": "LARGEST_OVERLAP"
}
```

## Complete Workflow Example

```python
import arcpy
import os

# Configuration constants
OUTPATH = r"C:\Analysis\Output"
GDB_NAME = "wildfire_analysis.gdb" 
SPATIAL_REF = arcpy.SpatialReference(3111)  # VicGrid94
GRID_EXTENT = arcpy.Extent(2125000, 2259000, 2940000, 2828000, 3111)
CELL_HA = 100  # 100 hectare cells
CELL_M = CELL_HA * 10000

# Data processing configuration
data_matrix = {
    "fire_stations": {
        "fieldname": "fire_station_count",
        "source": r"C:\Data\fire_stations.shp",
        "method": "count",
        "match_option": "WITHIN_A_DISTANCE",
        "distance": 5000  # 5km buffer
    },
    "population_density": {
        "fieldname": "pop_density", 
        "source": r"C:\Data\census_areas.shp",
        "method": "area_weighted_average",
        "source_field": "POP_DENSITY",
        "match_option": "INTERSECT"
    },
    "fire_history": {
        "fieldname": "fire_frequency",
        "source": r"C:\Data\fire_perimeters.shp", 
        "method": "count",
        "match_option": "INTERSECT",
        "filter": "FIRE_YEAR >= 2010"
    },
    "vegetation_type": {
        "fieldname": "veg_class",
        "source": r"C:\Data\vegetation.shp",
        "method": "classify", 
        "source_field": "VEG_CODE",
        "match_option": "LARGEST_OVERLAP"
    }
}

# Initialise fishnet
fishnet = Fishnet(
    workspace=OUTPATH,
    gdb_name=GDB_NAME,
    cell_size=CELL_M,
    shape="HEXAGON", 
    extent=GRID_EXTENT,
    spatial_reference=SPATIAL_REF
)

# Clip to study area 
study_area = r"C:\Data\region_boundary.shp"
fishnet.clip_to(study_area)

# Process all datasets
for dataset_name, config in data_matrix.items():    
    print(f"Processing {dataset_name}...")

    # validate configuration
    fishnet.check_config(config)

    try:
        fishnet.join_layer(dataset_name, config)
    except Exception as e:
        print(f"Error processing {dataset_name}: {e}")

# Normalise risk factors
fishnet.normalise("fire_frequency", output_field="fire_freq_norm")
fishnet.normalise("pop_density", output_field="pop_density_norm") 

print(f"Analysis complete. Results saved to: {fishnet.fishnet}")
```

## Best Practices

### Performance Optimisation

1. **Use appropriate cell sizes**: Larger cells = faster processing, but less spatial detail
2. **Clip early**: Use `clip_to()` to reduce processing area before joining data
3. **Filter datasets**: Use SQL filters to reduce feature counts before processing
4. **Choose efficient match options**: `INTERSECT` is typically faster than `LARGEST_OVERLAP`

### Data Quality

1. **Validate configurations**: The class includes built-in validation - check error messages
2. **Handle coordinate systems**: Ensure all datasets use compatible projections
3. **Check field types**: Verify numeric fields for mathematical operations
4. **Test with small areas**: Validate workflows on subsets before full processing
