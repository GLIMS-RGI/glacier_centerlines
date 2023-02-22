#  Glacier Centerlines

A re-implementation of the [Kienholz et al., 2014](https://tc.copernicus.org/articles/8/503/2014/) algorithm in python.

This algorithm has been implemented succesfully [in OGGM](https://docs.oggm.org/en/stable/flowlines.html) since several years already. This code, however, has a few drawbacks:
- it requires OGGM to be installed to run.
- the OGGM developpers have made certain choices with respect to algorithm design and grid handling which made sense for the model, a bit less for the centerline tool
- the agorithm code is a bit hidden in the large OGGM codebase, maybe impeding innovation and adaptations from the community.

The main goal of this project is to adress these shortcomings and develop a single tool with a simple purpose: given glacier outlines and a digital elevation domain, compute the centerlines of glaciers.

### Getting started with the implementation

Here are two files to test things on:
- https://cluster.klima.uni-bremen.de/~oggm/tutorials/Norway_Inventory_sel.zip
- https://cluster.klima.uni-bremen.de/~oggm/tutorials/Norway_DEM_sel.tif

Workflow:
- assume that the DEM and the outlines are in a cartesian projection (units: m)
- open the geotiff with rioxarray
- start with one outline, and crop the DEM to the outline + a few grid point
- compute a mask of the glacier.
- use the OGGM code to compute the heads, terminus, without simplifying geometries as done in OGGM. 

The tools you will need:
- rioxarray to read geotiff
- geopandas to read and write geometries
- shapely to do the geometrical stuff (as OGGM does)
- scipy for the routing algorithm (as OGGM does)

# UPDATE, version V1.0.2:
See documentation in `\docs`

# UPDATE pip installation
Now the tool is [pip installable](https://pypi.org/project/glacier-centerlines/). 
The tool has become an entity task that can be called using oggm already existing functions:

(e.g. from `snippet_run_rgi_centerlines.py`)
```
#import general execution for oggm taks
from oggm.workflow import execute_entity_task

#new package where to take the task from:
import glacier_centerlines as gc

# run
execute_entity_task(gc.centerlines_rgi.compute_centerlines_rgi, gdirs)

``` 

The new centerlines will be stored at the original oggm glacier directory as `centerlines_rgi.tar.gz
`
