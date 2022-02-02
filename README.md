#  Glacier Centerlines

A re-implementation of the [Kienholz et al., 2014](https://tc.copernicus.org/articles/8/503/2014/) algorithm in python.

This algorithm has been implemented succesfully [in OGGM](https://docs.oggm.org/en/stable/flowlines.html) since several years already. This code, however, has a few drawbacks:
- it requires OGGM to be installed to run
- the OGGM developpers have made certain choices with respect to algorithm design and grid handling which made sense for the model, a bit less for the centerline tool
- the agorithm code is a bit hidden in the large OGGM codebase, maybe impeding innovation and adaptations from the community.

The main goal of this project is to adress these shortcomings and develop a single tool with a simple purpose: given glacier outlines and a digital elevation domain, compute the centerlines of glaciers.

### Getting started with the implementation

Here are two files to test things on:
- https://cluster.klima.uni-bremen.de/~oggm/tutorials/Norway_Inventory_sel.zip
- https://cluster.klima.uni-bremen.de/~oggm/tutorials/Norway_DEM_sel.tif

There are no complex glaciers there, I'll provide a test case based on Baltoro glacier later. This is enough to get started though.

Workflow:
- assume that the DEM and the outlines are in a cartesian projection (units: m)
- open the geotiff with rioxarray
- start with one outline, and crop the DEM to the outline + a few grid point
- compute a mask of the glacier as in: https://github.com/OGGM/oggm/blob/447a49d7f936dae4870453d7c65bf2c6f861d0d8/oggm/core/gis.py#L798
- use the OGGM code to compute the heads, terminus, but: do not simplify geometries as done in OGGM. I would really try to see if its possible to work and compute glacier heads and terminus in the native geometry resolution. OGGM code: https://github.com/OGGM/oggm/blob/master/oggm/core/centerlines.py

I don't think there is a need for the OGGM Centerline object. All calculations should be doable with shapely only.

The tools you will need:

- rioxarray to read geotiff
- geopandas to read and write geometries
- shapely to do the geometrical stuff (as OGGM does)
- scipy for the routing algorithm (as OGGM does)

Good luck!

