# Configuration file for OGGM parameters

### not sure: add description
is_first_call = False

# declare general paths
out_path = "/home/francesc/results/glacier_centerlines/"

# Baltoro test data
data_path = "./test_data/"
dem_file = "dem_balto_tm.tif"
shape_file = "outlines.shp"

# Norway test data
#dem_file = "Norway_DEM_sel.tif"
#data_path = "/home/francesc/data/glacier_centerlines/"
#shape_file = "Norway_Inventory_sel/Norway_Inventory_sel.shp"

# RGI11
#data_path = "/home/francesc/data/OGGM/rgi/RGIV60/11_rgi60_CentralEurope/"
#dem_file = 'dummy'#'RGI60-11.00/RGI60-11.00002/NASADEM/dem.tif'
#shape_file = "11_rgi60_CentralEurope.shp"

#data_path = "/home/francesc/repositories/glacier_centerlines/test_data/"
#dem_file = 'RGI60-11.00/RGI60-11.00002/NASADEM/dem.tif'
#shape_file = "11_rgi60_CentralEurope.shp"


### Input/Output paths. Set to ~ to default to home directory

# Where OGGM will write its output
#working_dir =

# Users can specify their own topography file if they want to.
# This is useful for testing, or if you
# are simulating a single region with better data.
# the empty default is what most users should do
#dem_file =

# Use compression for the intermediate pickles? (might slow down I/O a bit)
# Both the performance loss (0% ?) and the space gain (-10%) seem to be low
use_compression = True

### CENTERLINE determination

# Decision on grid spatial resolution for each glacier
# 'fixed': dx (meters) = fixed_dx
# 'linear':  dx (meters) = d1 * AREA (km) + d2 ; clipped to dmax (e.g.: 5, 10, 200)
# 'square':  dx (meters) = d1 * sqrt(AREA) (km) + d2 ; clipped to dmax (e.g.: 20, 10, 200)

# Was default for a long time
# grid_dx_method = 'linear'
# d1 = 5.
# d2 = 10.
# dmax = 100.

# New default?
grid_dx_method = 'square'
d1 = 14.
d2 = 10.
dmax = 200.

# Ignored if grid_dx_method != 'fixed'
fixed_dx = 50.

# Which algorithm to use for interpolating the topography to the local grid
# 'bilinear' or 'cubic'
topo_interp = 'cubic'

# Grid border buffer around the glacier (in pixels)
# Make it large if you want to do past simulations.
border = 40

# For tidewater glaciers it doesn't make sense to have large maps
# if for some reason you still want this, set to false
clip_tidewater_border = True

# The glacier area, CenLon and CenLat are usually taken from the RGI
# shapefile, which is a good thing for default RGI files. If you use your
# own inventory, however, it might be a good idea to let OGGM compute these
# attributes at runtime: set to `False` in this case.
use_rgi_area = True

# Head determination: (approx) size in meters of the half-size window
# where to look for maxima
localmax_window = 500. #In units of dx

# DEM smoothing: (approx) size in meters of the smoothing window.
# Set to 0 for no smoothing
smooth_window = 251.


# Kienholz et al eq (1)
q1 = 2/10**6 # 1/m
q2 = 500 #m
rmax = 1000 #m

q1 = 2e-6
q2 = 500.
rmax = 1000.

# Kienholz et al eq (2) --> it is hardcoded in _make_costgrid() function!
f1 = 1000.
f2 = 3000.
a = 4.25 #4.25 in literature
b = 3.7

# Kienholz et al eq (8) but modified here
# Buffer in pixels where to cut the incoming centerlines
kbuffer = 1


# For water-terminating glaciers, use the percentile instead of minimum h?
# Set to zero if no special treatment for water terminating glaciers should be
# used, and to an integer > 0 to specify the percentile
terminus_search_percentile = 10
terminus_search_altitude_range = 20

### FLOWLINES definition parameters
# Whether the model should use the glacier intersects information
# given by the user
use_intersects = True
# Grid spacing of a flowline in pixel coordinates
flowline_dx = 2

# Number of pixels to arbitrarily remove at junctions
flowline_junction_pix = 3
# Gaussian smooth of the altitude along a flowline
# sigma, in pixel coordinates (sigma=1 -> smooth around a -4:+4 window)
flowline_height_smooth = 1

# Prevent too small slopes? (see also min_slope param below)
filter_min_slope = True

### Elevation band flowlines (or "collapsed flowlines") parameters
# Only used if using the alternative flowline definition
# The elevation binsize in m - it was 10m in Huss&Farinotti2012, 30m in Werder 2019
elevation_band_flowline_binsize = 30

### CATCHMENT WIDTHS computation parameters
# altitude range threshold for filtering
# This stuff has not been really optimized, it's also not very critical
width_alt_range_thres = 250.
# Minimum number of elements per bin for altitude-binsize definition
min_n_per_bin = 2
# Baseline binsize for the altitude-area distribution
base_binsize = 50.
# Smoothing of the widths after altitude-area matching? 0 means no smoothing,
# 1 means default (i.e. kernel size 9).
smooth_widths_window_size = 1

### INVERSION params
# Clip the flowline slope, in degrees
# This will crop the slope during the ice thickness inversion.
# This is quite a sensitive parameter!
min_slope = 1.5
min_slope_ice_caps = 1.5
# When converting the centerlines to flowlines, we prevent negative slopes.
# Ideally, this value should be set to `min_slope` for physical consistency,
# but it turns that many flat glaciers will have weird flowlines with this
# setting. Using zero works ok, and was the default in OGGM for long
min_slope_flowline_filter = 0

### File output options

# Whether to store the model geometry files during operational runs
# This can be useful for advanced applications needing the retrieval
# of glacier geometries after the run, but is not necessary if you
# are interested in diagnostics only (volume, length, etc.)
store_model_geometry = False

# Single flow: if True, a single, principal flowline is computed. 
# If False, tributaries are also computed.
single_fl = False
    
### Plots 
# Produce plots in the terminal? True = y; False = n
plot = False