#code snippet to try glaciercenterlines_rgi centerlines, that is now pip installable (pip install glacier_centelrines) 

#code available in https://github.com/GLIMS-RGI/glacier_centerlines/tree/dev-entity_task

#load libraries and oggm dependencies
import os
from oggm import cfg, utils
cfg.initialize(logging_level='WARNING')
from oggm import workflow

#set tmpdir
cfg.PATHS['working_dir'] = utils.gettempdir(dirname='OGGM-new-cls', reset=False)

# try hintereisferner, our beloved refference glacier in Tirol:
rgi_ids = ['RGI60-11.00897']

# bese url:
base_url = ('https://cluster.klima.uni-bremen.de/~oggm/gdirs/oggm_v1.6/L3-L5_files/centerlines/w5e5/qc0/pcpwin/match_geod_pergla/')

# Initialize glacier directories
gdirs = workflow.init_glacier_directories(rgi_ids, from_prepro_level=3, prepro_border=80, prepro_base_url=base_url)

#import general execution for oggm taks
from oggm.workflow import execute_entity_task

#new package where to take the task from:
import glacier_centerlines as gc

# run
execute_entity_task(gc.centerlines_rgi.compute_centerlines_rgi, gdirs)
