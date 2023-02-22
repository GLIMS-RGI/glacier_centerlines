# This is a dummy init.
__version__ = '1.0.1'

# API
# Some decorators used by many
from oggm.utils import entity_task, global_task

#add functions
from glacier_centerlines import functions_rgi

# Classes
from glacier_centerlines.centerlines_rgi import compute_centerlines_rgi

