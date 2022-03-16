"""Some useful functions that did not fit into the other modules.
"""
import netCDF4
import shapely
import shapely.geometry as shpg
import numpy as np
import copy

class ncDataset(netCDF4.Dataset):
    """Wrapper around netCDF4 setting auto_mask to False"""

    def __init__(self, *args, **kwargs):
        super(ncDataset, self).__init__(*args, **kwargs)
        self.set_auto_mask(False)

def nicenumber(number, binsize, lower=False):
    """Returns the next higher or lower "nice number", given by binsize.
    Examples:
    ---------
    >>> nicenumber(12, 10)
    20
    >>> nicenumber(19, 50)
    50
    >>> nicenumber(51, 50)
    100
    >>> nicenumber(51, 50, lower=True)
    50
    """

    e, _ = divmod(number, binsize)
    if lower:
        return e * binsize
    else:
        return (e + 1) * binsize
    
def clip_scalar(value, vmin, vmax):
    """A faster numpy.clip ON SCALARS ONLY.
    See https://github.com/numpy/numpy/issues/14281
    """
    return vmin if value < vmin else vmax if value > vmax else value

def _filter_heads(heads, heads_height, radius, polygon):
    """Filter the head candidates following Kienholz et al. (2014), Ch. 4.1.2
    Parameters
    ----------
    heads : list of shapely.geometry.Point instances
        The heads to filter out (in raster coordinates).
    heads_height : list
        The heads altitudes.
    radius : float
        The radius around each head to search for potential challengers
    polygon : shapely.geometry.Polygon class instance
        The glacier geometry (in raster coordinates).
    Returns
    -------
    a list of shapely.geometry.Point instances with the "bad ones" removed
    """

    heads = copy.copy(heads)
    heads_height = copy.copy(heads_height)

    i = 0
    # I think a "while" here is ok: we remove the heads forwards only
    while i < len(heads):
        head = heads[i]
        pbuffer = head.buffer(radius)
        inter_poly = pbuffer.intersection(polygon.exterior)
        if inter_poly.type in ['MultiPolygon',
                               'GeometryCollection',
                               'MultiLineString']:
            #  In the case of a junction point, we have to do a check
            # http://lists.gispython.org/pipermail/community/
            # 2015-January/003357.html
            if inter_poly.type == 'MultiLineString':
                inter_poly = shapely.ops.linemerge(inter_poly)

            if inter_poly.type != 'LineString':
                # keep the local polygon only
                for sub_poly in inter_poly.geoms:
                    if sub_poly.intersects(head):
                        inter_poly = sub_poly
                        break
        elif inter_poly.type == 'LineString': #### i have in treoduced "tuple()"in here
            inter_poly = shpg.Polygon(tuple(np.asarray(inter_poly.xy).T))
        elif inter_poly.type == 'Polygon':
            pass
        else:
            extext = ('Geometry collection not expected: '
                      '{}'.format(inter_poly.type))
            #raise InvalidGeometryError(extext)

        # Find other points in radius and in polygon
        _heads = [head]
        _z = [heads_height[i]]
        for op, z in zip(heads[i+1:], heads_height[i+1:]):
            if inter_poly.intersects(op):
                _heads.append(op)
                _z.append(z)

        # If alone, go to the next point
        if len(_heads) == 1:
            i += 1
            continue

        # If not, keep the highest
        _w = np.argmax(_z)

        for head in _heads:
            if not (head is _heads[_w]):
                heads_height = np.delete(heads_height, heads.index(head))
                heads.remove(head)

    return heads, heads_height

def find_nearest(array, value):
    array = np. asarray(array)
    idx = (np. abs(array - value)). argmin()    
    return array[idx]

class SuperclassMeta(type):
    """Metaclass for abstract base classes.
    http://stackoverflow.com/questions/40508492/python-sphinx-inherit-
    method-documentation-from-superclass
    """
    def __new__(mcls, classname, bases, cls_dict):
        cls = super().__new__(mcls, classname, bases, cls_dict)
        for name, member in cls_dict.items():
            if not getattr(member, '__doc__'):
                try:
                    member.__doc__ = getattr(bases[-1], name).__doc__
                except AttributeError:
                    pass
        return cls

# define classes
from utils import SuperclassMeta


class glacier_dir(object):
    def __init__(self, grid):
        self.grid = grid
    
class grid_inf(object):
    def __init__(self, grid):
        self.grid = grid

class Centerline(object, metaclass=SuperclassMeta):
    """Geometry (line and widths) and flow rooting properties, but no thickness
    """

    def __init__(self, line, dx=None, surface_h=None, orig_head=None,
                 rgi_id=None, map_dx=None):
        """ Initialize a Centerline
        Parameters
        ----------
        line : :py:class:`shapely.geometry.LineString`
            The geometrically calculated centerline
        dx : float
            Grid spacing of the initialised flowline in pixel coordinates
        surface_h :  :py:class:`numpy.ndarray`
            elevation [m] of the points on ``line``
        orig_head : :py:class:`shapely.geometry.Point`
            geometric point of the lines head
        rgi_id : str
            The glacier's RGI identifier
        map_dx : float
            the map's grid resolution. Centerline.dx_meter = dx * map_dx
        """

        self.line = None  # Shapely LineString
        self.head = None  # Shapely Point
        self.tail = None  # Shapely Point
        self.dis_on_line = None
        self.nx = None
        if line is not None:
            self.set_line(line)  # Init all previous properties
        else:
            self.nx = len(surface_h)
            self.dis_on_line = np.arange(self.nx) * dx

        self.order = None  # Hydrological flow level (~ Strahler number)

        # These are computed at run time by compute_centerlines
        self.flows_to = None  # pointer to a Centerline object (when available)
        self.flows_to_point = None  # point of the junction in flows_to
        self.inflows = []  # list of Centerline instances (when available)
        self.inflow_points = []  # junction points

        # Optional attrs
        self.dx = dx  # dx in pixels (assumes the line is on constant dx
        self.map_dx = map_dx  # the pixel spacing
        try:
            self.dx_meter = self.dx * self.map_dx
        except TypeError:
            # For backwards compatibility we allow this for now
            self.dx_meter = None
        self._surface_h = surface_h
        self._widths = None
        self.is_rectangular = None
        self.is_trapezoid = None
        self.orig_head = orig_head  # Useful for debugging and for filtering
        self.geometrical_widths = None  # these are kept for plotting and such
        self.apparent_mb = None  # Apparent MB, NOT weighted by width.
        self.mu_star = None  # the mu* associated with the apparent mb
        self.mu_star_is_valid = False  # if mu* leeds to good flux, keep it
        self.flux = None  # Flux (kg m-2)
        self.flux_needs_correction = False  # whether this branch was baaad
        self.rgi_id = rgi_id  # Useful if line is used with another glacier

    def set_flows_to(self, other, check_tail=True, to_head=False):
        """Find the closest point in "other" and sets all the corresponding
        attributes. Btw, it modifies the state of "other" too.
        Parameters
        ----------
        other : :py:class:`oggm.Centerline`
            another flowline where self should flow to
        """

        self.flows_to = other

        if check_tail:
            # Project the point and Check that its not too close
            prdis = other.line.project(self.tail, normalized=False)
            ind_closest = np.argmin(np.abs(other.dis_on_line - prdis)).item()
            n = len(other.dis_on_line)
            if n >= 9:
                ind_closest = utils.clip_scalar(ind_closest, 4, n-5)
            elif n >= 7:
                ind_closest = utils.clip_scalar(ind_closest, 3, n-4)
            elif n >= 5:
                ind_closest = utils.clip_scalar(ind_closest, 2, n-3)
            p = shpg.Point(other.line.coords[int(ind_closest)])
            self.flows_to_point = p
        elif to_head:
            self.flows_to_point = other.head
        else:
            # just the closest
            self.flows_to_point = _projection_point(other, self.tail)
        other.inflow_points.append(self.flows_to_point)
        other.inflows.append(self)
        
def line_interpol(line, dx):
    """Interpolates a shapely LineString to a regularly spaced one.
    Shapely's interpolate function does not guaranty equally
    spaced points in space. This is what this function is for.
    We construct new points on the line but at constant distance from the
    preceding one.
    Parameters
    ----------
    line: a shapely.geometry.LineString instance
    dx: the spacing
    Returns
    -------
    a list of equally distanced points
    """

    # First point is easy
    points = [line.interpolate(dx / 2.)]

    # Continue as long as line is not finished
    while True:
        pref = points[-1]
        pbs = pref.buffer(dx).boundary.intersection(line)
        if pbs.type == 'Point':
            pbs = [pbs]
        elif pbs.type == 'LineString':
            # This is rare
            pbs = [shpg.Point(c) for c in pbs.coords]
            assert len(pbs) == 2
        elif pbs.type == 'GeometryCollection':
            # This is rare
            opbs = []
            for p in pbs.geoms:
                if p.type == 'Point':
                    opbs.append(p)
                elif p.type == 'LineString':
                    opbs.extend([shpg.Point(c) for c in p.coords])
            pbs = opbs
        else:
            if pbs.type != 'MultiPoint':
                raise RuntimeError('line_interpol: we expect a MultiPoint '
                                   'but got a {}.'.format(pbs.type))

        try:
            # Shapely v2 compat
            pbs = pbs.geoms
        except AttributeError:
            pass

        # Out of the point(s) that we get, take the one farthest from the top
        refdis = line.project(pref)
        tdis = np.array([line.project(pb) for pb in pbs])
        p = np.where(tdis > refdis)[0]
        if len(p) == 0:
            break
        points.append(pbs[int(p[0])])

    return points


