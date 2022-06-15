"""Some useful functions that did not fit into the other modules.
"""
import shapely
import shapely.geometry as shpg
import numpy as np
import copy
from functools import (wraps)

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
        elif inter_poly.type == 'LineString':  # i have in treoduced "tuple()"
            inter_poly = shpg.Polygon(tuple(np.asarray(inter_poly.xy).T))
        elif inter_poly.type == 'Polygon':
            pass
        else:
            extext = ('Geometry collection not expected: '
                      '{}'.format(inter_poly.type))
            # raise InvalidGeometryError(extext)

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
    idx = (np. abs(array - value)).argmin()
    return array[idx]


class glacier_dir(object):
    def __init__(self, grid):
        self.grid = grid


class grid_inf(object):
    def __init__(self, grid):
        self.grid = grid


def lazy_property(fn):
    """Decorator that makes a property lazy-evaluated."""

    attr_name = '_lazy_' + fn.__name__

    @property
    @wraps(fn)
    def _lazy_property(self):
        if not hasattr(self, attr_name):
            setattr(self, attr_name, fn(self))
        return getattr(self, attr_name)

    return _lazy_property


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


# A faster numpy.clip when only one value is clipped (here: min).
clip_min = np.core.umath.maximum

# A faster numpy.clip when only one value is clipped (here: max).
clip_max = np.core.umath.minimum

def cls_to_geoline(cls):
    """
    list of Centerline object to list of shapely.geometry.lines

    Parameters
    ----------
    cls : list of centerlines
        list of centerlines instances.

    Returns
    -------
    list of shapely.geometry.lines

    """
    liness = []
    for cl in cls:
        liness.append(cl.line)
    lines = liness
    return lines
