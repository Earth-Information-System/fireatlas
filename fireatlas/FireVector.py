""" FireVector
This is the module used for vector related calculations
"""

import math
import geopandas as gpd
import shapely.geometry as geometry

from shapely.ops import unary_union, polygonize
from shapely.geometry import Polygon, MultiPoint, MultiLineString

from scipy.spatial import Delaunay, ConvexHull

from fireatlas import FireConsts
from fireatlas import settings


def doConcH(points, alpha):
    """
    Compute the alpha shape (concave hull) of a set
    of points.

    Parameters
    ----------
    locs : np.array (nx2)
        x, y values of all fire pixels
    alpha :
        alpha value to influence the
        gooeyness of the border. Smaller numbers
        don't fall inward as much as larger numbers.
        Too large, and you lose everything!
    source: http://blog.thehumangeo.com/2014/05/12/drawing-boundaries-in-python/
    """

    def add_edge(edges, edge_points, coords, i, j):
        """
        Add a line between the i-th and j-th points,
        if not in the list already
        """
        if (i, j) in edges or (j, i) in edges:
            # already added
            return
        edges.add((i, j))
        edge_points.append(coords[[i, j]])

    coords = points
    tri = Delaunay(coords)
    edges = set()
    edge_points = []

    # loop over triangles:
    for ia, ib, ic in tri.simplices:
        pa = coords[ia]
        pb = coords[ib]
        pc = coords[ic]
        # Lengths of sides of triangle
        a = math.sqrt((pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2)
        b = math.sqrt((pb[0] - pc[0]) ** 2 + (pb[1] - pc[1]) ** 2)
        c = math.sqrt((pc[0] - pa[0]) ** 2 + (pc[1] - pa[1]) ** 2)
        # Semiperimeter of triangle
        s = (a + b + c) / 2.0

        # Area of triangle by Heron's formula
        try:
            area = math.sqrt(s * (s - a) * (s - b) * (s - c))
        except:
            area = 0

        if area > 0:  # modified by YC
            circum_r = a * b * c / (4.0 * area)
        else:
            circum_r = 0
        # Here's the radius filter.
        # print circum_r
        # if circum_r < 1.0/alpha:
        if circum_r < alpha:
            add_edge(edges, edge_points, coords, ia, ib)
            add_edge(edges, edge_points, coords, ib, ic)
            add_edge(edges, edge_points, coords, ic, ia)
    m = geometry.MultiLineString(edge_points)
    triangles = list(polygonize(m))

    # return triangles
    return unary_union(triangles)


def doConvH(locs):
    """derive the convex hull given fire locations
    Parameters
    ----------
    locs : np.array (nx2)
        x, y values of all fire pixels
    Returns
    -------
    hull : Polygon object
        calculated hull shape
    """
    # calculate the convex hull using scipy.spatial.ConvexHull
    qhull = ConvexHull(locs)

    # derive qhull object vertices
    verts = locs[qhull.vertices]

    # convert vertices to polygon
    hull = Polygon(verts)

    return hull


def cal_hull(locs):
    """wrapper to calculate the hull given fire locations.
        the returned hull type depends on the pixel number
    Parameters
    ----------
    locs : np.array (nx2)
        x, y values of all fire pixels
    Returns
    -------
    hull : object
        calculated hull (a buffer of VIIRS half pixel size included)
    """
    # set buffer according to sensor
    if settings.FIRE_SENSOR == "viirs":
        buf = settings.VIIRSbuf
    elif settings.FIRE_SENSOR == "mcd64":
        buf = settings.MCD64buf

    # number of pixels
    nfp = len(locs)
    hull = None

    # if there are more than 3 points: try using concave hull
    if nfp > 3:
        hull = doConcH(locs, alpha=settings.valpha)

    # if you don't have a good hull yet and there are more
    # than 2 points: try using convex hull
    if nfp > 2 and (hull is None or hull.area == 0):
        hull = doConvH(locs)

    # if you don't have a good hull yet: make a MultiPoint
    if hull is None or hull.area == 0:
        hull = MultiPoint(locs)

    hull = hull.buffer(buf)

    return hull


def get_ext_pixels(pixels, hull):
    """calculate the exterior pixels around a hull
    Parameters
    ----------
    pixels : df
        dataframe containing x and y of all pixels
    hull : geometry, 'Polygon' | 'MultiPolygon'
        the existing hull for the fire
    Returns
    -------
    boolean np.array
        whether each pixel is part of the exterior
    """

    hts_buf = hull.buffer(-settings.extbuffer)
    pixel_arr = gpd.points_from_xy(pixels["x"], pixels["y"])

    return ~pixel_arr.intersects(hts_buf)


def get_fline_pixels(pixels, hull):
    """calculate the fline pixels around a hull
    Parameters
    ----------
    pixels : df
        dataframe containing x and y of all pixels
    hull : geometry, 'Polygon' | 'MultiPolygon'
        the existing hull for the fire
    Returns
    -------
    boolean np.array
        whether each pixel is part of the fline
    """
    pixel_arr = gpd.points_from_xy(pixels["x"], pixels["y"])

    # extract the pixels near the hull
    # if hull is a polygon, find new pixels near the hull
    if hull.geom_type == "Polygon":
        lr = hull.exterior.buffer(settings.fpbuffer)
        return pixel_arr.intersects(lr)

    # if hull is a multipolygon, find new pixels near the outer shell
    elif hull.geom_type == "MultiPolygon":
        mlr = MultiLineString([x.exterior for x in hull.geoms]).buffer(settings.fpbuffer)
        return pixel_arr.intersects(mlr)


def calConcHarea(hull):
    """calculate area given the concave hull (km2)
    Parameters
    ----------
    hull : geometry, 'Polygon' | 'MultiPoint'
        the hull for the fire
    Returns
    -------
    farea : float
        the area (km2) of the polygon enclosed by the vertices
    """
    farea = hull.area  # area in deg**2

    if farea > 0:
        lat = hull.bounds[1]  # latitude in deg

        # convert deg**2 to km2
        farea = (
            farea
            * settings.EARTH_RADIUS_KM**2
            * math.cos(math.radians(lat))
            * math.radians(1) ** 2
        )

    return farea
