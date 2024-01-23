""" FireVector
This is the module used for vector related calculations
"""


def doConcH(points, alpha):
    """
    Compute the alpha shape (concave hull) of a set
    of points.
    @param points: Iterable container of points.
    @param alpha: alpha value to influence the
        gooeyness of the border. Smaller numbers
        don't fall inward as much as larger numbers.
        Too large, and you lose everything!
    source: http://blog.thehumangeo.com/2014/05/12/drawing-boundaries-in-python/
    """
    from shapely.ops import unary_union, polygonize
    from scipy.spatial import Delaunay
    import numpy as np
    import math
    import shapely.geometry as geometry

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
    """ derive the convex hull given fire locations
    Parameters
    ----------
    locs : list (nx2)
        latitude and longitude values of all fire pixels
    Returns
    -------
    hull : Polygon object
        calculated hull shape
    """
    from scipy.spatial import ConvexHull
    from shapely.geometry import Polygon

    # calculate the convex hull using scipy.spatial.ConvexHull
    qhull = ConvexHull(locs)

    # derive qhull object vertices
    verts = locs[qhull.vertices]

    # convert vertices to polygon
    hull = Polygon(verts)

    return hull


def cal_hull(fp_locs, sensor="viirs"):
    """ wrapper to calculate the hull given fire locations.
        the returned hull type depends on the pixel number
    Parameters
    ----------
    fp_locs : np.array (nx2)
        x, y values of all fire pixels
    Returns
    -------
    hull : object
        calculated hull (a buffer of VIIRS half pixel size included)
    """
    from shapely.geometry import MultiPoint
    from FireConsts import valpha, VIIRSbuf

    # set buffer according to sensor
    if sensor == "viirs":
        buf = VIIRSbuf
    elif sensor == "mcd64":
        buf = MCD64buf
    else:
        print('please set sensor to viirs or mcd64')

    # number of pixels
    nfp = len(fp_locs)
    hull = None

    # if there are more than 3 points: try using concave hull
    if nfp > 3:
        hull = doConcH(fp_locs, alpha=valpha)
    
    # if you don't have a good hull yet and there are more 
    # than 2 points: try using convex hull
    if nfp > 2 and (hull is None or hull.area == 0):
        hull = doConvH(fp_locs)

    # if you don't have a good hull yet: make a MultiPoint
    if hull is None or hull.area == 0:
        hull = MultiPoint(fp_locs)
        
    hull = hull.buffer(buf)
    
    return hull


def cal_extpixels(fps, hull, alpha=100):
    """ calculate the exterior pixels around a hull
    Parameters
    ----------
    fps : list (nx3)
        list of all fire pixels
    hull : geometry, 'Polygon' | 'MultiPoint'
        the existing hull for the fire
    alpha : int
        alpha value (1/deg) to define the inward buffer
    Returns
    -------
    fps_ext : list (nx2)
        the exterior fire pixels for the fire
    """
    from shapely.geometry import Point
    from FireConsts import extbuffer

    # use alpha to define an inward buffer and an interior part of the hull
    # hts_buf = hull.buffer(-1/alpha)
    hts_buf = hull.buffer(-extbuffer)

    # loop over all pixel locations in fps to determine exterior pixels
    fps_ext = []
    for fp in fps:
        # point locations of each pixel (fp.oc[0]:lat, fp.loc[1]:lon)
        pt = Point(fp.loc[0], fp.loc[1])

        # exclude points in interior part of the hull
        if not hts_buf.contains(pt):
            fps_ext.append(fp)
    return fps_ext


def calConcHarea(hull):
    """ calculate area given the concave hull (km2)
    Parameters
    ----------
    hull : geometry, 'Polygon' | 'MultiPoint'
        the hull for the fire
    Returns
    -------
    farea : float
        the area (km2) of the polygon enclosed by the vertices
    """
    import math
    from FireConsts import EARTH_RADIUS_KM

    farea = hull.area  # area in deg**2

    if farea > 0:
        lat = hull.bounds[1]  # latitude in deg

        # convert deg**2 to km2
        farea = (
            farea
            * EARTH_RADIUS_KM ** 2
            * math.cos(math.radians(lat))
            * math.radians(1) ** 2
        )

    return farea