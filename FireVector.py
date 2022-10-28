""" FireVector
This is the module used for vector related calculations
"""


def alpha_shape(points, alpha):
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
    from shapely.ops import cascaded_union, polygonize
    from scipy.spatial import Delaunay
    import numpy as np
    import math
    import shapely.geometry as geometry

    if len(points) < 4:
        # When you have a triangle, there is no sense in computing an alpha shape.
        return geometry.MultiPoint(list(points)).convex_hull

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

    coords = np.array([point for point in points])
    tri = Delaunay(coords)
    edges = set()
    edge_points = []

    # loop over triangles:
    # ia, ib, ic = indices of corner points of the
    # triangle
    # for ia, ib, ic in tri.vertices:
    for ia, ib, ic in tri.vertices:
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

    # return triangles, edge_points
    return cascaded_union(triangles), edge_points


def addbuffer(geom_ll, vbuf):
    """ add a geometry in geographic projection with a buffer in meter
    Parameters
    ----------
    geom_ll : shapely geometry
        the geometry (in lat /lon)
    vbuf : float
        the buffer value (in meter)
    """
    # to speed up the calculation, use the following approach
    # the buffer in geographic projection is calculated using the centroid lat value
    # from FireConsts import EARTH_RADIUS_KM
    # import numpy as np
    # lat = geom_ll.centroid.y
    # ldeg = (EARTH_RADIUS_KM*np.cos(np.deg2rad(lat))*1000*2*np.pi/360)
    # vbufdeg = vbuf/ldeg    # convert vbuf in m to degs
    # geom_ll_buf = geom_ll.buffer(vbufdeg)

    return geom_ll.buffer(vbuf)


def doMultP(locs, buf):
    """ deirvve a MultipPolygon (bufferred MultiPoint) shape from given fire locations
    Parameters
    ----------
    locs : list (nx2)
        latitude and longitude values of all fire pixels
    Returns
    -------
    multP : MultiPolygon object
        calculated shape
    """
    from shapely.geometry import MultiPoint

    # MultiPoint shape
    geom = MultiPoint([(x, y) for x, y in locs])

    # Add buffer to form MultiPolygon shape
    geom = addbuffer(geom, buf)

    return geom


def doConvH(locs, buf):
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
    try:
        qhull = ConvexHull(locs)
    except:
        return None

    # derive qhull object vertices
    verts = [(locs[i][0], locs[i][1]) for i in qhull.vertices]

    # convert vertices to polygon
    hull = Polygon(verts)

    # Add buffer
    # hull = hull.buffer(VIIRSbuf)
    hull = addbuffer(hull, buf)
    return hull


def doConcH(locs, buf, alpha=100):
    """ derive the concave hull given fire locations
    Parameters
    ----------
    locs : list (nx2)
        latitude and longitude values of all fire pixels
    alpha : int
        alpha value (1/deg)
    Returns
    -------
    hull : Polygon or MultiPolygon object
        calculated hull shape
    """

    # calcluate alpha shape
    try:
        # switch lats and lons to use the alpha_shape
        # locs_lonlat = [(v[1],v[0]) for v in locs]
        # concave_hull, edge_points = alpha_shape(locs_lonlat,alpha=alpha)

        concave_hull, edge_points = alpha_shape(locs, alpha=alpha)
        if (
            concave_hull.area == 0
        ):  # sometimes the concave_hull algorithm returns a empty polygon
            concave_hull = doConvH(locs, buf)
            return concave_hull
        else:
            # Add buffer
            # hull = concave_hull.buffer(VIIRSbuf)
            hull = addbuffer(concave_hull, buf)
            return hull

    except:
        return None


def cal_hull(fp_locs, sensor="viirs"):
    """ wrapper to calculate the hull given fire locations.
        the returned hull type depends on the pixel number
    Parameters
    ----------
    fp_locs : list (nx2)
        latitude and longitude values of all fire pixels
    Returns
    -------
    hull : object
        calculated hull (a buffer of VIIRS half pixel size included)
    """
    from FireConsts import valpha, VIIRSbuf  # ,MCD64buf
    import numpy as np

    # set buffer according to sensor
    if sensor == "viirs":
        buf = VIIRSbuf
    elif sensor == "mcd64":
        buf = MCD64buf
    else:
        print('please set sensor to viirs or mcd64')

    # number of points
    nfp = len(fp_locs)

    # For cluster with 1-2 pixel, calculate hull using buffered points (MultiPolygon)
    if nfp < 3:
        hull = doMultP(fp_locs, buf)
    # For cluster with 3 pixels, calculate hull using convex hull
    elif nfp == 3:  # call doConvH to get the hull
        hull = doConvH(fp_locs, buf)
        if hull == None:  # in case where convex hull can't be determined, call doMultP
            hull = doMultP(fp_locs, buf)
        elif (
            hull.area == 0
        ):  # sometimes the fire pixels are in a traight line, in this case, also call doMultP
            hull = doMultP(fp_locs, buf)
    # For cluster with more than 3 pixels, calculate hull using alpha shape
    else:  # call doConcH to get the hull
        # derive the alpha value used in doConcH (in 1/deg)
        # x,y = zip(*fp_locs)
        # vdeg = sum(x)/len(x)
        # km1deg = 6371*np.cos(np.deg2rad(vdeg))*2*np.pi/360
        # valphadeg = 1/(valpha/1000/km1deg)   # in 1/deg
        hull = doConcH(fp_locs, buf, alpha=valpha)
        if hull == None:  # in case where convex hull can't be determined, call doMultP
            hull = doMultP(fp_locs, buf)
        elif hull.area == 0:
            hull = doMultP(fp_locs, buf)
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
    hts_buf = addbuffer(hull, -extbuffer)

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