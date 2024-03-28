""" FireClustering
This module include all functions used for doing fire clustering
"""
from .utils import timed
import rtree
import numpy as np
import math
import itertools

from sklearn.neighbors import BallTree


def build_rtree(geoms, fids=False):
    """Builds Rtree from a shapely multipolygon shape
    and optionally uses list of fids as identifier"""

    idx = rtree.index.Index()  # create new index
    for ind, geom in enumerate(geoms):
        if fids:
            idx.insert(ind, geom.bounds, fids[ind])
        else:
            idx.insert(ind, geom.bounds, ind)

    return idx


def idx_intersection(idx, bbox):
    """
    Finds all objects in an index where bounding boxes intersect with a geometry's bounding box
    """
    intersections = list(idx.intersection(bbox, objects=True))
    if len(intersections) > 0:
        fids, bbox = zip(*[(item.object, item.bbox) for item in intersections])
    else:
        fids = []
    return fids


def compute_all_spatial_distances(data, max_thresh_km):
    """ Derive neighbors for each point (with x,y)

    Parameters
    ----------
    data : pd DataFrame
        point location with 'x' and 'y' columns
    max_thresh_km : float
        maximum distance threshold (km) used for classifying neighbors

    Returns
    -------
    inds : np array of np array
        indices of neighbors (self excluded)
    """
    X = data[["x", "y"]].values

    bt = BallTree(X, leaf_size=20)

    inds = bt.query_radius(X, r=max_thresh_km * 1000, return_distance=False)
    
    # remove self from neighbors
    new_inds = []
    for i in range(len(inds)):
        pos = np.where(inds[i] == i)
        new_inds.append(np.delete(inds[i], pos))

    return np.array(new_inds, dtype=object)


@timed
def do_clustering(data, max_thresh_km):
    """ Do initial clustering for fire pixels

    Parameters
    ----------
    data : dataframe
        dataframe containing x and y values of all fire pixels
    max_thresh_km : float
        maximum distance threshold (km) used for classifying neighbors

    Returns
    -------
    point_to_cluster_id : list
        cluster id for each fire point
    """

    # copy the dataframe so the we don't modify it inplace
    data = data.copy()

    # value to fill in pixels without clustering
    NO_CLUSTER_VAL = -1

    # if number of points is 1 or 2, each point is one cluster
    num_points = len(data)
    if num_points < 3:
        data["initial_cid"] = list(range(num_points))
        return data

    # initialization
    cluster_id_counter = 0
    point_to_cluster_id = np.full(num_points, fill_value=NO_CLUSTER_VAL, dtype=np.int64)

    # compute neighbor pixels for each pixel
    neighbor_inds = compute_all_spatial_distances(data, max_thresh_km)

    # include all possible pixels in cluster
    to_check = np.full(num_points, fill_value=1, dtype=np.int8)
    while np.sum(to_check) > 0:
        start_ind = np.argmax(to_check == 1)  # catch first index to check

        neighbors_to_search = list(neighbor_inds[start_ind])
        all_neighbors = neighbors_to_search

        if (
            len(all_neighbors) == 0
        ):  # if no neighbor, record the current pixel as a separate cluster
            point_to_cluster_id[start_ind] = cluster_id_counter
            cluster_id_counter += 1
            to_check[start_ind] = 0

        else:  # if with neighbor, record all neighbors
            # find all neighbours of neighbours:
            searched_neighbours = [start_ind]
            while len(neighbors_to_search) > 0:
                # take the first of these
                px = neighbors_to_search[0]
                searched_neighbours.append(px)
                px_neighbors = list(neighbor_inds[px])
                all_neighbors = list(set(all_neighbors + px_neighbors))
                neighbors_to_search = list(
                    set(all_neighbors).difference(searched_neighbours)
                )
            # now we have all pixels in this cluster in all_neighbors
            point_to_cluster_id[all_neighbors] = cluster_id_counter
            cluster_id_counter += 1
            to_check[all_neighbors] = 0

    data["initial_cid"] = point_to_cluster_id

    return data


def cal_distance(loc1, loc2):
    """ Calculate the distance between two points

    Parameters
    ----------
    loc1 : list[lat,lon]
        position of first point
    loc2 : list[lat,lon]
        position of second point

    Returns
    -------
    distance : float
        distance (km) between two points
    """
    from .FireConsts import EARTH_RADIUS_KM

    lat1 = math.radians(loc1[0])
    lon1 = math.radians(loc1[1])
    lat2 = math.radians(loc2[0])
    lon2 = math.radians(loc2[1])

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    a = (
        math.sin(dlat / 2) ** 2
        + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2) ** 2
    )
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

    distance = EARTH_RADIUS_KM * c
    return distance


def cal_mindist(c1, c2):
    """ Calculate the minimum distance beween two clusters (may modify the algorithm to speed up this calculation)

    Parameters
    ----------
    c1 : list of [lat,lon]
        first cluster
    c2 : list of [lat,lon]
        second cluster

    Returns
    -------
    mindist : float
        the minimum distance (km) between c1 and c2
    """

    mindist = min([cal_distance(l1, l2) for l1, l2 in itertools.product(c1, c2)])

    return mindist
