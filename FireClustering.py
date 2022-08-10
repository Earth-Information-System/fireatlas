""" FireClustering
This module include all functions used for doing fire clustering
"""

def remove_self(inds, dist):
    ''' Remove self from the index and distance arrays

    Parameters
    ----------
    inds : np array of np array
        indices of neighbors for each point (self included)
    dist : np array of np array
        distance of neighbors for each point (self included)

    Returns
    -------
    new_inds : np array of np array
        indices of neighbors (self excluded)
    new_dist : np array of np array
        distance of neighbors (self excluded)
    '''
    import numpy as np

    new_inds = []
    new_dist = []
    for i in range(len(inds)):
       pos = np.where(inds[i] == i)
       new_inds.append(np.delete(inds[i], pos))
       new_dist.append(np.delete(dist[i], pos))

    return np.array(new_inds,dtype=object), np.array(new_dist,dtype=object)

def build_rtree(geoms, fids=False):
    '''Builds Rtree from a shapely multipolygon shape
    and optionally uses list of fids as identifier'''
    import rtree

    idx = rtree.index.Index() # create new index
    for ind, geom in enumerate(geoms):
        if fids:
            idx.insert(ind, geom.bounds, fids[ind])
        else:
            idx.insert(ind, geom.bounds, ind)

    return idx

def idx_intersection(idx, bbox):
    '''
    Finds all objects in an index whcih bounding boxes intersect with a geometry's bounding box
    '''
    intersections = list(idx.intersection(bbox, objects = True))
    if len(intersections)>0:
        fids, bbox = zip(*[(item.object, item.bbox) for item in intersections])
    else:
        fids = []
    return fids


def compute_all_spatial_distances(data, max_thresh_km):
    ''' Derive neighbors (and distances) for each point (with lat/lon)

    Parameters
    ----------
    data : pd DataFrame
        point location with 'latitude' and 'longitude' columns
    max_thresh_km : float
        maximum distance threshold (km) used for classifying neighbors

    Returns
    -------
    inds : np array of np array
        indices of neighbors (self exluded)
    dist : np array of np array
        distance (km) of neighbors (self excluded)
    '''
    import numpy as np
    import math
    from sklearn.neighbors import BallTree
    from FireConsts import EARTH_RADIUS_KM

    lats = data.latitude.values
    lons = data.longitude.values

    lats_rad = [math.radians(l) for l in lats]
    lons_rad = [math.radians(l) for l in lons]

    X = np.stack([lats_rad, lons_rad], axis=-1)

    bt = BallTree(X, metric='haversine', leaf_size=20)

    inds, dist = bt.query_radius(X, r = max_thresh_km / EARTH_RADIUS_KM, return_distance=True)
    inds, dist = remove_self(inds, dist)
    return inds, dist * EARTH_RADIUS_KM

def sort_neighbors(inds, dists):
    ''' Do neighbor sorting (based on distance) for all points

    Parameters
    ----------
    inds : np array of np array
        indices of neighbors (self exluded)
    dist : np array of np array
        distance (km) of neighbors (self excluded)

    Returns
    -------
    sorted_inds : np array of np array
        indices of neighbors (self exluded)
    sorted_dists : np array of np array
        distance (km) of neighbors (self excluded)
    '''
    import numpy as np

    sorted_inds = []
    sorted_dists = []

    for i in range(len(inds)):
        new_order = np.argsort(dists[i])

        sorted_inds.append(inds[i][new_order])
        sorted_dists.append(dists[i][new_order])

    return sorted_inds, sorted_dists

def do_clustering(data, max_thresh_km):
    ''' Do initial clustering for fire pixels

    Parameters
    ----------
    data : list (nx2)
        latitude and longitude values of all fire pixels
    max_thresh_km : float
        maximum distance threshold (km) used for classifying neighbors

    Returns
    -------
    point_to_cluster_id : list
        cluster id for each fire point
    '''
    import numpy as np
    import pandas as pd

    # value to fill in pixels without clustering
    NO_CLUSTER_VAL = -1

    # if number of points is 1 or 2, each point is one cluster
    num_points = len(data)
    if num_points < 3:
        cluster_id = list(range(num_points))
        return cluster_id

    # convert list to pd DataFrame
    dfdata = pd.DataFrame(data,columns=['latitude','longitude'])

    # initialization
    cluster_id_counter = 0
    point_to_cluster_id = np.full(num_points, fill_value=NO_CLUSTER_VAL, dtype=np.int64)

    # compute and sort neighbor pixels for each pixel
    neighbor_inds, neighbor_spatial_dists = compute_all_spatial_distances(dfdata, max_thresh_km)
    # neighbor_inds, neighbor_spatial_dists = sort_neighbors(neighbor_inds, neighbor_spatial_dists)

    # include all possible pixels in cluster
    to_check = np.full(num_points, fill_value=1, dtype=np.int8)
    while np.sum(to_check) > 0:
        start_ind = np.argmax(to_check == 1) # catch first index to check

        neighbors_to_search = list(neighbor_inds[start_ind])
        all_neighbors = neighbors_to_search

        if len(all_neighbors) == 0:  # if no neighbor, record the current pixel as a separate cluster
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
                neighbors_to_search = list(set(all_neighbors).difference(searched_neighbours))
            # now we have all pixels in this cluster in all_neighbors
            point_to_cluster_id[all_neighbors] = cluster_id_counter
            cluster_id_counter += 1
            to_check[all_neighbors] = 0

    # # try 1: retroactively change clusters (didnt work )
    # for i in range(num_points):
    #     neighbors = neighbor_inds[i]
    #     if len(neighbors) == 0:  # if no neighbor, record the current pixel as a separate cluster
    #         point_to_cluster_id[i] = cluster_id_counter
    #         cluster_id_counter += 1
    #     else:  # if with neighbor, record all neighbors
    #         # neighbor_cluster_ids = point_to_cluster_id[neighbors]   #  This won't work sometimes, e.g. 2017-11-20
    #         neighbor_cluster_ids = point_to_cluster_id[list(neighbors)]
    #         neighbor_cluster_ids = neighbor_cluster_ids[neighbor_cluster_ids != -1]

    #         if len(neighbor_cluster_ids) == 0:   # if no neighbor has been assigned to clusters, record current pixel as a separate cluster
    #             point_to_cluster_id[i] = cluster_id_counter
    #             cluster_id_counter += 1
    #         else:                               # if some neighbors have been assigned to clusters, assign current pixel to the nearest cluster
    #             point_to_cluster_id[i] = neighbor_cluster_ids[0]
    #             if len(neighbor_cluster_ids) > 1: # retroactively merge clusters
    #                 for clust in range(1, len(neighbor_cluster_ids)):
    #                     point_to_cluster_id[point_to_cluster_id == clust] = neighbor_cluster_ids[0]


    # # loop over all pixels and do clustering ORIGINAL
    # for i in range(num_points):
    #     neighbors = neighbor_inds[i]
    #     if len(neighbors) == 0:  # if no neighbor, record the current pixel as a separate cluster
    #         point_to_cluster_id[i] = cluster_id_counter
    #         cluster_id_counter += 1
    #     else:  # if with neighbor, record all neighbors
    #         # neighbor_cluster_ids = point_to_cluster_id[neighbors]   #  This won't work sometimes, e.g. 2017-11-20
    #         neighbor_cluster_ids = point_to_cluster_id[list(neighbors)]
    #         neighbor_cluster_ids = neighbor_cluster_ids[neighbor_cluster_ids != -1]

    #         if len(neighbor_cluster_ids) == 0:   # if no neighbor has been assigned to clusters, record current pixel as a separate cluster
    #             point_to_cluster_id[i] = cluster_id_counter
    #             cluster_id_counter += 1
    #         else:                               # if some neighbors have been assigned to clusters, assign current pixel to the nearest cluster
    #             point_to_cluster_id[i] = neighbor_cluster_ids[0]

    return point_to_cluster_id.tolist()

def cal_distance(loc1,loc2):
    ''' Calculate the distance between two points

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
    '''
    from FireConsts import EARTH_RADIUS_KM
    import math

    lat1 = math.radians(loc1[0])
    lon1 = math.radians(loc1[1])
    lat2 = math.radians(loc2[0])
    lon2 = math.radians(loc2[1])

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    a = math.sin(dlat / 2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

    distance = EARTH_RADIUS_KM * c
    return distance

def cal_mindist(c1,c2):
    ''' Calculate the minimum distance beween two clusters (may modify the algorithm to speed up this calculation)

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
    '''

    import itertools

    mindist = min([cal_distance(l1, l2) for l1,l2 in itertools.product(c1, c2)])

    return mindist

def filter_centroid(loc_tgt,locs,MAX_THRESH_CENT_KM=50):
    ''' Get cluster ids if centroid is close to the target cluster centroid

    Parameters
    ----------
    loc_tgt : list, [lat,lon]
        target centroid
    locs : list of [lat,lon]
        all existing fire centroid
    MAX_THRESH_CENT_KM : float
        maximum threshold to determine close clusters

    Returns
    -------
    closeindices : list
        ids of existing clusters that are close to the target cluster
    '''
    closeornot = [cal_distance(loc_tgt,loc) < MAX_THRESH_CENT_KM for loc in locs]
    closeindices = [i for i, x in enumerate(closeornot) if x == True]
    return closeindices

def cal_centroid(data):
    ''' Calculate the centroid of a list of points
    Parameters
    ----------
    data : list of [lat,lon]
        point location

    Returns
    -------
    xct : float
        lat of centroid
    yct : float
        lon of centroid
    '''
    x, y = zip(*(data))
    l = len(x)
    xct, yct =  sum(x) / l, sum(y) / l

    return xct,yct
