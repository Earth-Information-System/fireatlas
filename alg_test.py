# -*- coding: utf-8 -*-
"""
Compare vector layers from fire tracaking algorithm 
with reference datasets

Created on Mon Mar  7 12:58:48 2022

@author: Rebecca Scholten
"""

import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd
import sys
sys.path.insert(1, 'D:/fire_atlas/1_code/fireatlas')
import FireClustering

ddir = r'D:\fire_atlas\Data'
pdir = r'D:\fire_atlas\2_pipeline'

#%% Alaska 2015

# 1) load data
ref = gpd.read_file(ddir + '/AlaskaFireHistory/AlaskaFireHistory_Polygons.gdb')
ref = ref[ref['FIREYEAR'] == '2015']

fta = gpd.read_file(pdir + '/2015/Summary/final2015.gpkg').to_crs(ref.crs)
mcd = gpd.read_file(pdir + '/2015/Summary/final_mcd642015.gpkg').to_crs(ref.crs)
# lakes = gpd.read_file(pdir + '/2020/Summary/final_lakes2020.gpkg')

# 2) filter for larger fires (>25ha)
min(ref.Shape_Area)
ref_large = ref[ref['Shape_Area'] > 250000]
min(fta.area) > 250000

# 3) compare number of features and total area
len(ref_large)
len(fta)

# 3) identify overlapping fire scars and compute sizes


#%% NWT 2014

# 1) load data
ref = gpd.read_file(ddir + '/CanNFDB/NFDB_poly/NFDB_poly_20210707.shp')
ref = ref[(ref['YEAR'] == 2014) & 
          ((ref['SRC_AGENCY'] == 'NT') | (ref['SRC_AGENCY'] == 'PC-WB'))]
# filter out fires that are not in NWT (<60N)
ref_bds = ref.to_crs(4326).bounds
exclude = ref_bds[ref_bds.miny < 60].index
ref = ref[~ref.index.isin(exclude)]

fta = gpd.read_file(pdir + '/2014/Summary/final2014.gpkg').to_crs(ref.crs)
mcd = gpd.read_file(pdir + '/2014/Summary/final_mcd642014.gpkg').to_crs(ref.crs)
lakes = gpd.read_file(pdir + '/2014/Summary/final_lakes2014.gpkg')
lakes = lakes[lakes.geometry != None].to_crs(ref.crs)

# check number or fires
len(lakes)
len(ref)

# check total area
sum(ref.area)
sum(lakes.area)

# filter by area needed?
min(ref.area)
min(lakes.area)

# compare size distribution
plt.hist(ref.area, bins = 28) # these are not helpful
plt.hist(lakes.area, bins = 28)

# compare individual fires
lakes_geoms = [lakes.geometry[lakes.index == f][f] for f in lakes.index]
idx = FireClustering.build_rtree(lakes_geoms) # build index of fta fires
areas = []
for fire in ref.index:
    area_ref = ref.area[ref.index == fire][fire]
    geom_ref = ref.geometry[ref.index == fire][fire]
    id_cfs = FireClustering.idx_intersection(idx, geom_ref.bounds)
    if len(id_cfs) == 0:
        continue
    area_fta = 0
    if len(id_cfs) > 1:
        for id_cf in id_cfs:  # find the actual intersection
            if lakes_geoms[id_cf].intersects(geom_ref):
                area_fta += lakes_geoms[id_cf].area
    else:
        area_fta += lakes_geoms[id_cfs[0]].area
    areas.append((fire,area_ref,area_fta))
fid,x,y = zip(*areas)
df_area = pd.DataFrame({'refid':fid,'ref':x, 'fta':y})
df_area = df_area.assign(dif=(df_area.fta-df_area.ref)/1000000) # difference in km2
df_area = df_area.assign(dif_perc = df_area.dif/df_area.ref*1000000)
plt.scatter(df_area.ref/1000000, df_area.fta/1000000)
sum(df_area.dif) # total
df_area[df_area.dif == max(df_area.dif)]
plt.hist(df_area.dif, bins = 50)


