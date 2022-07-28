# -*- coding: utf-8 -*-
"""
Compare vector layers from fire tracaking algorithm 
with reference datasets

Created on Mon Mar  7 12:58:48 2022

@author: Rebecca Scholten
"""
import statistics
import glob
import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd
import sys
sys.path.insert(1, 'D:/fire_atlas/1_code/fireatlas')
import FireClustering

ddir = 'D:/fire_atlas/Data'
pdir = 'D:/fire_atlas/2_pipeline/'

#%% Alaska
ak = gpd.read_file(ddir + '/shapefiles/Alaska.gpkg')
ref = gpd.read_file(ddir + '/AlaskaFireHistory/AlaskaFireHistory_Polygons.gdb')
# ref = ref[ref['Shape_Area'] > 200000]
out = []

for year in range(2012,2022):
    # 1) load data
    ref_yr = ref[ref['FIREYEAR'] == str(year)]
    
    # fta = gpd.read_file(pdir+'AK/'+str(year)+'/Summary/final_viirs'+str(year)+'.gpkg').to_crs(ref.crs)
    fta = gpd.read_file(pdir+str(year)+'/Summary/final_viirs'+str(year)+'.gpkg')
    fta = fta[~fta.geometry.is_empty & fta.geometry.notna()& fta.geometry.is_valid]
    fta = gpd.clip(fta, ak.geometry[0]).to_crs(ref.crs)
    
    # 2) filter for larger fires (>20ha)
    if not min(fta.area) > 200000:
        print(min(fta.area))
    
    # add stats to list
    year_stats = [year,len(ref_yr),len(fta),#len(mcd),len(lakes),
                  sum(ref_yr.area)/1000000,sum(fta.area)/1000000]#,sum(mcd.area)/1000000,sum(lakes.area)/1000000]
    out.append(year_stats)

df = pd.DataFrame(out, columns = ['year','ref_no','viirs_no','ref_area','viirs_area'])
df = df.set_index('year')
df_no = df.filter(like='no',axis=1)
df_size = df.filter(like='area',axis=1)

# plot bar plots for comparison
df_no.plot.bar()
df_size.plot.bar()


# same for ignition points (only 2015)
ign_ref = gpd.read_file(ddir + '/AlaskaFireHistory/AlaskaFireHistory_Points.gdb')
ign_ref = ign_ref[ign_ref['FIRESEASON'] == '2015']
ign = gpd.read_file(pdir + '/2015/Summary/ign2015.gpkg').to_crs(ref.crs)

#%% NWT 2014
nt = gpd.read_file(ddir + '/shapefiles/nwt.gpkg')
ref = gpd.read_file(ddir + '/CanNFDB/NFDB_poly/NFDB_poly_20210707.shp')
ref = ref[((ref['SRC_AGENCY'] == 'NT') | (ref['SRC_AGENCY'] == 'PC-WB'))]
# filter out fires that are not in NWT (<60N)
ref_bds = ref.to_crs(4326).bounds
exclude = ref_bds[ref_bds.miny < 60].index
ref = ref[~ref.index.isin(exclude)]
out = []

for year in range(2012,2022):
    # 1) load data
    ref_yr = ref[ref['YEAR'] == year]
    
    fta = gpd.read_file(pdir+str(year)+'/Summary/final_viirs'+str(year)+'.gpkg')
    fta = fta[~fta.geometry.is_empty & fta.geometry.notna()& fta.geometry.is_valid]
    fta = gpd.clip(fta, nt.geometry[0]).to_crs(ref.crs)
    # mcd = gpd.read_file(pdir+'NT/'+str(year)+'/Summary/final_mcd64'+str(year)+'.gpkg').to_crs(ref.crs)
    # lakes = gpd.read_file(pdir+'NT/'+str(year)+'/Summary/final_lakes'+str(year)+'.gpkg')
    # lakes = lakes[lakes.geometry != None].to_crs(ref.crs)
    
    # 2) filter for larger fires (>20ha)
    if not min(fta.area) > 200000:
        print(min(fta.area))
    
    # add stats to list
    year_stats = [year,len(ref_yr),len(fta),#len(mcd),len(lakes),
                  sum(ref_yr.area)/1000000,sum(fta.area)/1000000]#,sum(mcd.area)/1000000,sum(lakes.area)/1000000]
    out.append(year_stats)

df = pd.DataFrame(out, columns = ['year','ref_no','viirs_no','ref_area','viirs_area'])
df = df.set_index('year')
df_no = df.filter(like='no',axis=1)
df_size = df.filter(like='area',axis=1)

# plot bar plots for comparison
df_no.plot.bar()
df_size.plot.bar()

#%% Canada
ca = gpd.read_file(ddir + '/shapefiles/canada.gpkg')
# ref = gpd.read_file(ddir + '/CanNFDB/NFDB_poly/NFDB_poly_20210707.shp')
filenames_ref = glob.glob(ddir + '/NBAC/*/*.shp')
out = []

for year in range(2012,2022):
    print(year)
    # print(filenames_ref[year-2012])
    # 1) load data
    ref_yr = gpd.read_file(filenames_ref[year-2012])
    fta = gpd.read_file(pdir+str(year)+'/Summary/final_viirs'+str(year)+'.gpkg')
    fta = fta[~fta.geometry.is_empty & fta.geometry.notna()& fta.geometry.is_valid]
    fta = gpd.clip(fta, ca.geometry[0]).to_crs(ref.crs)
    
    # 2) filter for larger fires (>20ha)
    if not min(fta.area) > 200000:
        print(min(fta.area))
    
    # add stats to list
    year_stats = [year,len(ref_yr),len(fta),#len(mcd),len(lakes),
                  sum(ref_yr.area)/1000000,sum(fta.area)/1000000]#,sum(mcd.area)/1000000,sum(lakes.area)/1000000]
    out.append(year_stats)

df = pd.DataFrame(out, columns = ['year','ref_no','viirs_no','ref_area','viirs_area'])
df = df.set_index('year')
df_no = df.filter(like='no',axis=1)
df_size = df.filter(like='area',axis=1)

# plot bar plots for comparison
df_no.plot.bar()
df_size.plot.bar()

#%% confusion matrix overlap

# canada 2021
ref = gpd.read_file(filenames_ref[-1]).to_crs('EPSG:3571')
fta = gpd.read_file(pdir+str(2021)+'/Summary/final_viirs'+str(2021)+'.gpkg')
fta = gpd.clip(fta, ca.geometry[0]).to_crs('EPSG:3571')

# NWT 2014
ref = gpd.read_file(filenames_ref[2]).to_crs('EPSG:3571')
ref = gpd.clip(ref, nt.to_crs('EPSG:3571').geometry[0])
fta = gpd.read_file(pdir+str(2014)+'/Summary/final_viirs'+str(2014)+'.gpkg')
fta = gpd.clip(fta, nt.geometry[0]).to_crs('EPSG:3571')

# AK 2015
ref = gpd.read_file(ddir + '/AlaskaFireHistory/AlaskaFireHistory_Polygons.gdb')
ref = ref[ref['FIREYEAR'] == str(2015)].to_crs('EPSG:3571')
fta = gpd.read_file(pdir+str(2015)+'/Summary/final_viirs'+str(2015)+'.gpkg')
fta = gpd.clip(fta, ak.geometry[0]).to_crs('EPSG:3571')


# compute overlap and difference
common = gpd.overlay(ref, fta, how = 'intersection')
fta_only = gpd.overlay(fta, ref, how = 'difference')
ref_only = gpd.overlay(ref, fta, how = 'difference')

sum(common.area)/1000000
sum(fta_only.area)/1000000
sum(fta_only.area)/sum(common.area)
sum(ref_only.area)/1000000
sum(ref_only.area)/sum(common.area)

#%% ### compare individual fires ###
# done per year

# canada 2021
ref = gpd.read_file(filenames_ref[-1])
test = gpd.read_file(pdir+str(2021)+'/Summary/final_viirs'+str(year)+'.gpkg')
test = gpd.clip(test, ca.geometry[0]).to_crs(ref.crs)

test_geoms = [test.geometry[test.index == f][f] for f in test.index]
idx = FireClustering.build_rtree(test_geoms) # build index of fta fires
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
            if test_geoms[id_cf].intersects(geom_ref):
                area_fta += test_geoms[id_cf].area
    else:
        area_fta += test_geoms[id_cfs[0]].area
    areas.append((fire,area_ref,area_fta))
fid,x,y = zip(*areas)
df_area = pd.DataFrame({'refid':fid,'ref':x, 'fta':y})
df_area = df_area.groupby('fta')[['ref']].sum()
df_area = df_area.reset_index()
df_area = df_area.assign(dif=(df_area.fta-df_area.ref)/1000000) # difference in km2
df_area = df_area.assign(dif_perc = df_area.dif/df_area.ref*1000000)

plt.scatter(df_area.ref/1000000, df_area.fta/1000000) # scatterplot
plt.xlim([0, 5500])
plt.ylim([0, 5500])
sum(df_area.dif) # total difference (km2)
statistics.mean(df_area.dif) # avg difference
statistics.median(df_area.dif) # median difference
statistics.median(df_area.dif_perc) # median perc difference
df_area[df_area.dif == max(df_area.dif)] # max difference
plt.hist(df_area.dif, bins = 50)
plt.hist(df_area.dif_perc, bins = 50)

#%% OLD stats

# check number or fires
len(test)
len(ref)

# check total area
sum(ref.area)/1000000
sum(test.area)/1000000

# filter by area needed?
min(ref.area)
min(test.area)

# compare size distribution
plt.hist(ref.area, bins = 28) # these are not helpful
plt.hist(test.area, bins = 28)

# compare number of ignitions
len(ign)
len(ign_ref)



