# module that interfaces the NIWA VCSN fetching module to source inputs
# Built for HBRC inhouse use, reuse should be notified
# Changes if any, made, should be uploaded to a branch to identify potential updates and changes

# Why python?
## the idea is that this module can be used irrespective of source and majority of tools are in python

#user Variables
thisStart = '2020-05-01'

#see imports for prerequisites

import matplotlib.pyplot as plt

import datetime

import numpy as np
import pandas as pd

import geopandas as gpd

from scipy.interpolate import interp2d

from statsmodels.tsa.seasonal import STL


import kFetchVCSN
#fetch VCSN data
kVc = kFetchVCSN.kFetchVCSN()
## the credentials line is hidden




#original file is here, remote file fetching fails some times, making a local copy
#read the catchment polygons
subCatchGDF = gpd.read_file('catchmentShapeFiles/SOURCECatchments_polys.shp')
#print(subCatchGDF)

#centroids
subCentroids = subCatchGDF.copy() #duplicate before calculations
subCentroids = subCentroids[['CatID','geometry']]
subCentroids.to_crs(epsg=4326,inplace=True)  #change it now to get wgs coordinates
subCentroids['arCenter'] = subCentroids.centroid 


#polygon to filter vcsn points
bndPolyGDF = subCatchGDF.dissolve()
#offset outside by 10km - vcsn grids are so
if bndPolyGDF.crs == 'epsg:4326':
    bndPolyGDF = bndPolyGDF.buffer(10/111.11) #km/km
    bndPolyGDF.to_crs(epsg=2193,inplace=True)
else: #everything else is assumed to be in meters
    bndPolyGDF = bndPolyGDF.buffer(10000)
bndPolyGDF = gpd.GeoDataFrame(geometry=gpd.GeoSeries(bndPolyGDF))

#vcsn table
tempDf = kVc._kFetchVCSN__gridTable
vcsnPtsDf = gpd.GeoDataFrame(tempDf,
                             geometry=gpd.points_from_xy(tempDf.LONGT, tempDf.LAT),
                             crs=4326)
#vcsnPtsDf = vcsnPtsDf.set_crs(4326, allow_override=True)
vcsnPtsDf.to_crs(epsg=2193,inplace=True)

#within_points = vcsnPtsDf.within(bndPolyGDF)
selectPoints = gpd.sjoin(vcsnPtsDf, bndPolyGDF, predicate = 'within')
selectPoints.to_crs(epsg=4326,inplace=True)
#print(selectPoints.head())

"""
#checks
temp = subCatchGDF.copy()
temp.to_crs(epsg=4326,inplace=True)
ax = temp.plot(color='grey')
selectPoints.plot(ax=ax)
plt.show()
"""

#get met data
precDf = pd.DataFrame()
petDf = pd.DataFrame()
for thisRow in selectPoints.iterrows():
    temp = thisRow[1]
    print(temp.LONGT,temp.LAT,temp.AGENT_NO)

    kVc.selectSite = (temp.LAT,temp.LONGT) #(-39,176.5)
    
    refPetDf = kVc.fetchData('Penman Evaporation', startTime=thisStart)
    refPetDf.rename(columns={"validityTime":"timestamp","Penman Evaporation":temp.AGENT_NO},inplace=True)
    refPetDf['timestamp'] = pd.to_datetime(refPetDf['timestamp']).apply(lambda dt: dt.replace(hour=0)).dt.tz_localize(None)
    refPetDf.set_index('timestamp', inplace=True)
    if petDf.empty:
        petDf = refPetDf.copy()
    else:
        petDf = pd.concat([petDf,refPetDf], axis=1, join="outer")
    
    refPrecDf = kVc.fetchData('Rain Accumulation', startTime=thisStart)
    refPrecDf.rename(columns={"validityTime":"timestamp","Rain Accumulation":temp.AGENT_NO},inplace=True)
    refPrecDf['timestamp']=pd.to_datetime(refPrecDf['timestamp']).apply(lambda dt: dt.replace(hour=0)).dt.tz_localize(None)
    refPrecDf.set_index('timestamp', inplace=True)
    if precDf.empty:
        precDf = refPrecDf.copy()
    else:
        precDf = pd.concat([precDf,refPrecDf], axis=1, join="outer")

#interpolation function
def myInterpFn(var,intrpFn,g):
    #print(g)
    thisMax = 250 if var == 'Rainfall' else 12
    temp = np.nan
    try:
        temp = intrpFn(g.x,g.y)[0]
        if not (0 <= temp <= thisMax):
            temp = np.nan
    except:
        pass
    return temp


#rainfall
for (thisItem,thisVar) in [(precDf,'Rainfall'), (petDf,'PET')]:
    thisSubCentroids = pd.DataFrame()
    for thisDate in thisItem.iterrows():
        temp = thisDate[1]
        print(thisDate[0])
        #print(temp.index)

        ptsDf = selectPoints.copy()
        ptsDf.set_index('AGENT_NO',inplace=True)
        #print(ptsDf.columns)
        ptsDf = ptsDf.merge(temp,left_index=True,right_index=True)
        #print(ptsDf.columns)
        #ptsDf.plot(column=thisDate[0], cmap='OrRd')
        #plt.show()

        intrpFn = interp2d(ptsDf['LONGT'], ptsDf['LAT'], ptsDf[thisDate[0]], kind='linear')
        if thisSubCentroids.empty:
            thisSubCentroids = subCentroids.copy()
        thisSubCentroids[thisDate[0]] = thisSubCentroids.arCenter.apply(lambda g: myInterpFn(thisVar,intrpFn,g)) #intrpFn(g.x,g.y)[0])
        #break
        thisSubCentroids = thisSubCentroids.copy()

    thisSubCentroids = thisSubCentroids[thisSubCentroids['CatID'].isin(subCentroids['CatID'].values)]
    thisSubCentroids = thisSubCentroids.T
    thisSubCentroids.columns = thisSubCentroids[thisSubCentroids.index == 'CatID'].values[0]
    thisSubCentroids.index.name = 'dTime'
    thisSubCentroids.reset_index(inplace=True)
    #print(thisSubCentroids.head())

    #remove non date time rows
    #thisSubCentroids = thisSubCentroids[isinstance(thisSubCentroids['oldIndex'],(datetime.datetime, np.datetime64, datetime.date))]
    thisSubCentroids = thisSubCentroids[pd.to_datetime(thisSubCentroids['dTime'],errors='coerce').notna()] 
    thisSubCentroids.set_index('dTime',inplace=True)

    #interpolate missing values
    for thisCol in thisSubCentroids.columns:
        thisSubCentroids[thisCol] = pd.to_numeric(thisSubCentroids[thisCol],errors = 'coerce') #objects to float series
    
    thisSubCentroids.interpolate(axis=0, inplace=True)
    
    thisSubCentroids.to_csv(thisVar+'.csv')

#ax = subCatchGDF.plot(color='yellow')
#within_points.plot(ax=ax)
#plt.show()

#bndPolyGDF.plot()
#plt.show()
