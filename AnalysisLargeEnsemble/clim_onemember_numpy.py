'''
for member in $(cat memberEM.txt);do
echo ${member}
python clim_onemember_numpy.py ${member}
done
'''
import numpy as np
import pandas as pd
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap,addcyclic
from mpl_toolkits.basemap import shiftgrid

import sys
##sys.path.insert(1, '/home/msantolaria/Documents/MyPythonLibrary/ClimAnag/')
#import climbasis as climb
#import domain as dom
#import dynad as dyn
#from dynad import *
import os

sourceData='/media/msantolaria/MARIAEXT2/CMIP6/IPSL/'
outputDir='/home/msantolaria/Documents/MyPythonLibrary/AnalysisLargeEnsemble/'
plotsDir=outputDir+'Plots/'
resultsDir=outputDir+'Results/'

##Parameters

iyr=1979
fyr=2014
nyr=fyr-iyr+1
model='IPSL-CM6A-LR'
member=sys.argv[1]
#member='r1i1p1f1'
season='MAM'
imon=3
fmon=5
domain='global'
variable='tas'

fileName='tas_Amon_IPSL-CM6A-LR_historical_'+member+'_gr_185001-201412.nc'

#dsY = xr.open_dataset(sourceData+fileNameY)[variableY]
#unitY=dsY.units
#fieldY=dom.field_dom(dsY,domain)
#valsY,anomsY=climb.monthly_selection(fieldY,rmon,iyrAmon,fyrAmon)
#vals,anoms=climb.seasonal_selection(field,season,imon,iyr,fmon,fyr)


fh = Dataset(sourceData+fileName, mode = 'r')
field01= fh.variables['tas'][:]  #when loading, field=field1; mfield=field
units=fh.variables['tas'].units
xlon01 = fh.variables['lon'][:]
#nlon=len(xlon)
ylat0 = fh.variables['lat'][:]
#nlat=len(ylat)
time=fh.variables['time']
fh.close()

##https://matplotlib.org/basemap/api/basemap_api.html
field0,xlon0 = shiftgrid(180., field01, xlon01, start=False)

##Domain--------------------------------------------------
xi=70;xf=89
xlon=xlon0[xi:xf]
nlon=len(xlon)
yi=92;yf=110
ylat=ylat0[yi:yf]
nlat=len(ylat)
field=field0[:,yi:yf,xi:xf]

#Selecting period----------------------------------

field=field[(12*(iyr-1850)):(12*(fyr+1-1850)),:,:]
#nlat=field_dom.shape[1]
#nlon=field_dom.shape[2]
#Monthly/season selection---------------------------

value=np.zeros((nyr,nlat,nlon))
#value=ma.masked_array(value1,mask=field[0:nyr].mask) # [0:nyr] to make data and mask compatible

#for i in range(nyr):
#    tmp=field[12*i+(rmon-1),:,:]
#    value[i,:,:]=tmp

for i in range(nyr):
    tmp=field[12*i+(imon-1):12*i+(fmon),:,:]
    value[i,:,:]=np.mean(tmp,axis=0)

clim=np.zeros((nlat,nlon))
clim=np.ma.mean(value,axis=0)   #Climatology of january over 1979-2005
std=np.ma.std(value,axis=0)

#Decomposing field---------------------------------
#
anom=np.zeros((nyr,nlat,nlon))
#anom=ma.masked_array(anom1,mask=field[0:nyr].mask)
#
for i in range(nyr):
    anom[i,:,:]=value[i,:,:]-clim

##Need to put detrending method---------------------


#Plotting clim+std,trend,anom
##-------------------------------
#CLIM+std
fig=plt.figure(figsize=(10,12))

#m = Basemap(projection='merc',llcrnrlat=-80,urcrnrlat=80,\
#            llcrnrlon=-180,urcrnrlon=180,lat_ts=20,resolution='c')

m = Basemap(projection='merc',llcrnrlat=28,urcrnrlat=48,\
            llcrnrlon=-5,urcrnrlon=40,lat_ts=20,resolution='c')


m.drawcoastlines()
#m.drawparallels(np.arange(-90,90,30),labels=[1,0,0,0])
#m.drawmeridians(np.arange(m.lonmin,m.lonmax+30,60),labels=[0,0,0,1])

m.drawparallels(np.arange(-90,90,5),labels=[1,0,0,0])
m.drawmeridians(np.arange(m.lonmin,m.lonmax+30,5),labels=[0,0,0,1])

lons, lats = np.meshgrid(xlon,ylat)
#clevs=np.arange(10,310,20)
clevs=np.arange(-50,55,5)
cmap=plt.cm.RdYlBu_r
CS1 = m.contourf(lons,lats,clim[:,:]-273.15,clevs,cmap=cmap,latlon=True,extend='both')
cb = fig.colorbar(CS1,orientation="horizontal", pad=0.05)
cb.set_label('$^{\circ}C$', fontsize=10)
#cb.set_label('%s'%(field_units), fontsize=12)
#plotname='Clim_%s_%s_%s_%s_%i-%i_paperI_let' %(variable,obs,domain,season,iyr_obs,fyr_obs)
plotname='clim_%s_%s_%s_%s_%s_%i-%i' %(variable,model,member,domain,season,iyr,fyr)
plt.title('%s'%(plotname), fontsize=12)
plt.tight_layout()
plt.savefig(plotsDir + plotname+ '.png',format='png')
plt.show()

##Save result in netCDF

ncout = Dataset(resultsDir+plotname+'.nc','w','NETCDF4'); # using netCDF3 for output format
ncout.createDimension('lon',nlon);
ncout.createDimension('lat',nlat);
lonvar = ncout.createVariable('lon','float32',('lon'));lonvar[:] = xlon;
latvar = ncout.createVariable('lat','float32',('lat'));latvar[:] = ylat;
myvar = ncout.createVariable('clim','float32',('lat','lon'));myvar.setncattr('units',units);myvar[:] = clim;
myvar1=ncout.createVariable('std','float32',('lat','lon'));myvar1.setncattr('units',units);myvar1[:] = std;
ncout.close();
