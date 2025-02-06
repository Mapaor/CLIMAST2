'''
for member in $(cat memberEM.txt);do
python spatialaverage_onemeber_numpy.py ${member}
done
'''
import xarray as xr
import numpy as np
import xarray as xr
import pandas as pd
from netCDF4 import Dataset
from scipy.stats import t
from scipy import stats
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap,addcyclic

import sys
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
domain='global'
season='annual'
imon=1
fmon=12
variable='tas'

fileName='tas_Amon_IPSL-CM6A-LR_historical_'+member+'_gr_185001-201412.nc'

#dsY = xr.open_dataset(sourceData+fileNameY)[variableY]
#unitY=dsY.units
#fieldY=dom.field_dom(dsY,domain)
#valsY,anomsY=climb.monthly_selection(fieldY,rmon,iyrAmon,fyrAmon)
#vals,anoms=climb.seasonal_selection(field,season,imon,iyr,fmon,fyr)


fh = Dataset(sourceData+fileName, mode = 'r')
field= fh.variables['tas'][:]  #when loading, field=field1; mfield=field
field_units=fh.variables['tas'].units
xlon = fh.variables['lon'][:]
nlon=len(xlon)
ylat = fh.variables['lat'][:]
nlat=len(ylat)
time=fh.variables['time']
fh.close()


##Domain--------------------------------------------------
#select
#field_dom
#Selecting period----------------------------------

field=field[(12*(iyr-1850)):(12*(fyr+1-1850)),:,:]
#nlat=field_dom.shape[1]
#nlon=field_dom.shape[2]

#Weighting latitude y_dom
#->as apply to latitude 1D no need to create a new axis (not a weighting matrix)
##
#y_dom=y[y_idom:y_fdom]
wgts = np.cos(np.deg2rad(ylat))
#wgts = np.sqrt(coslat)

tmp=np.ma.average(field,axis=1,weights=wgts) #average in latitude tmp.shape=(time,longitude)
spave_field_dom=np.ma.average(tmp,axis=1)  #average in longitude spave_field.shape=(time)
bigmatrix=spave_field_dom[0:field.shape[0]].reshape((nyr,12)) #reordering to have a matrix(row,col)=matrix(years,months)

ts_spatial=np.ma.mean(bigmatrix[:,imon:fmon+1],axis=1)
#unitY=dsY.units
#fieldY=dom.field_dom(dsY,domain)
plotname='timeseries_%s_%s_%s_%s_%s_%i-%i' %(variable,model,member,domain,season,iyr,fyr)
np.savetxt(resultsDir + plotname + '.txt',ts_spatial)
###
xd=np.array(range(nyr))
par = stats.linregress(xd,ts_spatial)
trend=par[0]
intercept=par[1]
rvalue=par[2]
pvalue=par[3]
stderr=par[4]

#calcular trend in time series y plotearla tambien

fig=plt.figure()
ax = fig.add_subplot(1, 1, 1)
xdplot=xd+1
my_ticks=np.arange(iyr,fyr+1,1)
frequency=5
#plt.ylim(0,100)
plt.ylabel('%s %s'%(variable,field_units))
plt.xlabel('Years')
plt.xticks(xdplot[::frequency],my_ticks[::frequency])
plt.title(plotname,fontsize=12)
plt.plot(xdplot,ts_spatial,linestyle='-',marker='o')
plt.plot(xdplot,trend*xdplot+intercept,label='%1.2f +- %1.2f'%(10*trend,10*stderr))
plt.legend()
plt.savefig(plotsDir + plotname+ '.png',format='png')
#plt.show()




