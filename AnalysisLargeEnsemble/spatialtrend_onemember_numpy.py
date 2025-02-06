'''
for member in $(cat memberEM.txt);do
echo ${member}
python spatialtrend_onemember_numpy.py ${member}
done
'''
import numpy as np
import pandas as pd
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap,addcyclic
from scipy.stats import t
from scipy import stats
import sys
import os

sourceData='/media/msantolaria/MARIAEXT2/CMIP6/IPSL/'
outputDir='/home/msantolaria/Documents/MyPythonLibrary/AnalysisLargeEnsemble/'
plotsDir=outputDir+'Plots/'
resultsDir=outputDir+'Results/'

##Parameters

iyr=1850
fyr=2014
nyr=fyr-iyr+1
model='IPSL-CM6A-LR'
member='EM'
#member='r1i1p1f1'
#member=sys.argv[1]
season='annual'
imon=1
fmon=12
domain='global'
variable='tas'
fileName='tas_Amon_IPSL-CM6A-LR_historical_'+member+'_gr_185001-201412.nc'


fh = Dataset(sourceData+fileName, mode = 'r')
field= fh.variables['tas'][:]  #when loading, field=field1; mfield=field
units=fh.variables['tas'].units
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
#Computing trends----------------------------------
trend=np.empty((nlat,nlon))
intercept=np.empty((nlat,nlon))
pvalue=np.empty((nlat,nlon))
rvalue=np.empty((nlat,nlon))
stderr=np.empty((nlat,nlon))
xd=np.array(range((nyr)))

for j in range(nlat):
    for i in range(nlon):
        yd=anom[:,j,i]
        par = stats.linregress(xd, yd)
        stderr[j,i]=par[4]
        pvalue[j,i]=par[3]
        rvalue[j,i]=par[2]
        intercept[j,i]=par[1]
        trend[j,i]=par[0]

#2)R as test of statistical significance
#
#T-Student test: t depends on confidence and degrees of freedom
#confidence=0.975
#df=nyr-2
#t=stats.t.ppf(confidence,df)
#rs=math.sqrt((t*t)/(t*t + df))

#pvalue_ma=ma.masked_where(pvalue>=0.10,pvalue)

#Plotting trend
##-------------------------------
#CLIM+std
fig=plt.figure(figsize=(10,12))

m = Basemap(projection='merc',llcrnrlat=-80,urcrnrlat=80,\
            llcrnrlon=-180,urcrnrlon=180,lat_ts=20,resolution='c')

m.drawcoastlines()
m.drawparallels(np.arange(-90,90,30),labels=[1,0,0,0])
m.drawmeridians(np.arange(m.lonmin,m.lonmax+30,60),labels=[0,0,0,1])
lons, lats = np.meshgrid(xlon,ylat)

clevs=np.arange(-1.4,1.5,0.1)
CS1=m.contourf(lons,lats,trend*10,clevs,cmap=plt.cm.RdBu_r,latlon=True)
#levels=[0,0.1,1.0]
#cs = m.contourf(lons,lats,pvalue,levels=levels,hatches=["+", ""], alpha=0.,latlon=True)
cb = fig.colorbar(CS1,orientation="horizontal",pad=0.05)
cb.set_label('%s/dec'%(units), fontsize=10)
plotname='spatialtrend_%s_%s_%s_%s_%s_%i-%i' %(variable,model,member,season,domain,iyr,fyr)
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
myvar0 = ncout.createVariable('trend','float32',('lat','lon'));myvar0.setncattr('units',units+'/yr');myvar0[:] = trend;
myvar1 = ncout.createVariable('intercept','float32',('lat','lon'));myvar1.setncattr('units',units);myvar1[:] = intercept;
myvar2 = ncout.createVariable('rvalue','float32',('lat','lon'));myvar2.setncattr('units',units);myvar2[:] = rvalue;
myvar3 = ncout.createVariable('pvalue','float32',('lat','lon'));myvar3.setncattr('units',units);myvar3[:] = pvalue;
myvar4 = ncout.createVariable('stderr','float32',('lat','lon'));myvar4.setncattr('units',units+'/yr');myvar4[:] = stderr;
ncout.close();

