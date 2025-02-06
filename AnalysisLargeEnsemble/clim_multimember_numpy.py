
import numpy as np
import pandas as pd
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap,addcyclic
import numpy.ma as ma
import sys
import os
import glob
sourceData='/media/msantolaria/MARIAEXT2/CMIP6/IPSL/'
resultsDir='/home/msantolaria/Documents/MyPythonLibrary/AnalysisLargeEnsemble/'
plotsDir=resultsDir+'Plots/'


##Parameters

iyr=1850
fyr=2014
nyr=fyr-iyr+1
model='IPSL-CM6A-LR'
season='annual'
imon=1
fmon=12
domain='global'
variable='tas'

##Reference
ref=sourceData+'tas_Amon_IPSL-CM6A-LR_historical_r1i1p1f1_gr_185001-201412.nc'
fhref= Dataset(ref, mode = 'r')
fieldref= fhref.variables['tas'][:]
xlon = fhref.variables['lon'][:]
nlon=len(xlon)
ylat = fhref.variables['lat'][:]
nlat=len(ylat)
time=fhref.variables['time']
nt=len(time)
fhref.close()

members=32
field=np.zeros((nt,nlat,nlon,members))
for k in np.arange(1,members+1,1):
    fileName='tas_Amon_IPSL-CM6A-LR_historical_r'+str(k)+'i1p1f1_gr_185001-201412.nc'
    print(fileName)
    fh=Dataset(sourceData+fileName)
    var=fh.variables['tas'][:]
    fh.close()
    field[:,:,:,k]=var[:,:,:]

'''#fileName='tas_Amon_IPSL-CM6A-LR_historical_'+member+'_gr_185001-201412.nc'
fileNameList=[]
fileNameList=glob.glob(sourceData+'tas_Amon_IPSL-CM6A-LR_historical_r*_gr_185001-201412.nc')
fileNameList=sorted(fileNameList,key=str.lower)

memberlist=[]
for elem in fileNameList:
    tmp=elem.split('_')[4]
    memberlist.append(tmp)
'''
fileArrayList=[]
for elem in fileNameList:
    fh = Dataset(elem, mode = 'r')
    field= fh.variables['tas'][:]
    field=field[(12*(iyr-1850)):(12*(fyr+1-1850)),:,:]
    fh.close()
    fileArrayList.append(field)
##Reference
ref=sourceData+'tas_Amon_IPSL-CM6A-LR_historical_r1i1p1f1_gr_185001-201412.nc'
fhref= Dataset(ref, mode = 'r')
fieldref= fhref.variables['tas'][:]
xlon = fhref.variables['lon'][:]
nlon=len(xlon)
ylat = fhref.variables['lat'][:]
nlat=len(ylat)
time=fhref.variables['time']
fhref.close()


##Domain--------------------------------------------------
#select
#field_dom
#Selecting period----------------------------------

#field=field[(12*(iyr-1850)):(12*(fyr+1-1850)),:,:]
#nlat=field_dom.shape[1]
#nlon=field_dom.shape[2]
#Monthly/season selection---------------------------
valueList=[]
for elem in fileArrayList:
    value=np.zeros((nyr,nlat,nlon))
    for i in range(nyr):
        tmp=elem[12*i+(imon-1):12*i+(fmon),:,:]
        value[i,:,:]=np.mean(tmp,axis=0)
    valueList.append(value)

climList=[]
stdList=[]
for val in valueList:
    clim=np.zeros((nlat,nlon))
    clim=np.ma.mean(val,axis=0)
    std=np.ma.std(val,axis=0)
    climList.append(clim)
    stdList.append(std)

#Decomposing field---------------------------------
#
anomList=[]
for e,elem in enumerate(valueList):
    anom=np.zeros((nyr,nlat,nlon))
    for i in range(nyr):
        anom[i,:,:]=elem[i,:,:]-climList[e]
    anomList.append(anom)

##Need to put detrending method---------------------


#Plotting clim+std,trend,anom
##-------------------------------
nrows=8
ncols=4
fig,axes=plt.subplots(nrows,ncols)
for i,ax in enumerate(axes.flatten()):
    m_ax = Basemap(ax=ax,projection='merc',llcrnrlat=-80,urcrnrlat=80,\
            llcrnrlon=-180,urcrnrlon=180,lat_ts=20,resolution='c')
    m_ax.drawcoastlines()
    #m_ax.drawparallels(np.arange(-90,90,30),labels=[1,0,0,0])
    #m_ax.drawmeridians(np.arange(m_ax.lonmin,m_ax.lonmax+30,60),labels=[0,0,0,1])
    lons, lats = np.meshgrid(xlon,ylat)
#clevs=np.arange(10,310,20)
    clevs=np.arange(-50,55,5)
    cmap=plt.cm.RdYlBu_r
    print(i)
    #clim=climList[i]
    CS1 = m_ax.contourf(lons,lats,clim[:,:,i]-273.15,clevs,cmap=cmap,latlon=True,extend='both')
    #ax.title("%s"%(memberlist[i]))
    ax.annotate("%i"%(i+1),xy=(0.07,0.08),bbox={'facecolor':'w','edgecolor':'white','alpha':0.08},fontsize=6)
plt.subplots_adjust(top=0.92,bottom=0.08,left=0.15,right=0.85,hspace=0.001,wspace=0.001)
cax=fig.add_axes([0.25,0.02,0.5,0.03])
cb = fig.colorbar(CS1,cax=cax,orientation="horizontal")
cb.set_label('$^{\circ}C$', fontsize=10)
#cb.set_label('%s'%(field_units), fontsize=12)
#plotname='Clim_%s_%s_%s_%s_%i-%i_paperI_let' %(variable,obs,domain,season,iyr_obs,fyr_obs)
plotname='Clim_%s_%s_multimember_%s_%s_%i-%i' %(variable,model,domain,season,iyr,fyr)
plt.suptitle('%s'%(plotname), fontsize=12)
plt.savefig(plotsDir + plotname+ '.pdf',format='pdf')
plt.show()


