
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
outputDir='/home/msantolaria/Documents/MyPythonLibrary/AnalysisLargeEnsemble/'
plotsDir=outputDir+'Plots/'
resultsDir=outputDir+'Results/'

##Parameters

iyr=1979
fyr=2014
nyr=fyr-iyr+1
model='IPSL-CM6A-LR'
season='annual'
domain='global'
variable='tas'

fileNameList=[]
fileNameList=glob.glob(resultsDir+'spatialtrend_tas_IPSL-CM6A-LR_r*_annual_global_1979-2014.nc')
fileNameList=sorted(fileNameList,key=str.lower)
print(fileNameList)

trendList=[]
#pvalueList=[]
for elem in fileNameList:
    fh = Dataset(elem, mode = 'r')
    trend= fh.variables['trend'][:]  #when loading, field=field1; mfield=field
    #pvalue=fh.variables['pvalue'][:]
    trendList.append(trend)
    #stdList.append(std)

ref=Dataset(fileNameList[0],mode='r')
units=ref.variables['trend'].units
xlon = ref.variables['lon'][:]
nlon=len(xlon)
ylat = ref.variables['lat'][:]
nlat=len(ylat)


memberlist=[]
for elem in fileNameList:
    tmp=elem.split('_')[3]
    memberlist.append(tmp)

#Decomposing field---------------------------------
#
##-------------------------------
nrows=4
ncols=8
fig=plt.figure(figsize=(18,16))
fig,axes=plt.subplots(nrows,ncols)
for i,ax in enumerate(axes.flatten()):
    m_ax = Basemap(ax=ax,projection='merc',llcrnrlat=-80,urcrnrlat=80,\
            llcrnrlon=-180,urcrnrlon=180,lat_ts=20,resolution='c')
    m_ax.drawcoastlines()
    #m_ax.drawparallels(np.arange(-90,90,30),labels=[1,0,0,0])
    #m_ax.drawmeridians(np.arange(m_ax.lonmin,m_ax.lonmax+30,60),labels=[0,0,0,1])
    lons, lats = np.meshgrid(xlon,ylat)
    clevs=np.arange(-4.0,4.25,0.25)
    cmap=plt.cm.seismic
    print(i)
    CS1 = m_ax.contourf(lons,lats,trendList[i]*10,clevs,cmap=cmap,latlon=True,extend='both')
    ax.annotate("%s"%(memberlist[i]),xy=(0.07,1.1),bbox={'facecolor':'w','edgecolor':'white','alpha':0.08},fontsize=8,xycoords='axes fraction')
plt.subplots_adjust(top=0.90, bottom=0.14, left=0.05, right=0.95, hspace=0.1,wspace=0.1)
#cax=fig.add_axes([left,bottom,width,height])
cax=fig.add_axes([0.20,0.1,0.5,0.03])
cb = fig.colorbar(CS1,cax=cax,orientation="horizontal")
cb.set_label('$^{\circ}C / dec$', fontsize=8)
#cb.set_label('%s'%(field_units), fontsize=12)
#plotname='Clim_%s_%s_%s_%s_%i-%i_paperI_let' %(variable,obs,domain,season,iyr_obs,fyr_obs)
plotname='spatialtrend_%s_%s_multimember_%s_%s_%i-%i' %(variable,model,domain,season,iyr,fyr)
plt.suptitle('%s'%(plotname), fontsize=12, y=0.95)
#plt.tight_layout()
plt.savefig(plotsDir + plotname+ '.png',format='png')
plt.show()


