import numpy as np
import pandas as pd
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap,addcyclic,shiftgrid
from scipy.stats import t
from scipy import stats
import glob
import sys
import os

sourceData='/media/msantolaria/MARIAEXT2/CMIP6/IPSL/'
outputDir='/home/msantolaria/Documents/MyPythonLibrary/AnalysisLargeEnsemble/'
plotsDir=outputDir+'Plots/'
resultsDir=outputDir+'Results/'

##Parameters

iyr=str(1850)
fyr=str(2014)
nyr=int(fyr)-int(iyr)+1
model='IPSL-CM6A-LR'
season='annual'
imon=1
fmon=12
domain='global'
variable='tas'
nbr=32
nrows=4
ncols=8
deg=1

##Reference
ref=sourceData+'tas_Amon_IPSL-CM6A-LR_historical_r1i1p1f1_gr_185001-201412.nc'
fhref= Dataset(ref, mode = 'r')
fieldref0= fhref.variables['tas'][:]
units=fhref.variables['tas'].units
xlon0 = fhref.variables['lon'][:]
##Shifting grid from 0 to 360 to -180 to 180
##https://matplotlib.org/basemap/api/basemap_api.html
fieldref,xlon = shiftgrid(180., fieldref0, xlon0, start=False)
nlon=len(xlon)
ylat = fhref.variables['lat'][:]
nlat=len(ylat)
time=fhref.variables['time']
nt=len(time)
fhref.close()


###Loading clim
fileNameList=[]
fileNameList=glob.glob(resultsDir+'clim'+'_tas_IPSL-CM6A-LR_r*_'+domain+'_'+season+'_'+iyr+'-'+fyr+'.nc')
fileNameList=sorted(fileNameList,key=str.lower)

climList=[]
stdList=[]
for elem in fileNameList:
    fh = Dataset(elem, mode = 'r')
    clim= fh.variables['clim'][:]  
    std=fh.variables['std'][:]
    climList.append(clim)
    stdList.append(std)

memberlist=[]
for elem in fileNameList:
    tmp=elem.split('_')[3]
    memberlist.append(tmp)

#Plot------------------------------------------------------------------------------------
fig=plt.figure(figsize=(18,16))
fig,axes=plt.subplots(nrows,ncols)
for i,ax in enumerate(axes.flatten()):
    m_ax = Basemap(ax=ax,projection='merc',llcrnrlat=-80,urcrnrlat=80,\
            llcrnrlon=-180,urcrnrlon=180,lat_ts=20,resolution='c')
    m_ax.drawcoastlines()
    #m_ax.drawparallels(np.arange(-90,90,30),labels=[1,0,0,0])
    #m_ax.drawmeridians(np.arange(m_ax.lonmin,m_ax.lonmax+30,60),labels=[0,0,0,1])
    lons, lats = np.meshgrid(xlon,ylat)
    clevs=np.arange(-50,55,5)
    cmap=plt.cm.RdYlBu_r
    CS1 = m_ax.contourf(lons,lats,climList[i]-273.15,clevs,cmap=cmap,latlon=True,extend='both')
    ax.annotate("%s"%(memberlist[i]),xy=(0.07,1.1),bbox={'facecolor':'w','edgecolor':'white','alpha':0.08},fontsize=8,xycoords='axes fraction')
plt.subplots_adjust(top=0.90, bottom=0.14, left=0.05, right=0.95, hspace=0.1,wspace=0.1)
#cax=fig.add_axes([left,bottom,width,height])
cax=fig.add_axes([0.20,0.1,0.5,0.03])
cb = fig.colorbar(CS1,cax=cax,orientation="horizontal")
cb.set_label('$^{\circ}C$', fontsize=8)
plotname='clim_%s_%s_multimember_%s_%s_%s-%s' %(variable,model,domain,season,iyr,fyr)
plt.suptitle('%s'%(plotname), fontsize=12, y=0.95)
#plt.tight_layout()
plt.savefig(plotsDir + plotname+ '.png',format='png')
plt.show()

fig=plt.figure(figsize=(18,16))
fig,axes=plt.subplots(nrows,ncols)
for i,ax in enumerate(axes.flatten()):
    m_ax = Basemap(ax=ax,projection='merc',llcrnrlat=-80,urcrnrlat=80,\
            llcrnrlon=-180,urcrnrlon=180,lat_ts=20,resolution='c')
    m_ax.drawcoastlines()
    lons, lats = np.meshgrid(xlon,ylat)
    clevs=np.arange(-3,3.25,0.25)
    cmap=plt.cm.rainbow
    CS1 = m_ax.contourf(lons,lats,stdList[i],clevs,cmap=cmap,latlon=True,extend='both')
    ax.annotate("%s"%(memberlist[i]),xy=(0.07,1.1),bbox={'facecolor':'w','edgecolor':'white','alpha':0.08},fontsize=8,xycoords='axes fraction')
plt.subplots_adjust(top=0.90, bottom=0.14, left=0.05, right=0.95, hspace=0.1,wspace=0.1)
#cax=fig.add_axes([left,bottom,width,height])
cax=fig.add_axes([0.20,0.1,0.5,0.03])
cb = fig.colorbar(CS1,cax=cax,orientation="horizontal")
cb.set_label('$^{\circ}C$', fontsize=8)
plotname='std_%s_%s_multimember_%s_%s_%s-%s' %(variable,model,domain,season,iyr,fyr)
#plt.tight_layout()
plt.suptitle('%s'%(plotname), fontsize=12, y=0.95)

plt.savefig(plotsDir + plotname+ '.png',format='png')
plt.show()

#Spatial trend-----------------------------------------------------------------------------------------
fileNameList=[]
fileNameList=glob.glob(resultsDir+'spatialtrend'+'_tas_IPSL-CM6A-LR_r*_'+season+'_'+domain+'_'+iyr+'-'+fyr+'.nc')
fileNameList=sorted(fileNameList,key=str.lower)
print(fileNameList)

trendList=[]
#pvalueList=[]
for elem in fileNameList:
    fh = Dataset(elem, mode = 'r')
    trend= fh.variables['trend'][:] 
    #pvalue=fh.variables['pvalue'][:]
    trendList.append(trend)
    #stdList.append(std)

fig=plt.figure(figsize=(18,16))
fig,axes=plt.subplots(nrows,ncols)
for i,ax in enumerate(axes.flatten()):
    m_ax = Basemap(ax=ax,projection='merc',llcrnrlat=-80,urcrnrlat=80,\
            llcrnrlon=-180,urcrnrlon=180,lat_ts=20,resolution='c')
    m_ax.drawcoastlines()
    #m_ax.drawparallels(np.arange(-90,90,30),labels=[1,0,0,0])
    #m_ax.drawmeridians(np.arange(m_ax.lonmin,m_ax.lonmax+30,60),labels=[0,0,0,1])
    lons, lats = np.meshgrid(xlon,ylat)
    clevs=np.arange(-1.0,1.1,0.1)
    cmap=plt.cm.seismic
    CS1 = m_ax.contourf(lons,lats,trendList[i]*10,clevs,cmap=cmap,latlon=True,extend='both')
    ax.annotate("%s"%(memberlist[i]),xy=(0.07,1.1),bbox={'facecolor':'w','edgecolor':'white','alpha':0.08},fontsize=8,xycoords='axes fraction')
plt.subplots_adjust(top=0.90, bottom=0.14, left=0.05, right=0.95, hspace=0.1,wspace=0.1)
#cax=fig.add_axes([left,bottom,width,height])
cax=fig.add_axes([0.20,0.1,0.5,0.03])
cb = fig.colorbar(CS1,cax=cax,orientation="horizontal")
cb.set_label('$^{\circ}C / dec$', fontsize=8)
plotname='spatialtrend_%s_%s_multimember_%s_%s_%s-%s' %(variable,model,domain,season,iyr,fyr)
plt.suptitle('%s'%(plotname), fontsize=12, y=0.95)
plt.savefig(plotsDir + plotname+ '.png',format='png')
#plt.show()
#------------------------------------
fileNameList=[]
fileNameList=glob.glob(resultsDir+'anomnotrend_deg'+str(deg)+'_std_'+'_tas_IPSL-CM6A-LR_r*_'+season+'_'+domain+'_'+iyr+'-'+fyr+'.nc')
fileNameList=sorted(fileNameList,key=str.lower)

stdnotrendList=[]
for elem in fileNameList:
    fh = Dataset(elem, mode = 'r')
    std=fh.variables['stdnotrend'][:]
    stdList.append(std)


fig,axes=plt.subplots(nrows,ncols)
for i,ax in enumerate(axes.flatten()):
    m_ax = Basemap(ax=ax,projection='merc',llcrnrlat=-80,urcrnrlat=80,\
            llcrnrlon=-180,urcrnrlon=180,lat_ts=20,resolution='c')
    m_ax.drawcoastlines()
    lons, lats = np.meshgrid(xlon,ylat)
    clevs=np.arange(-3,3.25,0.25)
    cmap=plt.cm.rainbow
    CS1 = m_ax.contourf(lons,lats,stdList[i],clevs,cmap=cmap,latlon=True,extend='both')
    ax.annotate("%s"%(memberlist[i]),xy=(0.07,1.1),bbox={'facecolor':'w','edgecolor':'white','alpha':0.08},fontsize=8,xycoords='axes fraction')
plt.subplots_adjust(top=0.90, bottom=0.14, left=0.05, right=0.95, hspace=0.1,wspace=0.1)
#cax=fig.add_axes([left,bottom,width,height])
cax=fig.add_axes([0.20,0.1,0.5,0.03])
cb = fig.colorbar(CS1,cax=cax,orientation="horizontal")
cb.set_label('$^{\circ}C$', fontsize=8)
plotname='std_notrend_deg%i_%s_%s_multimember_%s_%s_%s-%s' %(deg,variable,model,domain,season,iyr,fyr)
plt.suptitle('%s'%(plotname), fontsize=12, y=0.95)
#plt.tight_layout()
plt.savefig(plotsDir + plotname+ '.png',format='png')
#plt.show()


####Spatial average------------------------------------
domain='global'
anoms=[]
anoms=glob.glob(resultsDir+'timeseries_anoms_'+'tas_IPSL-CM6A-LR_r*_'+domain+'_'+season+'_'+iyr+'-'+fyr+'.txt')
anoms=sorted(anoms,key=str.lower)

anoms_notrend=[]
anoms_notrend=glob.glob(resultsDir+'timeseries_anoms_notrend_deg'+str(deg)+'_'+'tas_IPSL-CM6A-LR_r*_'+domain+'_'+season+'_'+iyr+'-'+fyr+'.txt')
anoms_notrend=sorted(anoms_notrend,key=str.lower)

plotname='timeseries_anoms_'+'tas_IPSL-CM6A-LR_multimember_'+domain+'_'+season+'_'+iyr+'-'+fyr

ts_anom=[]
ts_anom_notrend=[]
for i in range(len(anoms)):
    tmp=np.loadtxt(anoms[i])
    tmp1=np.loadtxt(anoms_notrend[i])
    ts_anom.append(tmp)
    ts_anom_notrend.append(tmp1)

ts_anom_EM=np.sum(np.asarray(ts_anom),axis=0)/nbr
ts_anom_notrend_EM=np.sum(np.asarray(ts_anom_notrend),axis=0)/nbr

trend_ts_anom=[]
trend_ts_anom_notrend=[]

xd=np.array(range(nyr))
for k in np.arange(0,nbr,1):
    tmp1= stats.linregress(xd,ts_anom[k])
    trend_ts_anom.append(tmp1)

trend_ts_anom_EM=stats.linregress(xd,ts_anom_EM)

fig, axs = plt.subplots(2)
xdplot=xd+1
my_ticks=np.arange(int(iyr),int(fyr)+1,1)
frequency=10
#plt.ylim(0,100)
axs[0].set_ylabel('%s %s'%(variable,units))
axs[0].set_xlabel('Years')
axs[0].set_xticks(xdplot[::frequency],my_ticks[::frequency])
axs[0].set_title('Anom %s %s %s-%s'%(domain,season,iyr,fyr))
for k in range(nbr):
    axs[0].plot(xdplot,ts_anom_EM,'red',linestyle='-')
    ##For deg=1 ( else use polyval function or change polynom)
    axs[0].plot(xdplot,trend_ts_anom_EM[0]*xdplot+trend_ts_anom_EM[1],'black')
    axs[0].plot(xdplot,ts_anom[k],'red',linestyle=':',alpha=0.5)
axs[1].set_ylabel('%s %s'%(variable,units))
axs[1].set_xlabel('Years')
axs[1].set_xticks(xdplot[::frequency],my_ticks[::frequency])
axs[1].set_title('Anom no trend deg %i %s %s %s-%s'%(1,domain,season,iyr,fyr))
for k in range(nbr):
    axs[1].plot(xdplot,ts_anom_notrend_EM,'red',linestyle='-')
    axs[1].plot(xdplot,ts_anom_notrend[k],'red',linestyle=':',alpha=0.5)

plt.savefig(plotsDir + plotname+ '.png',format='png')

plt.show()
###--------------------------------------------------------------

####Spatial average------------------------------------
domain='EuroMed'
anoms_dom=[]
anoms_dom=glob.glob(resultsDir+'timeseries_anoms_'+'tas_IPSL-CM6A-LR_r*_'+domain+'_'+season+'_'+iyr+'-'+fyr+'.txt')
anoms_dom=sorted(anoms_dom,key=str.lower)

anoms_notrend_dom=[]
anoms_notrend_dom=glob.glob(resultsDir+'timeseries_anoms_notrend_deg'+str(deg)+'_'+'tas_IPSL-CM6A-LR_r*_'+domain+'_'+season+'_'+iyr+'-'+fyr+'.txt')
anoms_notrend_dom=sorted(anoms_notrend_dom,key=str.lower)

plotname='timeseries_anoms_'+'tas_IPSL-CM6A-LR_multimember_'+domain+'_'+season+'_'+iyr+'-'+fyr

ts_anom_dom=[]
ts_anom_notrend_dom=[]
for i in range(len(anoms_dom)):
    tmp=np.loadtxt(anoms_dom[i])
    tmp1=np.loadtxt(anoms_notrend_dom[i])
    ts_anom_dom.append(tmp)
    ts_anom_notrend_dom.append(tmp1)

ts_anom_dom_EM=np.sum(np.asarray(ts_anom_dom),axis=0)/nbr
ts_anom_notrend_dom_EM=np.sum(np.asarray(ts_anom_notrend_dom),axis=0)/nbr

trend_ts_anom_dom=[]
trend_ts_anom_notrend_dom=[]

xd=np.array(range(nyr))
for k in np.arange(0,nbr,1):
    tmp1= stats.linregress(xd,ts_anom_dom[k])
    trend_ts_anom.append(tmp1)

trend_ts_anom_dom_EM=stats.linregress(xd,ts_anom_dom_EM)

fig, axs = plt.subplots(2)
xdplot=xd+1
my_ticks=np.arange(int(iyr),int(fyr)+1,1)
frequency=10
#plt.ylim(0,100)
axs[0].set_ylabel('%s %s'%(variable,units))
axs[0].set_xlabel('Years')
axs[0].set_xticks(xdplot[::frequency],my_ticks[::frequency])
axs[0].set_title('Anom %s %s %s-%s'%(domain,season,iyr,fyr))
for k in range(nbr):
    axs[0].plot(xdplot,ts_anom_dom_EM,'red',linestyle='-')
    ##For deg=1 ( else use polyval function or change polynom)
    axs[0].plot(xdplot,trend_ts_anom_dom_EM[0]*xdplot+trend_ts_anom_dom_EM[1],'black')
    axs[0].plot(xdplot,ts_anom[k],'red',linestyle=':',alpha=0.5)
    #axs[0].plot(xdplot,trend_ts_anom_dom[k][0]*xdplot+trend_ts_anom[k][1])
axs[1].set_ylabel('%s %s'%(variable,units))
axs[1].set_xlabel('Years')
axs[1].set_xticks(xdplot[::frequency],my_ticks[::frequency])
axs[1].set_title('Anom no trend deg %i %s %s %s-%s'%(1,domain,season,iyr,fyr))
for k in range(nbr):
    axs[1].plot(xdplot,ts_anom_notrend_EM,'red',linestyle='-')
    axs[1].plot(xdplot,ts_anom_notrend[k],'red',linestyle=':',alpha=0.5)

plt.savefig(plotsDir + plotname+ '.png',format='png')

plt.show()


###
##Scatter plot
###
domain='global'
subdomain='EuroMed'

trends=[]
trends_dom=[]
for k in range(nbr):
    tmp=trend_ts_anom[k][0]
    tmp1=trend_ts_anom_dom[k][0]
    trends.append(tmp)
    trends_dom.append(tmp1)

trends=np.asarray(trends)
trends_dom=np.asarray(trends_dom)
reg=stats.linregress(trends,trends_dom)
rvalue=reg[2]

fig=plt.figure()
ax = fig.add_subplot(1, 1, 1)
x=np.arange(0,nbr,1)
plt.plot(trends*10,trends_dom*10,'o')
plt.xlabel('%s'%(domain))
plt.ylabel('%s'%(subdomain))
plt.title("Anom trends %s %i-%i r=%1.2f"%(season,iyr,fyr,rvalue))
plotname='scatter_%s_%s_multimember_%s_%s_%i-%i' %(variable,model,subdomain,season,iyr,fyr)
plt.savefig(plotsDir+plotname+'.png',format='png')
plt.show()


