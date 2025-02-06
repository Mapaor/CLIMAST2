###
##Loading libraries
###
import numpy as np
import pandas as pd
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap,addcyclic,shiftgrid
from scipy.stats import t
from scipy import stats
import sys
import os

###
##Selecting output/input directories
###
sourceData='/media/msantolaria/MARIAEXT2/CMIP6/IPSL/'
outputDir='/home/msantolaria/Documents/MyPythonLibrary/AnalysisLargeEnsemble/'
plotsDir=outputDir+'Plots/'
resultsDir=outputDir+'Results/'

###
##Selecting parameters
###
#Period
iyr=1850 #initial year
fyr=2014 #final year 
nyr=fyr-iyr+1 # number of year
#Season/month
season='annual'
imon=1
fmon=12
#Region
domain='global'
#Model and variable
model='IPSL-CM6A-LR'
variable='tas'
##Number of members
nbr=32
##Adapt multiplot size depending of number of members
nrows=4
ncols=8

###
##Reference file for basic parameters
##
#Open
ref=sourceData+'tas_Amon_IPSL-CM6A-LR_historical_r1i1p1f1_gr_185001-201412.nc'
fhref= Dataset(ref, mode = 'r')
fieldref0= fhref.variables['tas'][:]
units=fhref.variables['tas'].units
xlon0 = fhref.variables['lon'][:]
##Shifting longitude from 0 to 360 to -180 to 180
##Info: https://matplotlib.org/basemap/api/basemap_api.html
fieldref,xlon = shiftgrid(180., fieldref0, xlon0, start=False)
nlon=len(xlon)
ylat = fhref.variables['lat'][:]
nlat=len(ylat)
time=fhref.variables['time']
nt=len(time)
fhref.close()

#Select period ref
fieldref=fieldref[(12*(iyr-1850)):(12*(fyr+1-1850)),:,:]
###
##Let's create a matrix containing all members selected
#We select period and season/month -> value array
###
value=np.zeros((nyr,nlat,nlon,nbr))
field=np.zeros((fieldref.shape[0],nlat,nlon))#,nbr))
for k in np.arange(0,nbr,1):
    fileName='tas_Amon_IPSL-CM6A-LR_historical_r'+str(k+1)+'i1p1f1_gr_185001-201412.nc'
    print(fileName)
    fh=Dataset(sourceData+fileName)
    field0=fh.variables['tas'][:]
    lon0=fh.variables['lon'][:]
    field,lon = shiftgrid(180., field0, lon0, start=False)
    field=field[(12*(iyr-1850)):(12*(fyr+1-1850)),:,:]
    for i in range(nyr):
        tmp=field[12*i+(imon-1):12*i+(fmon),:,:]
        value[i,:,:,k]=np.mean(tmp,axis=0)
    fh.close()
###
###
##Computing climatology and standard deviation
###
clim=np.zeros((nlat,nlon,nbr))
std=np.zeros((nlat,nlon,nbr))
for k in np.arange(0,nbr,1):
    clim[:,:,k]=np.mean(value[:,:,:,k],axis=0) 
    std[:,:,k]=np.std(value[:,:,:,k],axis=0)

###
#Plotting climatology and standard deviation
###
#fig=plt.figure(figsize=(18,16))
fig,axes=plt.subplots(nrows,ncols,figsize=(18,16))
for k,ax in enumerate(axes.flatten()):
    m_ax = Basemap(ax=ax,projection='merc',llcrnrlat=-80,urcrnrlat=80,\
            llcrnrlon=-180,urcrnrlon=180,lat_ts=20,resolution='c')
    m_ax.drawcoastlines()
    lons, lats = np.meshgrid(xlon,ylat)
    clevs=np.arange(-50,55,5)
    cmap=plt.cm.RdYlBu_r
    CS1 = m_ax.contourf(lons,lats,clim[:,:,k]-273.15,clevs,cmap=cmap,latlon=True,extend='both')
    ax.annotate("%i"%(k+1),xy=(0.02,0.95),bbox={'facecolor':'w','edgecolor':'white','alpha':0.08},fontsize=8,xycoords='axes fraction')
plt.subplots_adjust(top=0.90, bottom=0.14, left=0.05, right=0.95, hspace=0.1,wspace=0.1)
#cax=fig.add_axes([left,bottom,width,height])
cax=fig.add_axes([0.20,0.1,0.5,0.03])
cb = fig.colorbar(CS1,cax=cax,orientation="horizontal")
cb.set_label('$^{\circ}C$', fontsize=8)
plotname='clim_%s_%s_multimember_%s_%s_%i-%i' %(variable,model,domain,season,iyr,fyr)
plt.suptitle('%s'%(plotname), fontsize=12, y=0.95)
#plt.tight_layout()
plt.savefig(plotsDir + plotname+ '.png',format='png')
#plt.show()

fig,axes=plt.subplots(nrows,ncols,figsize=(18,16))
for k,ax in enumerate(axes.flatten()):
    m_ax = Basemap(ax=ax,projection='merc',llcrnrlat=-80,urcrnrlat=80,\
            llcrnrlon=-180,urcrnrlon=180,lat_ts=20,resolution='c')
    m_ax.drawcoastlines()
    #m_ax.drawparallels(np.arange(-90,90,30),labels=[1,0,0,0])
    #m_ax.drawmeridians(np.arange(m_ax.lonmin,m_ax.lonmax+30,60),labels=[0,0,0,1])
    lons, lats = np.meshgrid(xlon,ylat)
    clevs=np.arange(-2,2.1,0.1)
    cmap=plt.cm.rainbow
    CS1 = m_ax.contourf(lons,lats,std[:,:,k],clevs,cmap=cmap,latlon=True,extend='both')
    ax.annotate("%i"%(k+1),xy=(0.02,0.95),bbox={'facecolor':'w','edgecolor':'white','alpha':0.08},fontsize=8,xycoords='axes fraction')
plt.subplots_adjust(top=0.90, bottom=0.14, left=0.05, right=0.95, hspace=0.1,wspace=0.1)
#cax=fig.add_axes([left,bottom,width,height])
cax=fig.add_axes([0.20,0.1,0.5,0.03])
cb = fig.colorbar(CS1,cax=cax,orientation="horizontal")
cb.set_label('$^{\circ}C$', fontsize=8)
plotname='std_%s_%s_multimember_%s_%s_%i-%i' %(variable,model,domain,season,iyr,fyr)
plt.suptitle('%s'%(plotname), fontsize=12, y=0.95)
#plt.tight_layout()
plt.savefig(plotsDir + plotname+ '.png',format='png')
#plt.show()

###
##Computing anomalies
###
anom=np.zeros((nyr,nlat,nlon,nbr))
#anom=ma.masked_array(anom1,mask=field[0:nyr].mask)

for i in range(nyr):
    anom[i,:,:,:]=value[i,:,:,:]-clim[:,:,:]
###
##Computing detrended anomalies using different detrending methods ( linear, quadratic, cubic)
#Info: #---------------------
#https://numpy.org/doc/stable/reference/generated/numpy.polyfit.html
#-----------------------------
###
##
deg=1 #Choose the degree of the polynomy

anom_notrend=np.zeros((nyr,nlat,nlon,nbr))
xyears=np.arange(0,nyr,1)
for k in np.arange(0,nbr,1):
    for j in range(nlat):
        for i in range(nlon):
            if deg==1:
                poly=np.polyfit(xyears,anom[:,j,i,k],deg=1)
                anom_notrend[:,j,i,k]=anom[:,j,i,k]-(poly[0]*xyears+poly[1])
            elif deg==2:
                poly=np.polyfit(xyears,anom[:,j,i,k],deg=2)
                anom_notrend[:,j,i,k]=anom[:,j,i,k]-(poly[0]*xyears**2+poly[1]*xyears+poly[2])
            elif deg==3:
                poly=np.polyfit(xyears,anom[:,j,i,k],deg=3)
                anom_notrend[:,j,i,k]=anom[:,j,i,k]-(poly[0]*xyears**3+poly[1]*xyears**2+poly[2]*xyears+poly[3])

###
##Computing standard deviation of anomalies ( with or without trend)
###
stdnotrend=np.zeros((nlat,nlon,nbr))

for k in np.arange(0,nbr,1):
    stdnotrend[:,:,k]=np.std(anom_notrend[:,:,:,k],axis=0)

###
##Plotting standard deviation of anomalies
###
fig,axes=plt.subplots(nrows,ncols,figsize=(18,16))
for k,ax in enumerate(axes.flatten()):
    m_ax = Basemap(ax=ax,projection='merc',llcrnrlat=-80,urcrnrlat=80,\
            llcrnrlon=-180,urcrnrlon=180,lat_ts=20,resolution='c')
    m_ax.drawcoastlines()
    #m_ax.drawparallels(np.arange(-90,90,30),labels=[1,0,0,0])
    #m_ax.drawmeridians(np.arange(m_ax.lonmin,m_ax.lonmax+30,60),labels=[0,0,0,1])
    lons, lats = np.meshgrid(xlon,ylat)
    clevs=np.arange(-2,2.1,0.1)
    cmap=plt.cm.rainbow
    CS1 = m_ax.contourf(lons,lats,stdnotrend[:,:,k],clevs,cmap=cmap,latlon=True,extend='both')
    ax.annotate("%i"%(k+1),xy=(0.02,0.95),bbox={'facecolor':'w','edgecolor':'white','alpha':0.08},fontsize=8,xycoords='axes fraction')
plt.subplots_adjust(top=0.90, bottom=0.14, left=0.05, right=0.95, hspace=0.1,wspace=0.1)
#cax=fig.add_axes([left,bottom,width,height])
cax=fig.add_axes([0.20,0.1,0.5,0.03])
cb = fig.colorbar(CS1,cax=cax,orientation="horizontal")
cb.set_label('$^{\circ}C$', fontsize=8)
#cb.set_label('%s'%(field_units), fontsize=12)
plotname='stdnotrend_%s_%s_multimember_%s_%s_%i-%i' %(variable,model,domain,season,iyr,fyr)
plt.suptitle('%s'%(plotname), fontsize=12, y=0.95)
#plt.tight_layout()
plt.savefig(plotsDir + plotname+ '.png',format='png')
#plt.show()


###
## Computing spatial trend
###
trend=np.empty((nlat,nlon,nbr))
intercept=np.empty((nlat,nlon,nbr))
pvalue=np.empty((nlat,nlon,nbr))
rvalue=np.empty((nlat,nlon,nbr))
stderr=np.empty((nlat,nlon,nbr))
xd=np.array(range((nyr)))
for k in np.arange(0,nbr,1):
    for j in range(nlat):
        for i in range(nlon):
            yd=anom[:,j,i,k]
            par = stats.linregress(xd, yd)
            stderr[j,i,k]=par[4]
            pvalue[j,i,k]=par[3]
            rvalue[j,i,k]=par[2]
            intercept[j,i,k]=par[1]
            trend[j,i,k]=par[0]
###
##Plotting spatial trend
###
fig,axes=plt.subplots(nrows,ncols,figsize=(18,16))
for k,ax in enumerate(axes.flatten()):
    m_ax = Basemap(ax=ax,projection='merc',llcrnrlat=-80,urcrnrlat=80,\
            llcrnrlon=-180,urcrnrlon=180,lat_ts=20,resolution='c')
    m_ax.drawcoastlines()
    #m_ax.drawparallels(np.arange(-90,90,30),labels=[1,0,0,0])
    #m_ax.drawmeridians(np.arange(m_ax.lonmin,m_ax.lonmax+30,60),labels=[0,0,0,1])
    lons, lats = np.meshgrid(xlon,ylat)
    clevs=np.arange(-1.0,1.1,0.1)
    cmap=plt.cm.seismic
    CS1 = m_ax.contourf(lons,lats,trend[:,:,k]*10,clevs,cmap=cmap,latlon=True,extend='both')
    ax.annotate("%i"%(k+1),xy=(0.02,0.95),bbox={'facecolor':'w','edgecolor':'white','alpha':0.08},fontsize=8,xycoords='axes fraction')
plt.subplots_adjust(top=0.90, bottom=0.14, left=0.05, right=0.95, hspace=0.1,wspace=0.1)
#cax=fig.add_axes([left,bottom,width,height])
cax=fig.add_axes([0.20,0.1,0.5,0.03])
cb = fig.colorbar(CS1,cax=cax,orientation="horizontal")
cb.set_label('$^{\circ}C / dec$', fontsize=8)
plotname='spatialtrend_%s_%s_multimember_%s_%s_%i-%i' %(variable,model,domain,season,iyr,fyr)
plt.suptitle('%s'%(plotname), fontsize=12, y=0.95)
#plt.tight_layout()
plt.savefig(plotsDir + plotname+ '.png',format='png')
#plt.show()

###
##Spatial average
###
domain='global'

ts_value=np.zeros((nyr,nbr))
ts_anom=np.zeros((nyr,nbr))
ts_anom_notrend=np.zeros((nyr,nbr))
wgts= np.cos(np.deg2rad(ylat))

for k in np.arange(0,nbr,1):
    ts_value[:,k]=np.average(np.average(value[:,:,:,k],axis=1,weights=wgts),axis=1)
    ts_anom[:,k]=np.average(np.average(anom[:,:,:,k],axis=1,weights=wgts),axis=1)
    ts_anom_notrend[:,k]=np.average(np.average(anom_notrend[:,:,:,k],axis=1,weights=wgts),axis=1)

trend_ts_value=[]
trend_ts_anom=[]
trend_ts_anom_notrend=[]

xd=np.array(range(nyr))
for k in np.arange(0,nbr,1):
    tmp= stats.linregress(xd,ts_value[:,k])
    trend_ts_value.append(tmp)
    tmp1= stats.linregress(xd,ts_anom[:,k])
    trend_ts_anom.append(tmp1)

ts_value_EM=np.sum(ts_value,axis=1)/nbr
ts_anom_EM=np.sum(ts_anom,axis=1)/nbr
trend_ts_anom_EM=stats.linregress(xd,ts_anom_EM)
ts_anom_notrend_EM=np.sum(ts_anom_notrend,axis=1)/nbr


plotname='timeseries_%s_%s_multimember_%s_%s_%i-%i' %(variable,model,domain,season,iyr,fyr)
###
##Plotting time series
###
fig, axs = plt.subplots(3)
xdplot=xd+1
my_ticks=np.arange(iyr,fyr+1,1)
frequency=10
##Values
axs[0].set_ylabel('%s %s'%(variable,units))
axs[0].set_xlabel('Years')
axs[0].set_xticks(xdplot[::frequency],my_ticks[::frequency])
axs[0].set_title('Anom %s %s %i-%i'%(domain,season,iyr,fyr))
for k in range(nbr):
    axs[0].plot(xdplot,ts_value_EM,'red',linestyle='-')
    axs[0].plot(xdplot,ts_value[:,k],'red',linestyle=':',alpha=0.5)
##Anomalies
axs[1].set_ylabel('%s %s'%(variable,units))
axs[1].set_xlabel('Years')
axs[1].set_xticks(xdplot[::frequency],my_ticks[::frequency])
axs[1].set_title('Anom %s %s %i-%i'%(domain,season,iyr,fyr))
for k in range(nbr):
    axs[1].plot(xdplot,ts_anom_EM,'red',linestyle='-')
    axs[1].plot(xdplot,trend_ts_anom_EM[0]*xdplot+trend_ts_anom_EM[1],'black')
    axs[1].plot(xdplot,ts_anom[:,k],'red',linestyle=':',alpha=0.5)
##Anomalies detrended
axs[2].set_ylabel('%s %s'%(variable,units))
axs[2].set_xlabel('Years')
axs[2].set_xticks(xdplot[::frequency],my_ticks[::frequency])
axs[2].set_title('Anom no trend deg %i %s %s %i-%i'%(deg,domain,season,iyr,fyr))
for k in range(nbr):
    axs[2].plot(xdplot,ts_anom_notrend_EM,'red',linestyle='-')
    axs[2].plot(xdplot,ts_anom_notrend[:,k],'red',linestyle=':',alpha=0.5)

plt.savefig(plotsDir + plotname+ '.png',format='png')
#plt.show()
###
##Selecting domain
###
#domain:EuroMed
#longitude: 10W-40E; latitude: 30N-50N
#in Python: x[70:89]; y[92:110]
###
domain='EuroMed'
xi=70;xf=89
xlon_dom=xlon[xi:xf]
nlon_dom=len(xlon_dom)
yi=92;yf=110
ylat_dom=ylat[yi:yf]
nlat_dom=len(ylat_dom)

value_dom=value[:,yi:yf,xi:xf,:]
anom_dom=anom[:,yi:yf,xi:xf,:]
anom_notrend_dom=anom_notrend[:,yi:yf,xi:xf,:]

ts_value_dom=np.zeros((nyr,nbr))
ts_anom_dom=np.zeros((nyr,nbr))
ts_anom_notrend_dom=np.zeros((nyr,nbr))
wgts_dom= np.cos(np.deg2rad(ylat_dom))

for k in np.arange(0,nbr,1):
    ts_value_dom[:,k]=np.average(np.average(value_dom[:,:,:,k],axis=1,weights=wgts_dom),axis=1)
    ts_anom_dom[:,k]=np.average(np.average(anom_dom[:,:,:,k],axis=1,weights=wgts_dom),axis=1)
    ts_anom_notrend_dom[:,k]=np.average(np.average(anom_notrend_dom[:,:,:,k],axis=1,weights=wgts_dom),axis=1)


trend_ts_value_dom=[]
trend_ts_anom_dom=[]
trend_ts_anom_notrend_dom=[]
xd=np.array(range(nyr))
for k in np.arange(0,nbr,1):
    tmp= stats.linregress(xd,ts_value_dom[:,k])
    trend_ts_value_dom.append(tmp)
    tmp1= stats.linregress(xd,ts_anom_dom[:,k])
    trend_ts_anom_dom.append(tmp1)

ts_value_dom_EM=np.sum(ts_value_dom,axis=1)/nbr
ts_anom_dom_EM=np.sum(ts_anom_dom,axis=1)/nbr
trend_ts_anom_dom_EM=stats.linregress(xd,ts_anom_dom_EM)
ts_anom_notrend_dom_EM=np.sum(ts_anom_notrend_dom,axis=1)/nbr

plotname='timeseries_%s_%s_multimember_%s_%s_%i-%i' %(variable,model,domain,season,iyr,fyr)


fig, axs = plt.subplots(3)
xdplot=xd+1
my_ticks=np.arange(iyr,fyr+1,1)
frequency=10
#Value
axs[0].set_ylabel('%s %s'%(variable,units))
axs[0].set_xlabel('Years')
axs[0].set_xticks(xdplot[::frequency],my_ticks[::frequency])
axs[0].set_title('Anom %s %s %i-%i'%(domain,season,iyr,fyr))
for k in range(nbr):
    axs[0].plot(xdplot,ts_value_dom_EM,'red',linestyle='-')
    axs[0].plot(xdplot,ts_value_dom[:,k],'red',linestyle=':',alpha=0.5)
##Anomalies
axs[1].set_ylabel('%s %s'%(variable,units))
axs[1].set_xlabel('Years')
axs[1].set_xticks(xdplot[::frequency],my_ticks[::frequency])
axs[1].set_title('Anom %s %s %i-%i'%(domain,season,iyr,fyr))
for k in range(nbr):
    axs[1].plot(xdplot,ts_anom_dom_EM,'red',linestyle='-')
    axs[1].plot(xdplot,trend_ts_anom_dom_EM[0]*xdplot+trend_ts_anom_dom_EM[1],'black')
    axs[1].plot(xdplot,ts_anom_dom[:,k],'red',linestyle=':',alpha=0.5)
##Anomalies detrended
axs[2].set_ylabel('%s %s'%(variable,units))
axs[2].set_xlabel('Years')
axs[2].set_xticks(xdplot[::frequency],my_ticks[::frequency])
axs[2].set_title('Anom no trend deg %i %s %s %i-%i'%(deg,domain,season,iyr,fyr))
for k in range(nbr):
    axs[2].plot(xdplot,ts_anom_notrend_dom_EM,'red',linestyle='-')
    axs[2].plot(xdplot,ts_anom_notrend_dom[:,k],'red',linestyle=':',alpha=0.5)

plt.savefig(plotsDir + plotname+ '.png',format='png')
#plt.show()

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


