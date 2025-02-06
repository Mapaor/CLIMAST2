
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
units='K'
fileNameList=[]
fileNameList=glob.glob(resultsDir+'timeseries_tas_IPSL-CM6A-LR_r*_global_annual_1979-2014.txt')
fileNameList=sorted(fileNameList,key=str.lower)
print(fileNameList)

tsList=[]
#pvalueList=[]
for elem in fileNameList:
    tmp= np.loadtxt(elem)  #when loading, field=field1; mfield=field
    tsList.append(tmp)
print(tsList)
EM=np.loadtxt(resultsDir+'timeseries_tas_IPSL-CM6A-LR_EM_global_annual_1979-2014.txt')
memberlist=[]
for elem in fileNameList:
    tmp=elem.split('_')[3]
    memberlist.append(tmp)

#Decomposing field---------------------------------
#
plotname='timeseries_%s_%s_multimember_%s_%s_%i-%i' %(variable,model,domain,season,iyr,fyr)

xd=np.arange(nyr)
fig=plt.figure()
ax = fig.add_subplot(1, 1, 1)
xdplot=xd+1
my_ticks=np.arange(iyr,fyr+1,1)
frequency=5
#plt.ylim(0,100)
plt.ylabel('%s %s'%(variable,units))
plt.xlabel('Years')
plt.xticks(xdplot[::frequency],my_ticks[::frequency])
plt.title(plotname,fontsize=12)
for elem in tsList:
    plt.plot(xdplot,EM,linestyle='-',color='red')
    plt.plot(xdplot,elem,linestyle=':',color='red')#,label=memberList[i])
#plt.plot(xdplot,trend*xdplot+intercept,label='%1.2f +- %1.2f'%(10*trend,10*stderr))
plt.savefig(plotsDir + plotname+ '.png',format='png')
plt.show()


