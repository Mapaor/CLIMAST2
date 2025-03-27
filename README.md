#CLIMAST Program ('Interdisciplinary training in climate sciences for master students') 
#http://www.climast.fr/ 
#Contact: m.santolaria-otin@meteo.ub.edu

Python Practical Exercice 

Introduction of CMIP Architecture
Large Ensemble model
Basic statistical analysis

This repository contains:
1) The README.md file
2) /scripts/: some scrips .py and one notebook .pynb as the basic example using numpy for the practical exercise in CLIMAST.
3) /AnalysisLargeEnsemble/: some scripts to plot, compute basic statistics with Large Ensemble models that you can use or that they can inspire you to develop your own scripts.s

-----------------------------------------------------------

GET READY (PYTHON - NUMPY)

Use the terminal or an Anaconda Prompt for the following steps: To create an environment: conda create --name practicasenv

Let’s install some packages, so first let’s activate our environment, “to work inside it”:

conda activate practicasenv

and install some packages:

conda install -c anaconda basemap 
conda install -c anaconda scipy

Info: https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html

https://anaconda.org/anaconda/basemap https://anaconda.org/conda-forge/netcdf4 https://anaconda.org/anaconda/scipy

IF YOU WANT TO GO FURTHER...investigate xarray library
https://anaconda.org/anaconda/xarray 
https://docs.xarray.dev/en/stable/gallery.html


USEFUL INFO

About CMIP6: https://pcmdi.llnl.gov/CMIP6/Guide/dataUsers.html
