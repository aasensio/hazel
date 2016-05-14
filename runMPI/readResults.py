import numpy as np
import matplotlib.pyplot as pl
from scipy.io import netcdf

def readResults(rootFile):

# Read inverted profiles
	ff = netcdf.netcdf_file(rootFile+'.inversion', 'r')
	synthProf = ff.variables['map'][:]
	ff.close()
	
# Read inverted parameters
	ff = netcdf.netcdf_file(rootFile+'.parameters', 'r')
	pars = ff.variables['map'][:]
	ff.close()
	
# Read errors
	ff = netcdf.netcdf_file(rootFile+'.errors', 'r')
	errors = ff.variables['map'][:]
	ff.close()
	
	return synthProf, pars, errors

synthProf, pars, errors = readResults('test')