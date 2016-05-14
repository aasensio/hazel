import numpy as np
import matplotlib.pyplot as pl
import scipy.io as io

def readResults(rootFile):

# Read inverted profiles
	ff = io.netcdf_file(rootFile+'.inversion', 'r')
	synthProf = ff.variables['map'][:]
	ff.close()
	
# Read inverted parameters
	ff = io.netcdf_file(rootFile+'.parameters', 'r')
	pars = ff.variables['map'][:]
	ff.close()
	
# Read errors
	ff = io.netcdf_file(rootFile+'.errors', 'r')
	errors = ff.variables['map'][:]
	ff.close()
	
	return synthProf, pars, errors

synthProf, pars, errors = readResults('test')
f = io.netcdf_file('OBSERVATION/testEmission.nc')
obs = f.variables['map'][:]
l = f.variables['lambda'][:]

pl.plot(l, obs[0,:,0])
pl.plot(l, synthProf[0,:,0])