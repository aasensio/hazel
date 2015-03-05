import numpy as np
import matplotlib.pyplot as pl
from netCDF4 import Dataset as nf

def readResults(rootFile):

# Read inverted profiles
	ff = nf(rootFile+'.inversion', 'r')
	synthProf = ff.variables['map'][:]
	ff.close()
	
# Read inverted parameters
	ff = nf(rootFile+'.parameters', 'r')
	pars = ff.variables['map'][:]
	ff.close()
	
# Read errors
	ff = nf(rootFile+'.errors', 'r')
	errors = ff.variables['map'][:]
	ff.close()
	
	return synthProf, pars, errors

synthProf, pars, errors = readResults('test')