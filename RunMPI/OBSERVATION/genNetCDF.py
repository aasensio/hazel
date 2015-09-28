import numpy as np
from scipy.io import netcdf

def i0Allen(wavelength, muAngle):
	"""
	Return the solar intensity at a specific wavelength and heliocentric angle
	wavelength: wavelength in angstrom
	muAngle: cosine of the heliocentric angle
	"""
	C = 2.99792458e10
	
	lambdaIC = 1e4 * np.asarray([0.20,0.22,0.245,0.265,0.28,0.30,0.32,0.35,0.37,0.38,0.40,0.45,0.50,0.55,0.60,0.80,1.0,1.5,2.0,3.0,5.0,10.0])
	uData = np.asarray([0.12,-1.3,-0.1,-0.1,0.38,0.74,0.88,0.98,1.03,0.92,0.91,0.99,0.97,0.93,0.88,0.73,0.64,0.57,0.48,0.35,0.22,0.15])
	vData = np.asarray([0.33,1.6,0.85,0.90,0.57, 0.20, 0.03,-0.1,-0.16,-0.05,-0.05,-0.17,-0.22,-0.23,-0.23,-0.22,-0.20,-0.21,-0.18,-0.12,-0.07,-0.07])

	lambdaI0 = 1e4 * np.asarray([0.20,0.22,0.24,0.26,0.28,0.30,0.32,0.34,0.36,0.37,0.38,0.39,0.40,0.41,0.42,0.43,0.44,0.45,0.46,0.48,0.50,0.55,0.60,0.65,0.70,0.75,\
		0.80,0.90,1.00,1.10,1.20,1.40,1.60,1.80,2.00,2.50,3.00,4.00,5.00,6.00,8.00,10.0,12.0])
	I0 = np.asarray([0.06,0.21,0.29,0.60,1.30,2.45,3.25,3.77,4.13,4.23,4.63,4.95,5.15,5.26,5.28,5.24,5.19,5.10,5.00,4.79,4.55,4.02,3.52,3.06,2.69,2.28,2.03,\
		1.57,1.26,1.01,0.81,0.53,0.36,0.238,0.160,0.078,0.041,0.0142,0.0062,0.0032,0.00095,0.00035,0.00018])
	I0 *= 1e14 * (lambdaI0 * 1e-8)**2 / C

	u = np.interp(wavelength, lambdaIC, uData)
	v = np.interp(wavelength, lambdaIC, vData)
	i0 = np.interp(wavelength, lambdaI0, I0)
	
	return (1.0 - u - v + u * muAngle + v * muAngle**2)* i0

def genNetCDF(wavelength, stI, stQ, stU, stV, sigmaI, sigmaQ, sigmaU, sigmaV, boundary, height, obsTheta, obsGamma, mask, pars, outputFile):
	"""
	This routine generates a NetCDF file with the observations ready for Hazel-MPI
	
	Args:
	    wavelength (float): array of size [nlambda]
	    stI (float): array of size [npixel, nlambda] with Stokes I
	    stQ (float): array of size [npixel, nlambda] with Stokes Q
	    stU (float): array of size [npixel, nlambda] with Stokes U
	    stV (float): array of size [npixel, nlambda] with Stokes V
	    sigmaI (float): array of size [npixel, nlambda] with the noise in Stokes I
	    sigmaQ (float): array of size [npixel, nlambda] with the noise in Stokes Q
	    sigmaU (float): array of size [npixel, nlambda] with the noise in Stokes U
	    sigmaV (float): array of size [npixel, nlambda] with the noise in Stokes V
	    boundary (float): array of size [npixel, 4] with the boundary conditions [I0,Q0,U0,V0] for every pixel
	    height (float): array of size [npixel] indicating the height of the pixel over the surface in arcsec
	    obsTheta (float): array of size [npixel] indicating the angle of the observation in degrees 
	    obsGamma (float): array of size [npixel] the angle of the reference for Stokes Q
	    mask (float): array of the original dimensions of the observations that is used later to reconstruct the inverted maps [nx,ny]
	    pars (float): array of size [nparameters,npixel] that gives the initial value of the parameters
	    					The size depends on the radiative transfer option:
								* 1-component (vector of size 8): B, thetaB, chiB, tau, vdop, a, vmac, beta
								* 2-component 1+1 with same field (vector of size 10): B, thetaB, chiB, tau1, tau2, vdop, a, vmac1, vmac2, beta
								* 2-component 1+1 with different field (vector of size 14): B1, thetaB1, chiB1, B2, thetaB2, chiB2, tau1, tau2, vdop1, vdop2, a, vmac1, vmac2, beta
								* 2-component 2 with different field with ff (vector of size 14): B1, thetaB1, chiB1, B2, thetaB2, chiB2, tau1, tau2, vdop1, vdop2, a, vmac1, vmac2, ff
	    outputFile (float): output file
	"""
	
	nPixel, nLambda = stI.shape
	nCols = 8
	
	obsMap = np.zeros((8,nLambda,nPixel))
	
	obsMap[0,:,:] = stI.T
	obsMap[1,:,:] = stQ.T
	obsMap[2,:,:] = stU.T
	obsMap[3,:,:] = stV.T
	obsMap[4,:,:] = sigmaI.T
	obsMap[5,:,:] = sigmaQ.T
	obsMap[6,:,:] = sigmaU.T
	obsMap[7,:,:] = sigmaV.T
	
	obsMap = np.transpose(obsMap,axes=(2,1,0))
	
	dimMap = mask.shape
				
# Variable dimensions
	fileID = netcdf.netcdf_file(outputFile, 'w')
	nPixDim = fileID.createDimension('npixel', nPixel)
	nColDim = fileID.createDimension('ncolumns', nCols)
	nStokesParDim = fileID.createDimension('nstokes_par', 4)
	nParsDim = fileID.createDimension('nparameters', pars.shape[0])
	nLambdaDim = fileID.createDimension('nlambda', nLambda)
	nXDim = fileID.createDimension('nx', dimMap[0])
	nYDim = fileID.createDimension('ny', dimMap[1])

# Variable definition
# Remember that variables are written in C format, so that they are reversed with respect to Fortran
	lambdaID = fileID.createVariable('lambda','f8',('nlambda',))
	stokesID = fileID.createVariable('map','f8',('npixel','nlambda','ncolumns',))
	boundaryID = fileID.createVariable('boundary','f8',('npixel','nstokes_par',))
	heightID = fileID.createVariable('height','f8',('npixel',))
	obsThetaID = fileID.createVariable('obs_theta','f8',('npixel',))
	obsGammaID = fileID.createVariable('obs_gamma','f8',('npixel',))
	maskID = fileID.createVariable('mask','i2',('nx','ny',))
	parsInitID = fileID.createVariable('pars_initial','f8',('npixel','nparameters',))


	lambdaID[:] = wavelength
	stokesID[:] = obsMap
	boundaryID[:] = boundary
	heightID[:] = height
	obsThetaID[:] = obsTheta
	obsGammaID[:] = obsGamma
	maskID[:] = mask
	parsInitID[:] = pars.T
	
	fileID.close()