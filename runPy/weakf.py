import numpy as np
import pyhazel
import matplotlib.pyplot as pl
from ipdb import set_trace as stop

synModeInput = 5
nSlabsInput = 1
B1Input = np.asarray([300.0,20.0,41.0])
B2Input = np.asarray([0.0,0.0,0.0])
hInput = 3.e0
tau1Input = 1.e0
tau2Input = 0.e0
boundaryInput  = np.asarray([4.098e-5,0.0,0.0,0.0])
transInput = 1
atomicPolInput = 1
anglesInput = np.asarray([0.0,0.0,90.0])
lambdaAxisInput = np.asarray([-1.5e0,2.5e0])
nLambdaInput = 150
dopplerWidthInput = 6.e0
dopplerWidth2Input = 0.e0
dampingInput = 0.e0
dopplerVelocityInput = 0.e0
dopplerVelocity2Input = 0.e0
ffInput = 0.e0
nbarInput = np.asarray([0.1,0.1,0.1,0.1])
omegaInput = np.asarray([0.21,0.1,0.1,0.1])

pyhazel.init()

# Compute the Stokes parameters using many default parameters, using ad-hoc anisotropy and number of photons per mode
#[l, stokes, etaOutput, epsOutput] = hazel.hazel(B1Input, hInput, tau1Input, boundaryInput, anglesInput, lambdaAxisInput, nLambdaInput, 
						#dopplerWidthInput, dampingInput, dopplerVelocityInput, nbarinput=nbar, omegainput=omega)

# Compute the Stokes parameters using many default parameters, using Allen's data
B = np.linspace(0.0,1200.0,100)
inferred = np.zeros((2,100))
for i in range(100):
	B1Input = np.asarray([B[i],0.0,0.0])
	[l, stokes, etaOutput, epsOutput] = pyhazel.synth(synModeInput, nSlabsInput, B1Input, B2Input, hInput, 
	                        tau1Input, tau2Input, boundaryInput, transInput, atomicPolInput, anglesInput, 
	                        lambdaAxisInput, nLambdaInput, dopplerWidthInput, dopplerWidth2Input, dampingInput, 
	                        dopplerVelocityInput, dopplerVelocity2Input, ffInput, nbarInput, omegaInput)

	deltaL = l[1]-l[0]
	stokesI = stokes[0,80:]
	stokesV = stokes[3,80:]
	dIdl = np.gradient(stokesI, deltaL)
	geff = 1.42
	alpha = -4.67e-13 * geff * 10830.0**2
	inferred[0,i] = 1.0 / alpha * np.sum(stokesV * dIdl) / np.sum(dIdl**2)
	inferred[1,i] = B1Input[0] * np.cos(B1Input[1]*np.pi/180.)

pl.plot(inferred[1,:], inferred[0,:])
pl.plot(np.arange(1200),np.arange(1200), '--')

# Now plot the Stokes parameters
labels = ['I/Imax','Q/Imax','U/Imax','V/Imax']

