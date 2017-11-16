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
transInput = 1
atomicPolInput = 1
magoptInput = 1
anglesInput = np.asarray([0.0,0.0,90.0])
nLambdaInput = 150
lambdaAxisInput = np.linspace(-1.5e0,2.5e0,nLambdaInput)
dopplerWidthInput = 6.e0
dopplerWidth2Input = 0.e0
dampingInput = 0.e0
dopplerVelocityInput = 0.e0
dopplerVelocity2Input = 0.e0
ffInput = 0.e0
betaInput = 1.0
beta2Input = 1.0
nbarInput = np.asarray([0.1,0.1,0.1,0.1])
omegaInput = np.asarray([0.21,0.1,0.1,0.1])
normalization = 0
boundaryInput  = np.zeros((nLambdaInput,4))
boundaryInput[:,0] = 4.098e-5


pyhazel.init()

# Compute the Stokes parameters using many default parameters, using Allen's data
n = 20
B = np.linspace(0.0,1200.0,n)
inferred = np.zeros((2,n))
for i in range(n):
	B1Input = np.asarray([B[i],0.0,0.0])
	[l, stokes, etaOutput, epsOutput] = pyhazel.synth(synModeInput, nSlabsInput, B1Input, B2Input, hInput, 
                        tau1Input, tau2Input, boundaryInput, transInput, atomicPolInput, magoptInput, anglesInput, 
                        nLambdaInput, lambdaAxisInput, dopplerWidthInput, dopplerWidth2Input, dampingInput, 
                        dopplerVelocityInput, dopplerVelocity2Input, ffInput, betaInput, beta2Input, nbarInput, omegaInput, normalization)

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
pl.show()

# Now plot the Stokes parameters
labels = ['I/Imax','Q/Imax','U/Imax','V/Imax']

