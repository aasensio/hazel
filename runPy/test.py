import numpy as np
from hazel import hazel
import matplotlib.pyplot as pl
import time

out = hazel()

synModeInput = 5
nSlabsInput = 3
B1Input = np.asarray([10.0,70.0,0.0])
B2Input = np.asarray([10.0,70.0,0.0])
hInput = 3.e0
tau1Input = 1.e0
tau2Input = 1.e0
transInput = 1
atomicPolInput = 1
magoptInput = 1
anglesInput = np.asarray([0.0,0.0,0.0])
nLambdaInput = 228
lambdaAxisInput = np.linspace(-1.5e0,2.5e0,nLambdaInput)
dopplerWidthInput = 6.5e0
dopplerWidth2Input = 6.5e0
dampingInput = 0.e0
dopplerVelocityInput = -5.e0
dopplerVelocity2Input = 5.e0
ffInput = 0.e0
betaInput = 1.0
beta2Input = 1.0
nbarInput = np.asarray([0.0,0.0,0.0,0.0])
omegaInput = np.asarray([0.0,0.0,0.0,0.0])
normalization = 0
boundaryInput  = np.zeros((nLambdaInput,4))
boundaryInput[:,0] = 4.098e-5



# Compute the Stokes parameters using many default parameters, using Allen's data
[l, stokes, etaOutput, epsOutput] = out.synth(synModeInput, nSlabsInput, B1Input, B2Input, hInput, 
                        tau1Input, tau2Input, boundaryInput, transInput, atomicPolInput, magoptInput, anglesInput, 
                        nLambdaInput, lambdaAxisInput, dopplerWidthInput, dopplerWidth2Input, dampingInput, 
                        dopplerVelocityInput, dopplerVelocity2Input, ffInput, betaInput, beta2Input, nbarInput, omegaInput, normalization)


synModeInput = 5
nSlabsInput = -2
B1Input = np.asarray([10.0,70.0,0.0])
B2Input = np.asarray([10.0,70.0,0.0])
hInput = 3.e0
tau1Input = 1.e0
tau2Input = 1.e0
transInput = 1
atomicPolInput = 1
magoptInput = 1
anglesInput = np.asarray([0.0,0.0,0.0])
nLambdaInput = 128
lambdaAxisInput = np.linspace(-1.5e0,2.5e0,nLambdaInput)
dopplerWidthInput = 6.5e0
dopplerWidth2Input = 6.5e0
dampingInput = 0.e0
dopplerVelocityInput = -5.e0
dopplerVelocity2Input = 5.e0
ffInput = 0.e0
betaInput = 1.0
beta2Input = 1.0
nbarInput = np.asarray([0.0,0.0,0.0,0.0])
omegaInput = np.asarray([0.0,0.0,0.0,0.0])
normalization = 0
boundaryInput  = np.zeros((nLambdaInput,4))
boundaryInput[:,0] = 4.098e-5



# Compute the Stokes parameters using many default parameters, using Allen's data
[l2, stokes2, etaOutput2, epsOutput2] = out.synth(synModeInput, nSlabsInput, B1Input, B2Input, hInput, 
                        tau1Input, tau2Input, boundaryInput, transInput, atomicPolInput, magoptInput, anglesInput, 
                        nLambdaInput, lambdaAxisInput, dopplerWidthInput, dopplerWidth2Input, dampingInput, 
                        dopplerVelocityInput, dopplerVelocity2Input, ffInput, betaInput, beta2Input, nbarInput, omegaInput, normalization)


# Now plot the Stokes parameters
labels = ['I/Imax','Q/Imax','U/Imax','V/Imax']

pl.close('all')
f, ax = pl.subplots(ncols=2, nrows=2, figsize=(8,6))
ax = ax.flatten()

for i in range(4):
    ax[i].plot(l - 10829.0911, stokes[i,:], 'o')
    ax[i].plot(l2 - 10829.0911, stokes2[i,:])
    ax[i].set_xlabel('Wavelength [A]')
    ax[i].set_ylabel(labels[i])
	
pl.tight_layout()

pl.show()

l, stokes, RF = out.synth_RF(['v1'], 1e-3, synModeInput, nSlabsInput, B1Input, B2Input, hInput, 
                        tau1Input, tau2Input, boundaryInput, transInput, atomicPolInput, magoptInput, anglesInput, 
                        nLambdaInput, lambdaAxisInput, dopplerWidthInput, dopplerWidth2Input, dampingInput, 
                        dopplerVelocityInput, dopplerVelocity2Input, ffInput, betaInput, beta2Input, nbarInput, omegaInput, normalization)
