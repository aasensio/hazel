import numpy as np
import pyhazel
import matplotlib.pyplot as pl
import time

synModeInput = 5
nSlabsInput = 1
B1Input = np.asarray([500.0,90.0,0.0])
B2Input = np.asarray([0.0,0.0,0.0])
hInput = 3.e0
tau1Input = 1.e0
tau2Input = 0.e0
boundaryInput  = np.asarray([4.098e-5,0.0,0.0,0.0])
transInput = 1
atomicPolInput = 1
magoptInput = 1
anglesInput = np.asarray([0.0,0.0,0.0])
lambdaAxisInput = np.linspace(-1.5e0,2.5e0,150)
nLambdaInput = 128
dopplerWidthInput = 6.5e0
dopplerWidth2Input = 0.e0
dampingInput = 0.e0
dopplerVelocityInput = 0.e0
dopplerVelocity2Input = 0.e0
ffInput = 0.e0
betaInput = 1.0
beta2Input = 1.0
nbarInput = np.asarray([0.0,0.0,0.0,0.0])
omegaInput = np.asarray([0.0,0.0,0.0,0.0])
nbarInput = np.asarray([1.0,1.0,1.0,0.0])
omegaInput = np.asarray([1.0,1.0,1.0,1.0])
normalization = 0

pyhazel.init()

# Compute the Stokes parameters using many default parameters, using Allen's data
[l, stokes, etaOutput, epsOutput] = pyhazel.synth(synModeInput, nSlabsInput, B1Input, B2Input, hInput, 
                        tau1Input, tau2Input, boundaryInput, transInput, atomicPolInput, magoptInput, anglesInput, 
                        nLambdaInput, lambdaAxisInput, dopplerWidthInput, dopplerWidth2Input, dampingInput, 
                        dopplerVelocityInput, dopplerVelocity2Input, ffInput, betaInput, beta2Input, nbarInput, omegaInput, normalization)


# Now plot the Stokes parameters
labels = ['I/Imax','Q/Imax','U/Imax','V/Imax']

for i in range(4):
	pl.subplot(2,2,i+1)
	pl.plot(l - 10829.0911, stokes[i,:])
	pl.xlabel('Wavelength [A]')
	pl.ylabel(labels[i])
	
pl.tight_layout()

pl.show()

time0 = time.time()
for i in range(100):
    [l, stokes, etaOutput, epsOutput] = pyhazel.synth(synModeInput, nSlabsInput, B1Input, B2Input, hInput, 
                        tau1Input, tau2Input, boundaryInput, transInput, atomicPolInput, magoptInput, anglesInput, 
                        nLambdaInput, lambdaAxisInput, dopplerWidthInput, dopplerWidth2Input, dampingInput, 
                        dopplerVelocityInput, dopplerVelocity2Input, ffInput, betaInput, beta2Input, nbarInput, omegaInput, normalization)
delta = time.time() - time0
print("Time spent = {0} - Synthesis per second = {1}".format(delta, 100.0 / delta))