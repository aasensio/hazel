import numpy as np
import pyhazel
import matplotlib.pyplot as pl

synModeInput = 5
nSlabsInput = 1
B1Input = np.asarray([3.0,80.0,41.0])
B2Input = np.asarray([0.0,0.0,0.0])
hInput = 3.e0
tau1Input = 1.e0
tau2Input = 0.e0
boundaryInput  = np.asarray([0.0,0.0,0.0,0.0])
transInput = 1
atomicPolInput = 1
anglesInput = np.asarray([90.0,0.0,90.0])
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
[l, stokes, etaOutput, epsOutput] = pyhazel.synth(synModeInput, nSlabsInput, B1Input, B2Input, hInput, 
                        tau1Input, tau2Input, boundaryInput, transInput, atomicPolInput, anglesInput, 
                        lambdaAxisInput, nLambdaInput, dopplerWidthInput, dopplerWidth2Input, dampingInput, 
                        dopplerVelocityInput, dopplerVelocity2Input, ffInput, nbarInput, omegaInput)

# Now plot the Stokes parameters
labels = ['I/Imax','Q/Imax','U/Imax','V/Imax']

for i in range(4):
	pl.subplot(2,2,i+1)
	pl.plot(l - 10829.0911, stokes[i,:])
	pl.xlabel('Wavelength [A]')
	pl.ylabel(labels[i])
	
pl.tight_layout()

pl.show()
