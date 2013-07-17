import numpy
import hazel
import matplotlib.pyplot as p

synModeInput = 5
nSlabsInput = 1
B1Input = [3.0,80.0,41.0]
B2Input = [0.0,0.0,0.0]
hInput = 3.e0
tau1Input = 1.e0
tau2Input = 0.e0
boundaryInput  = [0.0,0.0,0.0,0.0]
transInput = 1
atomicPolInput = 1
anglesInput = [90.0,0.0,90.0]
lambdaAxisInput = [-1.5e0,2.5e0]
nLambdaInput = 150
dopplerWidthInput = 6.e0
dopplerWidth2Input = 0.e0
dampingInput = 0.e0
dopplerVelocityInput = 0.e0
dopplerVelocity2Input = 0.e0
ffInput = 0.e0

print hazel.hazel.__doc__
# Compute the Stokes parameters using many default parameters		
[l, stokes, eps, eta] = hazel.hazel(B1Input, hInput, tau1Input, boundaryInput, anglesInput, lambdaAxisInput, nLambdaInput, 
						dopplerWidthInput, dampingInput, dopplerVelocityInput)

# Now plot the Stokes parameters
labels = ['I/Imax','Q/Imax','U/Imax','V/Imax']

for i in range(4):
	p.subplot(2,2,i+1)
	p.plot(l - 10829.0911, stokes[i,:])
	p.xlabel('Wavelength [A]')
	p.ylabel(labels[i])
	
p.tight_layout()

p.show()
