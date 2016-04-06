pro test
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
	nbarInput = [0.1,0.1,0.1,0.1]
	omegaInput = [0.21,0.1,0.1,0.1]
	
	wavelengthOutput = dblarr(nLambdaInput)
	stokesOutput = dblarr(4,nLambdaInput)
	epsOutput = dblarr(4,nLambdaInput)
	etaOutput = dblarr(4,4,nLambdaInput)
	
	res = call_external('hazel.so','hazel_',synModeInput, nSlabsInput, B1Input, B2Input, hInput, tau1Input, tau2Input, boundaryInput, $
        transInput, atomicPolInput, anglesInput, lambdaAxisInput, nLambdaInput, dopplerWidthInput, dopplerWidth2Input, dampingInput, $
        dopplerVelocityInput, dopplerVelocity2Input, ffInput, nbarInput, omegaInput, wavelengthOutput, stokesOutput, epsOutput, etaOutput,/verbose)
        
	stop
end