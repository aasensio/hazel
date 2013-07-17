pro test
	synModeInput = 5
	nSlabsInput = 1
	B1Input = [3.d0,80.d0,41.d0]
	B2Input = [0.d0,0.d0,0.d0]
	hInput = 3.d0
	tau1Input = 1.d0
	tau2Input = 0.d0
	boundaryInput  = [0.d0,0.d0,0.d0,0.d0]
	transInput = 1
	atomicPolInput = 1
	anglesInput = [90.d0,0.d0,90.d0]
	lambdaAxisInput = [-3.d0,2.5e0]
	nLambdaInput = 150
	dopplerWidthInput = 6.d0
	dopplerWidth2Input = 0.d0
	dampingInput = 0.d0
	dopplerVelocityInput = 0.d0
	dopplerVelocity2Input = 0.d0
	ffInput = 0.d0
	
	res = call_external('hazel.so','hazel_',B1Input, hInput, tau1Input, boundaryInput, anglesInput, lambdaAxisInput, nLambdaInput, $
						dopplerWidthInput, dampingInput, dopplerVelocityInput,/verbose)
end