program test
implicit none
	integer :: synModeInput, nSlabsInput, transInput, atomicPolInput, nLambdaInput
	
	real(kind=8), dimension(2) :: lambdaAxisInput
	real(kind=8), dimension(3) :: B1Input, B2Input, anglesInput
	real(kind=8), dimension(4) :: boundaryInput
	real(kind=8) :: hInput, tau1Input, tau2Input, dopplerWidthInput, dopplerWidth2Input, dampingInput, dopplerVelocityInput, dopplerVelocity2Input, ffInput
	
	synModeInput = 5
	nSlabsInput = 1
	B1Input = (/3.0,80.0,41.0/)
	B2Input = (/0.0,0.0,0.0/)
	hInput = 3.d0
	tau1Input = 1.d0
	tau2Input = 0.d0
	boundaryInput  = (/0.0,0.0,0.0,0.0/)
	transInput = 1
	atomicPolInput = 1
	anglesInput = (/90.0,0.0,90.0/)
	lambdaAxisInput = (/-3.d0,2.5d0/)
	nLambdaInput = 150
	dopplerWidthInput = 6.d0
	dopplerWidth2Input = 0.d0
	dampingInput = 0.d0
	dopplerVelocityInput = 0.d0
	dopplerVelocity2Input = 0.d0
	ffInput = 0.d0
		
	call hazel(synModeInput, nSlabsInput, B1Input, B2Input, hInput, tau1Input, tau2Input, boundaryInput, &
		transInput, atomicPolInput, anglesInput, lambdaAxisInput, nLambdaInput, dopplerWidthInput, dopplerWidth2Input, dampingInput, &
		dopplerVelocityInput, dopplerVelocity2Input, ffInput)
end program test