from numpy cimport ndarray as ar
from numpy import empty

cdef extern from "pyhazel.h":
	void c_hazel(int* synModeInput, int* nSlabsInput, double* B1Input, double* B2Input, double* hInput, double* tau1Input, double* tau2Input, 
		double* boundaryInput, int* transInput, int* atomicPolInput, double* anglesInput, double* lambdaAxisInput, int* nLambdaInput, 
		double* dopplerWidthInput, double* dopplerWidth2Input, double* dampingInput, double* dopplerVelocityInput, 
		double* dopplerVelocity2Input, double* ffInput, double* nbarInput, double* omegaInput, 
		double* wavelengthOutput, double* stokesOutput, double* epsOutput, double* etaOutput)
		
	void c_init()

def synth(int synModeInput, int nSlabsInput, ar[double,ndim=1] B1Input, ar[double,ndim=1] B2Input, double hInput, 
	double tau1Input, double tau2Input, 
	ar[double,ndim=1] boundaryInput, int transInput, int atomicPolInput, ar[double,ndim=1] anglesInput, 
	ar[double,ndim=1] lambdaAxisInput, int nLambdaInput, 
	double dopplerWidthInput, double dopplerWidth2Input, double dampingInput, double dopplerVelocityInput, 
	double dopplerVelocity2Input, double ffInput, ar[double,ndim=1] nbarInput, ar[double,ndim=1] omegaInput):
	
	cdef:
		ar[double,ndim=1] wavelengthOutput = empty(nLambdaInput, order='F')
		ar[double,ndim=2] stokesOutput = empty((4,nLambdaInput), order='F')
		ar[double,ndim=2] epsOutput = empty((4,nLambdaInput), order='F')
		ar[double,ndim=3] etaOutput = empty((4,4,nLambdaInput), order='F')
   
	c_hazel(&synModeInput, &nSlabsInput, &B1Input[0], &B2Input[0], &hInput, &tau1Input, &tau2Input, 
		&boundaryInput[0], &transInput, &atomicPolInput, &anglesInput[0], &lambdaAxisInput[0], &nLambdaInput, 
		&dopplerWidthInput, &dopplerWidth2Input, &dampingInput, &dopplerVelocityInput, 
		&dopplerVelocity2Input, &ffInput, &nbarInput[0], &omegaInput[0], <double*> wavelengthOutput.data, 
		<double*> stokesOutput.data, <double*> epsOutput.data, <double*> etaOutput.data)
    
	return wavelengthOutput, stokesOutput, epsOutput, etaOutput
	
def init():
	c_init()
