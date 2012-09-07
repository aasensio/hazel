; This returns the reduced rotation matrix
function rot_matrix_reduced, K, M, N, beta

	max_t = min([K+M,K-N])
	min_t = max([0,M-N])

	n_t = abs(max_t-min_t+1)

	if (n_t gt 1) then begin
		t = findgen(n_t)/(n_t-1.d0)*(max_t-min_t)+min_t
	endif else begin
		t = findgen(1) + min_t
	endelse

	cos_beta = cos(beta/2.d0)
	sin_beta = sin(beta/2.d0)

	factor1 = sqrt(factorial(K+M)*factorial(K-M)*factorial(K+N)*factorial(K-N))
	factor2 = factorial(K+M-t)*factorial(K-N-t)*factorial(t)*factorial(t+N-M)

	suma = total((-1.d0)^t*cos_beta^(2*K+M-N-2*t)*sin_beta^(2*t-M+N) / factor2)
	result = factor1 * suma
	return, result
end
