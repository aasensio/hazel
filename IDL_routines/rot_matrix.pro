function rot_matrix, K, M, N, alpha, beta, gamm

	ii = complex(0.d0,1.d0)
	result = exp(-ii*(alpha*M+gamm*N)) * rot_matrix_reduced(K,M,N,beta)
	return, result
end
