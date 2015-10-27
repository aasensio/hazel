; This routine generates a NetCDF file with the
; observations ready for Hazel-MPI
; stI, stQ, stU, stV are arrays of size [npixel, nlambda]
; sigmaI, sigmaQ, sigmaU, sigmaV are arrays of size [npixel, nlambda]
; lambda is an array of size [nlambda]
; boundary is an array of size [npixel, 4] with the boundary conditions [I0,Q0,U0,V0] for every pixel
; height is an array of size [npixel] indicating the height of the pixel over the surface in arcsec
; obs_theta is an array of size [npixel] indicating the angle of the observation in degrees
; obs_gamma is an array of size [npixel] indicating the angle of the reference for Stokes Q
; mask is an array of the original dimensions of the observations that is used later to
;   reconstruct the inverted maps [nx,ny]
; normalization: 'continuum' or 'peak' indicating how the profiles are normalized
; pars is an array of size [nparameters,npixel] that gives the initial value of the parameters
pro gen_netcdf, lambda, stI, stQ, stU, stV, sigmaI, sigmaQ, sigmaU, sigmaV, boundary, height, $
	obs_theta, obs_gamma, mask, pars, normalization, outputfile

	npixel = n_elements(stI[*,0])
	nlambda = n_elements(lambda)
	ncol = 8

	map = dblarr(8,nlambda,npixel)

	map[0,*,*] = transpose(stI)
	map[1,*,*] = transpose(stQ)
	map[2,*,*] = transpose(stU)
	map[3,*,*] = transpose(stV)
	map[4,*,*] = transpose(sigmaI)
	map[5,*,*] = transpose(sigmaQ)
	map[6,*,*] = transpose(sigmaU)
	map[7,*,*] = transpose(sigmaV)
	boundary = transpose(boundary)
	
	dim_map = size(mask, /dimensions)
			
; Variable dimensions
	file_id = ncdf_create(outputfile, /clobber)
	npix_dim = ncdf_dimdef(file_id, 'npixel', npixel)
	ncol_dim = ncdf_dimdef(file_id, 'ncolumns', ncol)
	nstokespar_dim = ncdf_dimdef(file_id, 'nstokes_par', 4)
	npars_dim = ncdf_dimdef(file_id, 'nparameters', n_elements(pars[*,0]))
	nlambda_dim = ncdf_dimdef(file_id, 'nlambda', nlambda)
	nxmap_dim = ncdf_dimdef(file_id, 'nx', reform(dim_map[0]))
	nymap_dim = ncdf_dimdef(file_id, 'ny', reform(dim_map[1]))
	
; Variable definition
	normalization_id = ncdf_vardef(file_id, 'normalization', /double)
	lambda_id = ncdf_vardef(file_id, 'lambda', [nlambda_dim], /double)
	stI_id = ncdf_vardef(file_id, 'map', [ncol_dim,nlambda_dim,npix_dim], /double)
	boundary_id = ncdf_vardef(file_id, 'boundary', [nstokespar_dim, npix_dim], /double)
	height_id = ncdf_vardef(file_id, 'height', [npix_dim], /double)
	obstheta_id = ncdf_vardef(file_id, 'obs_theta', [npix_dim], /double)
	obsgamma_id = ncdf_vardef(file_id, 'obs_gamma', [npix_dim], /double)
	mask_id = ncdf_vardef(file_id, 'mask', [nxmap_dim, nymap_dim], /short)
	parsInit_id = ncdf_vardef(file_id, 'pars_initial', [npars_dim, npix_dim], /double)	
	ncdf_control, file_id, /endef
	
; Variable write	
	if (normalization eq 'continuum') then begin
		ncdf_varput, file_id, normalization_id, 0.d0
	endif
	if (normalization eq 'peak') then begin
		ncdf_varput, file_id, normalization_id, 1.d0
	endif
	ncdf_varput, file_id, lambda_id, lambda
	ncdf_varput, file_id, stI_id, map
	ncdf_varput, file_id, boundary_id, boundary
	ncdf_varput, file_id, height_id, height
	ncdf_varput, file_id, obstheta_id, obs_theta
	ncdf_varput, file_id, obsgamma_id, obs_gamma
	ncdf_varput, file_id, mask_id, mask
	ncdf_varput, file_id, parsInit_id, pars
	ncdf_close, file_id		
end