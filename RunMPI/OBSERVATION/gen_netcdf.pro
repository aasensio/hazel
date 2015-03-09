; Return the boundary condition given by the Allen CLV
function i0_allen, wl, mu
	ic = ddread('ic.dat',/noverb)
	cl = ddread('cl.dat',/noverb)
   PC = 2.99792458d10
   PH = 6.62606876d-27

; Wavelength in A
   ic(0,*) = 1.d4 * ic(0,*)
; I_lambda to I_nu
   ic(1,*) = 1.d14 * ic(1,*) * (ic(0,*)*1.d-8)^2 / PC

   cl(0,*) = 1.d4 * cl(0,*)

   u = interpol(cl(1,*),cl(0,*),wl)
   v = interpol(cl(2,*),cl(0,*),wl)
   i0 = interpol(ic(1,*),ic(0,*),wl)

   imu = 1.d0 - u - v + u * mu + v * mu^2
   
   return, i0*imu
end



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
; pars is an array of size [nparameters,npixel] that gives the initial value of the parameters
pro gen_netcdf, lambda, stI, stQ, stU, stV, sigmaI, sigmaQ, sigmaU, sigmaV, boundary, height, $
	obs_theta, obs_gamma, mask, pars, outputfile

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


pro test
	dat = ddread('prominence_ApJ_642_554.prof',offset=1,/count,/double)
	dat = congrid(dat,9,100)
	
	nlambda = n_elements(dat[0,*])
	npixel = 4
	map = dblarr(npixel,8,nlambda)
	lambda = reform(dat[0,*])
	for i = 0, npixel-1 do begin
		map[i,*,*] = dat[1:*,*]
	endfor
	
	mu = [1.d0, 0.3d0, 1.d0, 0.3d0]
	obs_theta = acos(mu) * 180.d0 / !DPI
 	boundary = fltarr(4, npixel)
 	boundary[0,*] = i0_allen(10830.d0, mu)
 	height = [3.d0, 5.d0, 3.d0, 5.d0]
 	obs_gamma = replicate(90.d0,npixel)
	mask = replicate(1.0,2,2)

; Initial parameters, depending on the radiative transfer option:
; 1-component (vector of size 8): B, thetaB, chiB, tau, vdop, a, vmac, beta
; 2-component 1+1 with same field (vector of size 10): B, thetaB, chiB, tau1, tau2, vdop, a, vmac1, vmac2, beta
; 2-component 1+1 with different field (vector of size 14): B1, thetaB1, chiB1, B2, thetaB2, chiB2, tau1, tau2, vdop1, vdop2, a, vmac1, vmac2, beta
; 2-component 2 with different field with ff (vector of size 14): B1, thetaB1, chiB1, B2, thetaB2, chiB2, tau1, tau2, vdop1, vdop2, a, vmac1, vmac2, ff

	pars = fltarr(8,npixel)
	
; B
	pars[0,*] = [10.d0, 10.d0, 20.d0, 30.d0]

; thB
	pars[1,*] = replicate(0.d0,npixel)

; chiB
	pars[2,*] = replicate(0.d0,npixel)

; tau
	pars[3,*] = replicate(1.d0,npixel)

; vdop
	pars[4,*] = replicate(6.d0,npixel)

; damp
	pars[5,*] = replicate(0.2d0,npixel)

; vmac
	pars[6,*] = replicate(0.d0,npixel)
	
; beta
	pars[7,*] = replicate(1.d0,npixel)
 	
	gen_netcdf, lambda, reform(map[*,0,*]), reform(map[*,1,*]), reform(map[*,2,*]), reform(map[*,3,*]), $
		reform(map[*,4,*]), reform(map[*,5,*]), reform(map[*,6,*]), reform(map[*,7,*]), boundary,$
		height, obs_theta, obs_gamma, mask, pars, 'test.nc'
	stop
end
