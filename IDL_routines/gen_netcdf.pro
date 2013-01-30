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
; pars is an array of size [npixel,nparameters] that gives the initial value of the parameters
pro gen_netcdf, lambda, stI, stQ, stU, stV, sigmaI, sigmaQ, sigmaU, sigmaV, boundary, height, $
	obs_theta, obs_gamma, mask, pars, outputfile

	npixel = n_elements(stI[*,0])
	nlambda = n_elements(lambda)
	ncol = 8

	map = dblarr(npixel,8,nlambda)

	map[*,0,*] = stI
	map[*,1,*] = stQ
	map[*,2,*] = stU
	map[*,3,*] = stV
	map[*,4,*] = sigmaI
	map[*,5,*] = sigmaQ
	map[*,6,*] = sigmaU
	map[*,7,*] = sigmaV

	dim_map = size(mask, /dimensions)

	file_id = ncdf_create(outputfile, /clobber)
	nx_dim = ncdf_dimdef(file_id, 'npixel', npixel)
	ncol_dim = ncdf_dimdef(file_id, 'ncolumns', ncol)
	nstokespar_dim = ncdf_dimdef(file_id, 'nstokes_par', 4)
	npars_dim = ncdf_dimdef(file_id, 'nparameters', 9)
	nlambda_dim = ncdf_dimdef(file_id, 'nlambda', nlambda)
	nxmap_dim = ncdf_dimdef(file_id, 'nx', reform(dim_map[0]))
	nymap_dim = ncdf_dimdef(file_id, 'ny', reform(dim_map[1]))
	lambda_id = ncdf_vardef(file_id, 'lambda', [nlambda_dim], /double)
	stI_id = ncdf_vardef(file_id, 'map', [nx_dim,ncol_dim,nlambda_dim], /double)
	boundary_id = ncdf_vardef(file_id, 'boundary', [nx_dim, nstokespar_dim], /double)
	height_id = ncdf_vardef(file_id, 'height', [nx_dim], /double)
	obstheta_id = ncdf_vardef(file_id, 'obs_theta', [nx_dim], /double)
	obsgamma_id = ncdf_vardef(file_id, 'obs_gamma', [nx_dim], /double)
	mask_id = ncdf_vardef(file_id, 'mask', [nxmap_dim, nymap_dim], /short)
	parsInit_id = ncdf_vardef(file_id, 'pars_initial', [nx_dim, npars_dim], /double)
	ncdf_control, file_id, /endef

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


; Generate the file with the observations but change the height to compensate for the
; fact that the object is not on the plane of the sky
; angle_away_from_plane -> it is zero if the object is in the plane of the sky. If not,
; this angle is added to 90deg
pro gen_prominence_not_in_sky, file, angle_away_from_plane
	print, 'Writing file ', file

	restore,file+'_red.sav'
	restore,file+'_cal.sav'

	ini_he = 50
	fin_he = 135

; Size of the maps
	nslit = 223
	nstep = n_elements(st[0,0,0,*])

	compute_mu_distance, nslit, nstep, 2*0.175d0, 6.001741d-1, ind_limbo, 0, 58, mu, dist

	save, mu, dist, filename=file+'_artificial_98deg.height.idl'

; Number of pixels
	npixel = nslit*nstep

; Wavelenght axis
	lambda = lambda[ini_he:fin_he] - 10829.0911d0

; Number of wavelengths
	nlambda = n_elements(st[ini_he:fin_he,0,0,0])

; Reform variables
	st = transpose(st[ini_he:fin_he,*,*,*],[1,3,2,0])
	st = reform(st, npixel, 4, nlambda)
	sigma = transpose(reform(sigma, 4, npixel))
	maxi = reform(maxi, npixel)
	mu = reform(mu, npixel)
	dist = reform(dist, npixel)

	save, maxi, filename=file+'_maxi_artificial_98deg.sav'

; Select only pixels with no ambiguity
	ind_signal = ind_senyal
	st = st[ind_signal,*,*]
	maxi = maxi[ind_signal]
	mu = mu[ind_signal]
	dist = dist[ind_signal]
	sigma = sigma[ind_signal,*]

	mask = replicate(1,nstep,nslit)
	mask[ind_signal] = 1

	npixel = n_elements(ind_signal)

; Map with the observations
	map = dblarr(npixel,8,nlambda)

; Fill the map array with the observations
	map[*,0:3,*] = st[*,0:3,*]

; Normalize to the maximum
	for i = 0, 3 do for j = 0, nlambda-1 do map[*,i,j] = map[*,i,j] / maxi

; Noises
	map[*,4,*] = reform(sigma[*,0] / maxi) # replicate(1.d0, nlambda)
	map[*,5,*] = reform(sigma[*,1] / maxi) # replicate(1.d0, nlambda)
	map[*,6,*] = reform(sigma[*,2] / maxi) # replicate(1.d0, nlambda)
	map[*,7,*] = reform(sigma[*,3] / maxi) # replicate(1.d0, nlambda)

; Use smaller noises close to the core of the Stokes Q and U to make
; the fit better there
	map[*,5,48:58] = map[*,5,48:58] * 0.3
	map[*,6,48:58] = map[*,6,48:58] * 0.3
	map[*,7,48:58] = map[*,7,48:58] * 0.3

; Observing angle obtained from mu
	obs_angle = acos(mu) * 180.d0 / !DPI
	obs_gamma = replicate(90.d0,npixel)

; *************************
; Artificially set the prominence at 8 deg away from the local limb	
	obs_angle = obs_angle + angle_prominence

; Boundary condition
 	boundary = fltarr(npixel,4)
 	boundary[*,0] = i0_allen(10830.d0, mu)
 	ind = where(mu eq 0.d0, count)
 	if (count ne 0) then boundary[ind,0] = 0.d0

; Height above the limb
	height = fltarr(npixel)
	ind = where(mu ne 0.d0, count)
	if (count ne 0) then height[ind] = 3.d0
	ind = where(mu eq 0.d0, count)
	if (count ne 0) then height[ind] = dist[ind]

; Modify the height to compensate for the position
	height = (960.d0+height)/cos(angle_prominence*!DPI/180.d0)-960.d0

; Verify that the height is not zero. Hazel will complain about it.
 	ind = where(height lt 0.01, count)
 	if (count ne 0) then height[ind] = 0.01

; Initial parameters. Read them from a previous inversion (only for the thermodynamical parameters)
	file_id = ncdf_open(file+'.pars_onlyI.nc')

   nx_dim = ncdf_dimid(file_id, 'npixel')
   ncol_dim = ncdf_dimid(file_id, 'ncolumns')

   ncdf_diminq, file_id, nx_dim, name, npixel
   ncdf_diminq, file_id, ncol_dim, name, ncol
   map_id = ncdf_varid(file_id, 'map')
	ncdf_varget, file_id, map_id, pars
   ncdf_close, file_id

; Set the field strength to 40 G (saturation)
	pars[*,0] = 40.d0

	gen_netcdf, lambda, reform(map[*,0,*]), reform(map[*,1,*]), reform(map[*,2,*]), reform(map[*,3,*]), $
		reform(map[*,4,*]), reform(map[*,5,*]), reform(map[*,6,*]), reform(map[*,7,*]), boundary,$
		height, obs_angle, obs_gamma, mask, pars, file+'_artificial_98deg.nc'

	save, ind_signal, nslit, nstep, filename=file+'_artificial_98deg.picked.idl'


end


pro prepare_inversions_all
	gen_prominence, '24apr11.006'
	gen_prominence, '24apr11.007'
	gen_prominence, '24apr11.010'
	gen_prominence, '24apr11.012'

	gen_prominence_artificial_increase_height,'24apr11.010'
end
