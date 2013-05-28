; Generate the appropriate files for running a final inversion of each
; ambiguous case, so that we find the really correct solution
pro generateInversionAmbiguous, thB_amb, chiB_amb, file_in, file_obs, thetaObs

; Read the full file with the observations
	print, 'Reading observed profiles'
	file = file_obs+'.nc'
	file_id = ncdf_open(file)

	nx_dim = ncdf_dimid(file_id, 'npixel')
	ncol_dim = ncdf_dimid(file_id, 'ncolumns')
	nlambda_dim = ncdf_dimid(file_id, 'nlambda')
	nstokespar_dim = ncdf_dimid(file_id, 'nstokes_par')
	npars_dim = ncdf_dimid(file_id, 'nparameters')
	nxmap_dim = ncdf_dimid(file_id, 'nx')
	nymap_dim = ncdf_dimid(file_id, 'ny')

	ncdf_diminq, file_id, nx_dim, name, npixel
	ncdf_diminq, file_id, ncol_dim, name, ncol
	ncdf_diminq, file_id, nlambda_dim, name, nlambda

	lambda_id = ncdf_varid(file_id, 'lambda')
	map_id = ncdf_varid(file_id, 'map')
	boundary_id = ncdf_varid(file_id, 'boundary')
	height_id = ncdf_varid(file_id, 'height')
	obstheta_id = ncdf_varid(file_id, 'obs_theta')
	obsgamma_id = ncdf_varid(file_id, 'obs_gamma')
	mask_id = ncdf_varid(file_id, 'mask')
	parsInit_id = ncdf_varid(file_id, 'pars_initial')

	ncdf_varget, file_id, lambda_id, lambda
	ncdf_varget, file_id, map_id, map
	ncdf_varget, file_id, boundary_id, boundary
	ncdf_varget, file_id, height_id, height
	ncdf_varget, file_id, obstheta_id, obs_theta
	ncdf_varget, file_id, obsgamma_id, obs_gamma
	ncdf_varget, file_id, mask_id, mask
	ncdf_varget, file_id, parsInit_id, pars_old
	
	ncdf_close, file_id

; Read the inverted profiles
	print, 'Reading inverted parameters'
	file = file_in+'.pars.nc'
   file_id = ncdf_open(file)

	nx_dim = ncdf_dimid(file_id, 'npixel')
	ncol_dim = ncdf_dimid(file_id, 'ncolumns')

	ncdf_diminq, file_id, nx_dim, name, npixel
	ncdf_diminq, file_id, ncol_dim, name, ncol

	map_id = ncdf_varid(file_id, 'map')

	ncdf_varget, file_id, map_id, pars
	ncdf_close, file_id

	restore,file_obs+'.picked.idl'

; Now generate three files with the observations with exactly the same structure as
; the original one but changing the values of the parameters
	for i = 1, 3 do begin
		print, FORMAT='(A,I1)', 'Writing ambiguous case ', i
		pars_new = pars
		
		thB = reform(thB_amb[i,*,*])
		chiB = reform(chiB_amb[i,*,*])
		
		pars_new[*,1] = thB[ind_signal]
		pars_new[*,2] = chiB[ind_signal]
		
		file_id = ncdf_create(file_obs+'_autoambig'+strtrim(string(i,FORMAT='(I1)'),2)+'.nc', /clobber)
		nx_dim = ncdf_dimdef(file_id, 'npixel', npixel)
		ncol_dim = ncdf_dimdef(file_id, 'ncolumns', ncol)
		nstokespar_dim = ncdf_dimdef(file_id, 'nstokes_par', 4)
		npars_dim = ncdf_dimdef(file_id, 'nparameters', 9)
		nlambda_dim = ncdf_dimdef(file_id, 'nlambda', nlambda)
		nxmap_dim = ncdf_dimdef(file_id, 'nx', nxmap_dim)
		nymap_dim = ncdf_dimdef(file_id, 'ny', nymap_dim)

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
		ncdf_varput, file_id, parsInit_id, pars_new
		ncdf_close, file_id

; And generate the appropriate configuration files
		n = file_lines('../config_inversion_24apr11.010_artificial_98deg_new.dat')
		str = strarr(n)
		openr,2,'../config_inversion_24apr11.010_artificial_98deg_new.dat'
		readf,2,str
		close,2

		str[41] = "'invert_parameters_ambig.dat'"
		str[32] = "'"+strmid(file_obs,3,strlen(file_obs)-3)+'_autoambig'+strtrim(string(i,FORMAT='(I1)'),2)+".nc'"
		str[35] = "'RESULTS/"+file_in+'_autoambig'+strtrim(string(i,FORMAT='(I1)'),2)+".prof.nc'"
		str[38] = "'RESULTS/"+file_in+'_autoambig'+strtrim(string(i,FORMAT='(I1)'),2)+".pars.nc'"
		
		openw,2,'../config_inversion_'+strtrim(file_in)+'_amb'+strtrim(string(i,FORMAT='(I1)'),2)+'.dat'
		for k = 0, n-1 do printf,2,str[k]
		close,2
		
	endfor

; Finally generate a file to run the inversions
	openw,2,'../run_ambiguities'
	for i = 1, 3 do begin
		printf,2,'mpirun -n 8 ./phazel config_inversion_'+strtrim(file_in)+'_amb'+strtrim(string(i,FORMAT='(I1)'),2)+'.dat'
	endfor
	close,2

	file_chmod, '../run_ambiguities', /u_exec

	stop
	
end