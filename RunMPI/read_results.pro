pro read_results

	file = 'test_inversion.nc'
   file_id = ncdf_open(file)

	nx_dim = ncdf_dimid(file_id, 'npixel')
	ncol_dim = ncdf_dimid(file_id, 'ncolumns')
	nlambda_dim = ncdf_dimid(file_id, 'nlambda')

	ncdf_diminq, file_id, nx_dim, name, npixel
	ncdf_diminq, file_id, ncol_dim, name, ncol
	ncdf_diminq, file_id, nlambda_dim, name, nlambda

	lambda_id = ncdf_varid(file_id, 'lambda')
	map_id = ncdf_varid(file_id, 'map')

	ncdf_varget, file_id, lambda_id, lambda_syn
	ncdf_varget, file_id, map_id, syn
	ncdf_close, file_id


	file = 'OBSERVATION/test.nc'
   file_id = ncdf_open(file)

	nx_dim = ncdf_dimid(file_id, 'npixel')
	ncol_dim = ncdf_dimid(file_id, 'ncolumns')
	nlambda_dim = ncdf_dimid(file_id, 'nlambda')

	ncdf_diminq, file_id, nx_dim, name, npixel
	ncdf_diminq, file_id, ncol_dim, name, ncol
	ncdf_diminq, file_id, nlambda_dim, name, nlambda

	lambda_id = ncdf_varid(file_id, 'lambda')
	map_id = ncdf_varid(file_id, 'map')

	ncdf_varget, file_id, lambda_id, lambda_obs
	ncdf_varget, file_id, map_id, obs
	ncdf_close, file_id

	file = 'test_parameters.nc'
   file_id = ncdf_open(file)

	nx_dim = ncdf_dimid(file_id, 'npixel')
	ncol_dim = ncdf_dimid(file_id, 'ncolumns')

	ncdf_diminq, file_id, nx_dim, name, npixel
	ncdf_diminq, file_id, ncol_dim, name, ncol

	map_id = ncdf_varid(file_id, 'map')

	ncdf_varget, file_id, map_id, pars
	ncdf_close, file_id

	stop

end