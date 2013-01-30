function index, ind_signal, nslit, slit, step
	pos = step * nslit + slit
	return, where(ind_signal eq pos)
end

pro map, file_obs, file_in, postcript=postcript

	print, 'Reading inverted profiles'
	file = file_in+'.prof.nc'
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

	print, 'Reading observed profiles'
	file = '../OBSERVATION/'+file_obs+'.nc'
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

	restore,'../OBSERVATION/'+file_obs+'.picked.idl'

	B = fltarr(nslit, nstep)
	thB = fltarr(nslit, nstep)
	chiB = fltarr(nslit, nstep)
	h = fltarr(nslit, nstep)
	tau = fltarr(nslit, nstep)
	vdop = fltarr(nslit, nstep)
	damp = fltarr(nslit, nstep)
	vmac = fltarr(nslit, nstep)

	for i = 0, 7 do begin
		ind = where(pars[*,i] gt 1.d10, n)
		if (n ne 0) then begin
			pars[ind,i] = 0.d0
			syn[ind,*,*] = 0.d0
		endif
	endfor
	
	B[ind_signal] = reform(pars[*,0])
	thB[ind_signal] = pars[*,1]
	chiB[ind_signal] = pars[*,2]
	h[ind_signal] = pars[*,3]
	tau[ind_signal] = pars[*,4]
	vdop[ind_signal] = pars[*,5]
	damp[ind_signal] = pars[*,6]
	vmac[ind_signal] = pars[*,7]
	
	
	obs_map = fltarr(nslit, nstep, 4, nlambda)
	syn_map = fltarr(nslit, nstep, 4, nlambda)
	temp = fltarr(nslit, nstep)
	for i = 0, 3 do begin
		for j = 0, nlambda-1 do begin
			temp2 = obs[*,i,j]
			temp[ind_signal] = reform(temp2)
			obs_map[*,*,i,j] = temp

			temp2 = syn[*,i,j]
			temp[ind_signal] = reform(temp2)
			syn_map[*,*,i,j] = temp
		endfor
	endfor


	chi2 = fltarr(4,npixel)
	for i = 0, npixel-1 do chi2[0,i] = total( (obs[i,0,*]-syn[i,0,*])^2 / obs[i,4,*]^2) / (1.d0*nlambda)
	for i = 0, npixel-1 do chi2[1,i] = total( (obs[i,1,*]-syn[i,1,*])^2 / obs[i,5,*]^2) / (1.d0*nlambda)
	for i = 0, npixel-1 do chi2[2,i] = total( (obs[i,2,*]-syn[i,2,*])^2 / obs[i,6,*]^2) / (1.d0*nlambda)
	for i = 0, npixel-1 do chi2[3,i] = total( (obs[i,3,*]-syn[i,3,*])^2 / obs[i,7,*]^2) / (1.d0*nlambda)

	chi2I_map = fltarr(nslit, nstep)
	chi2Q_map = fltarr(nslit, nstep)
	chi2U_map = fltarr(nslit, nstep)
	chi2V_map = fltarr(nslit, nstep)
	
	chi2I_map[ind_signal] = reform(chi2[0,*])
	chi2Q_map[ind_signal] = reform(chi2[1,*])
	chi2U_map[ind_signal] = reform(chi2[2,*])
	chi2V_map[ind_signal] = reform(chi2[3,*])

	loadct,4


	window,0,xsize=(size(tau))[1],ysize=(size(tau))[2]
	tvscl, thB

 	window,1,xsize=1200, ysize=600
 
 	wset,0

	!mouse.button = 0
	while (!MOUSE.button NE 4) do begin ;Repeat until the right button is pressed.
		cursor, x, y, /change, /device
		wset,1	  
!p.multi = [0,4,2]
		loadct,0,/silent
		plot, obs_map[x, y, 0, *], psym=1, tit=textoidl('\chi')+'!u2!n='+strtrim(string(chi2I_map[x, y]), 2)
		oplot, syn_map[x, y, 0, *]
		xyouts, 2, 0.9, 'B='+strtrim(string(B[x,y]),2)+' G'
		xyouts, 2, 0.8, 'thB='+strtrim(string(thB[x,y]),2)+' deg'
		xyouts, 2, 0.7, 'chiB='+strtrim(string(chiB[x,y]),2)+' deg'

		plot, obs_map[x, y, 1, *], psym=1, tit=textoidl('\chi')+'!u2!n='+strtrim(string(chi2Q_map[x, y]), 2)
		oplot, syn_map[x, y, 1, *]

		plot, obs_map[x, y, 2, *], psym=1, tit=textoidl('\chi')+'!u2!n='+strtrim(string(chi2U_map[x, y]), 2)
		oplot, syn_map[x, y, 2, *]

		plot, obs_map[x, y, 3, *], psym=1, tit=textoidl('\chi')+'!u2!n='+strtrim(string(chi2V_map[x, y]), 2)
		oplot, syn_map[x, y, 3, *]

		plot, obs_map[x,y,0,*] - syn_map[x,y,0,*], psym=1, yran=[-5.d-3, 5.d-3], $
			tit=textoidl('\sigma')+'='+strtrim(string(stddev(obs_map[x,y,0,*] - syn_map[x,y,0,*])), 2)
		very, 1.d-3, line=2
		very, 3.d-3, line=2
		very, -1.d-3, line=2
		very, -3.d-3, line=2
		
		plot, obs_map[x,y,1,*] - syn_map[x,y,1,*], psym=1, yran=[-5.d-3, 5.d-3], $
			tit=textoidl('\sigma')+'='+strtrim(string(stddev(obs_map[x,y,1,*] - syn_map[x,y,1,*])), 2)
		very, 1.d-3, line=2
		very, 3.d-3, line=2
		very, -1.d-3, line=2
		very, -3.d-3, line=2
		
		plot, obs_map[x,y,2,*] - syn_map[x,y,2,*], psym=1, yran=[-5.d-3, 5.d-3], $
			tit=textoidl('\sigma')+'='+strtrim(string(stddev(obs_map[x,y,2,*] - syn_map[x,y,2,*])), 2)
		very, 1.d-3, line=2
		very, 3.d-3, line=2
		very, -1.d-3, line=2
		very, -3.d-3, line=2
		
		plot, obs_map[x,y,3,*] - syn_map[x,y,3,*], psym=1, yran=[-5.d-3, 5.d-3], $
			tit=textoidl('\sigma')+'='+strtrim(string(stddev(obs_map[x,y,3,*] - syn_map[x,y,3,*])), 2)
		very, 1.d-3, line=2
		very, 3.d-3, line=2
		very, -1.d-3, line=2
		very, -3.d-3, line=2
		
		!p.multi = 0
		wset,0
	endwhile

	wdelete,0
	wdelete,1
	
end

pro test
	map, '24apr11.010'
end
