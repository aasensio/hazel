@graphical_disamb
@angle_transformation
@get_ambiguities_general
@polarity
@generateInversionAmbiguous

; Disambiguation of a map based on the saturation regime
; disamb, '24apr11.010_artificial_98deg_new', '../OBSERVATION/24apr11.010_artificial'
pro disamb, file_in, file_obs, thetaObs, load_previous=load_previous, verify_ambiguous=verify_ambiguous, gen_files_inversion=gen_files_inversion


; Read the observations and the results of the inversion
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
	file = file_obs+'.nc'
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

	restore,file_obs+'.picked.idl'

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
	noise_map = fltarr(nslit, nstep, 4)
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
	for i = 0, 3 do begin
		temp2 = obs[*,i+4,0]
		temp[ind_signal] = reform(temp2)
		noise_map[*,*,i] = temp
	endfor
	
; Calculate which points have measurable Stokes V
	yesStokesV = intarr(nslit,nstep)
	for i = 0, nslit-1 do begin
		for j = 0, nstep-1 do begin
			std = stddev(obs_map[i,j,3,*])
 			if (std gt 0.5*noise_map[i,j,3]) then begin
				yesStokesV[i,j] = 1
			endif
		endfor
	endfor
		
	thB_amb = fltarr(8,nslit,nstep)
	chiB_amb = fltarr(8,nslit,nstep)

	thSky_amb = fltarr(8,nslit,nstep)
	chiSky_amb = fltarr(8,nslit,nstep)

; Read the information about the polarity
; Marian computed the polarity incorrectly assuming that the line was in absorption
; Since the line is in emission, we have to multiply it by -1
	restore,'stokesV_ampl_pol.sav'
	pol = -pol
		
	for i = 0, nslit-2 do begin
		for j = 0, nstep-1 do begin
			if (B[i,j] ne 0) then begin

				getAmbiguitiesGeneral, thetaObs, thB[i,j], chiB[i,j], yesStokesV[i,j], thetaB_out, chiB_out, thetaSky_out, chiSky_out, StokesVpolarity=pol[i,j]


				nsol = n_elements(thetaB_out)
				thB_amb[0:nsol-1,i,j] = thetaB_out
				chiB_amb[0:nsol-1,i,j] = chiB_out

				thSky_amb[0:nsol-1,i,j] = thetaSky_out
				chiSky_amb[0:nsol-1,i,j] = chiSky_out

; Synthesize with Hazel to verify that the profiles are indeed ambiguous
				if (keyword_set(verify_ambiguous)) then begin
					boundary = fltarr(4)
					geometry = [thetaObs, 0.d0, 90.d0]

					pars = [5*B[i,j], thB_amb[0,i,j], chiB_amb[0,i,j], h[i,j], tau[i,j], vdop[i,j], damp[i,j], vmac[i,j]]
					synth_hazel, pars, boundary, geometry, Stokes

					if (thB_amb[2,i,j] ne 0.d0 and chiB_amb[2,i,j] ne 0.d0) then begin
						pars = [5*B[i,j], thB_amb[1,i,j], chiB_amb[1,i,j], h[i,j], tau[i,j], vdop[i,j], damp[i,j], vmac[i,j]]
						synth_hazel, pars, boundary, geometry, Stokes2
					endif else begin
						Stokes2 = Stokes * 0.d0
					endelse

					!p.multi = [0,2,2]
 					cgplot, lambda_obs, obs_map[i,j,0,*], psym=1
 					cgoplot, lambda_obs, syn_map[i,j,0,*], col='red'
 					cgoplot, Stokes[0,*], Stokes[1,*], col='blue', thick=4
 					cgoplot, Stokes[0,*], Stokes2[1,*], col='green', thick=2

 					cgplot, lambda_obs, obs_map[i,j,1,*], psym=1
 					cgoplot, lambda_obs, syn_map[i,j,1,*], col='red'
 					cgoplot, Stokes[0,*], Stokes[2,*], col='blue', thick=4
 					cgoplot, Stokes[0,*], Stokes2[2,*], col='green', thick=2

 					cgplot, lambda_obs, obs_map[i,j,2,*], psym=1
 					cgoplot, lambda_obs, syn_map[i,j,2,*], col='red'
 					cgoplot, Stokes[0,*], Stokes[3,*], col='blue', thick=4
 					cgoplot, Stokes[0,*], Stokes2[3,*], col='green', thick=2

 					cgplot, lambda_obs, obs_map[i,j,3,*], psym=1
 					cgoplot, lambda_obs, syn_map[i,j,3,*], col='red'
 					cgoplot, Stokes[0,*], Stokes[4,*], col='blue', thick=4
 					cgoplot, Stokes[0,*], Stokes2[4,*], col='green', thick=2

					stop
				endif
			endif
		endfor
	endfor

; Generate configuration files for the inversion of the ambiguous cases
	if (keyword_set(gen_files_inversion)) then begin
		generateInversionAmbiguous, thB_amb, chiB_amb, file_in, file_obs, thetaObs
	endif
	
	thB_amb[*,222,0] = 0
	thB_amb[*,222,1] = 180
	
	chiB_amb[*,222,0] = -180
	chiB_amb[*,222,1] = 180

	loadct,4
	window,0,xsize=2000,ysize=650
	!p.multi = [0,4,2]
	tvframe, thB_amb[0,*,*], /bar, /sam, tit='!7h!6!dB', charsize=2
	tvframe, thB_amb[1,*,*], /bar, /sam, tit='!7h!6!dB', charsize=2
	tvframe, thB_amb[2,*,*], /bar, /sam, tit='!7h!6!dB', charsize=2
	tvframe, thB_amb[3,*,*], /bar, /sam, tit='!7h!6!dB', charsize=2	

	tvframe, thSky_amb[0,*,*], /bar, /sam, tit='!7H!6!dB', charsize=2
	tvframe, thSky_amb[1,*,*], /bar, /sam, tit='!7H!6!dB', charsize=2
	tvframe, thSky_amb[2,*,*], /bar, /sam, tit='!7H!6!dB', charsize=2
 	tvframe, thSky_amb[3,*,*], /bar, /sam, tit='!7H!6!dB', charsize=2

	window,1,xsize=2000,ysize=650
	tvframe, chiB_amb[0,*,*], /bar, /sam, tit='!7v!6!dB', charsize=2
	tvframe, chiB_amb[1,*,*], /bar, /sam, tit='!7v!6!dB', charsize=2
	tvframe, chiB_amb[2,*,*], /bar, /sam, tit='!7v!6!dB', charsize=2
	tvframe, chiB_amb[3,*,*], /bar, /sam, tit='!7v!6!dB', charsize=2

	tvframe, chiSky_amb[0,*,*], /bar, /sam, tit='!7V!6!dB', charsize=2
	tvframe, chiSky_amb[1,*,*], /bar, /sam, tit='!7V!6!dB', charsize=2
	tvframe, chiSky_amb[2,*,*], /bar, /sam, tit='!7V!6!dB', charsize=2
	tvframe, chiSky_amb[3,*,*], /bar, /sam, tit='!7V!6!dB', charsize=2

	!p.multi = [0,2,2]
	window,2,xsize=850,ysize=650
	tvframe, B<50, /bar, /sam
	tvframe, B * cos(thSky_amb[0,*,*]*!DPI/180.d0)/cos(thSky_amb[1,*,*]*!DPI/180.d0)<50>0, /bar, /sam
	tvframe, B * cos(thSky_amb[0,*,*]*!DPI/180.d0)/cos(thSky_amb[2,*,*]*!DPI/180.d0)<50>0, /bar, /sam
	tvframe, B * cos(thSky_amb[0,*,*]*!DPI/180.d0)/cos(thSky_amb[3,*,*]*!DPI/180.d0)<50>0, /bar, /sam
	!p.multi = 0



; 	stop
; 	
; 	restore,'stokesV_ampl_all.sav'
; 	pol = congrid(mapav10, 223, 60)
; 	indYes = where(pol ne 999)
; 	indNo = where(pol eq 999)
; 	pol[indYes] = pol[indYes] / abs(pol[indYes])
; 	pol[indNo] = 0.d0
; 	
; 	save, thB_amb, chiB_amb, filename='ambiguous_maps.idl'
; 	
; 	graphical_disamb, thB_amb, chiB_amb, thetaB_final, chiB_final, pol, load_previous=load_previous
; 	
; 	save, thetaB_final, chiB_final, filename='final_maps.idl'
	
	
end
