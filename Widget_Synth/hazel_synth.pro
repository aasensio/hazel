;-----------------------------------------
; Call the F90 code which calculates the emerging profiles
;-----------------------------------------
pro synthesize, state, handler, plot_profiles=plot_profiles, texto=texto

	siz = 60
	file = strarr(siz)
	openr,2,'init_parameters.dat'
	tmp = ''
	for i = 0, siz-1 do begin
		readf,2,tmp
		file(i) = tmp
	endfor
	close,2

	Bfield_str = strtrim(string(state.Bfield),2)
	thetaBfield_str = strtrim(string(state.thetaBfield),2)
	chiBfield_str = strtrim(string(state.chiBfield),2)
	if (state.randomazimuth eq 1) then begin
		chiBfield_str = '999.d0'
	endif

	if (state.number_slabs eq 3 or state.number_slabs eq -2) then begin
		Bfield_str2 = strtrim(string(state.Bfield2),2)
		thetaBfield_str2 = strtrim(string(state.thetaBfield2),2)
		chiBfield_str2 = strtrim(string(state.chiBfield2),2)
		file[23] = Bfield_str + '  ' + thetaBfield_str + '  ' + chiBfield_str + ' '+ Bfield_str2 + '  ' + thetaBfield_str2 + '  ' + chiBfield_str2
	endif else begin
		file[23] = Bfield_str + '  ' + thetaBfield_str + '  ' + chiBfield_str
	endelse

	thetaObs_str = strtrim(string(state.thetaObs),2)
	chiObs_str = strtrim(string(state.chiObs),2)
	gammaObs_str = strtrim(string(state.gammaObs),2)
	file[44] = thetaObs_str + '  ' + chiObs_str + '  ' + gammaObs_str

	file[38] = strtrim(string(state.Multiplet,FORMAT='(I1)'),2)

	file[17] = strtrim(string(state.paschen),2)

	file[26] = strtrim(string(state.height),2)
	file[20] = strtrim(string(state.number_slabs),2)

; Add 3
	file[41] = strtrim(string(1-state.effects),2)
	if (state.number_slabs eq 1) then begin
		file[29] = strtrim(string(state.dtau_desired),2)
	endif
	if (state.number_slabs eq 2) then begin
		file[29] = strtrim(string(state.dtau_desired),2)+ ' '+strtrim(string(state.dtau_desired2),2)
	endif
	if (state.number_slabs eq -2 or state.number_slabs eq 3) then begin
		file[29] = strtrim(string(state.dtau_desired),2)+ ' '+strtrim(string(state.dtau_desired2),2)+ ' '+strtrim(string(state.ff),2)
	endif

	file[35] = strtrim(string(state.stokes0[0]),2)+ ' '+strtrim(string(state.stokes0[1]),2)+ ' '+$
		strtrim(string(state.stokes0[2]),2)+ ' '+strtrim(string(state.stokes0[3]),2)

	file[11] = '0'
	file[14] = strtrim(string(state.D2),2)
			
	if (state.D2 ne 0) then begin
		file[11] = '1'
	endif

; Multiterm code
	if (state.which_code eq 0) then begin

; He I
		if (state.which_atom eq 0) then begin
			file_with_atom = 'ATOMS/helium.mod'
		endif

; S I
		if (state.which_atom eq 1) then begin
			file_with_atom = 'ATOMS/sulfur.mod'
		endif

; Na I
		if (state.which_atom eq 2) then begin
			file_with_atom = 'ATOMS/sodium.mod'
		endif
	endif

; Multilevel with HFS
	if (state.which_code eq 1) then begin

; Na I HFS
		if (state.which_atom eq 0) then begin
			file_with_atom = 'ATOMS/sodium_hfs.mod'
		endif
	endif
			
	wl = wavelength_atom_multiplet(file_with_atom,state.Multiplet)

	wleft = -1.d8 * state.waveaxis[0] / wl^2
	wright = -1.d8 * state.waveaxis[1] / wl^2
	file(47) = strtrim(string(wright),2)+ ' '+strtrim(string(wleft),2)+ ' '+$
		strtrim(string(state.waveaxis[2],FORMAT='(I5)'),2)

	PC = 2.99792458d10
	PH = 6.62606876d-27

	nu = PC / (wl*1.d-8)

	wl_str = strtrim(string(wl),2)
	vdopp_str = strtrim(string(state.Doppler),2)
	damp_str = strtrim(string(state.damping),2)
	
	if (state.number_slabs le 2) then begin
		file(50) = wl_str + '  ' + vdopp_str + '  '+ damp_str
	endif
	if (state.number_slabs eq -2 or state.number_slabs eq 3) then begin
		vdopp2_str = strtrim(string(state.Doppler2),2)
		file(50) = wl_str + '  ' + vdopp_str + '  ' + vdopp2_str + '  '+ damp_str
	endif
		

	if (state.number_slabs eq 1) then begin
		file(53) = strtrim(string(state.vel),2)
	endif
	if (state.number_slabs eq 2 or state.number_slabs eq -2 or state.number_slabs eq 3) then begin
		file(53) = strtrim(string(state.vel),2)+ ' '+strtrim(string(state.vel2),2)
	endif

	file(59) = strtrim(string(state.stimulated),2)
	file(56) = strtrim(string(state.magneto_opt),2)

	openw,2,'experiment.dat'
	for i = 0, siz-1 do begin
		printf,2,file(i)
	endfor
	close,2

	file = strarr(54)
	openr,2,'config_inversion.idl'
	tmp = ''
	for i = 0, 53 do begin
		readf,2,tmp
		file(i) = tmp
	endfor
	close,2
	
	file(8) = "'experiment.dat'"
	file(50) = strtrim(string(state.observation),2)

; Multiterm code
	if (state.which_code eq 0) then begin

; He I
		if (state.which_atom eq 0) then begin
			file(5) = "'ATOMS/helium.mod'"
		endif

; S I
		if (state.which_atom eq 1) then begin
			file(5) = "'ATOMS/sulfur.mod'"
		endif

; Na I
		if (state.which_atom eq 2) then begin
			file(5) = "'ATOMS/sodium.mod'"
		endif
	endif

; Multilevel with HFS
	if (state.which_code eq 1) then begin

; Na I HFS
		if (state.which_atom eq 0) then begin
			file(5) = "'ATOMS/sodium_hfs.mod'"
		endif
	endif
					
	openw,2,'config_inversion.dat'
	for i = 0, 53 do begin
		printf,2,file(i)
	endfor
	close,2

; Read the atomic model and modify the value of the multiplication factor in the 10830 A line
	change_factors, file_with_atom, state.factor_10830_nbar, state.factor_10830_omega,$
		state.j10

	if (keyword_set(texto)) then begin
		erase
		xyouts,0.3,0.03,'Calculating '+texto,/normal,charsize=2
		wait,0.01
	endif else begin
		xyouts,0.4,0.03,'Calculating...',/normal,charsize=2
		wait,0.01
	endelse

	if (state.which_code eq 0) then begin
		spawn,'./hazel'
	endif else begin
		print, 'Not available'
	endelse

	if (keyword_set(plot_profiles)) then begin
; Read the emerging profiles
		stokes = ddread('test.inversion',/noverb,offset=1,/countall)

		plot_profiles, stokes, state

		if (state.postcript eq 1) then begin
			file = dialog_pickfile()
			abre_ps,file,/mitad
		endif

		if (state.postcript eq 2) then begin
			file = dialog_pickfile(PATH=state.path_save, GET_PATH=path_save)
			state.path_save = path_save
			widget_control, handler, SET_UVALUE=state
		endif

		plot_profiles, stokes, state
		
		if (state.postcript eq 1) then begin
			txt = 'B='+Bfield_str+' -- !7h!6!dB!n='+thetaBfield_str+$
				' -- !7v!6!dB!n='+chiBfield_str
			xyouts,0.20,0.96,txt,/normal
			txt = '!7h!6!dobs!n='+thetaObs_str+' -- !7v!6!dobs!n='+chiObs_str+$
				' -- !7c!6!dobs!n='+gammaObs_str
			xyouts,0.20,0.46,txt,/normal
			cierra_ps
		endif

		if (state.postcript eq 2) then begin
			openw,2,file,width=132
			printf,2,stokes
			close,2
		endif

	endif
	 
	 save, state, filename='state.idl'
	 
end