pro plot_observation_result, file_obs, file_result, widget

	widget_control, widget, GET_VALUE=wind
	
	wset, wind
	
 	openr,2,file_obs
 	n = 0
 	readf,2, n
 	close,2
; 	stokes_obs = dblarr(5,n)
; 	readf,2,stokes_obs
; 	close,2

	stokes_obs = ddread(file_obs,offset=1,/count, /noverbose)
	
	result = file_test(file_result, /ZERO_LENGTH)
	while (result eq 1) do begin
		result = file_test(file_result, /ZERO_LENGTH)
	endwhile
	result = file_lines(file_result)
	
	if (result eq n+1) then begin

		stokes = ddread(file_result, offset=1, /count, /noverbose)
		openr,2,file_result
		n = 0
		readf,2,n,chi2
		close,2
		
; 		stokes = dblarr(5,n)
; 		readf,2,stokes
; 		close,2
		
		angstrom = '!3' + STRING(197B) + '!X'
		
		change_symbol, 'circle', 0.3
		
		cwpal
		!p.multi=[0,2,2]
		plot,stokes_obs(0,*),stokes_obs(1,*), xtit='Wavelength shift ['+angstrom+']',ytit='I/I!dmax!n',psym=8,xsty=1
		oplot,stokes(0,*),stokes(1,*), thick=2, col=2
		plot,stokes_obs(0,*),stokes_obs(2,*), xtit='Wavelength shift ['+angstrom+']',ytit='Q/I!dmax!n',psym=8,xsty=1
		oplot,stokes(0,*),stokes(2,*), thick=2, col=2
		plot,stokes_obs(0,*),stokes_obs(3,*), xtit='Wavelength shift ['+angstrom+']',ytit='U/I!dmax!n',psym=8,xsty=1
		oplot,stokes(0,*),stokes(3,*), thick=2, col=2
		plot,stokes_obs(0,*),stokes_obs(4,*), xtit='Wavelength shift ['+angstrom+']',ytit='V/I!dmax!n',psym=8,xsty=1
		oplot,stokes(0,*),stokes(4,*), thick=2, col=2
		!p.multi=0
		loadct,0,/silent
	endif
	
end