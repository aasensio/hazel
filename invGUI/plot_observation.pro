pro plot_observation, file, widget

	widget_control, widget, GET_VALUE=wind
	
	wset, wind
	stokes = ddread(file,offset=1,/count, /noverbose)
; 	openr,2,file
; 	n = 0
; 	readf,2,n
; 	stokes = dblarr(5,n)
; 	readf,2,stokes
; 	close,2
	
	angstrom = '!3' + STRING(197B) + '!X'
	
	change_symbol, 'circle', 0.3
	
	!p.multi=[0,2,2]
	plot,stokes(0,*),stokes(1,*), xtit='Wavelength shift ['+angstrom+']',ytit='I/I!dmax!n',psym=8,xsty=1
	plot,stokes(0,*),stokes(2,*), xtit='Wavelength shift ['+angstrom+']',ytit='Q/I!dmax!n',psym=8,xsty=1
	plot,stokes(0,*),stokes(3,*), xtit='Wavelength shift ['+angstrom+']',ytit='U/I!dmax!n',psym=8,xsty=1
	plot,stokes(0,*),stokes(4,*), xtit='Wavelength shift ['+angstrom+']',ytit='V/I!dmax!n',psym=8,xsty=1
	!p.multi=0
	
end