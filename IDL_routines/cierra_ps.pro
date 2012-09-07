pro cierra_ps, FONT=font

; Reset the postcript window	

	device,/inches,xoffset=3./4,yoffset=5,xsize=7,ysize=5
	device,/close
	set_plot,'x'
	if (keyword_set(FONT)) then begin
	 	get_lun,u
	    openr,u,'~/idl/fuente.dat'
	 	readf,u,a
	 	close,u
		free_lun,u
	 	!p.charsize=a
	endif
end
