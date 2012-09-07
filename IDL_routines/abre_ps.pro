pro abre_ps, nombre, TODO=todo, ENCAPSULATED=encapsulated, COLOR=color, LANDSCAPE=landscape,$
	FONT=font, MITAD=mitad, SMALLFONT=smallfont, CENTER=center, XOFFSET=XOFFSET
	
	if (not keyword_set(XOFFSET)) then begin
		xoffset=0.5
	endif
	if (keyword_set(FONT)) then begin
	 	get_lun,u
	 	a = !p.charsize
	   openw,u,'~/idl/fuente.dat'
	 	printf,u,a
	 	close,u	
		free_lun,u
	 	!p.charsize=1.4
	endif
	if (keyword_set(SMALLFONT)) then begin
	 	get_lun,u
	 	a = !p.charsize
	   openw,u,'~/idl/fuente.dat'
	 	printf,u,a
	 	close,u	
		free_lun,u
	 	!p.charsize=1.0
	endif
	a = pswindow(/cm)
	set_plot,'ps',/copy
	encap = 0
	col = 0
	landsca = 0
	if (keyword_set(ENCAPSULATED)) then encap=1
	if (keyword_set(COLOR)) then col = 1
	if (keyword_set(LANDSCAPE)) then landsca = 1
	if (not keyword_set(TODO)) then begin
	 	  if (keyword_set(MITAD)) then begin
		   	device, filename=nombre, ENCAPSULATED=encap, COLOR=col, $
	 	  			LANDSCAPE=landsca, bits_per_pixel=8,xoffset=XOFFSET,yoffset=14.5,$
	 	  			xsize=20.5-XOFFSET,ysize=14	 	  			
		  endif else if (keyword_set(CENTER)) then begin		  
		   	device, filename=nombre, ENCAPSULATED=encap, COLOR=col, $
	 	  			LANDSCAPE=landsca, bits_per_pixel=8, $
					XSIZE=a.xsize, YSIZE=a.ysize, XOFFSET=a.xoffset, $
    	   		 YOFFSET=a.yoffset
		  endif else begin
		   	device, filename=nombre, ENCAPSULATED=encap, COLOR=col, $
	 	  			LANDSCAPE=landsca, bits_per_pixel=8
		  endelse
	 endif else begin
		device, filename=nombre,xoffset=XOFFSET,yoffset=1.5,xsize=19.5,ysize=26, ENCAPSULATED=encap,$
				COLOR=col, LANDSCAPE=landsca, bits_per_pixel=8
	 endelse
end
