;-----------------------------------------
; Call the F90 code to calculate the variation of the density matrix
;-----------------------------------------
pro calculate_field_variation, state
	
	state2 = state
	nfield = state.bfield_var[2]	
	ind = findgen(nfield)/(nfield-1.d0) * (alog10(state.bfield_var[1])-$
		alog10(state.bfield_var[0])) + alog10(state.bfield_var[0])
	field = 10.d0^ind
	
	if (state.which_code eq 0) then begin
		prefix = 'J='
	endif else begin
		prefix = 'F='
	endelse
	
	for i = 0, nfield-1 do begin
		state2.Bfield = field[i]
		state2.waveaxis[2] = 2
		texto = 'B='+strtrim(string(field[i]),2)
		synthesize, state2, handler, texto=texto
		
		rotate_refsystem, 'tanti.res'
		
		if (state.which_refframe eq 0) then begin
			temp = ddread('tanti_vertical.res')
		endif else begin
			temp = ddread('tanti_magnetic.res')
		endelse
		
		if (i eq 0) then begin
			n = n_elements(temp[0,*])
	 		rho = fltarr(nfield,8,n)
		endif
		rho[i,*,*] = temp
	endfor
	nterm = n_elements(uniq(temp[1,*]))
	!p.multi=[0,nterm,3]
	
; sigma20(F,F) for all terms
	for i = 0, nterm-1 do begin
		ind20 = where(rho[0,1,*] eq i+1 and rho[0,2,*] eq rho[0,3,*] and $
			rho[0,4,*] eq 2 and rho[0,5,*] eq 0)
		ind00 = where(rho[0,1,*] eq i+1 and rho[0,2,*] eq rho[0,3,*] and $
			rho[0,4,*] eq 0 and rho[0,5,*] eq 0)
					
		if (ind20[0] ne -1) then begin
			n = n_elements(ind20)
			sigma = rho[*,6,ind20]/rho[*,6,ind00]
			maxim = max(sigma)
			minim = min(sigma)						
			for j = 0, n-1 do begin					
				if (j eq 0) then begin
					plot,field,rho[*,6,ind20[j]]/rho[*,6,ind00[j]],/xlog,$
						xran=[state.bfield_var[0],state.bfield_var[1]],$
						yran=[minim,maxim],xsty=1						
				endif else begin
					oplot,field,rho[*,6,ind20[j]]/rho[*,6,ind00[j]],line=j
				endelse
				xyouts,1.d-3,rho[0,6,ind20[j]]/rho[0,6,ind00[j]]-(maxim-minim)*0.05,$
					prefix+strtrim(string(rho[0,2,ind20[j]]/2,FORMAT='(F3.1)'),2),/data,$
					charsize=1.0
			endfor			
		endif		
	endfor
	
; sigma20(F,F') for all terms
	cwpal
	tit = ['Real','Imaginary']
	if (state.which_rho_plot eq 0) then begin
		for imagin = 0, 1 do begin
			for i = 0, nterm-1 do begin
				ind20 = where(rho[0,1,*] eq i+1 and rho[0,2,*] ne rho[0,3,*] and $
					rho[0,4,*] eq 2 and rho[0,5,*] eq 0)
				if (ind20[0] ne -1) then begin
					n = n_elements(ind20)
									
					sigma = rho[*,6+imagin,ind20]
					maxim = max(sigma)
					minim = min(sigma)
				
					for j = 0, n-1 do begin
						if (j eq 0) then begin
							maxloc = fix(alog10(max(abs(rho[*,6+imagin,ind20]))))
							plot,field,rho[*,6+imagin,ind20[j]]/10.d0^maxloc,/xlog,$
								xran=[state.bfield_var[0],state.bfield_var[1]],$
								yran=[minim,maxim]/10.d0^maxloc,xsty=1,title=tit[imagin]+$
								'  x10!u'+strtrim(string(maxloc,FORMAT='(I2)'),2)
						endif else begin
							oplot,field,rho[*,6+imagin,ind20[j]]/10.d0^maxloc,line=j
						endelse
						xyouts,1.d-3,(rho[0,6+imagin,ind20[j]]-(maxim-minim)*0.05)/10.d0^maxloc,$
							prefix+'('+strtrim(string(rho[0,2,ind20[j]]/2,FORMAT='(F3.1)'),2)+','+$
							strtrim(string(rho[0,3,ind20[j]]/2,FORMAT='(F3.1)'),2)+')',/data,$
							charsize=1.0
					endfor			
				endif		
			endfor
		endfor
	endif
	
; sigma2Q(F,F) for all terms
	if (state.which_rho_plot eq 1) then begin
		for imagin = 0, 1 do begin
			for i = 0, nterm-1 do begin
				ind00 = where(rho[0,1,*] eq i+1 and rho[0,2,*] eq rho[0,3,*] and $
					rho[0,4,*] eq 0 and rho[0,5,*] eq 0)			
				ind21 = where(rho[0,1,*] eq i+1 and rho[0,2,*] eq rho[0,3,*] and $
					rho[0,4,*] eq 2 and rho[0,5,*] eq 1)
				ind22 = where(rho[0,1,*] eq i+1 and rho[0,2,*] eq rho[0,3,*] and $
					rho[0,4,*] eq 2 and rho[0,5,*] eq 2)
				if (ind21[0] ne -1) then begin
					n = n_elements(ind21)
				
					sigma21 = rho[*,6+imagin,ind21]/rho[*,6,ind00]
					sigma22 = rho[*,6+imagin,ind22]/rho[*,6,ind00]
					maxim = max([sigma21,sigma22])
					minim = min([sigma21,sigma22])
									
					for j = 0, n-1 do begin
						if (j eq 0) then begin							
							plot,field,rho[*,6+imagin,ind21[j]]/rho[*,6,ind00[j]],/xlog,$
								xran=[state.bfield_var[0],state.bfield_var[1]],$
								yran=[minim,maxim],xsty=1,title=tit[imagin]
							if (imagin eq 1 and i eq 0) then begin
								legend,['!7r!6!u2!d1!n','!7r!6!u2!d2!n'],psym=[88,88],col=[1,2],$
									charsize=1.3,box=0,spacing=0,/right,/bottom
							endif
						endif else begin
							oplot,field,rho[*,6+imagin,ind21[j]]/rho[*,6,ind00[j]],line=j
						endelse						
					endfor
					
					for j = 0, n-1 do begin						
						oplot,field,rho[*,6+imagin,ind22[j]]/rho[*,6,ind00[j]],line=j,col=2
					endfor
					
				endif		
			endfor
		endfor
	endif
	loadct,0
	
	!p.multi=0	
end