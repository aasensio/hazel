function set_imaginary, A, b
	return, complex(float(A),b)
end

function set_real, A, b
	return, complex(b,imaginary(A))
end

;-----------------------------------------
; Write two files containing the rho^K_Q on the vertical and the magnetic
; reference frame
;-----------------------------------------
pro rotate_refsystem, file

	thb = 10.d0 * !DPI / 180.d0
	chb = 90.d0 * !DPI / 180.d0
	n = 1
	
	openw,2,'tanti_vertical.res',width=132
	openw,3,'tanti_magnetic.res',width=132
	
	a = ddread(file,/noverb)
	ind_term = uniq(a[1,*])
	nterm = n_elements(ind_term)
	
	for term = 0, nterm-1 do begin
		ind = where(a[1,*] eq a[1,ind_term[term]])
		nrhos = n_elements(ind)
		J2max = max([a[2,ind],a[3,ind]])
		J2min = min([a[2,ind],a[3,ind]])		
		kmax = max(a[4,ind])
		kmin = min(a[4,ind])
		qmax = max(a[5,ind])
		rho = complexarr(J2max+1,J2max+1,kmax+1,2*qmax+1)
		for i = 0, nrhos-1 do begin
			loc = ind[i]
			J2tab = a[2,loc]
			J2ptab = a[3,loc]
			ktab = a[4,loc]
			qtab = qmax+a[5,loc]
			irtab = a[6,loc]
			if (irtab eq 1) then begin
				rho[J2tab,J2ptab,ktab,qtab] = set_real(rho[J2tab,J2ptab,ktab,qtab], a[7,loc])
			endif else begin
				rho[J2tab,J2ptab,ktab,qtab] = set_imaginary(rho[J2tab,J2ptab,ktab,qtab], a[7,loc])
			endelse
			if (J2tab eq J2ptab) then begin
				if (qtab eq 0) then begin
					rho[J2tab,J2ptab,ktab,qtab] = set_imaginary(rho[J2tab,J2ptab,ktab,qtab], 0.d0)
				endif else begin
					sign = 1.d0
					if (a[5,loc] mod 2 ne 0) then sign = -1.d0
					if (irtab eq 1) then begin
						rho[J2tab,J2ptab,ktab,qmax-a[5,loc]] = $
							set_real( rho[J2tab,J2ptab,ktab,qmax-a[5,loc]], sign*a[7,loc])
					endif else begin
						rho[J2tab,J2ptab,ktab,qmax-a[5,loc]] = $
							set_imaginary(rho[J2tab,J2ptab,ktab,qmax-a[5,loc]], -sign*a[7,loc])
					endelse
				endelse
			endif else begin
				sign = 1.d0
				if ((J2tab-J2ptab-2*a[5,loc]) mod 4 ne 0) then sign = -1.d0
				if (irtab eq 1) then begin
					rho[J2ptab,J2tab,ktab,qmax-a[5,loc]] = $
						set_real(rho[J2ptab,J2tab,ktab,qmax-a[5,loc]], sign*a[7,loc])
				endif else begin
					rho[J2ptab,J2tab,ktab,qmax-a[5,loc]] = $
						set_imaginary(rho[J2ptab,J2tab,ktab,qmax-a[5,loc]], -sign*a[7,loc])
				endelse
			endelse
		endfor
		
; Write rho^K_Q and also write those transformed to the magnetic field
; reference system		
		for j2 = J2min, J2max, 2 do begin
			for jp2 = J2min, J2max, 2 do begin
				kmin = abs(jp2-j2)/2
				kmax = (jp2+j2) / 2
				for k = kmin, kmax do begin
					for q = -k, k do begin
						suma = 0.d0
						for qp = -k, k do begin
							suma = suma + rho[j2,jp2,k,qp+qmax] * $
								rot_matrix(k, q, qp, 0.d0, -thb, -chb)
						endfor
						
						printf,2,FORMAT='(I4,2X,I2,2X,4(1X,I3),2(1x,e15.7))',n,term+1,j2, jp2, k, q,$
							float(rho[j2,jp2,k,q+qmax]), imaginary(rho[j2,jp2,k,q+qmax])
						printf,3,FORMAT='(I4,2X,I2,2X,4(1X,I3),2(1x,e15.7))',n,term+1,j2, jp2, k, q,$
							float(suma), imaginary(suma)
						n = n + 1
					endfor
				endfor
			endfor
		endfor 
		
	endfor
	
	close,2
	close,3	
end

;-----------------------------------------
; Return the wavelength of the multiplet
;-----------------------------------------
function wavelength_atom_multiplet, file_with_atom, multiplet	
	openr,2,file_with_atom
	readf,2,I2
	readf,2,nlev
	str = ''
	for i = 0, nlev-1 do begin
		readf,2,ind,J2
		fmax = (J2+I2)
		fmin = abs(J2-I2)
		nf = fix((fmax-fmin)/2) + 1		
		for j = 0, nf-1 do readf,2,str		
	endfor
	readf,2,ntran
	for i = 0, multiplet-2 do readf,2,str
	lambda = 0.d0
	readf,2,ind,up,low,aul,lambda,fac1,fac2,fac3
	close,2
	return,lambda
end

;-----------------------------------------
; Change the value of the factors nbar, omega and J10/J00 in a given file
;-----------------------------------------
pro change_factors, file_with_atom, fac_nbar, fac_omega, j10
	line = ''
	nlines = 0
	openr,2,file_with_atom
	while (not eof(2)) do begin
		readf,2,line
		nlines = nlines + 1
	endwhile
	close,2
	
	lines = strarr(nlines)
	openr,2,file_with_atom
	tmp = ''
	for i = 0, nlines-1 do begin	 	  
		readf,2,tmp
		lines[i] = tmp
	endfor
	close,2
	
	k = 1
	openr,2,file_with_atom
	readf,2,I2
	readf,2,nlev
	str = ''
	for i = 0, nlev-1 do begin
		readf,2,ind,J2
		k = k + 1
		fmax = (J2+I2)
		fmin = abs(J2-I2)
		nf = fix((fmax-fmin)/2) + 1		
		for j = 0, nf-1 do begin
			readf,2,str
			k = k + 1
		endfor
	endfor
	
	readf,2,ntran
	k = k + 1
		
	for i = 0, ntran-1 do begin
		line = strsplit(lines[k+i+1],/extract)
		line[5] = strtrim(string(fac_nbar),2)
		line[6] = strtrim(string(fac_omega),2)
		line[7] = strtrim(string(j10),2)		
		lines[k+i+1] = strjoin(line,'    ')
	endfor
	close,2
	
	openw,2,file_with_atom
	for i = 0, nlines-1 do begin
		printf,2,lines[i]
	endfor
	close,2		
end

;-----------------------------------------
; Plot the profiles
;-----------------------------------------
pro plot_profiles, stokes, state

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
	 	 	 	 
	 wavelength = wavelength_atom_multiplet(file_with_atom,state.Multiplet)
	 
	 plot_observation = 0
	 if (state.obs_file ne '') then begin
	 	  openr,2,state.obs_file
		  readf,2,nl
		  obs = fltarr(5,nl)
		  readf,2,obs
		  close,2
		  change_symbol, 'circle', 0.3
		  plot_observation = 1
		  cwpal
	 endif
	 
	 angstrom = '!3' + STRING(197B) + '!X'	 
	 !p.multi = [0,2,2]
	 
	 xmin = min(stokes(0,*))
	 xmax = max(stokes(0,*))
	 
	 wl = stokes(0,*)
	 
	 plot,wl,stokes(1,*),xtit='Wavelength shift ['+angstrom+']',ytit='I/I!dmax!n',$
	 	thick=2,xsty=1,tit='!7k!6='+strtrim(string(wavelength),2)
	 if (plot_observation eq 1) then begin
	 	  oplot,obs(0,*),obs(1,*),psym=8,col=2
	 endif
	 
	 abso = 1.d0
	 if (state.normaliz eq 0) then begin
	 	  abso = max(stokes(1,*)) - min(stokes(1,*))
	 endif
	 
	 maxq = max(stokes(2,*)/abso)
	 minq = min(stokes(2,*)/abso)
	 maxu = max(stokes(3,*)/abso)
	 minu = min(stokes(3,*)/abso)
	 max = max([maxq,maxu])
	 min = min([minq,minu])
	 	 	 
    plot,wl,stokes(2,*)/abso,xtit='Wavelength shift ['+angstrom+']',ytit='Q/I!dmax!n',thick=2, yran=[min,max],xsty=1
	 if (plot_observation eq 1) then begin
	 	  oplot,obs(0,*),obs(2,*),psym=8,col=2
	 endif
    
	 plot,wl,stokes(3,*)/abso,xtit='Wavelength shift ['+angstrom+']',ytit='U/I!dmax!n',thick=2, yran=[min,max],xsty=1
	 if (plot_observation eq 1) then begin
	 	  oplot,obs(0,*),obs(3,*),psym=8,col=2
	 endif
    
	 max = max(stokes(4,*)/abso)
	 min = min(stokes(4,*)/abso)
	 plot,wl,stokes(4,*)/abso,xtit='Wavelength shift ['+angstrom+']',ytit='V/I!dmax!n',thick=2, yran=[min,max],xsty=1
	 if (plot_observation eq 1) then begin
	 	  oplot,obs(0,*),obs(4,*),psym=8,col=2
	 endif
    
	 
	 !p.multi=0
end

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

;-----------------------------------------
; Call the F90 code which calculates the emerging profiles
;-----------------------------------------
pro synthesize, state, handler, plot_profiles=plot_profiles, texto=texto
	 siz = 57
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
	 file(20) = Bfield_str + '  ' + thetaBfield_str + '  ' + chiBfield_str
	 
	 thetaObs_str = strtrim(string(state.thetaObs),2)
	 chiObs_str = strtrim(string(state.chiObs),2)
	 gammaObs_str = strtrim(string(state.gammaObs),2)
	 file(41) = thetaObs_str + '  ' + chiObs_str + '  ' + gammaObs_str
	 
	 file(35) = strtrim(string(state.Multiplet,FORMAT='(I1)'),2)
	 
	 file(17) = strtrim(string(state.paschen),2)
	 
	 file(23) = strtrim(string(state.height),2)
	 
	 file(38) = strtrim(string(1-state.effects),2)	 
	 file(26) = strtrim(string(state.dtau_desired),2)
	 file(29) = strtrim(string(state.beta),2)
	 file(32) = strtrim(string(state.stokes0[0]),2)+ ' '+strtrim(string(state.stokes0[1]),2)+ ' '+$
	 	  strtrim(string(state.stokes0[2]),2)+ ' '+strtrim(string(state.stokes0[3]),2)

	 file(11) = '0'
	 file(14) = strtrim(string(state.D2),2)
	 	 		  
	 if (state.D2 ne 0) then begin
	 	  file(11) = '1'
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
; 	 file(44) = strtrim(string(state.waveaxis[0]),2)+ ' '+strtrim(string(state.waveaxis[1]),2)+ ' '+$
; 	 	  strtrim(string(state.waveaxis[2],FORMAT='(I5)'),2)
	 file(44) = strtrim(string(wright),2)+ ' '+strtrim(string(wleft),2)+ ' '+$
	 	  strtrim(string(state.waveaxis[2],FORMAT='(I5)'),2)

; ; He I
; 	 if (state.which_atom eq 0) then begin
; 	 	  case state.Multiplet of
; 	 	   	1 : wl = 10829.0911d0
; 		   	2 : wl = 3888.6046d0
; 		   	3 : wl = 7065.7085d0
; 		   	4 : wl = 5875.9663d0
; 	 	  endcase
; 	 endif
; 	 
; ; S I
; 	 if (state.which_atom eq 1) then begin
; 	 	  case state.Multiplet of
; 	 	   	1 : wl = 9237.538d0
; 	 	  endcase
; 	 endif
; 	 
; ; Na I	 
; 	if (state.which_atom eq 2) then begin
; 	 	  case state.Multiplet of
; 	 	   	1 : wl = 5889.95d0
; 	 	  endcase
; 	 endif
	 
	 PC = 2.99792458d10
	 PH = 6.62606876d-27
	 
	 nu = PC / (wl*1.d-8)
	 
	 wl_str = strtrim(string(wl),2)
	 vdopp_str = strtrim(string(state.Doppler),2)
	 damp_str = strtrim(string(state.damping),2)
	 file(47) = wl_str + '  ' + vdopp_str + '  '+ damp_str
	 file(50) = '0.d0'
	 file(56) = strtrim(string(state.stimulated),2)
	 file(53) = strtrim(string(state.magneto_opt),2)
	 
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
	 
; 	 file = strarr(44)	 
; 	 if (state.which_atom eq 0 and state.which_code eq 0) then begin
; 		  openr,2,'ATOMS/helium_original.mod'
; 		  tmp = ''
; 		  for i = 0, 43 do begin	 	  
; 	 			readf,2,tmp
; 				file(i) = tmp
; 		  endfor
; 		  close,2
; 		  line_10830 = strsplit(file(19),/extract)
; 		  line_10830(5) = strtrim(string(state.factor_10830_nbar),2)
; 		  line_10830(6) = strtrim(string(state.factor_10830_omega),2)
; 		  line_10830(7) = strtrim(string(state.j10),2)
; 		  file(19) = strjoin(line_10830,'    ')
; 		  openw,2,'ATOMS/helium.mod'
; 		  for i = 0, 43 do begin
; 	 			printf,2,file(i)
; 		  endfor
; 		  close,2		  
; 	 endif	 
	 
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


function return_i0_allen, state

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
	 
	 mu = cos(state.thetaObs * !DPI / 180.d0)
	 ic = ddread('CLV/ic.dat',/noverb)
	 cl = ddread('CLV/cl.dat',/noverb)
	 PC = 2.99792458d10
	 PH = 6.62606876d-27
	 
; Wavelength in A	 
	 ic(0,*) = 1.d4 * ic(0,*)	 
; I_lambda to I_nu	 
	 ic(1,*) = 1.d14 * ic(1,*) * (ic(0,*)*1.d-8)^2 / PC
	 
	 cl(0,*) = 1.d4 * cl(0,*)
	 
	 u = interpol(cl(1,*),cl(0,*),wl)
	 v = interpol(cl(2,*),cl(0,*),wl)
	 
	 imu = 1.d0 - u - v + u * mu + v * mu^2
	 i0 = interpol(ic(1,*),ic(0,*),wl)
	 	 
	 return, i0*imu
end

;-----------------------------------------
; Initialization routine. Draws the widget
;-----------------------------------------
function hazel_init

	 if (file_test('state.idl') eq 1) then begin
	 	  restore,'state.idl'
	 endif else begin
	 	  state = {baseWidget: 0L, plotWidget: 0L, Bfield: 40.d0, Bfieldmax: 500.d0, thetaBfield: 30.d0, chiBfield: 90.d0, $
	 	   	thetaObs: 90.d0, chiObs: 0.d0, gammaObs: 90.d0, Multiplet: 1L, Doppler: 8.d0, auto: 0L, BSlider: 0L, $
		   	paschen: 1L, effects: 0L, observation: 0L, stokes0: [1.d0,0.d0,0.d0,0.d0], dtau_desired: 3.d0, height: 20.d0, $
		   	postcript: 0L, dtaured: 3.d0, beta: 2.d0, damping: 0.d0, stimulated: 0, magneto_opt: 0, obs_file: '', $
				i0_allen: 0L, D2: 0.d0, factor_10830_nbar: 1.d0, factor_10830_omega: 1.d0, waveaxis: [-3.d0,2.5d0,200.d0],$
				which_atom: 0L, MultipletSlider: 0L, normaliz: 1L, j10: 0.d0, which_code: 0,$
				bfield_var: [1.d-4,1.d3,15], which_rho_plot: 0L, which_refframe: 0L,$
				path_save: '.', randomazimuth: 0L}
	 endelse
	 
	 if (state.which_code eq 0) then begin
	 	if (state.which_atom eq 0) then begin
	 	  	state.baseWidget = widget_base(TITLE='Hanle simulator : He I', MBAR=menuBar)	 	  	
	 	endif
	 
	 	if (state.which_atom eq 1) then begin
	 	  	state.baseWidget = widget_base(TITLE='Hanle simulator : S I', MBAR=menuBar)	 	  	
	 	endif
	 
	 	if (state.which_atom eq 2) then begin
	 	  	state.baseWidget = widget_base(TITLE='Hanle simulator : Na I', MBAR=menuBar)	 	  	
	 	endif
	 endif
	 
	 if (state.which_code eq 1) then begin
	 	if (state.which_atom eq 0) then begin
	 	  	state.baseWidget = widget_base(TITLE='Hanle simulator : Na I HFS', MBAR=menuBar)	 	  	
	 	endif
	 endif
	 
	 hanleBase = widget_base(state.baseWidget, /COLUMN, /BASE_ALIGN_CENTER)
	 
	 atomMenu = widget_button(menuBar, VALUE='Multiterm', /MENU)
	 heliumButton = widget_button(atomMenu, VALUE='He I', UVALUE='HELIUM')
	 sulfurButton = widget_button(atomMenu, VALUE='S I', UVALUE='SULFUR')
	 sodiumButton = widget_button(atomMenu, VALUE='Na I', UVALUE='SODIUM')
	 
; 	 atomMenu = widget_button(menuBar, VALUE='Multilevel HFS', /MENU)	 
; 	 sodiumButton = widget_button(atomMenu, VALUE='Na I', UVALUE='SODIUM_HFS')
	 
	 horizBase = widget_base(hanleBase, /ROW)
	 
	 state.plotWidget = widget_draw(horizBase, XSIZE=650, YSIZE=450, /FRAME)
	 	 
	 rightbase = widget_base(horizBase, /COLUMN, FRAME=1)
	 collislb = widget_label(rightBase, VALUE='Lower level D^(2):')
	 collis = widget_text(rightBase, VALUE=strtrim(string(state.D2),2),UVALUE='D2',/EDITABLE,XSIZE=8,YSIZE=1)
	 factor_10830lb = widget_label(rightBase, VALUE='nbar factor')
	 factor_10830 = widget_text(rightBase, VALUE=strtrim(string(state.factor_10830_nbar),2),UVALUE='fact_10830_nbar',/EDITABLE,XSIZE=8,YSIZE=1)
	 factor_10830_wlb = widget_label(rightBase, VALUE='w factor')
	 factor_10830_w = widget_text(rightBase, VALUE=strtrim(string(state.factor_10830_omega),2),UVALUE='fact_10830_w',/EDITABLE,XSIZE=8,YSIZE=1)
	 j10_wlb = widget_label(rightBase, VALUE='J10/J00')
	 j10 = widget_text(rightBase, VALUE=strtrim(string(state.j10),2),UVALUE='j10_tensor',/EDITABLE,XSIZE=8,YSIZE=1)
	 
; Emission or tangent observation
	 t4 = widget_base(rightBase, /COLUMN)
	 lab4 = widget_label(t4, VALUE='Observation')	 
	 obsBase = widget_base(t4, /COLUMN, /EXCLUSIVE)
	 emissionButton = widget_button(obsBase, VALUE='Slab (optically thin)', UVALUE='EMISSION')
	 simpleslabButton = widget_button(obsBase, VALUE='Simplified slab', UVALUE='SIMPLE_SLAB')
	 formalButton = widget_button(obsBase, VALUE='Slab (no MO)', UVALUE='FORMAL')
	 deloparButton = widget_button(obsBase, VALUE='Slab (DELO)', UVALUE='DELOPAR')
	 exactslabButton = widget_button(obsBase, VALUE='Slab (exact)', UVALUE='EXACT_SLAB')
	 milneButton = widget_button(obsBase, VALUE='Milne-Eddington', UVALUE='MILNE')
	 
	 case state.observation of
	 	0: begin
	 			widget_control, emissionButton, /SET_BUTTON   ; 0 -> pure emission
	 		end
	 	1: begin
	 			widget_control, formalButton, /SET_BUTTON   ; 1 -> slab (neglecting MO terms)
	 		end
	 	2: begin
	 			widget_control, milneButton, /SET_BUTTON   ; 2 -> Milne-Eddington
	 		end
	 	3: begin
	 			widget_control, deloparButton, /SET_BUTTON   ; 3 -> DELO
	 		end
	 	4: begin
				widget_control, simpleslabButton, /SET_BUTTON   ; 4 -> simplified slab (optically thin)
			end
	 	5: begin
	 			widget_control, exactslabButton, /SET_BUTTON   ; 5 -> exact slab
	 		end
	 endcase
	 
; Wavelength axis
	 waveaxis = widget_base(t4, /COLUMN)
	 wleftbase = widget_base(waveaxis, /ROW)
	 wleftlab = widget_label(wleftbase, VALUE='wl:')
	 wleft = widget_text(wleftbase, VALUE=strtrim(string(state.waveaxis(0)),2),$
	 	UVALUE='wleft',/EDITABLE,XSIZE=5,YSIZE=1)
	 ;wrightbase = widget_base(waveaxis, /ROW)
	 wrightlab = widget_label(wleftbase, VALUE='wr:')
	 wright = widget_text(wleftbase, VALUE=strtrim(string(state.waveaxis(1)),2),$
	 	UVALUE='wright',/EDITABLE,XSIZE=5,YSIZE=1)
	 wstepbase = widget_base(waveaxis, /ROW)
	 wsteplab = widget_label(wstepbase, VALUE='step:')
	 wstep = widget_text(wstepbase, VALUE=strtrim(string(state.waveaxis(2)),2),UVALUE='wstep',/EDITABLE,XSIZE=5,YSIZE=1)
	 
	 slidersbase = widget_base(horizBase, /COLUMN)
	 ;slidersBase = widget_base(hanleBase, /ROW)
	 
; Magnetic field information
	 fieldBase = widget_base(slidersBase, /COLUMN, FRAME=1)
	 
; Magnetic field strength	 
	 BSliderBase = widget_base(fieldBase, /COLUMN)
	 maxBBase = widget_base(BSliderBase, /ROW)
	 BSlidertext = widget_label(maxBBase, VALUE='B_max:')
	 BSliderMax = widget_text(maxBBase,VALUE=strtrim(string(state.Bfieldmax),2),$
	 	UVALUE='BSliderMax',/EDITABLE,XSIZE=8,YSIZE=1)	 
	 state.BSlider = cw_fslider(BSliderBase,TITLE='Magnetic field strength [G]',UVALUE='BSlider',$
	 	  XSIZE=255,MINIMUM=0.d0,MAXIMUM=state.Bfieldmax,VALUE=state.Bfield)	 

; Magnetic field inclination
	 thetaBSliderBase = widget_base(fieldBase, /ROW)
	 thetaBSlider = widget_slider(thetaBSliderBase,TITLE='Magnetic field inclination [deg]',UVALUE='thetaBSlider',$
	 	  XSIZE=255,MAXIMUM=180.d0,VALUE=state.thetaBfield)
	 
; Magnetic field azimuth
	 chiBSliderBase = widget_base(fieldBase, /ROW)
	 chiBSlider = widget_slider(chiBSliderBase,TITLE='Magnetic field azimuth [deg]',UVALUE='chiBSlider',$
	 	  XSIZE=255,MINIMUM=-180.d0,MAXIMUM=180.d0,VALUE=state.chiBfield)
	 chiBrandomBase = widget_base(chiBSliderBase, /COLUMN, /EXCLUSIVE)
	 chiBrandButton = widget_button(chiBrandomBase, VALUE='Random azi', UVALUE='RANDOMAZI_ON')
	 chiBnorandButton = widget_button(chiBrandomBase, VALUE='Normal azi', UVALUE='RANDOMAZI_OFF')
	 widget_control, (state.randomazimuth) ? chiBrandButton : chiBnorandButton, /SET_BUTTON
		  

	 observationBase = widget_base(slidersBase, /COLUMN, FRAME=1)
	 thetaOSlider = widget_slider(observationBase,TITLE='Observing theta angle [deg]',UVALUE='thetaOSlider',$
	 	  XSIZE=255,MAXIMUM=180.d0,VALUE=state.thetaObs)

	 chiOSlider = widget_slider(observationBase,TITLE='Observing chi angle [deg]',UVALUE='chiOSlider',$
	 	  XSIZE=255,MAXIMUM=360.d0,VALUE=state.chiObs)
		  
	 gammaOSlider = widget_slider(observationBase,TITLE='Observing gamma angle [deg]',UVALUE='gammaOSlider',$
	 	  XSIZE=255,MAXIMUM=180.d0,VALUE=state.gammaObs)
		  	 		  
; Synthesize
	 synth_globBase = widget_base(hanleBase, /ROW, /ALIGN_LEFT)
	 	 
	 synthBase = widget_base(synth_globBase, /COLUMN)
	 multiBase = widget_base(synthBase, /ROW, FRAME=1)
	 if (state.which_code eq 0) then begin
	 	if (state.which_atom eq 0) then begin
	 	  	state.MultipletSlider =$
	 	  		widget_slider(multiBase,TITLE='Multiplet',UVALUE='MultipletSlider',$
	 	   	XSIZE=120,MAXIMUM=4,MINIMUM=1,VALUE=state.multiplet)	 
	 	endif
	 	if (state.which_atom eq 1) then begin
	 	  	state.MultipletSlider =$
	 	  		widget_slider(multiBase,TITLE='Multiplet',UVALUE='MultipletSlider',$
	 	   	XSIZE=120,MAXIMUM=4,MINIMUM=1,VALUE=state.multiplet,sensitive=0)
	 	endif
	 	if (state.which_atom eq 2) then begin
	 	  	state.MultipletSlider =$
	 	  		widget_slider(multiBase,TITLE='Multiplet',UVALUE='MultipletSlider',$
	 	   	XSIZE=120,MAXIMUM=4,MINIMUM=1,VALUE=state.multiplet,sensitive=0)
	 	endif
	 endif
	 if (state.which_code eq 1) then begin
	 	if (state.which_atom eq 0) then begin
	 	  	state.MultipletSlider =$
	 	  		widget_slider(multiBase,TITLE='Multiplet',UVALUE='MultipletSlider',$
	 	   	XSIZE=120,MAXIMUM=2,MINIMUM=1,VALUE=state.multiplet,sensitive=1)
	 	endif
	 endif
	 
	 DopplerSlider = widget_slider(multiBase,TITLE='Doppler velocity [km/s]',UVALUE='DopplerSlider',$
	 	  XSIZE=200,MAXIMUM=25,MINIMUM=0.1,VALUE=state.Doppler)
	 heightSlider = widget_slider(multiBase,TITLE='Height (<0 if apparent) ["]',UVALUE='heightSlider',$
	 	  XSIZE=255,MINIMUM=-100.d0,MAXIMUM=900.d0,VALUE=state.height)
	 
	 
	 i0 = return_i0_allen(state)
	 state.i0_allen = widget_label(multiBase, VALUE='Allen: '+strtrim(string(i0),2))
	 	 	 	 
		  	 	 
	 formalBase = widget_base(synthBase, /ROW, FRAME=1)
	 I0lb = widget_label(formalBase, VALUE='I0:')
	 I0 = widget_text(formalBase, VALUE=strtrim(string(state.stokes0(0),FORMAT='(E10.4)'),2),UVALUE='I0',/EDITABLE,XSIZE=10,YSIZE=1)
	 Q0lb = widget_label(formalBase, VALUE='Q0:')
	 Q0 = widget_text(formalBase, VALUE=strtrim(string(state.stokes0(1),FORMAT='(E10.4)'),2),UVALUE='Q0',/EDITABLE,XSIZE=10,YSIZE=1)
	 U0lb = widget_label(formalBase, VALUE='U0:')
	 U0 = widget_text(formalBase, VALUE=strtrim(string(state.stokes0(2),FORMAT='(E10.4)'),2),UVALUE='U0',/EDITABLE,XSIZE=10,YSIZE=1)
	 V0lb = widget_label(formalBase, VALUE='V0:')
	 V0 = widget_text(formalBase, VALUE=strtrim(string(state.stokes0(3),FORMAT='(E10.4)'),2),UVALUE='V0',/EDITABLE,XSIZE=10,YSIZE=1)
	 dtauredlb = widget_label(formalBase, VALUE='Dtau(max):')
	 dtaured_wid = widget_text(formalBase, VALUE=strtrim(string(state.dtau_desired,FORMAT='(F6.3)'),2),UVALUE='DTAU',/EDITABLE,XSIZE=8,YSIZE=1)
	 dbetalb = widget_label(formalBase, VALUE='beta:')
	 beta_wid = widget_text(formalBase, VALUE=strtrim(string(state.beta,FORMAT='(F6.3)'),2),UVALUE='BETA',/EDITABLE,XSIZE=8,YSIZE=1)
	 ddamplb = widget_label(formalBase, VALUE='a:')
	 damp_wid = widget_text(formalBase, VALUE=strtrim(string(state.damping,FORMAT='(F6.3)'),2),UVALUE='DAMPING',/EDITABLE,XSIZE=8,YSIZE=1)
	 	 	 	 
	 buttonBase = widget_base(hanleBase, /ROW, /ALIGN_LEFT)

; Postcript output
	 state.postcript = 0
	 t0 = widget_base(buttonBase, /COLUMN)
	 lab0 = widget_label(t0, VALUE='Output')
	 postcriptBase = widget_base(t0, /COLUMN, /EXCLUSIVE)
	 postcriptnoButton = widget_button(postcriptBase, VALUE='PS', UVALUE='PS_ON')
	 postcriptyesButton = widget_button(postcriptBase, VALUE='Graphical', UVALUE='PS_OFF')
	 postcriptwriteButton = widget_button(postcriptBase, VALUE='Save', UVALUE='PS_WRITE')	 
	 widget_control, postcriptyesButton, /SET_BUTTON
	 
; Automatic synthesis	 
	 t1 = widget_base(buttonBase, /COLUMN)
	 lab1 = widget_label(t1, VALUE='Automatic')
	 autoBase = widget_base(t1, /COLUMN, /EXCLUSIVE)
	 autoyesButton = widget_button(autoBase, VALUE='Automatic', UVALUE='AUTO_ON')
	 autonoButton = widget_button(autoBase, VALUE='Manual', UVALUE='AUTO_OFF')
	 widget_control, (state.auto) ? autoyesButton : autonoButton, /SET_BUTTON

; Paschen-Back effect
	 t2 = widget_base(buttonBase, /COLUMN)
	 lab2 = widget_label(t2, VALUE='ZEEMAN')
	 paschenBase = widget_base(t2, /COLUMN, /EXCLUSIVE)
	 paschenyesButton = widget_button(paschenBase, VALUE='Paschen-Back', UVALUE='PASCHEN')
	 paschennoButton = widget_button(paschenBase, VALUE='Linear Zeeman', UVALUE='LINEAR')
	 widget_control, (state.paschen) ? paschenyesButton : paschennoButton, /SET_BUTTON
	 
; Only Zeeman effect or all
	 t3 = widget_base(buttonBase, /COLUMN)
	 lab3 = widget_label(t3, VALUE='Atom. pol.')
	 effectBase = widget_base(t3, /COLUMN, /EXCLUSIVE)
	 scat_zeemButton = widget_button(effectBase, VALUE='Yes', UVALUE='ALL')
	 zeemButton = widget_button(effectBase, VALUE='No', UVALUE='ZEEMAN')
	 widget_control, (state.effects) ? zeemButton : scat_zeemButton, /SET_BUTTON
	 
; Include stimulated emission in RT
	 t5 = widget_base(buttonBase, /COLUMN)
	 lab3 = widget_label(t5, VALUE='Stim. emis. in RT')
	 stimemiBase = widget_base(t5, /COLUMN, /EXCLUSIVE)
	 yesButton = widget_button(stimemiBase, VALUE='Yes', UVALUE='STIM')
	 NOButton = widget_button(stimemiBase, VALUE='No', UVALUE='NOSTIM')
	 widget_control, (state.stimulated) ? yesButton : noButton, /SET_BUTTON
	 
; Include magneto-optical effects
	 t6 = widget_base(buttonBase, /COLUMN)
	 lab3 = widget_label(t6, VALUE='Magneto-optical')
	 magnetoBase = widget_base(t6, /COLUMN, /EXCLUSIVE)
	 yesButton = widget_button(magnetoBase, VALUE='Yes', UVALUE='MAGNETOOPT')
	 noButton = widget_button(magnetoBase, VALUE='No', UVALUE='NOMAGNETOOPT')
	 widget_control, (state.magneto_opt) ? yesButton : noButton, /SET_BUTTON
	 
; Normalization
	 t7 = widget_base(buttonBase, /COLUMN)
	 lab3 = widget_label(t7, VALUE='Normalization')
	 absorBase = widget_base(t7, /COLUMN, /EXCLUSIVE)
	 maxButton = widget_button(absorBase, VALUE='Maximum', UVALUE='MAXIMUM')
	 absButton = widget_button(absorBase, VALUE='Absorption', UVALUE='ABSORPTION')
	 widget_control, maxButton, /SET_BUTTON
	 	 	 	 
; Rhos plot
	 t7 = widget_base(buttonBase, /COLUMN)
	 lab3 = widget_label(t7, VALUE='Plot rho^K_Q')
	 absorBase = widget_base(t7, /COLUMN, /EXCLUSIVE)
	 cohButton = widget_button(absorBase, VALUE='Non-diagonal', UVALUE='NONDIAGONAL')
	 diagButton = widget_button(absorBase, VALUE='Diagonal', UVALUE='DIAGONAL')
	 widget_control, (state.which_rho_plot) ? diagButton : cohButton, /SET_BUTTON
	 
; Reference frame
	 t7 = widget_base(buttonBase, /COLUMN)
	 lab3 = widget_label(t7, VALUE='Ref. frame rho^K_Q')
	 absorBase = widget_base(t7, /COLUMN, /EXCLUSIVE)
	 vertButton = widget_button(absorBase, VALUE='Vertical', UVALUE='VERTICAL_REFFRAME')
	 magnButton = widget_button(absorBase, VALUE='Magnetic', UVALUE='MAGNETIC_REFFRAME')
	 widget_control, (state.which_refframe) ? magnButton : vertButton, /SET_BUTTON
	 
	 buttonsBase = widget_base(synth_globBase, /COLUMN)
	 observ_include = widget_button(buttonsBase,VALUE='Load Observation',UVALUE='LOAD_OBSERVATION')
	 reset_observ = widget_button(buttonsBase,VALUE='Reset Observation',UVALUE='RESET_OBSERVATION')
	 
	 bfieldvarBase = widget_base(buttonsBase,/ROW)
	 bfield_var = widget_button(bfieldvarBase,$
	 	VALUE='rho(B)',UVALUE='FIELD_VARIATION')
	 dnfields = widget_label(bfieldvarBase, VALUE='N:')
	 nfields = widget_text(bfieldvarBase,$
	 	VALUE=strtrim(string(state.bfield_var[2],FORMAT='(I3)'),2),$
	 	UVALUE='NFIELDS_RHO',/EDITABLE,XSIZE=3,YSIZE=1)
	 dBmin = widget_label(bfieldvarBase, VALUE='Bmin:')
	 bmin = widget_text(bfieldvarBase,$
	 	VALUE=strtrim(string(state.bfield_var[0]),2),UVALUE='BMIN_RHO',$
	 	/EDITABLE,XSIZE=5,YSIZE=1)
	 dBmax = widget_label(bfieldvarBase, VALUE='Bmax:')
	 bmax = widget_text(bfieldvarBase,$
	 	VALUE=strtrim(string(state.bfield_var[1]),2),UVALUE='BMAX_RHO',$
	 	/EDITABLE,XSIZE=5,YSIZE=1)
	 	
	 synthButton = widget_button(buttonsBase,VALUE='Calculate',UVALUE='Calculate')
	 
	 widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
	 
	 return, state
end

;-----------------------------------------
; Event handler
;-----------------------------------------
pro hazel_Event, event
	 widget_control, Event.id, GET_UVALUE=Action
	 
	 hand = widget_info(Event.Handler, /CHILD)
	 widget_control, hand, GET_UVALUE=state
	 
	 case Action of
	 	  'Calculate' : 	begin
		   		 	   		 synthesize, state, hand, /plot_profiles		   		 	   		 
		   		 	   	end
		  'FIELD_VARIATION' : 	begin
		   		 	   		 calculate_field_variation, state
		   		 	   	end
	 	  'LOAD_OBSERVATION' : 	begin
		   		 	   		 state.obs_file = dialog_pickfile(default_extension='.prof')
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end
		  'RESET_OBSERVATION' : 	begin
		   		 	   		 state.obs_file = ''
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end						  
		  'BSlider' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.Bfield = value
									 widget_control, hand, SET_UVALUE=state
									 if (state.auto eq 1) then synthesize, state, hand, /plot_profiles
		   		 	   	end
	 	  'BSliderMax' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.Bfieldmax = value
									 widget_control, hand, SET_UVALUE=state
									 widget_control, state.Bslider, SET_VALUE=[state.Bfield, 0, state.Bfieldmax]
									 widget_control, state.Bslider, GET_VALUE=value
									 state.Bfield = value
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end
	 	  'thetaBSlider' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.thetaBfield = value
									 widget_control, hand, SET_UVALUE=state
									 if (state.auto eq 1) then synthesize, state, hand, /plot_profiles
		   		 	   	end
	 	  'chiBSlider' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.chiBfield = value
									 widget_control, hand, SET_UVALUE=state
									 if (state.auto eq 1) then synthesize, state, hand, /plot_profiles
		   		 	   	end
	 	  'thetaOSlider' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.thetaObs = value									 
									 widget_control, hand, SET_UVALUE=state
									 I0 = return_i0_allen(state)
									 widget_control, state.i0_allen, SET_VALUE='Allen: '+strtrim(string(i0),2)
									 if (state.auto eq 1) then synthesize, state, hand, /plot_profiles
		   		 	   	end
	 	  'chiOSlider' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.chiObs = value
									 widget_control, hand, SET_UVALUE=state
									 if (state.auto eq 1) then synthesize, state, hand, /plot_profiles
		   		 	   	end
	 	  'gammaOSlider' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.gammaObs = value
									 widget_control, hand, SET_UVALUE=state
									 if (state.auto eq 1) then synthesize, state, hand, /plot_profiles
		   		 	   	end
	 	  'MultipletSlider' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.Multiplet = value
									 widget_control, hand, SET_UVALUE=state
									 I0 = return_i0_allen(state)
									 widget_control, state.i0_allen, SET_VALUE='Allen: '+strtrim(string(i0),2)
									 if (state.auto eq 1) then synthesize, state, hand, /plot_profiles
		   		 	   	end
	 	  'DopplerSlider' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.Doppler = value
									 widget_control, hand, SET_UVALUE=state
									 if (state.auto eq 1) then synthesize, state, hand, /plot_profiles
		   		 	   	end
	 	  'heightSlider' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.height = value
									 widget_control, hand, SET_UVALUE=state
									 if (state.auto eq 1) then synthesize, state, hand, /plot_profiles
		   		 	   	end
	 	  'AUTO_ON' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.auto = 1
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end 
	 	  'AUTO_OFF' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.auto = 0
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end 
	 	  'PS_ON' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.postcript = 1
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end 
	 	  'PS_OFF' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.postcript = 0
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end 
		  'PS_WRITE' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.postcript = 2
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end 
		  'RANDOMAZI_ON' :  begin
		  								widget_control, Event.id, GET_VALUE=value
		   		 	   		   state.randomazimuth = 1
									   widget_control, hand, SET_UVALUE=state
		  						end
		  'RANDOMAZI_OFF' :  begin
		  								widget_control, Event.id, GET_VALUE=value
		   		 	   		   state.randomazimuth = 0
									   widget_control, hand, SET_UVALUE=state
		  						end
	 	  'PASCHEN' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.paschen = 1
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end 
	 	  'LINEAR' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.paschen = 0
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end 
	 	  'ALL' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.effects = 0
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end 
	 	  'ZEEMAN' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.effects = 1
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end
	 	  'EMISSION' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.observation = 0
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end 
	 	  'TANGENTIAL' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.observation = -1
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end 
	 	  'FORMAL' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.observation = 1
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end
		  'DELOPAR' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.observation = 3
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end
		  'EXACT_SLAB' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.observation = 5
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end
	 	  'MILNE' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.observation = 2
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end
	 	  'SIMPLE_SLAB' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.observation = 4
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end
	 	  'I0' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.stokes0[0] = value
									 widget_control, Event.id, SET_VALUE=strtrim(string(value,FORMAT='(E10.4)'),2)
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end
		  'wleft' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.waveaxis[0] = value
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end
		  'wright' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.waveaxis[1] = value
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end
		  'wstep' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.waveaxis[2] = value
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end
	 	  'Q0' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.stokes0[1] = value
									 widget_control, Event.id, SET_VALUE=strtrim(string(value,FORMAT='(E10.4)'),2)
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end 
	 	  'U0' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.stokes0[2] = value
									 widget_control, Event.id, SET_VALUE=strtrim(string(value,FORMAT='(E10.4)'),2)
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end 
	 	  'V0' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.stokes0[3] = value
									 widget_control, Event.id, SET_VALUE=strtrim(string(value,FORMAT='(E10.4)'),2)
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end 
	 	  'DTAU' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.dtau_desired = value
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end 
	 	  'BETA' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.beta = value
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end 
	 	  'DAMPING' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.damping = value
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end 
	 	  'STIM' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.stimulated = 1
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end 
	 	  'NOSTIM' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.stimulated = 0
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end 
	 	  'MAGNETOOPT' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.magneto_opt = 1
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end
		  'MAXIMUM' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.normaliz = 1
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end 
		  'ABSORPTION' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.normaliz = 0
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end 
	 	  'NOMAGNETOOPT' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.magneto_opt = 0
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end 
		  'D2' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.D2 = value
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end
		  'fact_10830_w' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.factor_10830_omega = value
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end
		  'fact_10830_nbar' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.factor_10830_nbar = value
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end
		  'j10_tensor' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.j10 = value
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end
		  'HELIUM' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.which_atom = 0
									 state.multiplet = 1
									 state.which_code = 0
									 widget_control, hand, SET_UVALUE=state
									 widget_control, state.MultipletSlider, SET_SLIDER_MAX=4
									 widget_control, state.MultipletSlider, SET_VALUE=1
									 widget_control, state.MultipletSlider, sensitive=1
									 widget_control, state.baseWidget, $
									 	BASE_SET_TITLE='Hanle Simulator : He I'
		   		 	   	end
		  'SULFUR' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.which_atom = 1
									 state.multiplet = 1
									 state.which_code = 0
									 widget_control, hand, SET_UVALUE=state
									 widget_control, state.MultipletSlider, SET_SLIDER_MAX=1
									 widget_control, state.MultipletSlider, SET_VALUE=1
									 widget_control, state.MultipletSlider, sensitive=0
									 widget_control, state.baseWidget, $
									 	BASE_SET_TITLE='Hanle Simulator : S I'
		   		 	   	end
		   'SODIUM' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.which_atom = 2
									 state.multiplet = 1
									 state.which_code = 0
									 widget_control, hand, SET_UVALUE=state
									 widget_control, state.MultipletSlider, SET_SLIDER_MAX=1
									 widget_control, state.MultipletSlider, SET_VALUE=1
									 widget_control, state.MultipletSlider, sensitive=0
									 widget_control, state.baseWidget, $
									 	BASE_SET_TITLE='Hanle Simulator : Na I'
		   		 	   	end
		   'SODIUM_HFS' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.which_atom = 0
									 state.multiplet = 1
									 state.which_code = 1
									 widget_control, hand, SET_UVALUE=state
									 widget_control, state.MultipletSlider, SET_SLIDER_MAX=2
									 widget_control, state.MultipletSlider, SET_VALUE=1
									 widget_control, state.MultipletSlider, sensitive=1
									 widget_control, state.baseWidget, $
									 	BASE_SET_TITLE='Hanle Simulator : Na I HFS'
		   		 	   	end
		   'NFIELDS_RHO' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.bfield_var[2] = value
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end
		   'BMIN_RHO' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.bfield_var[0] = value
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end 
		  	'BMAX_RHO' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.bfield_var[1] = value
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end 
		   'NONDIAGONAL' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.which_rho_plot = 0
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end
		   'DIAGONAL' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.which_rho_plot = 1
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end
		   'VERTICAL_REFFRAME' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.which_refframe = 0
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end
		   'MAGNETIC_REFFRAME' : 	begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 state.which_refframe = 1
									 widget_control, hand, SET_UVALUE=state
		   		 	   	end    
		   'tableSigma' : begin
		   		 	   		 widget_control, Event.id, GET_VALUE=value
		   		 	   		 
		   		 	   	end
	 endcase
	 save, state, filename='state.idl'

end


;-----------------------------------------
; Main routine
;-----------------------------------------
pro hazel
	 state = hazel_init()
	 
	 widget_control, state.baseWidget, /REALIZE
	 
	 xmanager, 'HAZEL', state.baseWidget, EVENT_HANDLER='hazel_Event'
end

;-----------------------------------------
; Add a new property to the state structure
;-----------------------------------------
pro add_to_state
	restore,'state.idl'
	temp = create_struct(state, 'randomazimuth', 0L)
	state = temp
	save, state, filename='state.idl'
end
