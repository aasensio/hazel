@angle_transformation
@polarity

; Filter the vector of solutions to avoid all complex solutions and
; all those that are very similar. Additionally, use only those whose
; value is >=0, because these refer to t=sin(thetaSky), which is always
; positive, given that thetaSky is in the range [0,Pi]
function filterSolutions, v, nsolutions
	ind = where(abs(imaginary(v)) lt 1.d-5 and real_part(v) ge 0.d0 and real_part(v) le 1.d0, nsolutions)
	if (nsolutions ne 0) then begin		
		roots = real_part(v[ind])

; Avoid repeated values
		roots_nonduplicate = roots[0]
		for i = 1, nsolutions-1 do begin
			ind = where(abs(roots[i]-roots_nonduplicate) lt 1.d-3, count)
			if (count eq 0) then begin
				roots_nonduplicate = [roots_nonduplicate, roots[i]]
			endif
		endfor
		nsolutions = n_elements(roots_nonduplicate)
	endif else begin
		roots_nonduplicate = 0.d0
		nsolutions = 0
	endelse
	return, roots_nonduplicate
end

; Given original thetaSky and chiSky, test which of the two possible solutions of sin(thetaSky_new)=t is the correct one
; by verifying that the equations for the saturation regime are fulfilled
; This is not enough for the thetaLOS=90º
function pickThetaSkyFromSolutions, thetaLOS, thetaSky, chiSky, thetaSky_new, chiSky_new
	costhetaB = cos(thetaLOS)*cos(thetaSky) - sin(thetaLOS)*sin(thetaSky)*cos(chiSky)
	K = (3.d0*costhetaB^2-1.d0) * sin(thetaSky)^2 * cos(2.d0*chiSky)

	costhetaB_new = cos(thetaLOS)*cos(thetaSky_new) - sin(thetaLOS)*sin(thetaSky_new)*cos(chiSky_new)
	K_new = (3.d0*costhetaB_new^2-1.d0) * sin(thetaSky_new)^2* cos(2.d0*chiSky_new)

	if (abs(K-K_new) lt 1.d-3) then begin
		out = thetaSky_new
	endif else begin
		out = !DPI - thetaSky_new
	endelse
	return, out
end

; Given all the possible ambiguous cases, reorder them so that the original one appears the first
pro reorderAngles, thetaB_all, chiB_all, thetaSky_all, chiSky_all, thetaB, chiB
	
	if (n_elements(thetaB_all) gt 1) then begin		
		ind = where(abs(thetaB_all-thetaB) lt 1.d-3 and abs(chiB_all-chiB) lt 1.d-3)
		
		old = thetaB_all[0]
		thetaB_all[0] = thetaB_all[ind]
		thetaB_all[ind] = old

		old = thetaSky_all[0]
		thetaSky_all[0] = thetaSky_all[ind]
		thetaSky_all[ind] = old

		old = chiB_all[0]
		chiB_all[0] = chiB_all[ind]
		chiB_all[ind] = old

		old = chiSky_all[0]
		chiSky_all[0] = chiSky_all[ind]
		chiSky_all[ind] = old
	endif
end


; Return the possible ambiguous inclination and azimuths for
; observation at an angle thetaObs and assuming saturation limit
; If StokesVpolarity is given, use this polarity. If not, we use the sign of cos(thetaSky) inferred from the input angles
; getAmbiguitiesGeneral,90.d0,30.,0.0,1.,t1,t2
pro getAmbiguitiesGeneral, thetaObs, thetaB_in, chiB_in, yesStokesV, thetaB_out, chiB_out, thetaSky_out, chiSky_out, verbose=verbose, StokesVpolarity=StokesVpolarity

	conv = !DPI / 180.d0

	thetaB = thetaB_in * conv
	chiB = chiB_in * conv
	thetaLOS = thetaObs * conv

	nsol = 0

; Maximum number of solutions
	thetaB_out = dblarr(8)
	chiB_out = dblarr(8)

	thetaSky_out = dblarr(8)
	chiSky_out = dblarr(8)
	
; Transform angles with respect to the vertical to angles with respect to the LOS
	temp = absolute_to_sky_general(thetaLOS,thetaB, chiB)

	thetaSky = temp[0]
	chiSky = temp[1]
	
	if (keyword_set(verbose)) then begin
		print, '(thetaSky,chiSky) = ', thetaSky / conv, chiSky / conv
		print, '(thetaB,chiB) = ', thetaB / conv, chiB / conv

		ct = cos(thetaLOS)*cos(thetaSky) - sin(thetaLOS)*sin(thetaSky)*cos(chiSky)
		print, '  Q_orig=', (3.d0*ct^2-1.d0) * sin(thetaSky)^2*cos(2.d0*(chiSky)), ' - U_orig=',(3.d0*ct^2-1.d0) * sin(thetaSky)^2*sin(2.d0*(chiSky))

		sinthSky2 = 1.d0 - (cos(thetaLOS)*cos(thetaB) + sin(thetaLOS)*sin(thetaB)*cos(chiB))^2
		cos2PChiSky = 2.d0 * (cos(thetaLOS)*sin(thetaB)*cos(chiB)-sin(thetaLOS)*cos(thetaB))^2 / sinthSky2 - 1.d0
		sin2PChiSky = 2.d0 * (sin(thetaB)*sin(chiB)*(cos(thetaLOS)*sin(thetaB)*cos(chiB)-sin(thetaLOS)*cos(thetaB))) / sinthSky2

		print, '  Q_orig=', (3.d0*cos(thetaB)^2-1.d0) * sinthSky2 * cos2PChiSky, ' - U_orig=',(3.d0*cos(thetaB)^2-1.d0) * sinthSky2 * sin2PChiSky
	endif

; Get polarity of original Stokes V
	if (keyword_set(StokesVpolarity)) then begin
		polarityOriginal = StokesVpolarity
	endif else begin
		polarityOriginal = polarity(thetaSky)
	endelse

;*************
; 1st solution -> PhiSky' = PhiSky
;*************
	A = -3.d0*cos(thetaLOS)^2 + 3.d0*sin(thetaLOS)^2 * cos(chiSky)^2
	B = 3.d0*cos(thetaLOS)^2 - 1.d0
	C = -6.d0*cos(thetaLOS)*sin(thetaLOS)*cos(chiSky)

; Original value of the angular dependence of Q and U
	costhetaB = cos(thetaLOS)*cos(thetaSky) - sin(thetaLOS)*sin(thetaSky)*cos(chiSky)
	K = (3.d0*costhetaB^2-1.d0) * sin(thetaSky)^2

	coef = dblarr(5)
	coef[0] = K^2
	coef[1] = -2.d0*K*B
	coef[2] = -2.d0*A*K+B^2
	coef[3] = -C^2+2.d0*B*A
	coef[4] = C^2+A^2

	roots = fz_roots(coef, eps=1.d-16, /double)

	roots_final = real_part(filterSolutions(sqrt(roots), nsolutions))
	

	if (keyword_set(verbose)) then print, 'PhiSky_new=PhySky'
	for i = 0, nsolutions-1 do begin

		thetaSky_new = asin(roots_final[i])	
		thetaSky_new = pickThetaSkyFromSolutions(thetaLOS, thetaSky, chiSky, thetaSky_new, chiSky)
		
		if ((polarity(thetaSky_new) eq polarityOriginal and yesStokesV eq 1) or (yesStokesV eq 0)) then begin
			temp = sky_to_absolute_general(thetaLOS, thetaSky_new, chiSky)
			thetaB_out[nsol] = temp[0]
			chiB_out[nsol] = temp[1]

			thetaSky_out[nsol] = thetaSky_new
			chiSky_out[nsol] = chiSky

			if (keyword_set(verbose)) then begin
				ct = cos(thetaLOS)*cos(thetaSky_new) - sin(thetaLOS)*sin(thetaSky_new)*cos(chiSky)
				print, 'Solution '+strtrim(string(nsol),2), thetaB_out[nsol]/conv, chiB_out[nsol]/conv, $
					'  Q=', (3.d0*ct^2-1.d0) * sin(thetaSky_new)^2*cos(2.d0*(chiSky)), ' - U=',(3.d0*ct^2-1.d0) * sin(thetaSky_new)^2*sin(2.d0*(chiSky)),$
					thetaSky_new/conv, chiSky/conv

				sinthSky2 = 1.d0 - (cos(thetaLOS)*cos(thetaB) + sin(thetaLOS)*sin(thetaB)*cos(chiB))^2
				cos2PChiSky = 2.d0 * (cos(thetaLOS)*sin(thetaB)*cos(chiB)-sin(thetaLOS)*cos(thetaB))^2 / sinthSky2 - 1.d0
				sin2PChiSky = 2.d0 * (sin(thetaB)*sin(chiB)*(cos(thetaLOS)*sin(thetaB)*cos(chiB)-sin(thetaLOS)*cos(thetaB))) / sinthSky2

; 				print, '  Q=', (3.d0*cos(thetaB)^2-1.d0) * sinthSky2 * cos2PChiSky, ' - U=',(3.d0*cos(thetaB)^2-1.d0) * sinthSky2 * sin2PChiSky
			endif
			
			nsol = nsol + 1
		endif
	endfor

;*************
; 2nd solution -> PhiSky' = PhiSky+pi
;*************
	A = -3.d0*cos(thetaLOS)^2 + 3.d0*sin(thetaLOS)^2 * cos(chiSky)^2
	B = 3.d0*cos(thetaLOS)^2 - 1.d0
	C = 6.d0*cos(thetaLOS)*sin(thetaLOS)*cos(chiSky)

; Original value of the angular dependence of Q and U
	costhetaB = cos(thetaLOS)*cos(thetaSky) - sin(thetaLOS)*sin(thetaSky)*cos(chiSky)
	K = (3.d0*costhetaB^2-1.d0) * sin(thetaSky)^2

	coef = dblarr(5)
	coef[0] = K^2
	coef[1] = -2.d0*K*B
	coef[2] = -2.d0*A*K+B^2
	coef[3] = -C^2+2.d0*B*A
	coef[4] = C^2+A^2

	roots = fz_roots(coef, eps=1.d-16, /double)

	roots_final = real_part(filterSolutions(sqrt(roots), nsolutions))

	theta = findgen(200) / 199.d0 * !DPI

	ct = cos(thetaLOS)*cos(theta) - sin(thetaLOS)*sin(theta)*cos(chiSky+!DPI)

	if (keyword_set(verbose)) then print, 'PhiSky_new=PhySky+pi'
	for i = 0, nsolutions-1 do begin

		thetaSky_new = asin(roots_final[i])
		thetaSky_new = pickThetaSkyFromSolutions(thetaLOS, thetaSky, chiSky, thetaSky_new, chiSky+!DPI)
		
		if ((polarity(thetaSky_new) eq polarityOriginal and yesStokesV eq 1) or (yesStokesV eq 0)) then begin
			temp = sky_to_absolute_general(thetaLOS, thetaSky_new, chiSky+!DPI)
			thetaB_out[nsol] = temp[0]
			chiB_out[nsol] = temp[1]

			thetaSky_out[nsol] = thetaSky_new
			chiSky_out[nsol] = chiSky+!DPI

			if (keyword_set(verbose)) then begin
				ct = cos(thetaLOS)*cos(thetaSky_new) - sin(thetaLOS)*sin(thetaSky_new)*cos(chiSky+!DPI)
				print, 'Solution '+strtrim(string(nsol),2), thetaB_out[nsol]/conv, chiB_out[nsol]/conv,$
					'  Q=', (3.d0*ct^2-1.d0) * sin(thetaSky_new)^2*cos(2.d0*(chiSky+!DPI)), ' - U=',(3.d0*ct^2-1.d0) * sin(thetaSky_new)^2*sin(2.d0*(chiSky+!DPI)),$
					thetaSky_new/conv, (chiSky+!DPI)/conv

				sinthSky2 = 1.d0 - (cos(thetaLOS)*cos(thetaB) + sin(thetaLOS)*sin(thetaB)*cos(chiB))^2
				cos2PChiSky = 2.d0 * (cos(thetaLOS)*sin(thetaB)*cos(chiB)-sin(thetaLOS)*cos(thetaB))^2 / sinthSky2 - 1.d0
				sin2PChiSky = 2.d0 * (sin(thetaB)*sin(chiB)*(cos(thetaLOS)*sin(thetaB)*cos(chiB)-sin(thetaLOS)*cos(thetaB))) / sinthSky2

; 				print, '  Q=', (3.d0*cos(thetaB)^2-1.d0) * sinthSky2 * cos2PChiSky, ' - U=',(3.d0*cos(thetaB)^2-1.d0) * sinthSky2 * sin2PChiSky
			endif
			
			nsol = nsol + 1			
		endif		
	endfor

;*************
; 3th solution -> PhiSky' = PhiSky+pi/2
;*************
	A = -3.d0*cos(thetaLOS)^2 + 3.d0*sin(thetaLOS)^2 * sin(chiSky)^2
	B = 3.d0*cos(thetaLOS)^2 - 1.d0
	C = 6.d0*cos(thetaLOS)*sin(thetaLOS)*sin(chiSky)

; Original value of the angular dependence of Q and U
	costhetaB = cos(thetaLOS)*cos(thetaSky) - sin(thetaLOS)*sin(thetaSky)*cos(chiSky)
	K = -(3.d0*costhetaB^2-1.d0) * sin(thetaSky)^2

	coef = dblarr(5)
	coef[0] = K^2
	coef[1] = -2.d0*K*B
	coef[2] = -2.d0*A*K+B^2
	coef[3] = -C^2+2.d0*B*A
	coef[4] = C^2+A^2

	roots = fz_roots(coef, eps=1.d-16, /double)

	roots_final = real_part(filterSolutions(sqrt(roots), nsolutions))

	if (keyword_set(verbose)) then print, 'PhiSky_new=PhySky+pi/2'
	for i = 0, nsolutions-1 do begin
		thetaSky_new = asin(roots_final[i])
		thetaSky_new = pickThetaSkyFromSolutions(thetaLOS, thetaSky, chiSky, thetaSky_new, chiSky+!DPI/2.d0)
		
		if ((polarity(thetaSky_new) eq polarityOriginal and yesStokesV eq 1) or (yesStokesV eq 0)) then begin
			temp = sky_to_absolute_general(thetaLOS, thetaSky_new, chiSky+!DPI/2.d0)
			thetaB_out[nsol] = temp[0]
			chiB_out[nsol] = temp[1]

			thetaSky_out[nsol] = thetaSky_new
			chiSky_out[nsol] = chiSky+!DPI/2.d0

			if (keyword_set(verbose)) then begin
				ct = cos(thetaLOS)*cos(thetaSky_new) - sin(thetaLOS)*sin(thetaSky_new)*cos(chiSky+!DPI/2.d0)
				print, 'Solution '+strtrim(string(nsol),2), thetaB_out[nsol]/conv, chiB_out[nsol]/conv,$
					'  Q=', (3.d0*ct^2-1.d0) * sin(thetaSky_new)^2*cos(2.d0*(chiSky+!DPI/2.d0)), ' - U=',(3.d0*ct^2-1.d0) * sin(thetaSky_new)^2*sin(2.d0*(chiSky+!DPI/2.d0)),$
					thetaSky_new/conv, (chiSky+!DPI/2.d0)/conv

				sinthSky2 = 1.d0 - (cos(thetaLOS)*cos(thetaB) + sin(thetaLOS)*sin(thetaB)*cos(chiB))^2
				cos2PChiSky = 2.d0 * (cos(thetaLOS)*sin(thetaB)*cos(chiB)-sin(thetaLOS)*cos(thetaB))^2 / sinthSky2 - 1.d0
				sin2PChiSky = 2.d0 * (sin(thetaB)*sin(chiB)*(cos(thetaLOS)*sin(thetaB)*cos(chiB)-sin(thetaLOS)*cos(thetaB))) / sinthSky2

; 				print, '  Q=', (3.d0*cos(thetaB)^2-1.d0) * sinthSky2 * cos2PChiSky, ' - U=',(3.d0*cos(thetaB)^2-1.d0) * sinthSky2 * sin2PChiSky
			endif
			
			nsol = nsol + 1
		endif
	endfor

;*************
; 4th solution -> PhiSky' = PhiSky-pi/2
;*************
	A = -3.d0*cos(thetaLOS)^2 + 3.d0*sin(thetaLOS)^2 * sin(chiSky)^2
	B = 3.d0*cos(thetaLOS)^2 - 1.d0
	C = -6.d0*cos(thetaLOS)*sin(thetaLOS)*sin(chiSky)

; Original value of the angular dependence of Q and U
	costhetaB = cos(thetaLOS)*cos(thetaSky) - sin(thetaLOS)*sin(thetaSky)*cos(chiSky)
	K = -(3.d0*costhetaB^2-1.d0) * sin(thetaSky)^2

	coef = dblarr(5)
	coef[0] = K^2
	coef[1] = -2.d0*K*B
	coef[2] = -2.d0*A*K+B^2
	coef[3] = -C^2+2.d0*B*A
	coef[4] = C^2+A^2

	roots = fz_roots(coef, eps=1.d-16, /double)

	roots_final = real_part(filterSolutions(sqrt(roots), nsolutions))

	if (keyword_set(verbose)) then print, 'PhiSky_new=PhySky-pi/2'
	for i = 0, nsolutions-1 do begin

		thetaSky_new = asin(roots_final[i])
		thetaSky_new = pickThetaSkyFromSolutions(thetaLOS, thetaSky, chiSky, thetaSky_new, chiSky-!DPI/2.d0)
		
		if ((polarity(thetaSky_new) eq polarityOriginal and yesStokesV eq 1) or (yesStokesV eq 0)) then begin
			temp = sky_to_absolute_general(thetaLOS, thetaSky_new, chiSky-!DPI/2.d0)
			thetaB_out[nsol] = temp[0]
			chiB_out[nsol] = temp[1]

			thetaSky_out[nsol] = thetaSky_new
			chiSky_out[nsol] = chiSky-!DPI/2.d0

			if (keyword_set(verbose)) then begin
				ct = cos(thetaLOS)*cos(thetaSky_new) - sin(thetaLOS)*sin(thetaSky_new)*cos(chiSky-!DPI/2.d0)
				print, 'Solution '+strtrim(string(nsol),2), thetaB_out[nsol]/conv, chiB_out[nsol]/conv, $
					'  Q=', (3.d0*ct^2-1.d0) * sin(thetaSky_new)^2*cos(2.d0*(chiSky-!DPI/2.d0)), ' - U=',(3.d0*ct^2-1.d0) * sin(thetaSky_new)^2*sin(2.d0*(chiSky-!DPI/2.d0)),$
					thetaSky_new/conv, (chiSky-!DPI/2.d0)/conv

				sinthSky2 = 1.d0 - (cos(thetaLOS)*cos(thetaB) + sin(thetaLOS)*sin(thetaB)*cos(chiB))^2
				cos2PChiSky = 2.d0 * (cos(thetaLOS)*sin(thetaB)*cos(chiB)-sin(thetaLOS)*cos(thetaB))^2 / sinthSky2 - 1.d0
				sin2PChiSky = 2.d0 * (sin(thetaB)*sin(chiB)*(cos(thetaLOS)*sin(thetaB)*cos(chiB)-sin(thetaLOS)*cos(thetaB))) / sinthSky2

; 				print, '  Q=', (3.d0*cos(thetaB)^2-1.d0) * sinthSky2 * cos2PChiSky, ' - U=',(3.d0*cos(thetaB)^2-1.d0) * sinthSky2 * sin2PChiSky
			endif
			
			nsol = nsol + 1
		endif
	endfor

; Cut the vector of angles to those that are compatible with the observations
	thetaB_out = thetaB_out[0:nsol-1]
	chiB_out = chiB_out[0:nsol-1]

	thetaSky_out = thetaSky_out[0:nsol-1]
	chiSky_out = chiSky_out[0:nsol-1]
	
	reorderAngles, thetaB_out, chiB_out, thetaSky_out, chiSky_out, thetaB, chiB
	
	thetaB_out = thetaB_out / conv
	chiB_out = chiB_out / conv

	thetaSky_out = thetaSky_out / conv
	chiSky_out = chiSky_out / conv

end
