; Return the angles in the plane of the sky given angles with respect
; to the vertical for observations on the limb
function absolute_to_sky, thetaB, chiB
	t1 = sin(thetaB) * sin(chiB)
	t2 = -cos(thetaB)
	t3 = sin(thetaB) * cos(chiB)

	thetaSky = acos(t3)
	sinthSky = sqrt(1.d0-t3^2)

	sinchiSky = t1 / sinthSky
	coschiSky = t2 / sinthSky

; Test for the quadrant
	chiSky_preliminary = acos(coschiSky)
	if (sign(sinchiSky) ge 0.d0) then begin
		chiSky = chiSky_preliminary
	endif else begin
		chiSky = -chiSky_preliminary
	endelse

	return, [thetaSky, chiSky]
end

; Return the angles in the vertical system given angles in the plane of the sky
; for observations on the limb
function sky_to_absolute, thetaSky, chiSky
	t1 = sin(thetaSky) * sin(chiSky)
	t2 = cos(thetaSky)
	t3 = -sin(thetaSky) * cos(chiSky)

	thetaB = acos(t3)
	sinthB = sqrt(1.d0-t3^2)

	sinchiB = t1 / sinthB
	coschiB = t2 / sinthB

; Test for the quadrant
	chiB_preliminary = acos(coschiB)
	if (sign(sinchiB) ge 0.d0) then begin
		chiB = chiB_preliminary
	endif else begin
		chiB = -chiB_preliminary
	endelse

	return, [thetaB, chiB]
end


; Return the angles in the plane of the sky given angles with respect
; to the vertical for observations at angle theta
function absolute_to_sky_general, theta, thetaB, chiB

	cosThetaSky = cos(theta) * cos(thetaB) + sin(theta) * sin(thetaB) * cos(chiB)
	sinThetaSky = sqrt(1.d0 - cosThetaSky^2)

	thetaSky = acos(cosThetaSky)

	cosChiSky = ( cos(theta) * sin(thetaB) * cos(chiB) - cos(thetaB) * sin(theta) ) / sinThetaSky
	sinChiSky = ( sin(thetaB) * sin(chiB) ) / sinThetaSky

; Test for the quadrant
	chiSky_preliminary = acos(cosChiSky < 1.d0 > (-1.d0))
	if (sign(sinchiSky) ge 0.d0) then begin
		chiSky = chiSky_preliminary
	endif else begin
		chiSky = -chiSky_preliminary
	endelse

	return, [thetaSky, chiSky]
end

; Return the angles in the plane of the sky given angles with respect
; to the vertical for observations at angle theta
function sky_to_absolute_general, theta, thetaSky, chiSky

	cosThetaB = cos(theta) * cos(thetaSky) - sin(theta) * sin(thetaSky) * cos(chiSky)
	sinThetaB = sqrt(1.d0 - cosThetaB^2)

	thetaB = acos(cosThetaB)

	cosChiB = ( cos(theta) * sin(thetaSky) * cos(chiSky) + cos(thetaSky) * sin(theta) ) / sinThetaB
	sinChiB = ( sin(thetaSky) * sin(chiSky) ) / sinThetaB

; Test for the quadrant
	chiB_preliminary = acos(cosChiB < 1.d0 > (-1.d0))
	if (sign(sinchiB) ge 0.d0) then begin
		chiB = chiB_preliminary
	endif else begin
		chiB = -chiB_preliminary
	endelse

	return, [thetaB, chiB]
end