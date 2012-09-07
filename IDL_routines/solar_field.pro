; Return the radiation field of the Sun using interpolation
; Input:
;     wl : wavelength in A
;     mu : cosine of the heliocentric angle
; Output: radiation field in cgs units

function solar_field, wl, mu
	 
	ic = ddread('ic.dat',/noverb)
	cl = ddread('cl.dat',/noverb)
	
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