function i0_allen, wl, mu
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
   i0 = interpol(ic(1,*),ic(0,*),wl)

   imu = 1.d0 - u - v + u * mu + v * mu^2
   
   return, i0*imu
end
