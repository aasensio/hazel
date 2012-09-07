function nbar, wl, h
	 Ic = ddread('ic.dat',/noverb)
; Unit conversion
	 wvl = Ic(0,*)

	 PC = 2.99792458d10
    PH = 6.62606876d-27
   
	 Ic(1,*) = Ic(1,*) * wvl^2 / (PC*1.d4)
	 Ic(0,*) = Ic(0,*) * 1.d4
	 
	 coeff = ddread('cl.dat',/noverb)
	 coeff(0,*) = coeff(0,*) * 1.d4
	 
	 Rsun = 976.6d0
	 
	 sg = Rsun / (h+Rsun)
	 cg = sqrt(1.d0-sg^2)
	 
	 a0 = 1.d0 - cg
	 a1 = cg - 0.5d0 - 0.5d0*cg^2/sg*alog((1.d0+sg)/cg)
	 a2 = (cg+2.d0)*(cg-1.d0) / (3.d0*(cg+1.d0))
	 
	 I0 = interpol(Ic(1,*),Ic(0,*),wl)
	 u1 = interpol(coeff(1,*),coeff(0,*),wl)
	 u2 = interpol(coeff(2,*),coeff(0,*),wl)
	 
	 J = 0.5d0 * I0 * (a0 + a1*u1 + a2*u2)
	 
	 return, 1d10*(wl^3*1d-24/(2d0*ph*pc)) * J
	 
end
