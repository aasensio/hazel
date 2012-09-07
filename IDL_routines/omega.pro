function omega, wl, h
	 Ic = ddread('ic.dat',/noverb)
; Unit conversion

	 PC = 2.99792458d10
	 Ic(1,*) = Ic(1,*) * Ic(0,*)^2 / (PC*1.d4)
	 Ic(0,*) = Ic(0,*) * 1.d4
	 
	 coeff = ddread('cl.dat',/noverb)
	 coeff(0,*) = coeff(0,*) * 1.d4
	 
	 Rsun = 976.6d0
	 
	 sg = Rsun / (h+Rsun)
	 cg = sqrt(1.d0-sg^2)
	 
	 a0 = 1.d0 - cg
	 a1 = cg - 0.5d0 - 0.5d0*cg^2/sg*alog((1.d0+sg)/cg)
	 a2 = (cg+2.d0)*(cg-1.d0) / (3.d0*(cg+1.d0))
	 
	 b0 = (1.d0-cg^3) / 3.d0
	 b1 = (8.d0*cg^3-3.d0*cg^2-2.d0) / 24.d0 - cg^4 / (8.d0*sg) * alog((1.d0+sg)/cg)
	 b2 = (cg-1.d0)*(3.d0*cg^3+6.d0*cg^2+4.d0*cg+2.d0) / (15.d0*(cg+1.d0))
	 
	 I0 = interpol(Ic(1,*),Ic(0,*),wl)
	 u1 = interpol(coeff(1,*),coeff(0,*),wl)
	 u2 = interpol(coeff(2,*),coeff(0,*),wl)

	 print, wl, I0, u1, u2
	 stop
	 
	 J = 0.5d0 * I0 * (a0 + a1*u1 + a2*u2)
	 K = 0.5d0 * I0 * (b0 + b1*u1 + b2*u2)
	 
	 return, (3.d0*K-J)/(2.d0*J)	 
	 
end
