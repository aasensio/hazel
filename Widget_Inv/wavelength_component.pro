; Return the wavelength of the multiplet component of the given atom
function wavelength_component, atom, multiplet
	case (atom) of
		0: file = 'ATOMS/helium.mod'
		1: file = 'ATOMS/sulfur.mod'
	endcase
	
	temp = ''
	openr,2,file
	readf,2,s2
	readf,2,n
	for i = 0, n-1 do begin
		readf,2, index, l2
		minim = abs(l2-s2)/2
		maxim = (l2+s2)/2
		for j = minim, maxim do begin
			readf,2,temp
		endfor
	endfor
	readf,2,n
	wl = 0.d0
	for i = 0, n-1 do begin
		readf,2,index, low, up, Aul, wl
		if (i eq multiplet) then begin
			close,2
			return, wl
		endif
	endfor
	close,2	
end
