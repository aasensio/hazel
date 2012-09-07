pro change_symbol, which, size, fill
	 if (which eq 'circle') then begin
	 	  symbo = FINDGEN(17) * (!PI*2/16.)
	 	  USERSYM, size*COS(symbo), size*SIN(symbo), FILL=fill
	 endif
	 
	 if (which eq 'rhomb') then begin
	 	  x = [-size,0,size,0,-size]
		  y = [0,size,0,-size,0]
	 	  USERSYM, x,y, FILL=fill
	 endif
end
