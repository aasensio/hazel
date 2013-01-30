; Get the polarity of Stokes V
function polarity, thetaSky
	if (thetaSky le !DPI/2.d0) then polarity = 1
	if (thetaSky gt !DPI/2.d0 and thetaSky le !DPI) then polarity = -1
	return, polarity
end