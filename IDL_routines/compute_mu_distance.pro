; Compute the distance to the limb and the value of mu for the observations
; of the prominence
pro compute_mu_distance, nslit, nstep, delta_slit, delta_step, ind_limbo, from_linfit, to_linfit, mu, dist

	dist = findgen(nslit,nstep)
	side = findgen(nslit,nstep)

; Define the straight line of the limb
	x1 = from_linfit * delta_step
	y1 = ind_limbo[from_linfit] * delta_slit

	x2 = to_linfit * delta_step
	y2 = ind_limbo[to_linfit] * delta_slit

; Calculate the distance to the straight line
	for i = 0, nslit-1 do begin
		for j = 0, nstep-1 do begin
			x0 = j * delta_step
			y0 = i * delta_slit
			den = sqrt((x2-x1)^2 + (y2-y1)^2)
			num = abs( (x2-x1)*(y1-y0) - (x1-x0)*(y2-y1) )
			sign = 1.d0
			dist[i,j] = num / den

			side[i,j] = (x2-x1)*(y1-y0) - (x1-x0)*(y2-y1)
		endfor
	endfor

	Rsun = 960.d0

	ind = where(side gt 0.d0, count)

; Compute mu
	mu = sqrt( 1.d0 - (Rsun - dist)^2 / Rsun^2)

; Set all values outside the limb to mu=0
	if (count ne 0) then begin
		mu[ind] = 0.d0
	endif
end