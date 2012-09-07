pro genprof
	restore,'Valignment+zeeman.sav'
	restore,'24apr11.004_cal.sav'

	openw,2,'orient_zeeman.prof',width=150
	printf,2,103
	for i = 47, 149 do begin
		printf,2,lambda[i]-10829.0911d0, int[i]/max(int), sq[i]/max(int), su[i]/max(int), $
			sv[i]/max(int), 0.01, 0.001d0, 0.001d0, 0.001d0
	endfor
	close,2

	stop
end