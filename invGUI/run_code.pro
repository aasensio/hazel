pro run_code, info, inversion=inversion
	if (keyword_set(inversion)) then begin
		file_delete, 'done.info', /quiet, /allow_nonexistent
		
		spawn,'./hazel &';, result_spawn
	
		exists = 0
		while (exists eq 0) do begin
		;widget_control, info.outputText, SET_VALUE=result_spawn
			wait,2
			plot_observation_result, info.obs_file,'temporal.prof', info.plotWidget
			result = file_info('done.info')
			exists = result.exists
		endwhile
	endif else begin
		spawn,'./hazel'
	endelse
end