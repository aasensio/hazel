pro draw_sphere2, state

	erase	
	t3d, /reset
	scale3, xran=[-3,3], yran=[-3,3], zran=[-3,3], ax=state.thetaObs, az=state.chiObs
	
	n = n_elements(state.x)
	
	for i = 0, n-1 do begin
		plots, (*state.x)[i], (*state.y)[i], (*state.z)[i],/t3d,/data,psym=4
	endfor
		
	for i = 0, state.quality-1 do begin		
		for j = 0, state.quality-1 do begin						
			k_from = (*state.links)[i,j mod state.quality]
			k_to = (*state.links)[i,(j+1) mod state.quality]			
			plots, (*state.x)[k_from],(*state.y)[k_from],(*state.z)[k_from],/t3d,/data
			plots,(*state.x)[k_to],(*state.y)[k_to],(*state.z)[k_to],/t3d,/data,/continue
		endfor
	endfor
	
	for i = 0, state.quality-1 do begin		
		for j = 0, state.quality-2 do begin						
			k_from = (*state.links)[j mod state.quality,i]
			k_to = (*state.links)[(j+1) mod state.quality,i]			
			plots, (*state.x)[k_from],(*state.y)[k_from],(*state.z)[k_from],/t3d,/data
			plots,(*state.x)[k_to],(*state.y)[k_to],(*state.z)[k_to],/t3d,/data,/continue
		endfor
	endfor
			
	thetaB = state.thetaB * !DPI / 180.d0
	chiB = state.chiB * !DPI / 180.d0
	
	x_from = [state.r, 0.d0, 0.d0]
	x_to = x_from + [cos(thetaB), sin(thetaB)*sin(chiB), sin(thetaB)*cos(chiB)]
		
	plots, x_from[0], x_from[1], x_from[2], /t3d, /data
	plots, x_to[0], x_to[1], x_to[2], /t3d, /data, thick=2, /continue
	
	mu = strtrim(string(cos(state.thetaObs * !DPI / 180.d0)),2)
	xyouts, -2.5d0, 1.5d0, mu
	
	
	
end

pro draw_sphere, state

; 	erase	
	t3d, /reset
	scale3, xran=[-3,3], yran=[-3,3], zran=[-3,3], ax=state.thetaObs, az=state.chiObs
	r = replicate(state.r,state.quality,state.quality)
	
	mesh_obj,4,v,p,r
	pshades = byte(255*randomu(4,n_elements(p[*,0]),n_elements(p[0,*])))
	
	tvscl, polyshade(v,p,/t3d,poly_shades=pshades)
		
	thetaB = state.thetaB * !DPI / 180.d0
	chiB = state.chiB * !DPI / 180.d0
	
	x_from = [state.r, 0.d0, 0.d0]
	x_to = x_from + [cos(thetaB), sin(thetaB)*sin(chiB), sin(thetaB)*cos(chiB)]
		
	plots, x_from[0], x_from[1], x_from[2], /t3d, /data
	plots, x_to[0], x_to[1], x_to[2], /t3d, /data, thick=2, /continue
	
	mu = strtrim(string(cos(state.thetaObs * !DPI / 180.d0)),2)
	xyouts, -2.5d0, 1.5d0, mu
	
end

pro define_sphere, state
	
	*state.x = findgen(state.quality*state.quality)
	*state.y = findgen(state.quality*state.quality)
	*state.z = findgen(state.quality*state.quality)
	*state.links = intarr(state.quality,state.quality)
		
	k = 0
	state.r = 3.d0
	
	for i = 0, state.quality-1 do begin
		theta = i * !DPI / (state.quality-1.d0)		
		for j = 0, state.quality-1 do begin
			phi = j * 2.d0 * !DPI / state.quality
			
			(*state.x)[k] = state.r*sin(theta)*cos(phi)
			(*state.y)[k] = state.r*sin(theta)*sin(phi)
			(*state.z)[k] = state.r*cos(theta)
			
			(*state.links)[i,j] = k
			
			k = k + 1
		endfor
	endfor		
	
end


function sphere_init
	
; Define the points of the sphere	
	quality = 20
; 	state = {baseWidget: 0L, thetaObs: 0.d0, chiObs: 0.d0, thetaB: 0.d0, chiB: 0.d0, $
; 		quality : 19, x : findgen(quality*quality), y : findgen(quality*quality), $
; 		z : findgen(quality*quality), r : 0.d0, links : intarr(quality,quality)}
	state = {baseWidget: 0L, thetaObs: 0.d0, chiObs: 0.d0, thetaB: 0.d0, chiB: 0.d0, $
		quality : 19, x : ptr_new(0), y : ptr_new(0), $
		z : ptr_new(0), r : 0.d0, links : ptr_new(0)}
			
	define_sphere, state
			
; Base widget	
	state.baseWidget = widget_base(TITLE='Sphere')
	
	sphereBase = widget_base(state.baseWidget, /COLUMN, /BASE_ALIGN_CENTER)
	
	plotWidget = widget_draw(sphereBase, XSIZE=650, YSIZE=450, /FRAME)
	
	slidersBase = widget_base(sphereBase, /ROW)
	
	precisionSlider = widget_slider(sphereBase,$
		TITLE='Quality of the plot',UVALUE='qualitySlider',$
		XSIZE=255,MAXIMUM=50.d0,VALUE=quality, /DRAG)
		
	obsAngleBase = widget_base(slidersBase, /COLUMN)
		
	temp = widget_label(obsAngleBase,VALUE='Observation angles')
	thetaOSlider = widget_slider(obsAngleBase,$
		TITLE='Observing theta angle [deg]',UVALUE='thetaObsSlider',$
		XSIZE=255,MAXIMUM=180.d0,VALUE=0.d0, /DRAG)

	chiOSlider = widget_slider(obsAngleBase,$
		TITLE='Observing chi angle [deg]',UVALUE='chiObsSlider',$
		XSIZE=255,MAXIMUM=360.d0,VALUE=0.d0, /DRAG)
		
	magAngleBase = widget_base(slidersBase, /COLUMN)
	
	temp = widget_label(magAngleBase,VALUE='Magnetic field angles')
	thetaBSlider = widget_slider(magAngleBase,$
		TITLE='Magnetic theta angle [deg]',UVALUE='thetaBSlider',$
		XSIZE=255,MAXIMUM=180.d0,VALUE=0.d0, /DRAG)

	chiBSlider = widget_slider(magAngleBase,$
		TITLE='Magnetic chi angle [deg]',UVALUE='chiBSlider',$
		XSIZE=255,MAXIMUM=360.d0,VALUE=0.d0, /DRAG)
		
	widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
			
	return, state
end

;-----------------------------------------
; Event handler
;-----------------------------------------
pro sphere_event, event
	widget_control, Event.id, GET_UVALUE=Action

	hand = widget_info(Event.Handler, /CHILD)
	widget_control, hand, GET_UVALUE=state
		
	case Action of
		'thetaBSlider' : 	begin
									widget_control, Event.id, GET_VALUE=value
									state.thetaB = value
									widget_control, hand, SET_UVALUE=state
									draw_sphere, state
      						end
     	'chiBSlider' : 	begin
									widget_control, Event.id, GET_VALUE=value
									state.chiB = value
									widget_control, hand, SET_UVALUE=state
									draw_sphere, state
      						end
      'thetaObsSlider' : 	begin
									widget_control, Event.id, GET_VALUE=value
									state.thetaObs = value
									widget_control, hand, SET_UVALUE=state
									draw_sphere, state
      						end
     	'chiObsSlider' : 	begin
									widget_control, Event.id, GET_VALUE=value
									state.chiObs = value
									widget_control, hand, SET_UVALUE=state
									draw_sphere, state
      						end
      'qualitySlider' : 	begin
									widget_control, Event.id, GET_VALUE=value
									state.quality = value
									widget_control, hand, SET_UVALUE=state
									define_sphere, state
									draw_sphere, state
      						end
	endcase
end

pro sphereB
	state = sphere_init()
	
	widget_control, state.baseWidget, /REALIZE
	
	draw_sphere, state

	xmanager, 'SphereB', state.baseWidget, EVENT_HANDLER='sphere_event'
end
