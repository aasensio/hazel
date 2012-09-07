@hazel_routines
@hazel_plot
@hazel_variation
@hazel_synth
@hazel_init
@hazel_event

;-----------------------------------------
; Main routine
;-----------------------------------------
pro hazel, reset_state=reset_state
	 state = hazel_init(reset_state=reset_state)
	 
	 widget_control, state.baseWidget, /REALIZE
	 
	 xmanager, 'HAZEL', state.baseWidget, EVENT_HANDLER='hazel_Event'
end

;-----------------------------------------
; Add a new property to the state structure
;-----------------------------------------
pro add_to_state
	restore,'state.idl'
	temp = create_struct(state, 'randomazimuth', 0L)
	state = temp
	save, state, filename='state.idl'
end
