@plot_observation
@plot_result
@plot_observation_result
@generate_conf_files
@wavelength_component
@plot_result
@hazel_inv_init
@hazel_inv_event
@run_code

pro hazel_inv, reset_state=reset_state   
	
	info = inv_init(reset_state=reset_state)

   widget_control, info.baseWidget, /REALIZE
   xmanager, 'Inversion', info.baseWidget, EVENT_HANDLER='Inv_Event'
end
