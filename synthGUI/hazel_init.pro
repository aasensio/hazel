;-----------------------------------------
; Initialization routine. Draws the widget
;-----------------------------------------
function hazel_init, reset_state=reset_state
	
	 sliderSize = 128
	 sliderSize2 = 80
	 if (keyword_set(reset_state) or file_test('state.idl') eq 0) then begin	 	  
	 	  state = {baseWidget: 0L, plotWidget: 0L, Bfield: 40.d0, Bfield2: 40.d0, Bfieldmax: 500.d0, Bfieldmax2: 500.d0, $
	 	  		thetaBfield: 30.d0, thetaBfield2: 30.d0, chiBfield: 90.d0, chiBfield2: 90.d0, $
	 	   	thetaObs: 90.d0, chiObs: 0.d0, gammaObs: 90.d0, Multiplet: 1L, Doppler: 8.d0, Doppler2: 8.d0, auto: 0L, BSlider: 0L, BSlider2: 0L, $
		   	paschen: 1L, effects: 0L, observation: 0L, stokes0: [1.d0,0.d0,0.d0,0.d0], dtau_desired: 3.d0, height: 20.d0, $
		   	postcript: 0L, dtaured: 3.d0, beta: 1.d0, beta2: 1.0, damping: 0.d0, stimulated: 0, magneto_opt: 0, obs_file: '', $
				i0_allen: 0L, D2: 0.d0, factor_10830_nbar: 1.d0, factor_10830_omega: 1.d0, waveaxis: [-3.d0,2.5d0,200.d0],$
				which_atom: 0L, MultipletSlider: 0L, normaliz: 1L, j10: 0.d0, which_code: 0,$
				bfield_var: [1.d-4,1.d3,15], which_rho_plot: 0L, which_refframe: 0L,$
				path_save: '.', randomazimuth: 0L, number_slabs: 1, dtau_desired2: 0.d0, vel: 0.d0, vel2: 0.d0, ff: 0.d0, useAllen: 1}
	endif else begin
		restore,'state.idl'
	 endelse
	 
	 if (state.which_code eq 0) then begin
	 	if (state.which_atom eq 0) then begin
	 	  	state.baseWidget = widget_base(TITLE='Hazel Synthesis : He I', MBAR=menuBar)	 	  	
	 	endif
	 
	 	if (state.which_atom eq 1) then begin
	 	  	state.baseWidget = widget_base(TITLE='Hazel Synthesis : Ca II', MBAR=menuBar)	 	  	
	 	endif
	 
	 	if (state.which_atom eq 2) then begin
	 	  	state.baseWidget = widget_base(TITLE='Hazel Synthesis : Na I', MBAR=menuBar)	 	  	
	 	endif
	 endif
	 
	 if (state.which_code eq 1) then begin
	 	if (state.which_atom eq 0) then begin
	 	  	state.baseWidget = widget_base(TITLE='Hazel Synthesis : Na I HFS', MBAR=menuBar)	 	  	
	 	endif
	 endif
	 
	 hanleBase = widget_base(state.baseWidget, /COLUMN, /BASE_ALIGN_CENTER)
	 
	 atomMenu = widget_button(menuBar, VALUE='Multiterm', /MENU)
	 heliumButton = widget_button(atomMenu, VALUE='He I', UVALUE='HELIUM')
	 sulfurButton = widget_button(atomMenu, VALUE='Ca II', UVALUE='CALCIUM')
	 sodiumButton = widget_button(atomMenu, VALUE='Na I', UVALUE='SODIUM')
	 
; 	 atomMenu = widget_button(menuBar, VALUE='Multilevel HFS', /MENU)	 
; 	 sodiumButton = widget_button(atomMenu, VALUE='Na I', UVALUE='SODIUM_HFS')
	 
	 horizBase = widget_base(hanleBase, /ROW)
	 
	 state.plotWidget = widget_draw(horizBase, XSIZE=650, YSIZE=450, /FRAME)
	 	 
	 rightbase = widget_base(horizBase, /COLUMN, FRAME=1)
	 ; collislb = widget_label(rightBase, VALUE='Lower level D^(2):')
	 ; collis = widget_text(rightBase, VALUE=strtrim(string(state.D2),2),UVALUE='D2',/EDITABLE,XSIZE=8,YSIZE=1)
	 ; factor_10830lb = widget_label(rightBase, VALUE='nbar factor')
	 ; factor_10830 = widget_text(rightBase, VALUE=strtrim(string(state.factor_10830_nbar),2),UVALUE='fact_10830_nbar',/EDITABLE,XSIZE=8,YSIZE=1)
	 ; factor_10830_wlb = widget_label(rightBase, VALUE='w factor')
	 ; factor_10830_w = widget_text(rightBase, VALUE=strtrim(string(state.factor_10830_omega),2),UVALUE='fact_10830_w',/EDITABLE,XSIZE=8,YSIZE=1)
	 ; j10_wlb = widget_label(rightBase, VALUE='J10/J00')
	 ; j10 = widget_text(rightBase, VALUE=strtrim(string(state.j10),2),UVALUE='j10_tensor',/EDITABLE,XSIZE=8,YSIZE=1)
	 
; Emission or tangent observation
	 t4 = widget_base(rightBase, /COLUMN, FRAME=1)
	 lab4 = widget_label(t4, VALUE='Observation')	 
	 obsBase = widget_base(t4, /COLUMN, /EXCLUSIVE)
	 emissionButton = widget_button(obsBase, VALUE='Optically thin', UVALUE='EMISSION')
	 ; simpleslabButton = widget_button(obsBase, VALUE='Simplified', UVALUE='SIMPLE_SLAB')
	 ; formalButton = widget_button(obsBase, VALUE='No MO', UVALUE='FORMAL')
	 ; deloparButton = widget_button(obsBase, VALUE='DELO', UVALUE='DELOPAR')
	 exactslabButton = widget_button(obsBase, VALUE='Exact', UVALUE='EXACT_SLAB')
; 	 milneButton = widget_button(obsBase, VALUE='Milne-Eddington', UVALUE='MILNE')
	 
	 case state.observation of
	 	0: begin
	 			widget_control, emissionButton, /SET_BUTTON   ; 0 -> pure emission
	 		end
	 	1: begin
	 			widget_control, formalButton, /SET_BUTTON   ; 1 -> slab (neglecting MO terms)
	 		end
	 	2: begin
	 			widget_control, milneButton, /SET_BUTTON   ; 2 -> Milne-Eddington
	 		end
	 	3: begin
	 			widget_control, deloparButton, /SET_BUTTON   ; 3 -> DELO
	 		end
	 	4: begin
				widget_control, simpleslabButton, /SET_BUTTON   ; 4 -> simplified slab (optically thin)
			end
	 	5: begin
	 			widget_control, exactslabButton, /SET_BUTTON   ; 5 -> exact slab
	 		end
	 endcase
	 
; Wavelength axis
	 waveaxis = widget_base(t4, /COLUMN, FRAME=1)
	 lab0 = widget_label(waveaxis, VALUE='Wavelength axis')
	 wleftbase = widget_base(waveaxis, /ROW)
	 wleftlab = widget_label(wleftbase, VALUE='wl:')
	 wleft = widget_text(wleftbase, VALUE=strtrim(string(state.waveaxis(0)),2),$
	 	UVALUE='wleft',/EDITABLE,XSIZE=5,YSIZE=1)
	 ;wrightbase = widget_base(waveaxis, /ROW)
	 wrightlab = widget_label(wleftbase, VALUE='wr:')
	 wright = widget_text(wleftbase, VALUE=strtrim(string(state.waveaxis(1)),2),$
	 	UVALUE='wright',/EDITABLE,XSIZE=5,YSIZE=1)
	 wstepbase = widget_base(waveaxis, /ROW)
	 wsteplab = widget_label(wstepbase, VALUE='step:')
	 wstep = widget_text(wstepbase, VALUE=strtrim(string(state.waveaxis(2)),2),UVALUE='wstep',/EDITABLE,XSIZE=5,YSIZE=1)

; Postcript output
	 state.postcript = 0
	 t0 = widget_base(t4, /COLUMN, FRAME=1)
	 lab0 = widget_label(t0, VALUE='Output')
	 postcriptBase = widget_base(t0, /COLUMN, /EXCLUSIVE)
	 postcriptnoButton = widget_button(postcriptBase, VALUE='PS', UVALUE='PS_ON')
	 postcriptyesButton = widget_button(postcriptBase, VALUE='Graphical', UVALUE='PS_OFF')
	 postcriptwriteButton = widget_button(postcriptBase, VALUE='Save', UVALUE='PS_WRITE')	 
	 widget_control, postcriptyesButton, /SET_BUTTON
	 
; Automatic synthesis	 
	 t1 = widget_base(t4, /COLUMN, FRAME=1)
	 lab1 = widget_label(t1, VALUE='Automatic')
	 autoBase = widget_base(t1, /COLUMN, /EXCLUSIVE)
	 autoyesButton = widget_button(autoBase, VALUE='Automatic', UVALUE='AUTO_ON')
	 autonoButton = widget_button(autoBase, VALUE='Manual', UVALUE='AUTO_OFF')
	 widget_control, (state.auto) ? autoyesButton : autonoButton, /SET_BUTTON

	 
	 slidersbase = widget_base(horizBase, /COLUMN)
	 ;slidersBase = widget_base(hanleBase, /ROW)
	 
; Magnetic field information
	 fieldBase_all = widget_base(slidersBase, /ROW, FRAME=1)
	 fieldBase = widget_base(fieldBase_all, /COLUMN, FRAME=1)
	 
; FIRST COMPONENT
; Magnetic field strength
	 BSliderBase = widget_base(fieldBase, /COLUMN)
	 maxBBase = widget_base(BSliderBase, /ROW)
	 BSlidertext = widget_label(maxBBase, VALUE='B_max:')
	 BSliderMax = widget_text(maxBBase,VALUE=strtrim(string(state.Bfieldmax),2),$
	 	UVALUE='BSliderMax',/EDITABLE,XSIZE=8,YSIZE=1)	 
	 state.BSlider = cw_fslider(BSliderBase,TITLE='B1 strength [G]',UVALUE='BSlider',$
	 	  XSIZE=sliderSize,MINIMUM=0.d0,MAXIMUM=state.Bfieldmax,VALUE=state.Bfield, EDIT=1)	 

; Magnetic field inclination
	 thetaBSliderBase = widget_base(fieldBase, /ROW)
	 thetaBSlider = widget_slider(thetaBSliderBase,TITLE='B1 incl. [deg]',UVALUE='thetaBSlider',$
	 	  XSIZE=sliderSize,MAXIMUM=180.d0,VALUE=state.thetaBfield)
	 
; Magnetic field azimuth	 
	 chiBSliderBase = widget_base(fieldBase, /ROW)
	 chiBSlider = widget_slider(chiBSliderBase,TITLE='B1 azimuth [deg]',UVALUE='chiBSlider',$
	 	  XSIZE=sliderSize,MINIMUM=-180.d0,MAXIMUM=180.d0,VALUE=state.chiBfield)
	 	  
; 	 chiBrandomBase = widget_base(chiBSliderBase, /COLUMN, /EXCLUSIVE)
	 chiBrandomBase = widget_base(fieldBase, /COLUMN, /EXCLUSIVE)
	 chiBrandButton = widget_button(chiBrandomBase, VALUE='Random azi', UVALUE='RANDOMAZI_ON')
	 chiBnorandButton = widget_button(chiBrandomBase, VALUE='Normal azi', UVALUE='RANDOMAZI_OFF')
	 widget_control, (state.randomazimuth) ? chiBrandButton : chiBnorandButton, /SET_BUTTON

; SECOND COMPONENT
; Magnetic field strength
	 fieldBase = widget_base(fieldBase_all, /COLUMN, FRAME=1)
	 BSliderBase = widget_base(fieldBase, /COLUMN)
	 maxBBase = widget_base(BSliderBase, /ROW)
	 BSlidertext = widget_label(maxBBase, VALUE='B_max:')
	 BSliderMax = widget_text(maxBBase,VALUE=strtrim(string(state.Bfieldmax2),2),$
	 	UVALUE='BSliderMax2',/EDITABLE,XSIZE=8,YSIZE=1)	 
	 state.BSlider2 = cw_fslider(BSliderBase,TITLE='B2 strength [G]',UVALUE='BSlider2',$
	 	  XSIZE=sliderSize,MINIMUM=0.d0,MAXIMUM=state.Bfieldmax2,VALUE=state.Bfield2, EDIT=1)

; Magnetic field inclination
	 thetaBSliderBase = widget_base(fieldBase, /ROW)
	 thetaBSlider = widget_slider(thetaBSliderBase,TITLE='B2 inclin [deg]',UVALUE='thetaBSlider2',$
	 	  XSIZE=sliderSize,MAXIMUM=180.d0,VALUE=state.thetaBfield2)
	 
; Magnetic field azimuth	 
	 chiBSliderBase = widget_base(fieldBase, /ROW)
	 chiBSlider = widget_slider(chiBSliderBase,TITLE='B2 azimuth [deg]',UVALUE='chiBSlider2',$
	 	  XSIZE=sliderSize,MINIMUM=-180.d0,MAXIMUM=180.d0,VALUE=state.chiBfield2)
		  

	 observationBase = widget_base(slidersBase, /ROW, FRAME=1)
	 thetaOSlider = widget_slider(observationBase,TITLE='Obs theta [deg]',UVALUE='thetaOSlider',$
	 	  XSIZE=sliderSize,MAXIMUM=180.d0,VALUE=state.thetaObs)

	 chiOSlider = widget_slider(observationBase,TITLE='Obs. chi [deg]',UVALUE='chiOSlider',$
	 	  XSIZE=sliderSize,MAXIMUM=360.d0,VALUE=state.chiObs)
		  
	 gammaOSlider = widget_slider(observationBase,TITLE='Obs. gamma [deg]',UVALUE='gammaOSlider',$
	 	  XSIZE=sliderSize,MAXIMUM=180.d0,VALUE=state.gammaObs)
		  	 		  
; Synthesize
	 synth_globBase = widget_base(hanleBase, /ROW, /ALIGN_LEFT)
	 	 
	 synthBase = widget_base(synth_globBase, /COLUMN)
	 multiBase = widget_base(synthBase, /ROW, FRAME=1)
	 if (state.which_code eq 0) then begin
	 	if (state.which_atom eq 0) then begin
	 	  	state.MultipletSlider =$
	 	  		widget_slider(multiBase,TITLE='Multiplet',UVALUE='MultipletSlider',$
	 	   	XSIZE=0.5*sliderSize,MAXIMUM=4,MINIMUM=1,VALUE=state.multiplet)	 
	 	endif
	 	if (state.which_atom eq 1) then begin
	 	  	state.MultipletSlider =$
	 	  		widget_slider(multiBase,TITLE='Multiplet',UVALUE='MultipletSlider',$
	 	   	XSIZE=0.5*sliderSize,MAXIMUM=4,MINIMUM=1,VALUE=state.multiplet,sensitive=0)
	 	endif
	 	if (state.which_atom eq 2) then begin
	 	  	state.MultipletSlider =$
	 	  		widget_slider(multiBase,TITLE='Multiplet',UVALUE='MultipletSlider',$
	 	   	XSIZE=0.5*sliderSize,MAXIMUM=4,MINIMUM=1,VALUE=state.multiplet,sensitive=0)
	 	endif
	 endif
	 if (state.which_code eq 1) then begin
	 	if (state.which_atom eq 0) then begin
	 	  	state.MultipletSlider =$
	 	  		widget_slider(multiBase,TITLE='Multiplet',UVALUE='MultipletSlider',$
	 	   	XSIZE=0.5*sliderSize,MAXIMUM=2,MINIMUM=1,VALUE=state.multiplet,sensitive=1)
	 	endif
	 endif
	 	 
	 heightSlider = widget_slider(multiBase,TITLE='Height ["]',UVALUE='heightSlider',$
	 	  XSIZE=sliderSize,MINIMUM=-100.d0,MAXIMUM=100.d0,VALUE=state.height)	 	  	 	 
	 
	 ; i0 = return_i0_allen(state)
	 ; state.i0_allen = widget_label(multiBase, VALUE='Allen: '+strtrim(string(i0),2))

	 effectBase = widget_base(multiBase, /COLUMN, /NONEXCLUSIVE)
	 t0 = widget_button(effectBase, VALUE='Use Allen', UVALUE='USEALLEN')
	 if (state.useAllen eq 1) then begin
	 	widget_control, t0, /SET_BUTTON
	 endif else begin
	 	widget_control, t0, SET_BUTTON=0
	 endelse


	 t0 = widget_base(multiBase, /COLUMN)
	 I0lb = widget_label(t0, VALUE='I0:')
	 state.i0_allen = widget_text(t0, VALUE=strtrim(string(state.stokes0(0),FORMAT='(E10.4)'),2),UVALUE='I0',/EDITABLE,XSIZE=10,YSIZE=1)

	 

	 ; Q0lb = widget_label(multiBase, VALUE='Q0:')
	 ; Q0 = widget_text(multiBase, VALUE=strtrim(string(state.stokes0(1),FORMAT='(E10.4)'),2),UVALUE='Q0',/EDITABLE,XSIZE=10,YSIZE=1)
	 ; U0lb = widget_label(multiBase, VALUE='U0:')
	 ; U0 = widget_text(multiBase, VALUE=strtrim(string(state.stokes0(2),FORMAT='(E10.4)'),2),UVALUE='U0',/EDITABLE,XSIZE=10,YSIZE=1)
	 ; V0lb = widget_label(multiBase, VALUE='V0:')
	 ; V0 = widget_text(multiBase, VALUE=strtrim(string(state.stokes0(3),FORMAT='(E10.4)'),2),UVALUE='V0',/EDITABLE,XSIZE=10,YSIZE=1)

	 
	 buttonsBase = widget_base(multiBase, /ROW)
	 observ_include = widget_button(buttonsBase,VALUE='Load Observation',UVALUE='LOAD_OBSERVATION')
	 reset_observ = widget_button(buttonsBase,VALUE='Reset Observation',UVALUE='RESET_OBSERVATION')
	 synthButton = widget_button(buttonsBase,VALUE='Calculate',UVALUE='Calculate')
	 	 	 	 
		  	 	 
	 formalBase = widget_base(synthBase, /ROW, FRAME=1)
	 
	 base = widget_base(formalBase, /COLUMN)
	 DopplerSlider = cw_fslider(base,TITLE='Doppler v1 [km/s]',UVALUE='DopplerSlider',$
	 	  XSIZE=sliderSize,MAXIMUM=50,MINIMUM=0.1,VALUE=state.Doppler, EDIT=1)
	 DopplerSlider = cw_fslider(base,TITLE='Doppler v2 [km/s]',UVALUE='DopplerSlider2',$
	 	  XSIZE=sliderSize,MAXIMUM=25,MINIMUM=0.1,VALUE=state.Doppler2, EDIT=1)

	 base = widget_base(formalBase, /COLUMN)
	 dtauredlb = cw_fslider(base,TITLE='Dtau1',UVALUE='DTAU1',$
	 	  XSIZE=sliderSize2,MINIMUM=0.d0,MAXIMUM=10.0,VALUE=state.dtau_desired, FORMAT='(F8.2)', EDIT=1)

	 dtauredlb = cw_fslider(base,TITLE='Dtau2',UVALUE='DTAU2',$
	 	  XSIZE=sliderSize2,MINIMUM=0.d0,MAXIMUM=10.0,VALUE=state.dtau_desired2, FORMAT='(F8.2)', EDIT=1)

	 base = widget_base(formalBase, /COLUMN)
	 dtauredlb = cw_fslider(base,TITLE='v1 [km/s]',UVALUE='VEL1',$
	 	  XSIZE=sliderSize2,MINIMUM=-20.d0,MAXIMUM=20.0,VALUE=state.vel, FORMAT='(F8.2)', EDIT=1)

	 dtauredlb = cw_fslider(base,TITLE='v2 [km/s]',UVALUE='VEL2',$
	 	  XSIZE=sliderSize2,MINIMUM=-20.d0,MAXIMUM=20.0,VALUE=state.vel2, FORMAT='(F8.2)', EDIT=1)

	 base = widget_base(formalBase, /COLUMN)
	 dtauredlb = cw_fslider(base,TITLE='a',UVALUE='DAMPING',$
	 	  XSIZE=sliderSize2,MINIMUM=0.d0,MAXIMUM=2.0,VALUE=state.damping, FORMAT='(F8.2)', EDIT=1)

	 dtauredlb = cw_fslider(base,TITLE='ff',UVALUE='FF',$
	 	  XSIZE=sliderSize2,MINIMUM=0.d0,MAXIMUM=1.0,VALUE=state.ff, FORMAT='(F8.2)', EDIT=1)

	 base = widget_base(formalBase, /COLUMN)
	 dtauredlb = cw_fslider(base,TITLE='beta1',UVALUE='SE1',$
	 	  XSIZE=sliderSize2,MINIMUM=1.d0,MAXIMUM=10.0,VALUE=state.beta, FORMAT='(F8.2)', EDIT=1)

	 dtauredlb = cw_fslider(base,TITLE='beta2',UVALUE='SE2',$
	 	  XSIZE=sliderSize2,MINIMUM=1.d0,MAXIMUM=10.0,VALUE=state.beta2, FORMAT='(F8.2)', EDIT=1)

	 ; dtauredlb = widget_label(formalBase, VALUE='Dtau1:')
	 ; dtaured_wid = widget_text(formalBase, VALUE=strtrim(string(state.dtau_desired,FORMAT='(F6.3)'),2),UVALUE='DTAU1',/EDITABLE,XSIZE=8,YSIZE=1)
	 ; dtauredlb = widget_label(formalBase, VALUE='Dtau2:')
	 ; dtaured_wid = widget_text(formalBase, VALUE=strtrim(string(state.dtau_desired2,FORMAT='(F6.3)'),2),UVALUE='DTAU2',/EDITABLE,XSIZE=8,YSIZE=1)
	 ; dbetalb = widget_label(formalBase, VALUE='v1 [km/s]:')
	 ; beta_wid = widget_text(formalBase, VALUE=strtrim(string(state.vel,FORMAT='(F6.3)'),2),UVALUE='VEL1',/EDITABLE,XSIZE=8,YSIZE=1)
	 ; dbetalb = widget_label(formalBase, VALUE='v2 [km/s]:')
	 ; beta_wid = widget_text(formalBase, VALUE=strtrim(string(state.vel2,FORMAT='(F6.3)'),2),UVALUE='VEL2',/EDITABLE,XSIZE=8,YSIZE=1)
	 ; ddamplb = widget_label(formalBase, VALUE='a:')
	 ; damp_wid = widget_text(formalBase, VALUE=strtrim(string(state.damping,FORMAT='(F6.3)'),2),UVALUE='DAMPING',/EDITABLE,XSIZE=8,YSIZE=1)
	 ; ddamplb = widget_label(formalBase, VALUE='ff(1):')
	 ; ff_wid = widget_text(formalBase, VALUE=strtrim(string(state.ff,FORMAT='(F6.3)'),2),UVALUE='FF',/EDITABLE,XSIZE=8,YSIZE=1)
	 ; ddamplb = widget_label(formalBase, VALUE='Se(1):')
	 ; ff_wid = widget_text(formalBase, VALUE=strtrim(string(state.beta,FORMAT='(F6.3)'),2),UVALUE='SE1',/EDITABLE,XSIZE=8,YSIZE=1)
	 ; ddamplb = widget_label(formalBase, VALUE='Se(2):')
	 ; ff_wid = widget_text(formalBase, VALUE=strtrim(string(state.beta2,FORMAT='(F6.3)'),2),UVALUE='SE2',/EDITABLE,XSIZE=8,YSIZE=1)
	 	 	 	 
	 buttonBase = widget_base(hanleBase, /ROW, /ALIGN_LEFT)


; Paschen-Back effect
	 ; t2 = widget_base(buttonBase, /COLUMN)
	 ; lab2 = widget_label(t2, VALUE='ZEEMAN')
	 ; paschenBase = widget_base(t2, /COLUMN, /EXCLUSIVE)
	 ; paschenyesButton = widget_button(paschenBase, VALUE='Paschen-Back', UVALUE='PASCHEN')
	 ; paschennoButton = widget_button(paschenBase, VALUE='Linear Zeeman', UVALUE='LINEAR')
	 ; widget_control, (state.paschen) ? paschenyesButton : paschennoButton, /SET_BUTTON
	 
; Only Zeeman effect or all
	 t3 = widget_base(synth_globBase, /COLUMN)
	 lab3 = widget_label(t3, VALUE='Atom. pol.')
	 effectBase = widget_base(t3, /COLUMN, /EXCLUSIVE)
	 scat_zeemButton = widget_button(effectBase, VALUE='Yes', UVALUE='ALL')
	 zeemButton = widget_button(effectBase, VALUE='No', UVALUE='ZEEMAN')
	 zeemButton = widget_button(effectBase, VALUE='No Anisotropy', UVALUE='NOANISOT')
	 widget_control, (state.effects) ? zeemButton : scat_zeemButton, /SET_BUTTON
	 
; Include stimulated emission in RT
	 ; t5 = widget_base(buttonBase, /COLUMN)
	 ; lab3 = widget_label(t5, VALUE='Stim. emis. in RT')
	 ; stimemiBase = widget_base(t5, /COLUMN, /EXCLUSIVE)
	 ; yesButton = widget_button(stimemiBase, VALUE='Yes', UVALUE='STIM')
	 ; NOButton = widget_button(stimemiBase, VALUE='No', UVALUE='NOSTIM')
	 ; widget_control, (state.stimulated) ? yesButton : noButton, /SET_BUTTON
	 
; Include magneto-optical effects
	 ; t6 = widget_base(buttonBase, /COLUMN)
	 ; lab3 = widget_label(t6, VALUE='Magneto-optical')
	 ; magnetoBase = widget_base(t6, /COLUMN, /EXCLUSIVE)
	 ; yesButton = widget_button(magnetoBase, VALUE='Yes', UVALUE='MAGNETOOPT')
	 ; noButton = widget_button(magnetoBase, VALUE='No', UVALUE='NOMAGNETOOPT')
	 ; widget_control, (state.magneto_opt) ? yesButton : noButton, /SET_BUTTON
	 
; Normalization
	 t7 = widget_base(synth_globBase, /COLUMN)
	 lab3 = widget_label(t7, VALUE='Normalization')
	 absorBase = widget_base(t7, /COLUMN, /EXCLUSIVE)
	 maxButton = widget_button(absorBase, VALUE='Maximum', UVALUE='MAXIMUM')
	 absButton = widget_button(absorBase, VALUE='Absorption', UVALUE='ABSORPTION')
	 contButton = widget_button(absorBase, VALUE='Continuum', UVALUE='CONTINUUM')
	 widget_control, maxButton, /SET_BUTTON
	 	 	 	 
; Rhos plot
	 ; t7 = widget_base(buttonBase, /COLUMN)
	 ; lab3 = widget_label(t7, VALUE='Plot rho^K_Q')
	 ; absorBase = widget_base(t7, /COLUMN, /EXCLUSIVE)
	 ; cohButton = widget_button(absorBase, VALUE='Non-diagonal', UVALUE='NONDIAGONAL')
	 ; diagButton = widget_button(absorBase, VALUE='Diagonal', UVALUE='DIAGONAL')
	 ; widget_control, (state.which_rho_plot) ? diagButton : cohButton, /SET_BUTTON
	 
; Reference frame
	 ; t7 = widget_base(buttonBase, /COLUMN)
	 ; lab3 = widget_label(t7, VALUE='Ref. frame rho^K_Q')
	 ; absorBase = widget_base(t7, /COLUMN, /EXCLUSIVE)
	 ; vertButton = widget_button(absorBase, VALUE='Vertical', UVALUE='VERTICAL_REFFRAME')
	 ; magnButton = widget_button(absorBase, VALUE='Magnetic', UVALUE='MAGNETIC_REFFRAME')
	 ; widget_control, (state.which_refframe) ? magnButton : vertButton, /SET_BUTTON

; Number of slabs
	 t7 = widget_base(synth_globBase, /COLUMN)
	 lab3 = widget_label(t7, VALUE='Number of slabs')
	 absorBase = widget_base(t7, /COLUMN, /EXCLUSIVE)
	 oneButton = widget_button(absorBase, VALUE='One', UVALUE='ONE_SLAB')
	 twoButton = widget_button(absorBase, VALUE='1+1 (same B)', UVALUE='TWO_SLABS')
	 twodiffButton = widget_button(absorBase, VALUE='1+1 (diff. B)', UVALUE='TWO_SLABS_DIFFIELD')
	 twocompButton = widget_button(absorBase, VALUE='Two (diff. B)', UVALUE='TWO_SLABS_COMPO')
	 case(state.number_slabs) of
	 	1: widget_control, oneButton, /SET_BUTTON
	 	2: widget_control, twoButton, /SET_BUTTON
	 	3: widget_control, twodiffButton, /SET_BUTTON
	 	-2: widget_control, twocompButton, /SET_BUTTON
	 endcase
	 
; 	 buttonsBase = widget_base(synth_globBase, /COLUMN)
; 	 observ_include = widget_button(buttonsBase,VALUE='Load Observation',UVALUE='LOAD_OBSERVATION')
; 	 reset_observ = widget_button(buttonsBase,VALUE='Reset Observation',UVALUE='RESET_OBSERVATION')
	 
; 	 bfieldvarBase = widget_base(buttonsBase,/ROW)
; 	 bfield_var = widget_button(bfieldvarBase,$
; 	 	VALUE='rho(B)',UVALUE='FIELD_VARIATION')
; 	 dnfields = widget_label(bfieldvarBase, VALUE='N:')
; 	 nfields = widget_text(bfieldvarBase,$
; 	 	VALUE=strtrim(string(state.bfield_var[2],FORMAT='(I3)'),2),$
; 	 	UVALUE='NFIELDS_RHO',/EDITABLE,XSIZE=3,YSIZE=1)
; 	 dBmin = widget_label(bfieldvarBase, VALUE='Bmin:')
; 	 bmin = widget_text(bfieldvarBase,$
; 	 	VALUE=strtrim(string(state.bfield_var[0]),2),UVALUE='BMIN_RHO',$
; 	 	/EDITABLE,XSIZE=5,YSIZE=1)
; 	 dBmax = widget_label(bfieldvarBase, VALUE='Bmax:')
; 	 bmax = widget_text(bfieldvarBase,$
; 	 	VALUE=strtrim(string(state.bfield_var[1]),2),UVALUE='BMAX_RHO',$
; 	 	/EDITABLE,XSIZE=5,YSIZE=1)
	 	
; 	 synthButton = widget_button(buttonsBase,VALUE='Calculate',UVALUE='Calculate')
	 
	 widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
	 
	 return, state
end