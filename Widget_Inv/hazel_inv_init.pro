function inv_init, reset_state=reset_state
	
   if (keyword_set(reset_state) or file_test('state.idl') eq 0) then begin
		info = {baseWidget: 0L, output_file_widget: 0L, obs_file_widget: 0L, obsplotWidget: 0L, $
			multipletWidget: 0L, plotWidget: 0L,$
			zeeman_pb: 1, atom_pol: 1, stimulated: 1, mag_opt: 1, rt_mode: 3, boundary: fltarr(4),$
			ncycles: 4, inv_bfield: [0,0,1,1], inv_thetaB: [0,0,1,1], inv_chiB: [0,0,1,1],$
			inv_vdopp: [1,1,0,0], inv_tau: [1,1,0,0], inv_tau2: [1,1,0,0], inv_depol: [0,0,0,0], $
			inv_vmacro: [1,1,0,0], inv_vmacro2: [1,1,0,0],$
			inv_damping: [1,1,0,0], inv_source: [0,0,0,0], inv_height: [0,0,0,0],$
			inv_Bfield2: [0,0,1,1], inv_thetaB2: [0,0,1,1], inv_chiB2: [0,0,1,1], inv_vdopp2: [1,1,0,0],$
			inv_ff: [0,0,1,1],$
			Bfield: 1200.d0, thetaB: 0.d0, chiB: 0.d0, vdopp: 0.d0, tau: 0.d0, depol: 0.d0, vmacro: 0.d0,$
			damping: 0.d0, source: 0.d0, height: 0.d0, tau2: 0.d0, vmacro2: 0.d0,$
			Bfield2: 0.d0, thetaB2: 0.d0, chiB2: 0.d0, vdopp2: 0.d0, ff: 0.d0, $
			inversion_mode: [2,1,2,1], $
			stI_weight: [1.d0,1.d0,1.d0,1.d0], stQ_weight: [0.d0,0.d0,1.d0,1.d0], $
			stU_weight: [0.d0,0.d0,1.d0,1.d0], stV_weight: [0.d0,0.d0,1.d0,1.d0], $
			theta_obs: 90.d0, chi_obs: 90.d0, gamm_obs: 0.d0,$
			dir_output_file: './direct.location',$
			dir_stopping: 1, dir_feval: -1, dir_redvol: 0.001, $
			dir_range_bfield: [0.d0,1000.d0], dir_range_thetaB: [0.d0, 180.d0],$
			dir_range_chiB: [0.d0,180.d0], dir_range_vdopp: [0.d0, 20.d0],$
			dir_range_tau: [0.d0, 1.d0], dir_range_depol: [0.d0, 18.d0], $
			dir_range_vmacro: [-10.d0,10.d0], dir_range_damping: [0.d0, 4.d0], $
			dir_range_beta: [0.d0, 1.d0], dir_range_height: [0.d0, 4.d0],$
			dir_range_tau2: [0.d0, 1.d0], dir_range_vmacro2: [-10.d0, 10.d0],$
			dir_range_bfield2: [0.d0,1000.d0], dir_range_thetaB2: [0.d0, 180.d0],$
			dir_range_chiB2: [0.d0,180.d0], dir_range_vdopp2: [0.d0, 20.d0],$
			dir_range_ff: [0.d0, 1.d0], verbose: 0, linear_solver: 0,$
			output_file: 'out.dat', obs_file: '', multiplet: 1, atom: 0, itermax: 20,$
			Btext: 0L, thetaBtext: 0L, chiBtext: 0L, Dopplertext: 0L, Tautext: 0L, Depoltext: 0L, $
   		Macrotext: 0L, Dampingtext: 0L, Betatext: 0L, Heighttext: 0L, label_multiplet: 0L, $
			outputText: 0L, path_obs: '.', number_slabs: 1, dtau_desired2: 0.d0, vel: 0.d0, vel2: 0.d0,$
			Tautext2: 0L, Macrotext2: 0L, i0_allen: 0L, Btext2: 0L, thetaBtext2: 0L, chiBtext2: 0L,$
			Dopplertext2: 0L, buttonsB2: 0L, buttonsthetaB2: 0L, buttonschiB2: 0L,$
			buttonsDoppler2: 0L, buttonsMacro2: 0L, buttonstau2: 0L, j10: '0.0',$
			fftext: 0L, buttonsff: 0L}
	endif else begin	
		restore,'state.idl'
	endelse

	scrsz=get_screen_size()
	xs=scrsz[0]
	ys=scrsz[1]

	if (xs ge 750 and ys ge 950) then begin
		info.baseWidget = widget_base(TITLE='Hazel Inversion', x_scroll_size=750, y_scroll_size=950,/scroll)
	endif else begin
		if xs lt 750 then xwind = xs else xwind = 750
		if ys lt 950 then ywind = ys else ywind = 950
		info.baseWidget = widget_base(TITLE='Hazel Inversion', x_scroll_size=xwind,y_scroll_size=ywind*.9,/scroll)
	endelse
		
;    info.baseWidget = widget_base(TITLE='Hazel Inversion')

   invBase = widget_base(info.baseWidget, /COLUMN, /BASE_ALIGN_CENTER)
   tabs = widget_tab(info.baseWidget, VALUE='hola', UVALUE='TABS')

   tab1 = widget_base(tabs, TITLE='Run', /COLUMN, UVALUE='TAB1')
   tab2 = widget_base(tabs, TITLE='Observations', /COLUMN, UVALUE='TAB2')
   tab3 = widget_base(tabs, TITLE='DIRECT', /COLUMN, UVALUE='TAB3')
   tab4 = widget_base(tabs, TITLE='Parameters', /COLUMN, UVALUE='TAB4')

;-----------------------------
; TAB 1 - RUN
;-----------------------------
   horizBase = widget_base(tab1, /ROW)
   info.plotWidget = widget_draw(horizBase, XSIZE=650, YSIZE=450, /FRAME)

   t1 = widget_base(tab1, /ROW)

; Verbose mode
   baseVerb = widget_base(t1, /COLUMN, FRAME=1)
   label = widget_label(baseVerb, VALUE='Verbose')
   tempBase = widget_base(baseVerb, /COLUMN, /EXCLUSIVE)
   yesButton = widget_button(tempBase, VALUE='On', UVALUE='VERBOSE')
   noButton = widget_button(tempBase, VALUE='Off', UVALUE='NOVERBOSE')
	if (info.verbose eq 0) then begin
		widget_control, noButton, /SET_BUTTON
	endif else begin
		widget_control, yesButton, /SET_BUTTON
	endelse

; Linear solver
   baseLinear = widget_base(t1, /COLUMN, FRAME=1)
   label = widget_label(baseLinear, VALUE='Linear solver')
   tempBase = widget_base(baseLinear, /COLUMN, /EXCLUSIVE)
   yesButton = widget_button(tempBase, VALUE='LU decomposition', UVALUE='LU_SOLVER')
   noButton = widget_button(tempBase, VALUE='Conjugate gradient', UVALUE='CG_SOLVER')   
	if (info.linear_solver eq 0) then begin
		widget_control, yesButton, /SET_BUTTON
	endif else begin
		widget_control, noButton, /SET_BUTTON
	endelse
	
	baseAtom = widget_base(t1, /ROW, FRAME=1)
; Atom
	buttons = cw_bgroup(baseAtom, ['He I','S I'], /COLUMN, /EXCLUSIVE, LABEL_TOP='Atom', SET_VALUE=info.atom, UVALUE='ATOM')
	
	baseAtom2 = widget_base(baseAtom, /COLUMN)
; Multiplet	
	case (info.atom) of
; He I
		0: info.MultipletWidget = widget_slider(baseAtom2,TITLE='Multiplet',UVALUE='MULTIPLET',XSIZE=120,MAXIMUM=4,MINIMUM=1,VALUE=info.multiplet)
; S I
		1: info.MultipletWidget = widget_slider(baseAtom2,TITLE='Multiplet',UVALUE='MULTIPLET',XSIZE=120,MAXIMUM=3,MINIMUM=1,VALUE=info.multiplet)
	endcase
	
	wav = wavelength_component(info.atom,info.multiplet-1)
	info.label_multiplet = widget_label(baseAtom2, VALUE='lambda: '+strtrim(string(wav),2))

; J10
	baseJ10 = widget_base(t1, /ROW, FRAME=1)
	j10Label = widget_label(baseJ10, VALUE='J10/J00 10830 A : ')
   j10widget = widget_text(baseJ10, VALUE=info.j10, /EDITABLE, XSIZE=20,YSIZE=1, UVALUE='J10_VALUE')

; Output file
   outBase = widget_base(tab1, /ROW)
   outLabel = widget_label(outBase, VALUE='Output file : ')
   info.output_file_widget = widget_text(outBase, VALUE=info.output_file, /EDITABLE, XSIZE=70,YSIZE=1, UVALUE='OUTPUT_FILE')
   outButton = widget_button(outBase, VALUE='Select file', UVALUE='OUTPUT_FILE_DIALOG')
; Run inversion button
   runBase = widget_base(tab1, /ROW)
   runButton = widget_button(runBase, VALUE='Run inversion', XSIZE=100,YSIZE=100, UVALUE='RUN_INV')
	runButton = widget_button(runBase, VALUE='Run synthesis', XSIZE=100,YSIZE=100, UVALUE='RUN_SYNTH')
	
	info.outputText = widget_text(tab1, VALUE='', XSIZE=120, YSIZE=10, UVALUE='OUTPUT_TEXT')
;-----------------------------
; TAB 2 - OBSERVATIONS
;-----------------------------

; Plot region
   horizBase = widget_base(tab2, /ROW)
   info.obsplotWidget = widget_draw(horizBase, XSIZE=650, YSIZE=450, /FRAME)

; Observation
   filesBase = widget_base(tab2, /COL)
   obsBase = widget_base(filesBase, /ROW)
   obsLabel = widget_label(obsBase, VALUE='File with observations : ')
   info.obs_file_widget = widget_text(obsBase, VALUE=info.obs_file, /EDITABLE, XSIZE=70,YSIZE=1, UVALUE='OBS_FILE')
   obsButton = widget_button(obsBase, VALUE='Load observation',  UVALUE='OBS_FILE_DIALOG')
   plotBase = widget_base(filesBase, /ROW)
   plotObsButton = widget_button(plotBase, VALUE='Plot observation', UVALUE='PLOT_OBSERVATION')

;-----------------------------
; TAB 3 - DIRECT
;-----------------------------
   t1 = widget_base(tab3, /COLUMN)
   outputFile = cw_field(t1, VALUE=info.dir_output_file, TITLE='Output file for DIRECT :',/RETURN_EVENTS,UVALUE='DIR_OUTPUTFILE')
   stopBase = widget_base(t1, /ROW, FRAME=1)
   stopButton = cw_bgroup(stopBase, ['Number of function evaluations ','Reduction in volume'], /COL, $
                          LABEL_LEFT='Stopping criterium :',/EXCLUSIVE, IDS=stopButtonIDS)
	if (info.dir_stopping eq 0) then begin
   	widget_control, stopButtonIDS(0), /SET_BUTTON
	endif else begin
		widget_control, stopButtonIDS(1), /SET_BUTTON
	endelse
	
   stop2Base = widget_base(stopBase, /COLUMN)
   funcEval = cw_field(stop2Base, VALUE=strtrim(string(info.dir_feval,FORMAT='(I6)'),2), TITLE='Function evaluations :', XSIZE=8, $
		/RETURN_EVENTS, UVALUE='DIR_FUNCEVAL')
   funcRedu = cw_field(stop2Base, VALUE=strtrim(string(info.dir_redvol,FORMAT='(F7.5)'),2), TITLE='Reduction in volume :', XSIZE=8, $
		/RETURN_EVENTS, UVALUE='DIR_REDVOL')

   t1globglob = widget_base(t1, /ROW, FRAME=1)
   
   t1glob = widget_base(t1globglob, /COLUMN, FRAME=1)
   label = widget_label(t1glob, VALUE='Ranges of variation')
   textsize = 100
; Bfield
   t2 = widget_base(t1glob, /ROW)
   label = widget_label(t2, VALUE='B [G]',XSIZE=textsize)
   label = widget_label(t2, VALUE='min: ')
   Bmin = widget_text(t2, VALUE=strtrim(string(info.dir_range_bfield(0),FORMAT='(F7.2)'),2),UVALUE='B_DIR_MIN',/EDITABLE,XSIZE=8,YSIZE=1)
   label = widget_label(t2, VALUE='    max: ')
   Bmax = widget_text(t2, VALUE=strtrim(string(info.dir_range_bfield(1),FORMAT='(F7.2)'),2),UVALUE='B_DIR_MAX',/EDITABLE,XSIZE=8,YSIZE=1)

; thetaB
   t2 = widget_base(t1glob, /ROW)
   label = widget_label(t2, VALUE='thetaB [deg]',XSIZE=textsize)
   label = widget_label(t2, VALUE='min: ')
   thetamin = widget_text(t2, VALUE=strtrim(string(info.dir_range_thetaB(0),FORMAT='(F7.2)'),2),UVALUE='THETAB_DIR_MIN',/EDITABLE,XSIZE=8,YSIZE=1)
   label = widget_label(t2, VALUE='    max: ')
   thetamax = widget_text(t2, VALUE=strtrim(string(info.dir_range_thetaB(1),FORMAT='(F7.2)'),2),UVALUE='THETAB_DIR_MAX',/EDITABLE,XSIZE=8,YSIZE=1)

; chiB
   t2 = widget_base(t1glob, /ROW)
   label = widget_label(t2, VALUE='chiB [deg]',XSIZE=textsize)
   label = widget_label(t2, VALUE='min: ')
   chimin = widget_text(t2, VALUE=strtrim(string(info.dir_range_chiB(0),FORMAT='(F7.2)'),2),UVALUE='CHIB_DIR_MIN',/EDITABLE,XSIZE=8,YSIZE=1)
   label = widget_label(t2, VALUE='    max: ')
   chimax = widget_text(t2, VALUE=strtrim(string(info.dir_range_chiB(1),FORMAT='(F7.2)'),2),UVALUE='CHIB_DIR_MAX',/EDITABLE,XSIZE=8,YSIZE=1)

; vdopp
   t2 = widget_base(t1glob, /ROW)
   label = widget_label(t2, VALUE='vdop [km/s]',XSIZE=textsize)
   label = widget_label(t2, VALUE='min: ')
   vdopmin = widget_text(t2, VALUE=strtrim(string(info.dir_range_vdopp(0),FORMAT='(F7.2)'),2),UVALUE='VDOPP_DIR_MIN',/EDITABLE,XSIZE=8,YSIZE=1)
   label = widget_label(t2, VALUE='    max: ')
   vdopmax = widget_text(t2, VALUE=strtrim(string(info.dir_range_vdopp(1),FORMAT='(F7.2)'),2),UVALUE='VDOPP_DIR_MAX',/EDITABLE,XSIZE=8,YSIZE=1)

; dtau
   t2 = widget_base(t1glob, /ROW)
   label = widget_label(t2, VALUE='Dtau ',XSIZE=textsize)
   label = widget_label(t2, VALUE='min: ')
   dtaumin = widget_text(t2, VALUE=strtrim(string(info.dir_range_tau(0),FORMAT='(F7.2)'),2),UVALUE='DTAU_DIR_MIN',/EDITABLE,XSIZE=8,YSIZE=1)
   label = widget_label(t2, VALUE='    max: ')
   dtaumax = widget_text(t2, VALUE=strtrim(string(info.dir_range_tau(1),FORMAT='(F7.2)'),2),UVALUE='DTAU_DIR_MAX',/EDITABLE,XSIZE=8,YSIZE=1)

; delta
   t2 = widget_base(t1glob, /ROW)
   label = widget_label(t2, VALUE='D^(2) ',XSIZE=textsize)
   label = widget_label(t2, VALUE='min: ')
   d2min = widget_text(t2, VALUE=strtrim(string(info.dir_range_depol(0),FORMAT='(F7.2)'),2),UVALUE='D2_DIR_MIN',/EDITABLE,XSIZE=8,YSIZE=1)
   label = widget_label(t2, VALUE='    max: ')
   d2max = widget_text(t2, VALUE=strtrim(string(info.dir_range_depol(1),FORMAT='(F7.2)'),2),UVALUE='D2_DIR_MAX',/EDITABLE,XSIZE=8,YSIZE=1)

; vmacro
   t2 = widget_base(t1glob, /ROW)
   label = widget_label(t2, VALUE='vmac [km/s] ',XSIZE=textsize)
   label = widget_label(t2, VALUE='min: ')
   vmacmin = widget_text(t2, VALUE=strtrim(string(info.dir_range_vmacro(0),FORMAT='(F7.2)'),2),UVALUE='VMAC_DIR_MIN',/EDITABLE,XSIZE=8,YSIZE=1)
   label = widget_label(t2, VALUE='    max: ')
   vmacmax = widget_text(t2, VALUE=strtrim(string(info.dir_range_vmacro(1),FORMAT='(F7.2)'),2),UVALUE='VMAC_DIR_MAX',/EDITABLE,XSIZE=8,YSIZE=1)

; damping
   t2 = widget_base(t1glob, /ROW)
   label = widget_label(t2, VALUE='damping ',XSIZE=textsize)
   label = widget_label(t2, VALUE='min: ')
   dampmin = widget_text(t2, VALUE=strtrim(string(info.dir_range_damping(0),FORMAT='(F7.2)'),2),UVALUE='DAMP_DIR_MIN',/EDITABLE,XSIZE=8,YSIZE=1)
   label = widget_label(t2, VALUE='    max: ')
   dampmax = widget_text(t2, VALUE=strtrim(string(info.dir_range_damping(1),FORMAT='(F7.2)'),2),UVALUE='DAMP_DIR_MAX',/EDITABLE,XSIZE=8,YSIZE=1)

; beta
   t2 = widget_base(t1glob, /ROW)
   label = widget_label(t2, VALUE='beta ',XSIZE=textsize)
   label = widget_label(t2, VALUE='min: ')
   betamin = widget_text(t2, VALUE=strtrim(string(info.dir_range_beta(0),FORMAT='(F7.2)'),2),UVALUE='BETA_DIR_MIN',/EDITABLE,XSIZE=8,YSIZE=1)
   label = widget_label(t2, VALUE='    max: ')
   betamax = widget_text(t2, VALUE=strtrim(string(info.dir_range_beta(1),FORMAT='(F7.2)'),2),UVALUE='BETA_DIR_MAX',/EDITABLE,XSIZE=8,YSIZE=1)

; Height
   t2 = widget_base(t1glob, /ROW)
   label = widget_label(t2, VALUE='Height ["] ',XSIZE=textsize)
   label = widget_label(t2, VALUE='min: ')
   hmin = widget_text(t2, VALUE=strtrim(string(info.dir_range_height(0),FORMAT='(F7.2)'),2),UVALUE='H_DIR_MIN',/EDITABLE,XSIZE=8,YSIZE=1)
   label = widget_label(t2, VALUE='    max: ')
   hmax = widget_text(t2, VALUE=strtrim(string(info.dir_range_height(1),FORMAT='(F7.2)'),2),UVALUE='H_DIR_MAX',/EDITABLE,XSIZE=8,YSIZE=1)


t1glob = widget_base(t1globglob, /COLUMN, FRAME=1)
   label = widget_label(t1glob, VALUE='Ranges of variation component 2')
   textsize = 100
; Bfield2
   t2 = widget_base(t1glob, /ROW)
   label = widget_label(t2, VALUE='B [G]',XSIZE=textsize)
   label = widget_label(t2, VALUE='min: ')
   Bmin = widget_text(t2, VALUE=strtrim(string(info.dir_range_bfield2(0),FORMAT='(F7.2)'),2),UVALUE='B2_DIR_MIN',$
   	/EDITABLE,XSIZE=8,YSIZE=1)
   label = widget_label(t2, VALUE='    max: ')
   Bmax = widget_text(t2, VALUE=strtrim(string(info.dir_range_bfield2(1),FORMAT='(F7.2)'),2),UVALUE='B2_DIR_MAX',$
   	/EDITABLE,XSIZE=8,YSIZE=1)

; thetaB
   t2 = widget_base(t1glob, /ROW)
   label = widget_label(t2, VALUE='thetaB [deg]',XSIZE=textsize)
   label = widget_label(t2, VALUE='min: ')
   thetamin = widget_text(t2, VALUE=strtrim(string(info.dir_range_thetaB2(0),FORMAT='(F7.2)'),2),UVALUE='THETAB2_DIR_MIN',$
   	/EDITABLE,XSIZE=8,YSIZE=1)
   label = widget_label(t2, VALUE='    max: ')
   thetamax = widget_text(t2, VALUE=strtrim(string(info.dir_range_thetaB2(1),FORMAT='(F7.2)'),2),UVALUE='THETAB2_DIR_MAX',$
   	/EDITABLE,XSIZE=8,YSIZE=1)

; chiB
   t2 = widget_base(t1glob, /ROW)
   label = widget_label(t2, VALUE='chiB [deg]',XSIZE=textsize)
   label = widget_label(t2, VALUE='min: ')
   chimin = widget_text(t2, VALUE=strtrim(string(info.dir_range_chiB2(0),FORMAT='(F7.2)'),2),UVALUE='CHIB2_DIR_MIN',$
   	/EDITABLE,XSIZE=8,YSIZE=1)
   label = widget_label(t2, VALUE='    max: ')
   chimax = widget_text(t2, VALUE=strtrim(string(info.dir_range_chiB2(1),FORMAT='(F7.2)'),2),UVALUE='CHIB2_DIR_MAX',$
   	/EDITABLE,XSIZE=8,YSIZE=1)

; vdopp
   t2 = widget_base(t1glob, /ROW)
   label = widget_label(t2, VALUE='vdop [km/s]',XSIZE=textsize)
   label = widget_label(t2, VALUE='min: ')
   vdopmin = widget_text(t2, VALUE=strtrim(string(info.dir_range_vdopp2(0),FORMAT='(F7.2)'),2),UVALUE='VDOPP2_DIR_MIN',$
   	/EDITABLE,XSIZE=8,YSIZE=1)
   label = widget_label(t2, VALUE='    max: ')
   vdopmax = widget_text(t2, VALUE=strtrim(string(info.dir_range_vdopp2(1),FORMAT='(F7.2)'),2),UVALUE='VDOPP2_DIR_MAX',$
   	/EDITABLE,XSIZE=8,YSIZE=1)

; dtau2
   t2 = widget_base(t1glob, /ROW)
   label = widget_label(t2, VALUE='Dtau2 ',XSIZE=textsize)
   label = widget_label(t2, VALUE='min: ')
   betamin = widget_text(t2, VALUE=strtrim(string(info.dir_range_tau2(0),FORMAT='(F7.2)'),2),UVALUE='DTAU2_DIR_MIN',$
   	/EDITABLE,XSIZE=8,YSIZE=1)
   label = widget_label(t2, VALUE='    max: ')
   betamax = widget_text(t2, VALUE=strtrim(string(info.dir_range_tau2(1),FORMAT='(F7.2)'),2),UVALUE='DTAU2_DIR_MAX',$
   	/EDITABLE,XSIZE=8,YSIZE=1)

; vmac2
   t2 = widget_base(t1glob, /ROW)
   label = widget_label(t2, VALUE=' vmac2 [km/s] ',XSIZE=textsize)
   label = widget_label(t2, VALUE='min: ')
   hmin = widget_text(t2, VALUE=strtrim(string(info.dir_range_vmacro2(0),FORMAT='(F7.2)'),2),UVALUE='VMAC2_DIR_MIN',$
   	/EDITABLE,XSIZE=8,YSIZE=1)
   label = widget_label(t2, VALUE='    max: ')
   hmax = widget_text(t2, VALUE=strtrim(string(info.dir_range_vmacro2(1),FORMAT='(F7.2)'),2),UVALUE='VMAC2_DIR_MAX',$
   	/EDITABLE,XSIZE=8,YSIZE=1)

; ff
   t2 = widget_base(t1glob, /ROW)
   label = widget_label(t2, VALUE=' ff1 ',XSIZE=textsize)
   label = widget_label(t2, VALUE='min: ')
   hmin = widget_text(t2, VALUE=strtrim(string(info.dir_range_ff(0),FORMAT='(F7.2)'),2),UVALUE='FF_DIR_MIN',$
   	/EDITABLE,XSIZE=8,YSIZE=1)
   label = widget_label(t2, VALUE='    max: ')
   hmax = widget_text(t2, VALUE=strtrim(string(info.dir_range_ff(1),FORMAT='(F7.2)'),2),UVALUE='FF_DIR_MAX',$
   	/EDITABLE,XSIZE=8,YSIZE=1)

;-----------------------------
; TAB 4 - PARAMETERS
;-----------------------------
   t1 = widget_base(tab4, /ROW)
; Zeeman or PB
   basePB = widget_base(t1, /COLUMN, FRAME=1, YSIZE=10)
   label = widget_label(basePB, VALUE='ZEEMAN/PB')
   tempBase = widget_base(basePB, /COLUMN, /EXCLUSIVE)
   yesButton = widget_button(tempBase, VALUE='Paschen-Back', UVALUE='PASCHEN')
   noButton = widget_button(tempBase, VALUE='Linear Zeeman', UVALUE='LINEAR')
	if (info.zeeman_pb eq 1) then begin
   	widget_control, yesButton, /SET_BUTTON
	endif else begin
		widget_control, noButton, /SET_BUTTON
	endelse
   
; Atomic polarization
   baseAtomPol = widget_base(t1, /COLUMN, FRAME=1)
   label = widget_label(baseAtomPol, VALUE='Atom. pol.')
   tempBase = widget_base(baseAtomPol, /COLUMN, /EXCLUSIVE)
   yesButton = widget_button(tempBase, VALUE='Yes', UVALUE='ATOMPOL')
   noButton = widget_button(tempBase, VALUE='No', UVALUE='NOATOMPOL')
	if (info.atom_pol eq 1) then begin
   	widget_control, yesButton, /SET_BUTTON
	endif else begin
		widget_control, noButton, /SET_BUTTON
	endelse

; Stimulated emission
   baseStimEmis = widget_base(t1, /COLUMN, FRAME=1)
   label = widget_label(baseStimEmis, VALUE='Stimulated emission')
   tempBase = widget_base(baseStimEmis, /COLUMN, /EXCLUSIVE)
   yesButton = widget_button(tempBase, VALUE='Yes', UVALUE='STIM')
   noButton = widget_button(tempBase, VALUE='No', UVALUE='NOSTIM')
	if (info.stimulated eq 1) then begin
		widget_control, yesButton, /SET_BUTTON
	endif else begin
   	widget_control, noButton, /SET_BUTTON
	endelse

; Magneto-optical
   baseMagOpt = widget_base(t1, /COLUMN, FRAME=1)
   label = widget_label(baseMagOpt, VALUE='Magneto-optical')
   tempBase = widget_base(baseMagOpt, /COLUMN, /EXCLUSIVE)
   yesButton = widget_button(tempBase, VALUE='Yes', UVALUE='MAGOPT')
   noButton = widget_button(tempBase, VALUE='No', UVALUE='NOMAGOPT')
	if (info.mag_opt eq 1) then begin
   	widget_control, yesButton, /SET_BUTTON
	endif else begin
		widget_control, noButton, /SET_BUTTON
	endelse

; Radiative transfer mode
   baseRT = widget_base(t1, /COLUMN, FRAME=1)
   label = widget_label(baseRT, VALUE='RT mode')
   tempBase = widget_base(baseRT, /COLUMN, /EXCLUSIVE)
   emissionButton = widget_button(tempBase, VALUE='Slab (optically thin)', UVALUE='EMISSION')
   formalButton = widget_button(tempBase, VALUE='Slab (no MO)', UVALUE='FORMAL')
   deloparButton = widget_button(tempBase, VALUE='Slab (DELO)', UVALUE='DELOPAR')
   milneButton = widget_button(tempBase, VALUE='Milne-Eddington', UVALUE='MILNE')
   exactslabButton = widget_button(tempBase, VALUE='Slab (exact)', UVALUE='EXACT_SLAB')
   simpleslabButton = widget_button(tempBase, VALUE='Simplified slab', UVALUE='SIMPLE_SLAB')
	case (info.rt_mode) of
		0: widget_control, emissionButton, /SET_BUTTON  ; 0 -> pure emission
		1: widget_control, formalButton, /SET_BUTTON  ; 1 -> slab (neglecting MO terms)
		2: widget_control, milneButton, /SET_BUTTON  ; 2 -> Milne-Eddington
		3: widget_control, deloparButton, /SET_BUTTON  ; 3 -> DELO
		4: widget_control, simpleslabButton, /SET_BUTTON    ; 4 -> simplified slab (optically thin)
		5: widget_control, exactslabButton, /SET_BUTTON   ; 5 -> exact slab
	endcase

; Number of slabs
	 t7 = widget_base(t1, /COLUMN, FRAME=1)
	 lab3 = widget_label(t7, VALUE='Number of slabs')
	 absorBase = widget_base(t7, /COLUMN, /EXCLUSIVE)
	 oneButton = widget_button(absorBase, VALUE='One', UVALUE='ONE_SLAB')
	 twoButton = widget_button(absorBase, VALUE='1+1', UVALUE='TWO_SLABS')
	 twodifButton = widget_button(absorBase, VALUE='1+1 (diff. field)', UVALUE='TWO_SLABS_DIFFIELD')
	 twocompButton = widget_button(absorBase, VALUE='2 (diff. field)', UVALUE='TWO_SLABS_COMPO')

	trightglob = widget_base(t1, /COLUMN)
   t3glob = widget_base(trightglob, /COLUMN, FRAME=1)
   label = widget_label(t3glob, VALUE='Boundary conditions')
   t3 = widget_base(t3glob, /ROW)
   label = widget_label(t3, VALUE='I0:')
   I0 = widget_text(t3, VALUE=strtrim(string(info.boundary(0),FORMAT='(E8.2)')),UVALUE='I0',/EDITABLE,XSIZE=8,YSIZE=1)
   label = widget_label(t3, VALUE='Q0:')
   Q0 = widget_text(t3, VALUE=strtrim(string(info.boundary(1),FORMAT='(E8.2)')),UVALUE='Q0',/EDITABLE,XSIZE=8,YSIZE=1)
   label = widget_label(t3, VALUE='U0:')
   U0 = widget_text(t3, VALUE=strtrim(string(info.boundary(2),FORMAT='(E8.2)')),UVALUE='U0',/EDITABLE,XSIZE=8,YSIZE=1)
   label = widget_label(t3, VALUE='V0:')
   V0 = widget_text(t3, VALUE=strtrim(string(info.boundary(3),FORMAT='(E8.2)')),UVALUE='V0',/EDITABLE,XSIZE=8,YSIZE=1)
   t4 = widget_base(t2, /ROW, FRAME=1)

	t3glob = widget_base(trightglob, /COLUMN, FRAME=1)
   ncycles = cw_bgroup(t3glob, ['1','2','3','4'], /ROW, /EXCLUSIVE, LABEL_LEFT='Number of cycles', UVALUE='NCYCLES', SET_VALUE=info.ncycles-1)
   
;    t2 = widget_base(tab4, /ROW)
;    t3glob = widget_base(t2, /COLUMN, FRAME=1)
;    label = widget_label(t3glob, VALUE='Boundary conditions')
;    t3 = widget_base(t3glob, /ROW)
;    label = widget_label(t3, VALUE='I0:')
;    I0 = widget_text(t3, VALUE=strtrim(string(info.boundary(0),FORMAT='(E8.2)')),UVALUE='I0',/EDITABLE,XSIZE=8,YSIZE=1)
;    label = widget_label(t3, VALUE='Q0:')
;    Q0 = widget_text(t3, VALUE=strtrim(string(info.boundary(1),FORMAT='(E8.2)')),UVALUE='Q0',/EDITABLE,XSIZE=8,YSIZE=1)
;    label = widget_label(t3, VALUE='U0:')
;    U0 = widget_text(t3, VALUE=strtrim(string(info.boundary(2),FORMAT='(E8.2)')),UVALUE='U0',/EDITABLE,XSIZE=8,YSIZE=1)
;    label = widget_label(t3, VALUE='V0:')
;    V0 = widget_text(t3, VALUE=strtrim(string(info.boundary(3),FORMAT='(E8.2)')),UVALUE='V0',/EDITABLE,XSIZE=8,YSIZE=1)
;    t4 = widget_base(t2, /ROW, FRAME=1)
; 	  
;    ncycles = cw_bgroup(t4, ['1','2','3','4'], /ROW, /EXCLUSIVE, LABEL_LEFT='Number of cycles', UVALUE='NCYCLES', SET_VALUE=info.ncycles-1)
	
   tglob = widget_base(tab4, /ROW)
   
   t4glob = widget_base(tglob, /COLUMN)
   t4 = widget_base(t4glob, /COLUMN, FRAME=1)
   label = widget_label(t4, VALUE='Parameters to invert')

   textsize = 100

; Bfield
   tBfield = widget_base(t4, /ROW)
   label = widget_label(tBfield, VALUE='B [G]', XSIZE=textsize)
   info.Btext = widget_text(tBfield, VALUE=strtrim(string(info.Bfield,FORMAT='(F8.2)'),2), UVALUE='BVALUE', /EDITABLE, XSIZE=8)
   buttons = cw_bgroup(tBfield, ['1','2','3','4'], /ROW, /NONEXCLUSIVE, LABEL_LEFT='Cycles', SET_VALUE=info.inv_bfield, UVALUE='CYC_BVALUE')
	ind = where(info.inv_bfield eq 0)
	if (n_elements(ind) eq 4) then begin
		widget_control, info.Btext, sensitive=0
	endif

; Inclination
   tIncli = widget_base(t4, /ROW)
   label = widget_label(tIncli, VALUE='thetaB [deg]', XSIZE=textsize)
   info.thetaBtext = widget_text(tIncli, VALUE=strtrim(string(info.thetaB,FORMAT='(F6.2)'),2), UVALUE='THETAVALUE', /EDITABLE, XSIZE=8, YSIZE=1)
   buttons = cw_bgroup(tIncli, ['1','2','3','4'], /ROW, /NONEXCLUSIVE, LABEL_LEFT='Cycles', SET_VALUE=info.inv_thetaB, UVALUE='CYC_THETAB')   
	ind = where(info.inv_thetaB eq 0)
	if (n_elements(ind) eq 4) then begin
		widget_control, info.thetaBtext, sensitive=0
	endif
	
; Azimuth
   tAzi = widget_base(t4, /ROW)
   label = widget_label(tAzi, VALUE='chiB [deg]', XSIZE=textsize)
   info.chiBtext = widget_text(tAzi, VALUE=strtrim(string(info.chiB,FORMAT='(F6.2)'),2), UVALUE='CHIVALUE', /EDITABLE, XSIZE=8, YSIZE=1)
   buttons = cw_bgroup(tAzi, ['1','2','3','4'], /ROW, /NONEXCLUSIVE, LABEL_LEFT='Cycles', SET_VALUE=info.inv_chiB, UVALUE='CYC_CHIB')
	ind = where(info.inv_chiB eq 0)
	if (n_elements(ind) eq 4) then begin
		widget_control, info.chiBtext, sensitive=0
	endif

; Doppler width
   tDoppler = widget_base(t4, /ROW)
   label = widget_label(tDoppler, VALUE='Doppler [km/s]', XSIZE=textsize)
   info.Dopplertext = widget_text(tDoppler, VALUE=strtrim(string(info.vdopp,FORMAT='(F6.2)'),2), UVALUE='DOPPLERVALUE', /EDITABLE, XSIZE=8, YSIZE=1)
   buttons = cw_bgroup(tDoppler, ['1','2','3','4'], /ROW, /NONEXCLUSIVE, LABEL_LEFT='Cycles', SET_VALUE=info.inv_vdopp, UVALUE='CYC_VDOPP')
	ind = where(info.inv_vdopp eq 0)
	if (n_elements(ind) eq 4) then begin
		widget_control, info.Dopplertext, sensitive=0
	endif

; Optical depth
   tTau = widget_base(t4, /ROW)
   label = widget_label(tTau, VALUE='Dtau', XSIZE=textsize)
   info.Tautext = widget_text(tTau, VALUE=strtrim(string(info.tau,FORMAT='(F6.2)'),2), UVALUE='DTAUVALUE', /EDITABLE, XSIZE=8, YSIZE=1)
   buttons = cw_bgroup(tTau, ['1','2','3','4'], /ROW, /NONEXCLUSIVE, LABEL_LEFT='Cycles', SET_VALUE=info.inv_tau, UVALUE='CYC_TAU')
	ind = where(info.inv_tau eq 0)
	if (n_elements(ind) eq 4) then begin
		widget_control, info.Tautext, sensitive=0
	endif

; Depolarization rate
   tDepol = widget_base(t4, /ROW)
   label = widget_label(tDepol, VALUE='D^(2)', XSIZE=textsize)
   info.Depoltext = widget_text(tDepol, VALUE=strtrim(string(info.depol,FORMAT='(E8.2)'),2), UVALUE='DEPOLVALUE', /EDITABLE, XSIZE=8, YSIZE=1)
   buttons = cw_bgroup(tDepol, ['1','2','3','4'], /ROW, /NONEXCLUSIVE, LABEL_LEFT='Cycles', SET_VALUE=info.inv_depol, UVALUE='CYC_DEPOL')
	ind = where(info.inv_depol eq 0)
	if (n_elements(ind) eq 4) then begin
		widget_control, info.Depoltext, sensitive=0
	endif

; Macroscopic velocity
   tMacro = widget_base(t4, /ROW)
   label = widget_label(tMacro, VALUE='v_mac [km/s]', XSIZE=textsize)
   info.Macrotext = widget_text(tMacro, VALUE=strtrim(string(info.vmacro,FORMAT='(F6.2)'),2), UVALUE='VMACVALUE', /EDITABLE, XSIZE=8, YSIZE=1)
   buttons = cw_bgroup(tMacro, ['1','2','3','4'], /ROW, /NONEXCLUSIVE, LABEL_LEFT='Cycles', SET_VALUE=info.inv_vmacro, UVALUE='CYC_VMACRO')
	ind = where(info.inv_vmacro eq 0)
	if (n_elements(ind) eq 4) then begin
		widget_control, info.Macrotext, sensitive=0
	endif

; Damping
   tDamping = widget_base(t4, /ROW)
   label = widget_label(tDamping, VALUE='Damping', XSIZE=textsize)
   info.Dampingtext = widget_text(tDamping, VALUE=strtrim(string(info.damping,FORMAT='(F6.2)'),2), UVALUE='DAMPVALUE', /EDITABLE, XSIZE=8, YSIZE=1)
   buttons = cw_bgroup(tDamping, ['1','2','3','4'], /ROW, /NONEXCLUSIVE, LABEL_LEFT='Cycles', SET_VALUE=info.inv_damping, UVALUE='CYC_DAMP')
	ind = where(info.inv_damping eq 0)
	if (n_elements(ind) eq 4) then begin
		widget_control, info.Dampingtext, sensitive=0
	endif

; Source function gradient
   tBeta = widget_base(t4, /ROW)
   label = widget_label(tBeta, VALUE='Beta', XSIZE=textsize)
   info.Betatext = widget_text(tBeta, VALUE=strtrim(string(info.source,FORMAT='(F6.2)'),2), UVALUE='BETAVALUE', /EDITABLE, XSIZE=8, YSIZE=1)
   buttons = cw_bgroup(tBeta, ['1','2','3','4'], /ROW, /NONEXCLUSIVE, LABEL_LEFT='Cycles', SET_VALUE=info.inv_source, UVALUE='CYC_BETA')
	ind = where(info.inv_source eq 0)
	if (n_elements(ind) eq 4) then begin
		widget_control, info.Betatext, sensitive=0
	endif

; Height
   tHeight = widget_base(t4, /ROW)
   label = widget_label(tHeight, VALUE='Height', XSIZE=textsize)
   info.Heighttext = widget_text(tHeight, VALUE=strtrim(string(info.height,FORMAT='(F6.2)'),2), UVALUE='HEIGHTVALUE', /EDITABLE, XSIZE=8, YSIZE=1)
   buttons = cw_bgroup(tHeight, ['1','2','3','4'], /ROW, /NONEXCLUSIVE, LABEL_LEFT='Cycles', SET_VALUE=info.inv_height, UVALUE='CYC_HEIGHT')
	ind = where(info.inv_height eq 0)
	if (n_elements(ind) eq 4) then begin
		widget_control, info.Heighttext, sensitive=0
	endif

; COMPONENT 2   
   t4glob = widget_base(tglob, /COLUMN)
   t4 = widget_base(t4glob, /COLUMN, FRAME=1)
   label = widget_label(t4, VALUE='Parameters to invert component 2')

; Bfield
   tBfield = widget_base(t4, /ROW)
   label = widget_label(tBfield, VALUE='B [G]', XSIZE=textsize)
   info.Btext2 = widget_text(tBfield, VALUE=strtrim(string(info.Bfield2,FORMAT='(F8.2)'),2), UVALUE='BVALUE2', /EDITABLE, XSIZE=8)
   info.buttonsB2 = cw_bgroup(tBfield, ['1','2','3','4'], /ROW, /NONEXCLUSIVE, LABEL_LEFT='Cycles', SET_VALUE=info.inv_bfield2, UVALUE='CYC_BVALUE2')
	ind = where(info.inv_bfield2 eq 0)
	if (n_elements(ind) eq 4) then begin
		widget_control, info.Btext2, sensitive=0
	endif

; Inclination
   tIncli = widget_base(t4, /ROW)
   label = widget_label(tIncli, VALUE='thetaB [deg]', XSIZE=textsize)
   info.thetaBtext2 = widget_text(tIncli, VALUE=strtrim(string(info.thetaB2,FORMAT='(F6.2)'),2), UVALUE='THETAVALUE2', /EDITABLE, XSIZE=8, YSIZE=1)
   info.buttonsthetaB2 = cw_bgroup(tIncli, ['1','2','3','4'], /ROW, /NONEXCLUSIVE, LABEL_LEFT='Cycles', SET_VALUE=info.inv_thetaB2, UVALUE='CYC_THETAB2')
	ind = where(info.inv_thetaB2 eq 0)
	if (n_elements(ind) eq 4) then begin
		widget_control, info.thetaBtext2, sensitive=0
	endif
	
; Azimuth
   tAzi = widget_base(t4, /ROW)
   label = widget_label(tAzi, VALUE='chiB [deg]', XSIZE=textsize)
   info.chiBtext2 = widget_text(tAzi, VALUE=strtrim(string(info.chiB2,FORMAT='(F6.2)'),2), UVALUE='CHIVALUE2', /EDITABLE, XSIZE=8, YSIZE=1)
   info.buttonschiB2 = cw_bgroup(tAzi, ['1','2','3','4'], /ROW, /NONEXCLUSIVE, LABEL_LEFT='Cycles', SET_VALUE=info.inv_chiB2, UVALUE='CYC_CHIB2')
	ind = where(info.inv_chiB2 eq 0)
	if (n_elements(ind) eq 4) then begin
		widget_control, info.chiBtext2, sensitive=0
	endif

; Doppler width
   tDoppler = widget_base(t4, /ROW)
   label = widget_label(tDoppler, VALUE='Doppler [km/s]', XSIZE=textsize)
   info.Dopplertext2 = widget_text(tDoppler, VALUE=strtrim(string(info.vdopp2,FORMAT='(F6.2)'),2), UVALUE='DOPPLERVALUE2', /EDITABLE, XSIZE=8, YSIZE=1)
   info.buttonsDoppler2 = cw_bgroup(tDoppler, ['1','2','3','4'], /ROW, /NONEXCLUSIVE, LABEL_LEFT='Cycles', SET_VALUE=info.inv_vdopp2, UVALUE='CYC_VDOPP2')
	ind = where(info.inv_vdopp2 eq 0)
	if (n_elements(ind) eq 4) then begin
		widget_control, info.Dopplertext2, sensitive=0
	endif

; Optical depth2
   tTau = widget_base(t4, /ROW)
   label = widget_label(tTau, VALUE='Dtau2', XSIZE=textsize)
   info.Tautext2 = widget_text(tTau, VALUE=strtrim(string(info.tau2,FORMAT='(F6.2)'),2), UVALUE='DTAU2VALUE', /EDITABLE, XSIZE=8, YSIZE=1)
   info.buttonstau2 = cw_bgroup(tTau, ['1','2','3','4'], /ROW, /NONEXCLUSIVE, LABEL_LEFT='Cycles', SET_VALUE=info.inv_tau2, UVALUE='CYC_TAU2')
	ind = where(info.inv_tau2 eq 0)
	if (n_elements(ind) eq 4) then begin
		widget_control, info.Tautext2, sensitive=0
	endif

; Macroscopic velocity
   tMacro = widget_base(t4, /ROW)
   label = widget_label(tMacro, VALUE='v_mac2 [km/s]', XSIZE=textsize)
   info.Macrotext2 = widget_text(tMacro, VALUE=strtrim(string(info.vmacro2,FORMAT='(F6.2)'),2), UVALUE='VMAC2VALUE', /EDITABLE, XSIZE=8, YSIZE=1)
   info.buttonsMacro2 = cw_bgroup(tMacro, ['1','2','3','4'], /ROW, /NONEXCLUSIVE, LABEL_LEFT='Cycles', SET_VALUE=info.inv_vmacro2, UVALUE='CYC_VMACRO2')
	ind = where(info.inv_vmacro2 eq 0)
	if (n_elements(ind) eq 4) then begin
		widget_control, info.Macrotext2, sensitive=0
	endif

; Filling factor
   tMacro = widget_base(t4, /ROW)
   label = widget_label(tMacro, VALUE='ff1 [km/s]', XSIZE=textsize)
   info.fftext = widget_text(tMacro, VALUE=strtrim(string(info.vmacro2,FORMAT='(F6.2)'),2), UVALUE='FFVALUE', /EDITABLE, XSIZE=8, YSIZE=1)
   info.buttonsff = cw_bgroup(tMacro, ['1','2','3','4'], /ROW, /NONEXCLUSIVE, LABEL_LEFT='Cycles', SET_VALUE=info.inv_ff, UVALUE='CYC_FF')
	ind = where(info.inv_ff eq 0)
	if (n_elements(ind) eq 4) then begin
		widget_control, info.fftext, sensitive=0
	endif

   
; Weights
   trightglob = widget_base(tglob, /COLUMN)
	tnitersmax = widget_base(trightglob, /COLUMN, FRAME=1)
	itermax = cw_field(tnitersmax, VALUE=strtrim(string(info.itermax,FORMAT='(I3)'),2), TITLE='Max. LM iters. per cycle :', $
		XSIZE=8, UVALUE='MAX_ITER', /RETURN_EVENTS)
	
   tright = widget_base(trightglob, /COLUMN, FRAME=1)
   label = widget_label(tright, VALUE='Weigths in each cycle')
   tI = widget_base(tright, /ROW)
   label = widget_label(tI, VALUE='I weights :')
   in1 = widget_text(tI, VALUE=strtrim(string(info.stI_weight(0),FORMAT='(F6.2)'),2), XSIZE=5, /EDITABLE, UVALUE='WI1')
   in2 = widget_text(tI, VALUE=strtrim(string(info.stI_weight(1),FORMAT='(F6.2)'),2), XSIZE=5, /EDITABLE, UVALUE='WI2')
   in3 = widget_text(tI, VALUE=strtrim(string(info.stI_weight(2),FORMAT='(F6.2)'),2), XSIZE=5, /EDITABLE, UVALUE='WI3')
   in4 = widget_text(tI, VALUE=strtrim(string(info.stI_weight(3),FORMAT='(F6.2)'),2), XSIZE=5, /EDITABLE, UVALUE='WI4')
   tQ = widget_base(tright, /ROW)
   label = widget_label(tQ, VALUE='Q weights :')
   Q1 = widget_text(tQ, VALUE=strtrim(string(info.stQ_weight(0),FORMAT='(F6.2)'),2), XSIZE=5, /EDITABLE, UVALUE='WQ1')
   Q2 = widget_text(tQ, VALUE=strtrim(string(info.stQ_weight(1),FORMAT='(F6.2)'),2), XSIZE=5, /EDITABLE, UVALUE='WQ2')
   Q3 = widget_text(tQ, VALUE=strtrim(string(info.stQ_weight(2),FORMAT='(F6.2)'),2), XSIZE=5, /EDITABLE, UVALUE='WQ3')
   Q4 = widget_text(tQ, VALUE=strtrim(string(info.stQ_weight(3),FORMAT='(F6.2)'),2), XSIZE=5, /EDITABLE, UVALUE='WQ4')
   tU = widget_base(tright, /ROW)
   label = widget_label(tU, VALUE='U weights :')
   U1 = widget_text(tU, VALUE=strtrim(string(info.stU_weight(0),FORMAT='(F6.2)'),2), XSIZE=5, /EDITABLE, UVALUE='WU1')
   U2 = widget_text(tU, VALUE=strtrim(string(info.stU_weight(1),FORMAT='(F6.2)'),2), XSIZE=5, /EDITABLE, UVALUE='WU2')
   U3 = widget_text(tU, VALUE=strtrim(string(info.stU_weight(2),FORMAT='(F6.2)'),2), XSIZE=5, /EDITABLE, UVALUE='WU3')
   U4 = widget_text(tU, VALUE=strtrim(string(info.stU_weight(3),FORMAT='(F6.2)'),2), XSIZE=5, /EDITABLE, UVALUE='WU4')
   tV = widget_base(tright, /ROW)
   label = widget_label(tV, VALUE='V weights :')
   V1 = widget_text(tV, VALUE=strtrim(string(info.stV_weight(0),FORMAT='(F6.2)'),2), XSIZE=5, /EDITABLE, UVALUE='WV1')
   V2 = widget_text(tV, VALUE=strtrim(string(info.stV_weight(1),FORMAT='(F6.2)'),2), XSIZE=5, /EDITABLE, UVALUE='WV2')
   V3 = widget_text(tV, VALUE=strtrim(string(info.stV_weight(2),FORMAT='(F6.2)'),2), XSIZE=5, /EDITABLE, UVALUE='WV3')
   V4 = widget_text(tV, VALUE=strtrim(string(info.stV_weight(3),FORMAT='(F6.2)'),2), XSIZE=5, /EDITABLE, UVALUE='WV4')

; Observation angles
   tobsang = widget_base(trightglob, /COLUMN, FRAME=1)
   label = widget_label(tobsang, VALUE='Observation angles')
   theta = cw_field(tobsang, VALUE=strtrim(string(info.theta_obs,FORMAT='(F6.2)'),2), TITLE='theta [deg] :', XSIZE=8, UVALUE='OBS_THETA',$
		/RETURN_EVENTS)
   chi = cw_field(tobsang, VALUE=strtrim(string(info.chi_obs,FORMAT='(F6.2)'),2), TITLE='chi [deg]  :', XSIZE=8, UVALUE='OBS_CHI', /RETURN_EVENTS)
   gamm = cw_field(tobsang, VALUE=strtrim(string(info.gamm_obs,FORMAT='(F6.2)'),2), TITLE='gamma [deg] :', XSIZE=8, UVALUE='OBS_GAMMA',$
		/RETURN_EVENTS)
		
	i0 = i0_allen(wavelength_component(info.atom,info.multiplet-1),cos(info.theta_obs*!DPI/180.d0))
	info.i0_allen = widget_label(tobsang, VALUE='Allen: '+strtrim(string(i0),2))
		
	
; Inversion modes
   tModeLM = widget_base(trightglob, /ROW, FRAME=1)
   tCyc1 = cw_bgroup(tModeLM, ['LM','DIRECT','SIMPLE'], /COLUMN, /EXCLUSIVE, LABEL_TOP='Cycle 1', SET_VALUE=info.inversion_mode(0)-1, UVALUE='INVMODE_CYC1')
   tCyc2 = cw_bgroup(tModeLM, ['LM','DIRECT','SIMPLE'], /COLUMN, /EXCLUSIVE, LABEL_TOP='Cycle 2', SET_VALUE=info.inversion_mode(1)-1, UVALUE='INVMODE_CYC2')
   tCyc3 = cw_bgroup(tModeLM, ['LM','DIRECT','SIMPLE'], /COLUMN, /EXCLUSIVE, LABEL_TOP='Cycle 3', SET_VALUE=info.inversion_mode(2)-1, UVALUE='INVMODE_CYC3')
   tCyc4 = cw_bgroup(tModeLM, ['LM','DIRECT','SIMPLE'], /COLUMN, /EXCLUSIVE, LABEL_TOP='Cycle 4', SET_VALUE=info.inversion_mode(3)-1, UVALUE='INVMODE_CYC4')


	case(info.number_slabs) of
	 	1: begin
	 			widget_control, oneButton, /SET_BUTTON				
				widget_control, info.Btext2, sensitive=0
				widget_control, info.thetaBtext2, sensitive=0
				widget_control, info.chiBtext2, sensitive=0
				widget_control, info.Dopplertext2, sensitive=0
				widget_control, info.Tautext2, sensitive=0
				widget_control, info.Macrotext2, sensitive=0
				widget_control, info.buttonsB2, sensitive=0
				widget_control, info.buttonsthetaB2, sensitive=0
				widget_control, info.buttonschiB2, sensitive=0
				widget_control, info.buttonsDoppler2, sensitive=0
				widget_control, info.buttonsMacro2, sensitive=0
				widget_control, info.buttonstau2, sensitive=0
				widget_control, info.buttonsff, sensitive=0
			end
	 	2: begin
	 			widget_control, twoButton, /SET_BUTTON				
				widget_control, info.Btext2, sensitive=0
				widget_control, info.thetaBtext2, sensitive=0
				widget_control, info.chiBtext2, sensitive=0
				widget_control, info.Dopplertext2, sensitive=0
				widget_control, info.Tautext2, sensitive=1
				widget_control, info.Macrotext2, sensitive=1
				widget_control, info.buttonsB2, sensitive=0
				widget_control, info.buttonsthetaB2, sensitive=0
				widget_control, info.buttonschiB2, sensitive=0
				widget_control, info.buttonsDoppler2, sensitive=0
				widget_control, info.buttonsMacro2, sensitive=1
				widget_control, info.buttonstau2, sensitive=1
				widget_control, info.buttonsff, sensitive=0
			end
	 	3: begin
	 			widget_control, twodifButton, /SET_BUTTON				
				widget_control, info.Btext2, sensitive=1
				widget_control, info.thetaBtext2, sensitive=1
				widget_control, info.chiBtext2, sensitive=1
				widget_control, info.Dopplertext2, sensitive=1
				widget_control, info.Tautext2, sensitive=1
				widget_control, info.Macrotext2, sensitive=1
				widget_control, info.buttonsB2, sensitive=1
				widget_control, info.buttonsthetaB2, sensitive=1
				widget_control, info.buttonschiB2, sensitive=1
				widget_control, info.buttonsDoppler2, sensitive=1
				widget_control, info.buttonsMacro2, sensitive=1
				widget_control, info.buttonstau2, sensitive=0
				widget_control, info.buttonsff, sensitive=0
			end
		-2: begin
	 			widget_control, twocompButton, /SET_BUTTON
				widget_control, info.Btext2, sensitive=1
				widget_control, info.thetaBtext2, sensitive=1
				widget_control, info.chiBtext2, sensitive=1
				widget_control, info.Dopplertext2, sensitive=1
				widget_control, info.Tautext2, sensitive=1
				widget_control, info.Macrotext2, sensitive=1
				widget_control, info.buttonsB2, sensitive=1
				widget_control, info.buttonsthetaB2, sensitive=1
				widget_control, info.buttonschiB2, sensitive=1
				widget_control, info.buttonsDoppler2, sensitive=1
				widget_control, info.buttonsMacro2, sensitive=1
				widget_control, info.buttonstau2, sensitive=1
				widget_control, info.buttonsff, sensitive=1
			end
	 endcase
   
   widget_control, info.baseWidget, SET_UVALUE = info	
   return, info
end