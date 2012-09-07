state = {zeeman_pb: 1, atom_pol: 1, stimulated: 1, mag_opt: 1, rt_mode: 3, boundary: fltarr(4),$
	ncycles: 4, inv_bfield: [0,0,1,1], inv_inclination: [0,0,1,1], inv_azimuth: [0,0,1,1],$
	inv_doppler: [1,1,0,0], inv_tau: [1,1,0,0], inv_depol: [0,0,0,0], inv_macro: [1,1,0,0],$
	inv_damping: [1,1,0,0], inv_source: [0,0,0,0], inv_height: [0,0,0,0],$
	Bfield: 0.d0, thetaB: 0.d0, chiB: 0.d0, vdopp: 0.d0, tau: 0.d0, depol: 0.d0, vmacro: 0.d0,$
	damping: 0.d0, source: 0.d0, height: 0.d0,$
	inversion_mode: [2,1,2,1], $
	stI_weight: [1.d0,1.d0,1.d0,1.d0], stQ_weight: [1.d0,1.d0,1.d0,1.d0], $
	stU_weight: [1.d0,1.d0,1.d0,1.d0], stV_weight: [1.d0,1.d0,1.d0,1.d0], $
	theta_obs: 90.d0, chi_obs: 90.d0, gamm_obs: 0.d0,$
	dir_output_file: './direct.location',$
	dir_stopping: 1, dir_feval: -1, dir_redvol: 0.001, $
	dir_range_bfield: [0.d0,1000.d0], dir_range_thetaB: [0.d0, 180.d0],$
	dir_range_chiB: [0.d0,180.d0], dir_range_vdopp: [0.d0, 20.d0],$
	dir_range_tau: [0.d0, 1.d0], dir_range_depol: [0.d0, 18.d0], $
	dir_range_vmacro: [-10.d0,10.d0], dir_range_damping: [0.d0, 4.d0], $
	dir_range_beta: [0.d0, 1.d0], dir_range_height: [0.d0, 4.d0],$
	verbose: 0, linear_solver: 0
	}