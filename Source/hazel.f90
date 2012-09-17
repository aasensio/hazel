program multiterm
use vars
use maths
use io
use SEE
use rt_coef
use marquardt
use synth
use allen
implicit none

	integer :: iter, nworst_chisq, loop_cycle, successful
	real(kind=8) :: chisq_relative_change, params_relative_change
	logical :: correct
	character(len=120) :: temporal_file
	
	temporal_file = 'temporal.prof'
	
! Initialize the random number generator
	call random_seed
	
! Initialize Allen's data
	call read_allen_data
		
! Fill the factorial array
	call factrl
	
! Fill the 3j symbol array
!	call init_regge(10)
			
! Read the configuration file	
	call read_config
	
! Read the file with the experiment
	call read_experiment(input_experiment, params, fixed)
	
! Read the atomic model	
	call read_model_file(input_model_file)
	
	inversion%min_chisq = 1.d10
	successful = 0

	
!*********************************
!** SYNTHESIS MODE
!*********************************	
	if (working_mode == 0) then
		if (verbose_mode == 1) then
			print *, 'Working in SYNTHESIS mode'
		endif
		observation%n = fixed%no
		allocate(observation%wl(observation%n))
		if (.not.associated(inversion%stokes_unperturbed)) allocate(inversion%stokes_unperturbed(0:3,fixed%no))		
		call do_synthesis(params, fixed, observation, inversion%stokes_unperturbed)
		call write_final_profiles(output_inverted_profiles,observation,inversion)
	endif
	
!*********************************
!** INVERSION MODE
!*********************************
	if (working_mode >= 1) then
		if (verbose_mode == 1) then
			print *, 'Working in INVERSION mode'
		endif

! Read the parameters to be inverted
 		call read_parameters_to_invert(input_inverted_parameters, params, inversion, fixed)
		
! Read the observed profile	
		call read_observation(input_observed_profiles, fixed, observation)
		
		if (.not.associated(inversion%stokes_unperturbed)) allocate(inversion%stokes_unperturbed(0:3,fixed%no))
		if (.not.associated(inversion%stokes_perturbed)) allocate(inversion%stokes_perturbed(0:3,fixed%no))
		if (.not.associated(inversion%dydx)) allocate(inversion%dydx(0:3,fixed%no,params%n_total))
	
! If the first cycle is LM, carry out a first synthesis
! with the original values of the parameters
		if (inversion%algorithm(1) == 1) then
			inversion%loop_cycle = 1		
			call do_synthesis(params, fixed, observation, inversion%stokes_unperturbed)
			inversion%chisq = compute_chisq(observation,inversion)
		endif		
		
	! Loop over the number of cycles
		do loop_cycle = 1, inversion%n_cycles
			
			inversion%loop_cycle = loop_cycle
				
	! Select the parameters to invert in this cycle
			params%inverted = inversion%cycles(:,loop_cycle)			
			
			write(*,FMT='(A)') '*******************************'
			write(*,FMT='(A,I2,A1,I2)') 'Starting cycle ', inversion%loop_cycle, '/', inversion%n_cycles
			write(*,FMT='(A,I2,A)') ' Inverting ', sum(params%inverted),' parameters'
			write(*,FMT='(A)') '*******************************'
				
			write(*,FMT='(A,4(2X,F5.3))') 'Stokes parameters weights : ', inversion%stokes_weights(0:3,loop_cycle)
	
	! LEVENBERG-MARQUARDT		
			if (inversion%algorithm(loop_cycle) == 1) then
	
	! If the weights of the previous step were different, recalculate the value of chi^2
				if (loop_cycle > 1) then
					if (inversion%stokes_weights(0,loop_cycle-1) /= inversion%stokes_weights(0,loop_cycle) .or. &
						inversion%stokes_weights(1,loop_cycle-1) /= inversion%stokes_weights(1,loop_cycle) .or. &
						inversion%stokes_weights(2,loop_cycle-1) /= inversion%stokes_weights(2,loop_cycle) .or. &
						inversion%stokes_weights(3,loop_cycle-1) /= inversion%stokes_weights(3,loop_cycle)) then
					
						call do_synthesis(params, fixed, observation, inversion%stokes_unperturbed)
						inversion%chisq = compute_chisq(observation,inversion)
					endif
				endif
				
				print *, 'LEVENBERG-MARQUARDT MODE'
	
				inversion%lambda = 0.01d0
				nworst_chisq = 0
				inversion%min_chisq = inversion%chisq
				inversion%chisq = 1.d10
				inversion%chisq_old = inversion%chisq			
				chisq_relative_change = 1.d10
				iter = 1
				successful = 1

				call do_synthesis(params, fixed, observation, inversion%stokes_unperturbed)
	
	! Main inversion loop	
				do while (iter < inversion%iter_max .and. nworst_chisq < 10 .and. abs(chisq_relative_change) > 1.d-4 )
	
					write(*,FMT='(A,I3,A1,I3,A,F11.6)') 'Iteration ', iter, '/', inversion%iter_max, ' -- Lambda : ', inversion%lambda

! Only recalculate derivatives if the step has been successful (we located a point with a smaller chi^2)
					if (successful /= 0) call compute_dydx(params,fixed,inversion,observation)
					
					call compute_trial_params(params,fixed,inversion,observation,trial)
	
					call check_boundaries(params,trial,correct)
					
					call print_parameters(params,'  -Old parameters : ',.TRUE.)
					call print_parameters(trial,'  -New parameters : ',.FALSE.)
	
					if (correct) then
						call do_synthesis(trial, fixed, observation, inversion%stokes_unperturbed)
						
						inversion%chisq_old = inversion%chisq
						inversion%chisq = compute_chisq(observation,inversion)
	
! Verify if the chisq is smaller or larger than the previous fit
! Better model
						if (inversion%chisq < inversion%min_chisq) then
							successful = 1							

! Change lambda
							if (inversion%lambda >= 1.d4) then
								inversion%lambda = inversion%lambda / 100.d0
							else if (inversion%lambda >= 1.d-4 .and. inversion%lambda < 1.d4) then
								inversion%lambda = inversion%lambda / 10.d0
							else if (inversion%lambda < 1.d-4) then
								inversion%lambda = inversion%lambda / 5.d0
							endif
							
							if (inversion%lambda < 1.d-6) inversion%lambda = 1.d-6
							!params_relative_change = compute_params_relative_change(params,trial)
							params = trial
							print *, 'New chi^2 : ', inversion%chisq, inversion%lambda
							chisq_relative_change = (inversion%chisq-inversion%chisq_old) / inversion%chisq * 100.d0						
							print *, 'chi^2 relative change [%] : ', chisq_relative_change
							print *, 'Relative change in the parameters [%] : ', params_relative_change
							print *, 'Successful step'
							nworst_chisq = 0
							call write_final_profiles(temporal_file,observation,inversion)
							inversion%min_chisq = inversion%chisq
							!successful = successful + 1
						else
! Worse model
							successful = 0

							if (inversion%lambda >= 1.d4) then
								inversion%lambda = inversion%lambda * 100.d0
							else if (inversion%lambda >= 1.d-4 .and. inversion%lambda < 1.d4) then
								inversion%lambda = inversion%lambda * 10.d0
							else if (inversion%lambda < 1.d-4) then
								inversion%lambda = inversion%lambda * 5.d0
							endif
														
							print *, 'Larger chi^2 : ', inversion%chisq, ' --- Increasing Lambda', inversion%lambda
							nworst_chisq = nworst_chisq + 1
							write(*,FMT='(A,I1,A1,I1)') '   Step : ', nworst_chisq, '/', 5
	
						endif			
					else
						print *, 'Unphysical values. Trying again...'
						inversion%lambda = inversion%lambda * 10.d0
					endif
	
					print *, '---------------------'
					print *
	
					iter = iter + 1				
	
				enddo  ! Iterations
				
				print *, 'Number of calls to forward modeling routines : ', fixed%total_forward_modeling
	
				call do_synthesis(params, fixed, observation, inversion%stokes_unperturbed)
				inversion%chisq = compute_chisq(observation,inversion)
				
	! Write the final profiles
				call write_final_profiles(output_inverted_profiles,observation,inversion)
	
	! Write the final parameters in a file so that it can be used for restarting the inversion code
				call write_experiment(params, fixed)
				
				call print_parameters(params,'-Final Parameters : ',.TRUE.)
				print *, 'Final chi^2 : ', inversion%chisq
			endif
				
					
	!*********************************
	!** INVERSION MODE WITH DIRECT
	!*********************************	
			if (inversion%algorithm(loop_cycle) == 2 .or. inversion%algorithm(loop_cycle) == 3) then
				print *, 'DIRECT MODE'
					
				trial = params
				call invert_with_direct(trial, inversion%algorithm(loop_cycle))
				params = trial
				
				call do_synthesis(params, fixed, observation, inversion%stokes_unperturbed)
				inversion%chisq = compute_chisq(observation,inversion)
			
	! Write the final profiles
				call write_final_profiles(output_inverted_profiles,observation,inversion)
				
	! Write the final parameters in a file so that it can be used for restarting the inversion code
				call write_experiment(params, fixed)
				
				call print_parameters(params,'-Final Parameters : ',.TRUE.)
				print *, 'Final chi^2 : ', inversion%chisq
							
			endif
		
	!*********************************
	!** INVERSION MODE WITH PIKAIA
	!*********************************	
! 			if (inversion%algorithm(loop_cycle) == 3) then			
! 				print *, 'PIKAIA MODE'
! 	
! 				trial = params
! 				call invert_with_pikaia(trial)
! 				params = trial
! 				
! 				call do_synthesis(params, fixed, observation, inversion%stokes_unperturbed)
! 				inversion%chisq = compute_chisq(observation,inversion)
! 			
! 	! Write the final profiles
! 				call write_final_profiles(output_inverted_profiles,observation,inversion)
! 				
! 	! Write the final parameters in a file so that it can be used for restarting the inversion code
! 				call write_experiment(params, fixed)
! 				
! 				call print_parameters(params,'-Final Parameters : ',.TRUE.)
! 				print *, 'Final chi^2 : ', inversion%chisq
! 			endif
			
			print *
			write(*,FMT='(A)') '*******************************'
			write(*,FMT='(A,I2,A1,I2)') 'End of cycle ', inversion%loop_cycle, '/', inversion%n_cycles
			write(*,FMT='(A)') '*******************************'
			print *
				
		enddo  ! Cycles
		
	endif
	
	open(unit=12,file='done.info',action='write',status='replace')
	write(12,*) 'Done.'
	close(12)
	

! Clearn all the allocated memory
!	call clean
	
end program multiterm
