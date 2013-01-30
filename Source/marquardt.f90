module marquardt
use vars
use SEE
use rt_coef
use synth
use svd
use maths
use io
use l_bfgs_b, only : setulb
implicit none
contains

!*************************************************************
!*************************************************************
! MARQUARDT METHOD
!*************************************************************
!*************************************************************

!------------------------------------------------------------
! Calculate the chi^2 function
!------------------------------------------------------------
	function compute_chisq(in_observation,in_inversion)
	type(type_inversion) :: in_inversion
	type(type_observation) :: in_observation
	real(kind=8) :: compute_chisq, weight
	integer :: nf, i, j
		
		compute_chisq = 0.d0

! Compute the chisq using the appropriate weights for each cycle
		do i = 0, 3
			weight = in_inversion%stokes_weights(i,in_inversion%loop_cycle)			
			compute_chisq = compute_chisq + weight * &
				sum((in_observation%stokes(i,:)-in_inversion%stokes_unperturbed(i,:))**2/in_observation%sigma(i,:)**2) / &
				(4.d0*in_observation%n)
		enddo
		
	end function compute_chisq
	
!------------------------------------------------------------
! Calculate the maximum relative change in the parameters
!------------------------------------------------------------
	function compute_params_relative_change(params,trial)
	real(kind=8) :: compute_params_relative_change
	type(variable_parameters) :: params, trial
	real(kind=8), allocatable :: rel_change(:)
	integer :: i, j
	
		allocate(rel_change(params%n_total))

		if (params%bgauss /= 0.d0) then
			rel_change(1) = abs((params%bgauss-trial%bgauss) / params%bgauss)
		else
			rel_change(1) = 0.d0
		endif
		
		if (params%thetabd /= 0.d0) then
			rel_change(2) = abs((params%thetabd-trial%thetabd) / params%thetabd)
		else
			rel_change(2) = 0.d0
		endif
		
		if (params%chibd /= 0.d0) then
			rel_change(3) = abs((params%chibd-trial%chibd) / params%chibd)
		else
			rel_change(3) = 0.d0
		endif
		
		if (params%vdopp /= 0.d0) then
			rel_change(4) = abs((params%vdopp-trial%vdopp) / params%vdopp)
		else
			rel_change(4) = 0.d0
		endif
		
		if (params%dtau /= 0.d0) then
			rel_change(5) = abs((params%dtau-trial%dtau) / params%dtau)
		else
			rel_change(5) = 0.d0
		endif
		
		if (params%delta_collision /= 0.d0) then
			rel_change(6) = abs((params%delta_collision-trial%delta_collision) / params%delta_collision)
		else
			rel_change(6) = 0.d0
		endif
		
		if (params%vmacro /= 0.d0) then
			rel_change(7) = abs((params%vmacro-trial%vmacro) / params%vmacro)
		else
			rel_change(7) = 0.d0
		endif
		
		if (params%damping /= 0.d0) then
			rel_change(8) = abs((params%damping-trial%damping) / params%damping)
		else
			rel_change(8) = 0.d0
		endif
		
		if (params%beta /= 0.d0) then
			rel_change(9) = abs((params%beta-trial%beta) / params%beta)
		else
			rel_change(9) = 0.d0
		endif
		
		if (params%height /= 0.d0) then
			rel_change(10) = abs((params%height-trial%height) / params%height)
		else
			rel_change(10) = 0.d0
		endif

		if (params%dtau2 /= 0.d0) then
			rel_change(11) = abs((params%dtau2-trial%dtau2) / params%dtau2)
		else
			rel_change(11) = 0.d0
		endif

		if (params%vmacro2 /= 0.d0) then
			rel_change(12) = abs((params%vmacro2-trial%vmacro2) / params%vmacro2)
		else
			rel_change(12) = 0.d0
		endif

		if (params%bgauss2 /= 0.d0) then
			rel_change(13) = abs((params%bgauss2-trial%bgauss2) / params%bgauss2)
		else
			rel_change(13) = 0.d0
		endif
		
		if (params%thetabd2 /= 0.d0) then
			rel_change(14) = abs((params%thetabd2-trial%thetabd2) / params%thetabd2)
		else
			rel_change(14) = 0.d0
		endif
		
		if (params%chibd2 /= 0.d0) then
			rel_change(15) = abs((params%chibd2-trial%chibd2) / params%chibd2)
		else
			rel_change(15) = 0.d0
		endif

		if (params%vdopp2 /= 0.d0) then
			rel_change(16) = abs((params%vdopp2-trial%vdopp2) / params%vdopp2)
		else
			rel_change(16) = 0.d0
		endif
		
		if (params%ff /= 0.d0) then
			rel_change(17) = abs((params%ff-trial%ff) / params%ff)
		else
			rel_change(17) = 0.d0
		endif
		
		compute_params_relative_change = maxval(rel_change)*100.d0		

		deallocate(rel_change)
		
		
	end function compute_params_relative_change
	
	
!------------------------------------------------------------
! Check if the values are inside physical boundaries
! If not inside the boundaries, use the previous value
!------------------------------------------------------------
	subroutine check_boundaries(original,trial,correct)
	type(variable_parameters) :: original, trial
	logical :: correct
	real(kind=8) :: upper(17), lower(17)
	integer :: i
	
! Read ranges
		open(unit=24,file=direct_ranges,action='read',status='old')
		call lb(24,12)
		do i = 1, 17
			call lb(24,2)
			read(24,*) lower(i), upper(i)
		enddo
		close(24)
	
! Verify if the parameters are inside the boundaries. If any of them is outside, use the previous value
! of this parameter.
	
		if (trial%bgauss < lower(1) .or. trial%bgauss > upper(1)) then
			trial%bgauss = original%bgauss
		endif
		if (trial%thetabd < lower(2) .or. trial%thetabd > upper(2)) then
			trial%thetabd = original%thetabd
		endif
		if (trial%chibd < lower(3) .or. trial%chibd > upper(3)) then
			trial%chibd = original%chibd
		endif
		if (trial%vdopp < lower(4) .or. trial%vdopp > upper(4)) then
			trial%vdopp = original%vdopp
		endif
		if (trial%dtau < lower(5) .or. trial%dtau > upper(5)) then
			trial%dtau = original%dtau
		endif
		if (original%nslabs /= 1) then
			if (trial%dtau2 < lower(11) .or. trial%dtau2 > upper(11)) then
				trial%dtau2 = original%dtau2
			endif
		endif
		if (trial%delta_collision < lower(6) .or. trial%delta_collision > upper(6)) then
			trial%delta_collision = original%delta_collision
		endif		
		if (trial%vmacro < lower(7) .or. trial%vmacro > upper(7)) then
			trial%vmacro = original%vmacro
		endif
		if (original%nslabs /= 1) then
			if (trial%vmacro2 < lower(12) .or. trial%vmacro2 > upper(12)) then
				trial%vmacro2 = original%vmacro2
			endif
		endif
		if (trial%damping < lower(8) .or. trial%damping > upper(8)) then
			trial%damping = original%damping
		endif
		if (trial%beta < lower(9) .or. trial%beta > upper(9)) then
			trial%beta = original%beta
		endif
		if (trial%height < lower(10) .or. trial%height > upper(10)) then
			trial%height = original%height
		endif

! If two components with two different fields
		if (original%nslabs == 3 .or. original%nslabs == -2) then
			if (trial%bgauss2 < lower(13) .or. trial%bgauss2 > upper(13)) then
				trial%bgauss2 = original%bgauss2
			endif
			if (trial%thetabd2 < lower(14) .or. trial%thetabd2 > upper(14)) then
				trial%thetabd2 = original%thetabd2
			endif
			if (trial%chibd2 < lower(15) .or. trial%chibd2 > upper(15)) then
				trial%chibd2 = original%chibd2
			endif
			if (trial%vdopp2 < lower(16) .or. trial%vdopp2 > upper(16)) then
				trial%vdopp2 = original%vdopp2
			endif
		endif
		
! Two slabs in the same pixel
		if (original%nslabs == -2) then			
			if (trial%ff < lower(17) .or. trial%ff > upper(17)) then
				trial%ff = original%ff
			endif
		endif
		
		correct = .TRUE.
		
	end subroutine check_boundaries
	

!------------------------------------------------------------
! Perturb a given parameter
!------------------------------------------------------------
	subroutine perturb_parameter(in_params,out_params,which,perturbation,delta)
	type(variable_parameters) :: in_params, out_params
	integer :: which
	real(kind=8) :: perturbation, delta

		select case(which)
			case(1) 
				if (in_params%bgauss > 1.d-2) then
					out_params%bgauss = in_params%bgauss * (1.d0+perturbation)
				else
					out_params%bgauss = in_params%bgauss + perturbation
				endif
				delta = out_params%bgauss - in_params%bgauss
				
			case(2)
				if (abs(in_params%thetabd) > 1.d-2) then
					out_params%thetabd = in_params%thetabd * (1.d0+perturbation)
				else
					out_params%thetabd = in_params%thetabd + perturbation
				endif
				delta = out_params%thetabd - in_params%thetabd
				
			case(3) 
				if (abs(in_params%chibd) > 1.d-2) then				
					out_params%chibd = in_params%chibd * (1.d0+perturbation)
				else
					out_params%chibd = in_params%chibd + perturbation
				endif
				delta = out_params%chibd - in_params%chibd
				
			case(4) 
				if (in_params%vdopp > 1.d-2) then
					out_params%vdopp = in_params%vdopp * (1.d0+perturbation)
				else
					out_params%vdopp = in_params%vdopp + perturbation
				endif
				delta = out_params%vdopp - in_params%vdopp
				
			case(5) 
				if (in_params%dtau > 1.d-2) then
					out_params%dtau = in_params%dtau * (1.d0+perturbation)
				else
					out_params%dtau = in_params%dtau + perturbation
				endif
				delta = out_params%dtau - in_params%dtau
				
			case(6) 
				if (in_params%delta_collision > 1.d-2) then
					out_params%delta_collision = in_params%delta_collision * (1.d0+perturbation)
				else
					out_params%delta_collision = in_params%delta_collision + perturbation
				endif
				delta = out_params%delta_collision - in_params%delta_collision
				
			case(7)
				if (abs(in_params%vmacro) > 1.d-2) then				
					out_params%vmacro = in_params%vmacro * (1.d0+perturbation)
				else
					out_params%vmacro = in_params%vmacro + perturbation
				endif
				delta = out_params%vmacro - in_params%vmacro
				
			case(8)
				if (in_params%damping > 1.d-2) then
					out_params%damping = in_params%damping * (1.d0+perturbation)
				else
					out_params%damping = in_params%damping + perturbation
				endif
				delta = out_params%damping - in_params%damping
				
			case(9)
				if (abs(in_params%beta) > 1.d-2) then
					out_params%beta = in_params%beta * (1.d0+perturbation)
				else
					out_params%beta = in_params%beta + perturbation
				endif
				delta = out_params%beta - in_params%beta
				
			case(10)
				if (in_params%height > 1.d-2) then
					out_params%height = in_params%height * (1.d0+perturbation)
				else
					out_params%height = in_params%height + perturbation
				endif
				delta = out_params%height - in_params%height

			case(11) 
				if (in_params%dtau2 > 1.d-2) then
					out_params%dtau2 = in_params%dtau2 * (1.d0+perturbation)
				else
					out_params%dtau2 = in_params%dtau2 + perturbation
				endif
				delta = out_params%dtau2 - in_params%dtau2

			case(12)
				if (abs(in_params%vmacro2) > 1.d-2) then				
					out_params%vmacro2 = in_params%vmacro2 * (1.d0+perturbation)
				else
					out_params%vmacro2 = in_params%vmacro2 + perturbation
				endif
				delta = out_params%vmacro2 - in_params%vmacro2

			case(13) 
				if (in_params%bgauss2 > 1.d-2) then
					out_params%bgauss2 = in_params%bgauss2 * (1.d0+perturbation)
				else
					out_params%bgauss2 = in_params%bgauss2 + perturbation
				endif
				delta = out_params%bgauss2 - in_params%bgauss2
				
			case(14)
				if (abs(in_params%thetabd2) > 1.d-2) then
					out_params%thetabd2 = in_params%thetabd2 * (1.d0+perturbation)
				else
					out_params%thetabd2 = in_params%thetabd2 + perturbation
				endif
				delta = out_params%thetabd2 - in_params%thetabd2
				
			case(15) 
				if (abs(in_params%chibd2) > 1.d-2) then				
					out_params%chibd2 = in_params%chibd2 * (1.d0+perturbation)
				else
					out_params%chibd2 = in_params%chibd2 + perturbation
				endif
				delta = out_params%chibd2 - in_params%chibd2

			case(16) 
				if (in_params%vdopp2 > 1.d-2) then
					out_params%vdopp2 = in_params%vdopp2 * (1.d0+perturbation)
				else
					out_params%vdopp2 = in_params%vdopp2 + perturbation
				endif
				delta = out_params%vdopp2 - in_params%vdopp2
				
			case(17) 
				if (in_params%ff > 1.d-2) then
					out_params%ff = in_params%ff * (1.d0+perturbation)
				else
					out_params%ff = in_params%ff + perturbation
				endif
				delta = out_params%ff - in_params%ff
			
		end select
		
	end subroutine perturb_parameter
	
!------------------------------------------------------------
! Perturb a given parameter
!------------------------------------------------------------
	subroutine add_to_parameter(in_params,out_params,which,new_value)
	type(variable_parameters) :: in_params, out_params
	integer :: which
	real(kind=8) :: new_value, rel_change

! Limit the maximum relative change to 50%
		if (new_value > 0.5d0) then
			write(*,FMT='(A,A)') 'Clipping ', trim(adjustl(parameters_name(which)))
			new_value = 0.5d0
		endif
		select case(which)
			case(1)				
				out_params%bgauss = in_params%bgauss + new_value
			case(2) 
				out_params%thetabd = in_params%thetabd + new_value
			case(3) 
				out_params%chibd = in_params%chibd + new_value
			case(4) 
				out_params%vdopp = in_params%vdopp + new_value
			case(5) 
				out_params%dtau = in_params%dtau + new_value
			case(6) 
				out_params%delta_collision = in_params%delta_collision + new_value
			case(7) 
				out_params%vmacro = in_params%vmacro + new_value	
			case(8) 
				out_params%damping = in_params%damping + new_value
			case(9) 
				out_params%beta = in_params%beta + new_value
			case(10) 
				out_params%height = in_params%height + new_value
			case(11) 
				out_params%dtau2 = in_params%dtau2 + new_value
			case(12) 
				out_params%vmacro2 = in_params%vmacro2 + new_value
			case(13)				
				out_params%bgauss2 = in_params%bgauss2 + new_value
			case(14) 
				out_params%thetabd2 = in_params%thetabd2 + new_value
			case(15) 
				out_params%chibd2 = in_params%chibd2 + new_value
			case(16) 
				out_params%vdopp2 = in_params%vdopp2 + new_value
			case(17) 
				out_params%ff = in_params%ff + new_value
		end select
		
	end subroutine add_to_parameter	
	
!------------------------------------------------------------
! Calculate the derivatives of the Stokes profiles
!------------------------------------------------------------
	subroutine compute_dydx(in_params,in_fixed,in_inversion,in_observation)
	type(variable_parameters) :: in_params, params_scaled, params_scaled_modified, params_modified
	type(fixed_parameters) :: in_fixed
	type(type_inversion) :: in_inversion
	type(type_observation) :: in_observation
	integer :: i
	real(kind=8) :: delta
								
! Do the synthesis without perturbing the magnetic field		
! 		call do_synthesis(in_params, in_fixed, in_observation, in_inversion%stokes_unperturbed)
		
		params_scaled = in_params		
		params_modified = in_params
		
		call compress_parameters(in_params, params_scaled)		
		params_scaled_modified = params_scaled

		do i = 1, in_params%n_total
		
			if (in_params%inverted(i) == 1) then
				write(*,FMT='(A,I2,A,A)') 'Perturbing parameter ', i, '  --- ', trim(adjustl(parameters_name(i)))
		
! Perturb the given parameter and do another synthesis														
				call perturb_parameter(params_scaled,params_scaled_modified,i,0.0001d0,delta)
				
				call expand_parameters(params_scaled_modified, params_modified)
				
				call do_synthesis(params_modified, in_fixed, in_observation, in_inversion%stokes_perturbed)
		
				in_inversion%dydx(:,:,i) = (in_inversion%stokes_perturbed - in_inversion%stokes_unperturbed) / delta
			endif
		enddo
	
	end subroutine compute_dydx
	
!------------------------------------------------------------
! Calculate the derivatives of the Stokes profiles
!------------------------------------------------------------
	subroutine compute_trial_params(in_params,in_fixed,in_inversion,in_observation,in_trial)
	type(variable_parameters) :: in_params, in_trial, in_scaled, in_temp
	type(fixed_parameters) :: in_fixed
	type(type_inversion) :: in_inversion
	type(type_observation) :: in_observation
	real(kind=8), allocatable :: alpha(:,:), beta(:), w(:), v(:,:), x(:)
	real(kind=8) :: obs, syn, sig, weight, wmax
	integer :: i, j, k, l, np
			
		in_trial = in_params		
		call compress_parameters(in_params,in_scaled)
		in_temp = in_scaled
				
		np = in_params%n_total
		
		allocate(alpha(np,np))
		allocate(beta(np))
		allocate(v(np,np))
		allocate(w(np))
		allocate(x(np))
		
		alpha = 0.d0
		beta = 0.d0
						
		do i = 1, in_fixed%no
			do j = 0, 3
				weight = in_inversion%stokes_weights(j,in_inversion%loop_cycle)
				do k = 1, np
					
					obs = in_observation%stokes(j,i)
					syn = in_inversion%stokes_unperturbed(j,i)
					sig = in_observation%sigma(j,i)
										
					if (in_params%inverted(k) == 1) then
						beta(k) = beta(k) - 2.d0 * weight * (obs-syn) / sig**2 * in_inversion%dydx(j,i,k) / (4.d0*in_fixed%no)
						do l = 1, np
							if (in_params%inverted(l) == 1) then
								alpha(k,l) = alpha(k,l) + 2.d0 * weight * &
									in_inversion%dydx(j,i,k) * in_inversion%dydx(j,i,l) / sig**2 / (4.d0*in_fixed%no)
							endif
						enddo
					endif
				enddo
			enddo
		enddo
		
		beta = -0.5d0 * beta
		alpha = 0.5d0 * alpha
		
! Apply the Lambda parameter of the Levenverb-Marquardt		
		do i = 1, np
			alpha(i,i) = (1.d0+inversion%lambda) * alpha(i,i)
		enddo
		
		np = in_params%n_total
		call svdcmp(alpha,np,np,np,np,w,v)

! Regularize the linear system by setting all singular values below a certain
! threshold equal to zero
		wmax = maxval(w)
		do i = 1, np
			if (w(i) < wmax * 1.d-6) then
				w(i) = 0.d0
			endif
		enddo
		
		call svbksb(alpha,w,v,np,np,np,np,beta,x)				

		do i = 1, np
			if (in_params%inverted(i) == 1) then				
				call add_to_parameter(in_scaled,in_temp,i,x(i)/1.d0)				
			endif				
		enddo

		call expand_parameters(in_temp,in_trial)
		
		deallocate(alpha)
		deallocate(beta)
		deallocate(v)
		deallocate(w)
		deallocate(x)
	
	end subroutine compute_trial_params
	
!------------------------------------------------------------
! Re-scale the parameters to make them of the order of unity
!------------------------------------------------------------
	subroutine compress_parameters(in_params,in_params_scaled)
	type(variable_parameters) :: in_params, in_params_scaled
		in_params_scaled%bgauss = in_params%bgauss / 1000.d0
		in_params_scaled%thetabd = in_params%thetabd / 10.d0
		in_params_scaled%chibd = in_params%chibd / 10.d0		
		in_params_scaled%vdopp = in_params%vdopp / 5.d0
		in_params_scaled%dtau = in_params%dtau / 1.5d0
		in_params_scaled%delta_collision = in_params%delta_collision
		in_params_scaled%vmacro = in_params%vmacro / 5.d0
		in_params_scaled%damping = in_params%damping
		in_params_scaled%beta = in_params%beta / 3.d0
		in_params_scaled%height = in_params%height / 30.d0

! Two components with same field and Doppler width
		if (in_params%nslabs /= 1) then
			in_params_scaled%dtau2 = in_params%dtau2 / 1.5d0
			in_params_scaled%vmacro2 = in_params%vmacro2 / 5.d0
		endif

! Two components with different field and Doppler width
		if (in_params%nslabs == 3 .or. in_params%nslabs == -2) then
			in_params_scaled%bgauss2 = in_params%bgauss2 / 1000.d0
			in_params_scaled%thetabd2 = in_params%thetabd2 / 10.d0
			in_params_scaled%chibd2 = in_params%chibd2 / 10.d0
			in_params_scaled%vdopp2 = in_params%vdopp2 / 5.d0
		endif
		
! Two components in the same pixel
		if (in_params%nslabs == -2) then
			in_params_scaled%ff = in_params%ff
		endif

		
	end subroutine compress_parameters	
	
!------------------------------------------------------------
! Re-scale the parameters to their original values
!------------------------------------------------------------
	subroutine expand_parameters(in_params_scaled,in_params)
	type(variable_parameters) :: in_params, in_params_scaled
		in_params%bgauss = in_params_scaled%bgauss * 1000.d0
		in_params%thetabd = in_params_scaled%thetabd * 10.d0
		in_params%chibd = in_params_scaled%chibd * 10.d0
		in_params%vdopp = in_params_scaled%vdopp * 5.d0
		in_params%dtau = in_params_scaled%dtau * 1.5d0
		in_params%delta_collision = in_params_scaled%delta_collision
		in_params%vmacro = in_params_scaled%vmacro * 5.d0
		in_params%damping = in_params_scaled%damping
		in_params%beta = in_params_scaled%beta * 3.d0
		in_params%height = in_params_scaled%height * 30.d0

! Two components with same field and Doppler width
		if (in_params%nslabs /= 1) then
			in_params%dtau2 = in_params_scaled%dtau2 * 1.5d0
			in_params%vmacro2 = in_params_scaled%vmacro2 * 5.d0
		endif

! Two components with different field and Doppler width
		if (in_params%nslabs == 3 .or. in_params%nslabs == -2) then
			in_params%bgauss2 = in_params_scaled%bgauss2 * 1000.d0
			in_params%thetabd2 = in_params_scaled%thetabd2 * 10.d0
			in_params%chibd2 = in_params_scaled%chibd2 * 10.d0
			in_params%vdopp2 = in_params_scaled%vdopp2 * 5.d0
		endif
		
! Two components in the same pixel
		if (in_params%nslabs == -2) then
			in_params%ff = in_params_scaled%ff
		endif
			
	end subroutine expand_parameters	



!*************************************************************
!*************************************************************
! DIRECT METHOD
!*************************************************************
!*************************************************************

!------------------------------------------------------------
! Invert some parameters with the DIRECT method
!------------------------------------------------------------
	subroutine invert_with_direct(out_params, synth_option)
	real(kind=8) :: DIReps, DIRf
	integer :: ndim, DIRmaxf, DIRmaxT, synth_option
	integer :: DIRalg
	integer :: IError, logfile
	real(kind=8) :: fglobal, fglper, volper, sigmaper
	real(kind=8), allocatable :: u(:), l(:), upper(:), lower(:)
	integer :: n, i, j
	real(kind=8), allocatable :: DIRx(:)
	
	integer, parameter :: iisize=300, idsize=300, icsize=30
	integer :: iidata(iisize)
	real(kind=8) :: ddata(idsize)
	character(len=40) :: cdata(icsize)
	character(len=120) :: location_file_direct
	
	type(variable_parameters) :: out_params
						
		open(unit=23,file='logfile',action='write',status='replace')
		
		logfile = 23
		fglobal = -1.d100
		fglper = 0.d0
		sigmaper  = -1.d0
		
		DIReps = 1.d-6
		DIRmaxT = 6000
		
		ndim = sum(params%inverted)
		allocate(u(ndim))
		allocate(l(ndim))
		allocate(DIRx(ndim))
		
		
		allocate(upper(17))
		allocate(lower(17))
		open(unit=24,file=direct_ranges,action='read',status='old')
		call lb(24,3)
		
		call lb(24,2)
		read(24,*) location_file_direct
		
		open(unit=25,file=location_file_direct,action='write',status='replace')
		
		call lb(24,2)
		read(24,*) DIRmaxf
		if (DIRmaxf < 0) DIRmaxf = 10000
		print *, 'Maximum number of function evaluations : ', DIRmaxf
		
		call lb(24,2)
		read(24,*) volper
		if (volper < 0) volper = 1.d-20
		print *, 'Stopping when the hypervolume with respect to the original is less than ', volper
		
		do i = 1, 17
			call lb(24,2)
			read(24,*) lower(i), upper(i)
		enddo
		close(24)
		
		j = 1
		do i = 1, params%n_total
			if (params%inverted(i) == 1) then
				l(j) = lower(i)
				u(j) = upper(i)
				print *, 'Inverting ', trim(adjustl(parameters_name(i)))
				j = j + 1
			endif
		enddo
		
		deallocate(upper)
		deallocate(lower)

		DIRalg = 1

		if (synth_option == 2) then
			call DIRECT(fcn,DIRx,ndim,DIReps,DIRmaxf,DIRmaxT, DIRf, l, u, DIRalg, Ierror, logfile, &
				fglobal, fglper, volper, sigmaper, iidata, iisize, ddata, idsize, cdata, icsize)
		else
			call DIRECT(fcn_simplified_StokesI,DIRx,ndim,DIReps,DIRmaxf,DIRmaxT, DIRf, l, u, DIRalg, Ierror, logfile, &
				fglobal, fglper, volper, sigmaper, iidata, iisize, ddata, idsize, cdata, icsize)
		endif
			
		j = 1
		do i = 1, params%n_total
			if (params%inverted(i) == 1) then
				select case(i)
					case(1) 
						out_params%bgauss = DIRx(j)
					case(2) 
						out_params%thetabd = DIRx(j)
					case(3) 
						out_params%chibd = DIRx(j)
					case(4) 
						out_params%vdopp = DIRx(j)
					case(5) 
						out_params%dtau = DIRx(j)
					case(6) 
						out_params%delta_collision = DIRx(j)
					case(7)
						out_params%vmacro = DIRx(j)
					case(8)
						out_params%damping = DIRx(j)
					case(9)
						out_params%beta = DIRx(j)
					case(10)
						out_params%height = DIRx(j)
					case(11)
						out_params%dtau2 = DIRx(j)
					case(12)
						out_params%vmacro2 = DIRx(j)
					case(13) 
						out_params%bgauss2 = DIRx(j)
					case(14) 
						out_params%thetabd2 = DIRx(j)
					case(15) 
						out_params%chibd2 = DIRx(j)
					case(16) 
						out_params%vdopp2 = DIRx(j)
					case(17) 
						out_params%ff = DIRx(j)
				end select
				j = j + 1
			endif
		enddo
			
		close(23)
		close(25)

! K-means clustering
! Read the sampling carried out by DIRECT
! 		open(unit=25,file=location_file_direct,action='read',status='old')
! 		
! 		close(25)
! 		call cluster_initialize_5 ( 2, point_num, cluster_num, point, &
! 			cluster_center )
! 
! 		call kmeans_01 ( dim_num, point_num, cluster_num, it_max, it, point, &
! 			cluster, cluster_center, cluster_population, cluster_energy )
			
	end subroutine invert_with_direct
	
!------------------------------------------------------------
! Function that returns the chi^2 for the DIRECT problem. This solves the full RT problem
!------------------------------------------------------------
	subroutine fcn(n, x, f, flag, iidata, iisize, ddata, idsize, cdata, icsize)
	integer :: n,flag
   real(kind=8) :: x(n)
	real(kind=8) :: f
	integer :: iisize, idsize, icsize, i, j
	integer :: iidata(iisize)
	real(kind=8) :: ddata(idsize)
	character(len=40) :: cdata(icsize)
	type(variable_parameters) :: trial
	character(len=120) :: temporal_file

		temporal_file = 'temporal.prof'
		
		j = 1
		do i = 1, params%n_total
			if (params%inverted(i) == 1) then
				select case(i)
					case(1) 
						params%bgauss = x(j)
					case(2) 
						params%thetabd = x(j)
					case(3) 
						params%chibd = x(j)
					case(4) 
						params%vdopp = x(j)
					case(5) 
						params%dtau = x(j)
					case(6) 
						params%delta_collision = x(j)
					case(7)
						params%vmacro = x(j)
					case(8)
						params%damping = x(j)
					case(9)
						params%beta = x(j)
					case(10)
						params%height = x(j)
					case(11)
						params%dtau2 = x(j)
					case(12)
						params%vmacro2 = x(j)
					case(13) 
						params%bgauss2 = x(j)
					case(14) 
						params%thetabd2 = x(j)
					case(15) 
						params%chibd2 = x(j)
					case(16) 
						params%vdopp2 = x(j)
					case(17) 
						params%ff = x(j)
				end select
				j = j + 1
			endif
		enddo
			
		call do_synthesis(params, fixed, observation, inversion%stokes_unperturbed)
		call write_final_profiles(temporal_file, observation, inversion)
		
		f = compute_chisq(observation,inversion)
		
		call print_parameters(params,'      -Parameters : ',.TRUE.)
		
 		print *, 'chi^2 : ', f
		write(25,FMT='(11(E15.5,2X))') x, f
		
		flag = 0

	end subroutine fcn


!------------------------------------------------------------
! Function that returns the chi^2 for the DIRECT problem. This solves an approximate RT problem for Stokes I
!------------------------------------------------------------
	subroutine fcn_simplified_StokesI(n, x, f, flag, iidata, iisize, ddata, idsize, cdata, icsize)
	integer :: n,flag
   real(kind=8) :: x(n)
	real(kind=8) :: f
	integer :: iisize, idsize, icsize, i, j
	integer :: iidata(iisize)
	real(kind=8) :: ddata(idsize)
	character(len=40) :: cdata(icsize)
	type(variable_parameters) :: trial
	character(len=120) :: temporal_file
	real(kind=8), allocatable :: onum(:), prof(:), prof2(:)
	real(kind=8) :: va, dnum, adamp, onum0, onum1, onum2

		temporal_file = 'temporal.prof'

		j = 1
		do i = 1, params%n_total
			if (params%inverted(i) == 1) then
				select case(i)
					case(1)
						params%bgauss = x(j)
					case(2)
						params%thetabd = x(j)
					case(3)
						params%chibd = x(j)
					case(4)
						params%vdopp = x(j)
					case(5)
						params%dtau = x(j)
					case(6)
						params%delta_collision = x(j)
					case(7)
						params%vmacro = x(j)
					case(8)
						params%damping = x(j)
					case(9)
						params%beta = x(j)
					case(10)
						params%height = x(j)
					case(11)
						params%dtau2 = x(j)
					case(12)
						params%vmacro2 = x(j)
					case(13)
						params%bgauss2 = x(j)
					case(14)
						params%thetabd2 = x(j)
					case(15)
						params%chibd2 = x(j)
					case(16)
						params%vdopp2 = x(j)
					case(17)
						params%ff = x(j)
				end select
				j = j + 1
			endif
		enddo

		allocate(onum(fixed%no))
		allocate(prof(fixed%no))

		onum = -1.d8 * observation%wl / fixed%wl**2

		inversion%stokes_unperturbed = 0.d0

! Only one slab
		if (params%nslabs == 1) then
			dnum = params%vdopp*1.d5 / (fixed%wl*1.d-8*PC)
			va = params%vmacro*1.d5 / (fixed%wl*1.d-8*PC)
			adamp = params%damping

 			if (fixed%damping_treatment == 1) then
 				adamp = fixed%wl * 1.d-8 / (params%vdopp*1.d5) * aesto(fixed%nemiss) * abs(params%damping)
 			endif

			onum0 = -1.26d0*1.d-8 / (fixed%wl*1.d-8)**2
			onum1 = (-1.26d0+0.09d0)*1.d-8 / (fixed%wl*1.d-8)**2
			onum2 = 0.d0
			inversion%stokes_unperturbed(0,:) = profile(adamp,(onum0-onum-va)/dnum) + &
				0.6*profile(adamp,(onum1-onum-va)/dnum) + &
				0.2*profile(adamp,(onum2-onum-va)/dnum)

			if (fixed%Stokes_incident(0) == 0) then
				inversion%stokes_unperturbed(0,:) = 1.d0 - exp(-params%dtau / 1.56d0 * inversion%stokes_unperturbed(0,:))
			else
				inversion%stokes_unperturbed(0,:) = exp(-params%dtau / 1.56d0 * inversion%stokes_unperturbed(0,:))
			endif
		endif

! Two slabs
 		if (params%nslabs == 2 .or. params%nslabs == 3) then
! First slab
			dnum = params%vdopp*1.d5 / (fixed%wl*1.d-8*PC)
			va = params%vmacro*1.d5 / (fixed%wl*1.d-8*PC)
			adamp = params%damping

 			if (fixed%damping_treatment == 1) then
 				adamp = fixed%wl * 1.d-8 / (params%vdopp*1.d5) * aesto(fixed%nemiss) * abs(params%damping)
 			endif

			onum0 = -1.26d0*1.d-8 / (fixed%wl*1.d-8)**2
			onum1 = (-1.26d0+0.09d0)*1.d-8 / (fixed%wl*1.d-8)**2
			onum2 = 0.d0
			prof = profile(adamp,(onum0-onum-va)/dnum) + &
				0.6*profile(adamp,(onum1-onum-va)/dnum) + &
				0.2*profile(adamp,(onum2-onum-va)/dnum)

			if (fixed%Stokes_incident(0) == 0) then
				inversion%stokes_unperturbed(0,:) = 1.d0 - exp(-params%dtau / 1.56d0 * prof)
			else
				inversion%stokes_unperturbed(0,:) = exp(-params%dtau / 1.56d0 * prof)
			endif

! Second slab
! If nslabs=3, then we use a different Doppler width. If not, we use the same for both components
			if (params%nslabs == 3) then
				dnum = params%vdopp2*1.d5 / (fixed%wl*1.d-8*PC)
				if (fixed%damping_treatment == 1) then
 				adamp = fixed%wl * 1.d-8 / (params%vdopp2*1.d5) * aesto(fixed%nemiss) * abs(params%damping)
 			endif
 			
			endif
			va = params%vmacro2*1.d5 / (fixed%wl*1.d-8*PC)
			adamp = params%damping

			onum0 = -1.26d0*1.d-8 / (fixed%wl*1.d-8)**2
			onum1 = (-1.26d0+0.09d0)*1.d-8 / (fixed%wl*1.d-8)**2
			onum2 = 0.d0
			prof = profile(adamp,(onum0-onum-va)/dnum) + &
				0.6*profile(adamp,(onum1-onum-va)/dnum) + &
				0.2*profile(adamp,(onum2-onum-va)/dnum)

			inversion%stokes_unperturbed(0,:) = exp(-params%dtau2 / 1.56d0 * prof) * inversion%stokes_unperturbed(0,:)

 		endif

! Two slabs on the same pixel
		if (params%nslabs == -2) then

			allocate(prof2(fixed%no))

! First slab
			dnum = params%vdopp*1.d5 / (fixed%wl*1.d-8*PC)
			va = params%vmacro*1.d5 / (fixed%wl*1.d-8*PC)
			adamp = params%damping

 			if (fixed%damping_treatment == 1) then
 				adamp = fixed%wl * 1.d-8 / (params%vdopp*1.d5) * aesto(fixed%nemiss) * abs(params%damping)
 			endif

			onum0 = -1.26d0*1.d-8 / (fixed%wl*1.d-8)**2
			onum1 = (-1.26d0+0.09d0)*1.d-8 / (fixed%wl*1.d-8)**2
			onum2 = 0.d0
			prof = profile(adamp,(onum0-onum-va)/dnum) + &
				0.6*profile(adamp,(onum1-onum-va)/dnum) + &
				0.2*profile(adamp,(onum2-onum-va)/dnum)

! Second slab
			dnum = params%vdopp2*1.d5 / (fixed%wl*1.d-8*PC)
			va = params%vmacro2*1.d5 / (fixed%wl*1.d-8*PC)
			adamp = params%damping

 			if (fixed%damping_treatment == 1) then
 				adamp = fixed%wl * 1.d-8 / (params%vdopp2*1.d5) * aesto(fixed%nemiss) * abs(params%damping)
 			endif

			onum0 = -1.26d0*1.d-8 / (fixed%wl*1.d-8)**2
			onum1 = (-1.26d0+0.09d0)*1.d-8 / (fixed%wl*1.d-8)**2
			onum2 = 0.d0
			prof2 = profile(adamp,(onum0-onum-va)/dnum) + &
				0.6*profile(adamp,(onum1-onum-va)/dnum) + &
				0.2*profile(adamp,(onum2-onum-va)/dnum)

			if (fixed%Stokes_incident(0) == 0) then
				inversion%stokes_unperturbed(0,:) = params%ff * (1.d0 - exp(-params%dtau / 1.56d0 * prof)) + &
					(1.d0-params%ff) * (1.d0 - exp(-params%dtau2 / 1.56d0 * prof2))
			else
				inversion%stokes_unperturbed(0,:) = params%ff * exp(-params%dtau / 1.56d0 * prof) + &
					(1.d0-params%ff) * exp(-params%dtau2 / 1.56d0 * prof2)
			endif

			deallocate(prof2)
		endif


		f = compute_chisq(observation,inversion)

		deallocate(onum)
		deallocate(prof)

 		call write_final_profiles(temporal_file, observation, inversion)

		call print_parameters(params,'      -Parameters : ',.TRUE.)

 		print *, 'chi^2 : ', f
		write(25,FMT='(11(E15.5,2X))') x, f

		flag = 0

	end subroutine fcn_simplified_StokesI

!*************************************************************
!*************************************************************
! PIKAIA METHOD (GENETIC ALGORITHM)
!*************************************************************
!*************************************************************

!------------------------------------------------------------
! Invert some parameters with the GENETIC ALGORITHM (PIKAIA)
!------------------------------------------------------------
	subroutine invert_with_pikaia(out_params)
	integer :: n, seed, status, i, j
	real(kind=4) :: ctrl(12), f
	real(kind=4), allocatable :: x(:)
	type(variable_parameters) :: out_params
			
		n = sum(params%inverted)
		
! Use default values
		ctrl(1) = 100
		ctrl(2) = 5
		ctrl(3) = 5
		ctrl(4) = 0.85
		ctrl(5) = 2
		ctrl(6) = 0.005
		ctrl(7) = 0.0005
		ctrl(8) = 0.25
		ctrl(9) = 1.0
		ctrl(10) = 1
		ctrl(11) = 1
		ctrl(12) = 1
		
		seed = 1234
		call rninit(seed)
		
		allocate(x(n))
		
		call pikaia(chisq_pikaia,n,ctrl,x,f,status)
		
		j = 1
		do i = 1, params%n_total
			if (params%inverted(i) == 1) then
				select case(i)
					case(1) 
						out_params%bgauss = x(j)
					case(2) 
						out_params%thetabd = x(j)
					case(3) 
						out_params%chibd = x(j)
					case(4) 
						out_params%vdopp = x(j)
					case(5) 
						out_params%dtau = x(j)
					case(6) 
						out_params%delta_collision = x(j)
					case(7)
						out_params%vmacro = x(j)
					case(8)
						out_params%damping = x(j)
					case(9)
						out_params%beta = x(j)
					case(10)
						out_params%height = x(j)
					case(11)
						out_params%dtau2 = x(j)
					case(12)
						out_params%vmacro2 = x(j)
				end select
				j = j + 1
			endif
		enddo
		
		deallocate(x)
	
	end subroutine invert_with_pikaia
	
!------------------------------------------------------------
! Function that returns the chisq for the PIKAIA algorithm
!------------------------------------------------------------
	function chisq_pikaia(n, in_x)
	real(kind=4) :: chisq_pikaia
	integer :: n, i, j
   real(kind=4) :: x(n), in_x(n), f
					
		j = 1		

		do i = 1, params%n_total
			if (params%inverted(i) == 1) then

! Undo the change of variables that puts everything in the range [0,1]
				x(j) = minim_pikaia(i) - (minim_pikaia(i) - maxim_pikaia(i)) * in_x(i)
				
				select case(i)
					case(1)
						params%bgauss = x(j)
					case(2) 
						params%thetabd = x(j)
					case(3) 
						params%chibd = x(j)
					case(4) 
						params%vdopp = x(j)
					case(5) 
						params%dtau = x(j)
					case(6) 
						params%delta_collision = x(j)
					case(7)
						params%vmacro = x(j)
					case(8)
						params%damping = x(j)
					case(9)
						params%beta = x(j)
					case(10)
						params%height = x(j)
					case(11)
						params%dtau2 = x(j)
					case(12)
						params%vmacro2 = x(j)
				end select
				j = j + 1
			endif
		enddo
		
		call do_synthesis(params, fixed, observation, inversion%stokes_unperturbed)
		
		f = compute_chisq(observation,inversion)

		call print_parameters(params,'      -Parameters : ',.TRUE.)
						
		chisq_pikaia = 1.d0 / f		
		print *, 'chi^2 : ', f, chisq_pikaia
		
	end function chisq_pikaia


!*************************************************************
!*************************************************************
! SIMPLIFIED RADIATIVE TRANSFER FOR STOKES I
!*************************************************************
!*************************************************************

!------------------------------------------------------------
! Invert some parameters with the DIRECT method
!------------------------------------------------------------
	subroutine invert_with_simplified_transfer(out_params)
	real(kind=8), allocatable :: u(:), l(:), upper(:), lower(:)
	integer :: n, i, j, ndim, iprint, m
	real(kind=8), allocatable :: x(:), g(:), xold(:), wa(:)
	real(kind=8) :: f, factr, pgtol, fnew
	type(variable_parameters) :: out_params
	character(len=60)      :: task, csave
	logical                :: lsave(4)
	integer                :: isave(44)
	real(kind=8)               :: dsave(29)
	integer,  allocatable  :: nbd(:), iwa(:)


		factr = 1.d7
		pgtol = 1.d-5
		iprint = 0
		m = 5
		
		ndim = sum(params%inverted)
		allocate(u(ndim))
		allocate(l(ndim))
		allocate(x(ndim))
		allocate(xold(ndim))

		allocate(nbd(ndim))
		allocate(iwa(3*ndim))
		allocate(g(ndim))
		allocate(wa(2*m*n + 5*n + 11*m*m + 8*m))

		nbd = 2

		allocate(upper(17))
		allocate(lower(17))
		open(unit=24,file=direct_ranges,action='read',status='old')
		call lb(24,3)

		call lb(24,2)
		read(24,*)

		call lb(24,2)
		read(24,*)
		call lb(24,2)
		read(24,*)

		do i = 1, 17
			call lb(24,2)
			read(24,*) lower(i), upper(i)
		enddo
		close(24)

		j = 1
		do i = 1, params%n_total
			if (params%inverted(i) == 1) then
				l(j) = lower(i)
				u(j) = upper(i)
				print *, 'Inverting ', trim(adjustl(parameters_name(i)))
				j = j + 1
			endif
		enddo

		deallocate(upper)
		deallocate(lower)

! Invert the Stokes I profile with the simplified RT method

		task = 'START'

		do while(task(1:2) == 'FG' .or. task == 'NEW_X' .or. task == 'START')

!     This is the call to the L-BFGS-B code.
         call setulb ( ndim, m, x, l, u, nbd, f, g, factr, pgtol, &
                       wa, iwa, task, iprint,&
                       csave, lsave, isave, dsave )

! Evaluate the forward model and the gradient
         if (task(1:2) == 'FG') then

				call simplified_transfer(ndim, x, f)

				xold = x

! Calculate gradient numerically
				do i = 1, ndim
					x(i) = x(i) + 1.d-3
					call simplified_transfer(ndim, x, fnew)
					g(i) = (fnew - f) / 1.d-3
					x(i) = xold(i)					
				enddo

				print *, 'f,g= ', f, g

         endif
      enddo

      stop

		j = 1
		do i = 1, params%n_total
			if (params%inverted(i) == 1) then
				select case(i)
					case(1)
						out_params%bgauss = x(j)
					case(2)
						out_params%thetabd = x(j)
					case(3)
						out_params%chibd = x(j)
					case(4)
						out_params%vdopp = x(j)
					case(5)
						out_params%dtau = x(j)
					case(6)
						out_params%delta_collision = x(j)
					case(7)
						out_params%vmacro = x(j)
					case(8)
						out_params%damping = x(j)
					case(9)
						out_params%beta = x(j)
					case(10)
						out_params%height = x(j)
					case(11)
						out_params%dtau2 = x(j)
					case(12)
						out_params%vmacro2 = x(j)
					case(13)
						out_params%bgauss2 = x(j)
					case(14)
						out_params%thetabd2 = x(j)
					case(15)
						out_params%chibd2 = x(j)
					case(16)
						out_params%vdopp2 = x(j)
					case(17)
						out_params%ff = x(j)
				end select
				j = j + 1
			endif
		enddo

	end subroutine invert_with_simplified_transfer

!------------------------------------------------------------
! Function that returns the the DIRECT
!------------------------------------------------------------
	subroutine simplified_transfer(n, x, f)
	integer :: n, i, j
   real(kind=8) :: x(n)
	real(kind=8) :: f, g(n)
	type(variable_parameters) :: trial
	character(len=120) :: temporal_file
	real(kind=8), allocatable :: onum(:)
	real(kind=8) :: va, dnum, adamp, onum0, onum1, onum2

		j = 1
		do i = 1, params%n_total
			if (params%inverted(i) == 1) then
				select case(i)
					case(1)
						params%bgauss = x(j)
					case(2)
						params%thetabd = x(j)
					case(3)
						params%chibd = x(j)
					case(4)
						params%vdopp = x(j)
					case(5)
						params%dtau = x(j)
					case(6)
						params%delta_collision = x(j)
					case(7)
						params%vmacro = x(j)
					case(8)
						params%damping = x(j)
					case(9)
						params%beta = x(j)
					case(10)
						params%height = x(j)
					case(11)
						params%dtau2 = x(j)
					case(12)
						params%vmacro2 = x(j)
					case(13)
						params%bgauss2 = x(j)
					case(14)
						params%thetabd2 = x(j)
					case(15)
						params%chibd2 = x(j)
					case(16)
						params%vdopp2 = x(j)
					case(17)
						params%ff = x(j)
				end select
				j = j + 1
			endif
		enddo

		allocate(onum(fixed%no))
		
		onum = -1.d8 * observation%wl / fixed%wl**2
		
		if (params%nslabs == 1) then
			dnum = params%vdopp*1.d5 / (fixed%wl*1.d-8*PC)
			va = params%vmacro*1.d5 / (fixed%wl*1.d-8*PC)
			adamp = params%damping

			onum0 = 0.d0
			onum1 = -0.09d0*1.d-8 / (fixed%wl*1.d-8)**2
			onum2 = -1.26d0*1.d-8 / (fixed%wl*1.d-8)**2			
			inversion%stokes_unperturbed(0,:) = profile(adamp,(onum0-onum-va)/dnum) + &
				0.6*profile(adamp,(onum1-onum-va)/dnum) + &
				0.2*profile(adamp,(onum2-onum-va)/dnum)

			if (fixed%Stokes_incident(0) == 0) then
				inversion%stokes_unperturbed(0,:) = 1.d0 - exp(-params%dtau * inversion%stokes_unperturbed(0,:))
			else
				inversion%stokes_unperturbed(0,:) = exp(-params%dtau * inversion%stokes_unperturbed(0,:))
			endif

			print *, 'x,f = ' ,x, f
			
			f = compute_chisq(observation,inversion)

		endif

		deallocate(onum)

		call print_parameters(params,'      -Parameters : ',.TRUE.)

 		print *, 'chi^2 : ', f
		write(25,FMT='(11(E15.5,2X))') x, f

	end subroutine simplified_transfer

end module marquardt
