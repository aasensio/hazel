module marquardt
use vars
use SEE
use rt_coef
use synth
use svd
use maths
use io
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
	
! Verify if the parameters are inside the boundaries and return false if any of them is outside

!		correct = .FALSE.
		
!		if (trial%bgauss >= 0.d0 .and. trial%bgauss < 4000.d0 .and. &
!			 trial%thetabd >= 0.d0 .and. trial%thetabd < 180.d0 .and. &
!			 trial%chibd >= 0.d0 .and. trial%chibd < 360.d0 .and. &
!			 trial%vdopp >= 0.d0 .and. trial%vdopp < 30.d0) then
			 
!			 correct = .TRUE.
!		endif
	
! Verify if the parameters are inside the boundaries. If any of them is outside, use the previous value
! of this parameter.
	
		if (trial%bgauss < 0.d0 .or. trial%bgauss > 4000.d0) then
			trial%bgauss = original%bgauss
		endif
		if (trial%thetabd < 0.d0 .or. trial%thetabd > 180.d0) then
			trial%thetabd = original%thetabd
		endif
		if (trial%chibd < -180.d0 .or. trial%chibd > 180.d0) then
			trial%chibd = original%chibd
		endif
		if (trial%vdopp < 0.d0 .or. trial%vdopp > 30.d0) then
			trial%vdopp = original%vdopp
		endif
		if (trial%dtau < 0.d0 .or. trial%dtau > 300.d0) then
			trial%dtau = original%dtau
		endif
		if (trial%delta_collision < 0.d0 .or. trial%delta_collision > 18.d0) then
			trial%delta_collision = original%delta_collision
		endif		
		if (trial%vmacro < -15.d0 .or. trial%vmacro > 40.d0) then
			trial%vmacro = original%vmacro
		endif
		if (trial%damping < 0.d0 .or. trial%damping > 10.d0) then
			trial%damping = original%damping
		endif
		if (trial%beta < -10.d0 .or. trial%beta > 10.d0) then
			trial%beta = original%beta
		endif
		if (trial%height < 0.d0 .or. trial%height > 100.d0) then
			trial%height = original%height
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
			write(*,FMT='(A,A)') 'Clipping ', parameters_name(which)
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
		call do_synthesis(in_params, in_fixed, in_observation, in_inversion%stokes_unperturbed)
		params_scaled = in_params		
		params_modified = in_params
		
		call compress_parameters(in_params, params_scaled)		
		params_scaled_modified = params_scaled

		do i = 1, in_params%n_total
		
			if (in_params%inverted(i) == 1) then
				write(*,FMT='(A,I2,A,A)') 'Perturbing parameter ', i, '  --- ', parameters_name(i)
		
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
	real(kind=8) :: obs, syn, sig, weight
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
	end subroutine expand_parameters	



!*************************************************************
!*************************************************************
! DIRECT METHOD
!*************************************************************
!*************************************************************

!------------------------------------------------------------
! Invert some parameters with the DIRECT method
!------------------------------------------------------------
	subroutine invert_with_direct(out_params)
	real(kind=8) :: DIReps, DIRf
	integer :: ndim, DIRmaxf, DIRmaxT
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
		
		
		allocate(upper(10))
		allocate(lower(10))
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
		
		do i = 1, 10
			call lb(24,2)
			read(24,*) lower(i), upper(i)
		enddo
		close(24)
		
		j = 1
		do i = 1, params%n_total
			if (params%inverted(i) == 1) then
				l(j) = lower(i)
				u(j) = upper(i)
				print *, 'Inverting ', parameters_name(i)
				j = j + 1
			endif
		enddo
		
		deallocate(upper)
		deallocate(lower)

		DIRalg = 1
		
		call DIRECT(fcn,DIRx,ndim,DIReps,DIRmaxf,DIRmaxT, DIRf, l, u, DIRalg, Ierror, logfile, &
			fglobal, fglper, volper, sigmaper, iidata, iisize, ddata, idsize, cdata, icsize)
			
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
				end select
				j = j + 1
			endif
		enddo
			
		close(23)
		close(25)
			
	end subroutine invert_with_direct
	
!------------------------------------------------------------
! Function that returns the the DIRECT
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
				end select
				j = j + 1
			endif
		enddo
			
		call do_synthesis(params, fixed, observation, inversion%stokes_unperturbed)
		call write_final_profiles(temporal_file, observation, inversion)
		
		f = compute_chisq(observation,inversion)
		
		call print_parameters(params,'     -Parameters : ',.TRUE.)
		
		print *, 'chi^2 : ', f
		write(25,FMT='(11(E15.5,2X))') x, f
		
		flag = 0

	end subroutine fcn

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
				end select
				j = j + 1
			endif
		enddo
		
		call do_synthesis(params, fixed, observation, inversion%stokes_unperturbed)
		
		f = compute_chisq(observation,inversion)

		call print_parameters(params,'     -Parameters : ',.TRUE.)
						
		chisq_pikaia = 1.d0 / f		
		print *, 'chi^2 : ', f, chisq_pikaia
		
	end function chisq_pikaia

end module marquardt
