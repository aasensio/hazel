module io
use vars
use maths
implicit none
contains

!------------------------------------------------------------
! Read the main configuration file
!------------------------------------------------------------
	subroutine read_config
	
		open(unit=12,file='config_inversion.dat',action='read',status='old')
		call lb(12,5)
		read(12,*) input_model_file
		call lb(12,2)
		read(12,*) input_experiment
		call lb(12,2)
		read(12,*) direct_ranges
		call lb(12,2)
		read(12,*) output_rho_vertical_upper
		call lb(12,2)
		read(12,*) output_rho_vertical_lower
		call lb(12,2)
		read(12,*) output_rho_magnetic_upper
		call lb(12,2)
		read(12,*) output_rho_magnetic_lower
		call lb(12,2)
		read(12,*) output_rtcoef
		call lb(12,2)
		read(12,*) output_rtcoef_zeeman
		call lb(12,2)
		read(12,*) input_observed_profiles
		call lb(12,2)
		read(12,*) output_inverted_profiles
		call lb(12,2)
		read(12,*) output_final_parameters
		call lb(12,2)
		read(12,*) input_inverted_parameters
		call lb(12,2)
		read(12,*) verbose_mode
		call lb(12,2)
		read(12,*) linear_solver
		call lb(12,2)
		read(12,*) synthesis_mode
		call lb(12,2)
		read(12,*) working_mode		
		close(12)
		
	end subroutine read_config
	

!------------------------------------------------------------
! Read the file with the data for the experiment
!------------------------------------------------------------
	subroutine read_experiment(file,in_params,in_fixed)
	character(len=120) :: file
	type(variable_parameters) :: in_params
	type(fixed_parameters) :: in_fixed
	real(kind=8) :: delta, height
	integer :: j
	
		open(unit=12,file=file,action='read',status='old')

! Read if we want to include the stimulated emission
		call lb(12,5)
		read(12,*) isti
		if (isti == 1 .and. verbose_mode == 1) then
			print *, 'Including stimulated emission...'
		endif
	
! Read if we want to include the magnetic field
		call lb(12,2)
		read(12,*) imag
		if (imag == 1 .and. verbose_mode == 1) then
			print *, 'Including magnetic field ...'
		endif
		
! Read if we want to include the depolarizing collisions
		call lb(12,2)
		read(12,*) idep
		if (idep == 1 .and. verbose_mode == 1) then
			print *, 'Including depolarizing collisions ...'
		endif

! Value of the delta parameter for the depolarizing collisions		
		call lb(12,2)
		read(12,*) in_params%delta_collision
		
! Use Paschen-Back calculation or linear Zeeman effect
		call lb(12,2)
		read(12,*) use_paschen_back
		if (verbose_mode == 1) then
			if (use_paschen_back == 1) then
				print *, 'Taking into account the Paschen-Back effect...'
			else
				print *, 'Not taking into account the Paschen-Back effect...'
			endif
		endif
		
! Values characterizing the magnetic field
		call lb(12,2)
		read(12,*) in_params%bgauss, in_params%thetabd, in_params%chibd
		if (imag == 1 .and. verbose_mode == 1) then
			print *, 'Field properties'
			write(*,FMT='(4X,A,F8.3,3X,A,F6.2,3X,A,F6.2)') 'B : ', in_params%bgauss, 'thetab : ', in_params%thetabd, 'chib : ',&
				in_params%chibd
		endif
		
! Values of the height for the calculation of the radiation anisotropy
		call lb(12,2)
		read(12,*) in_params%height
				
! Number of slabs
		call lb(12,2)
		read(12,*) in_params%nslabs
		
! Value of the optical depth in the maximum (SLAB) or the strength of the line (ME)
		call lb(12,2)
		read(12,*) in_params%dtau
		
! Value of the gradient of the source function (ME)
		call lb(12,2)
		read(12,*) in_params%beta
		
! Incident Stokes parameters
		call lb(12,2)
		read(12,*) (in_fixed%Stokes_incident(j),j=0,3)
		
! Transition where to compute the emission
		call lb(12,2)
		read(12,*) in_fixed%nemiss
		if (verbose_mode == 1) then			
			write(*,FMT='(A,3X,I3)') ' Calculating emission from line : ', in_fixed%nemiss
		endif
		
! Use atomic polarization
		call lb(12,2)
		read(12,*) in_fixed%use_atomic_pol
		if (verbose_mode == 1) then
			if (in_fixed%use_atomic_pol == 1) then
				write(*,FMT='(A)') ' Using atomic polarization'
			else
				write(*,FMT='(A)') ' Not using atomic polarization'
			endif
		endif		
		
! Observation angle with respect to the vertical
		call lb(12,2)
		read(12,*) in_fixed%thetad, in_fixed%chid, in_fixed%gammad
		if (verbose_mode == 1) then
			print *, 'LOS reference system'
			write(*,FMT='(4X,A,F6.2,3X,A,F6.2,3X,A,F6.2)') 'thetad : ', in_fixed%thetad, 'chid : ', in_fixed%chid, 'gammad : ', &
				in_fixed%gammad
		endif
		
! Wavelength axis
		call lb(12,2)
		read(12,*) in_fixed%omin, in_fixed%omax, in_fixed%no
		if (verbose_mode == 1) then
			print *, 'Wavelength axis'
			write(*,FMT='(4X,A,F8.3,3X,A,F8.3,3X,A,I5)') 'Min : ', in_fixed%omin, 'Max : ', in_fixed%omax, 'N. steps : ', &
				in_fixed%no
		endif
		
! Line wavelength and Doppler width
		call lb(12,2)
		read(12,*) in_fixed%wl, in_params%vdopp, in_params%damping
		if (verbose_mode == 1) then
			write(*,FMT='(4X,A,F10.4,3X,A,F6.3)') 'Wavelength [A]: ', in_fixed%wl, 'Doppler vel. [km/s] : ', &
				in_params%vdopp, 'Damping : ', in_params%damping
		endif
		
! Macroscopic velocity
		call lb(12,2)
		read(12,*) in_params%vmacro
		
! Use magneto-optical terms
		call lb(12,2)
		read(12,*) use_mag_opt_RT
		
! Use stimulated emission terms in the absorption coefficients
		call lb(12,2)
		read(12,*) use_stim_emission_RT
		
	
		close(12)
		
		in_fixed%total_forward_modeling = 0
		
! Take into account the observation angle (scattering angle) to calculate the height from the apparent height
		if (in_params%height < 0.d0) then
			delta = abs(90.d0 - in_fixed%thetad)
			height = (RSUN + in_params%height) / (cos(delta*PI/180.d0)) - RSUN
			write(*,FMT='(A,F7.3)') 'The height taking into account the scattering angle is : ', height
		endif
		
	
	end subroutine read_experiment
	
!------------------------------------------------------------
! Read the main configuration file
!------------------------------------------------------------
	subroutine read_model_file(file)
	character(len=120) :: file
	integer :: i, j, nJ, n, i1, i2, ir, jmax2, jmin2, jp2, j2, kmin, kmax, k, q
	real(kind=8) :: rnujjp

		open(unit=12,file=file,action='read',status='old')
		
		read(12,*) is2
		read(12,*) n_terms		
		
		allocate(lsto2(n_terms))
		
		nrhos = 0
		
		open(unit=13,file='mterm.tab',action='write',status='replace')		
		
		jlimit2 = 0
		do i = 1, n_terms
			read(12,*) n, lsto2(i)
			jmin2 = abs(is2-lsto2(i))
			jmax2 = is2 + lsto2(i)
			do j2 = jmin2, jmax2, 2
				read(12,*) 
			enddo
			if (verbose_mode == 1) then
				print *, 'Level ', i
				print *, '         Jmin : ', jmin2/2.d0
				print *, '         Jmax : ', jmax2/2.d0
			endif
			jlimit2 = max(jmax2,jlimit2)
		enddo
		
		allocate(energy(n_terms,0:jlimit2))
		
		rewind(12)
		call lb(12,2)
		file_pointer = 2
		
		do i = 1, n_terms
			read(12,*) n, lsto2(i)
			file_pointer = file_pointer + 1
			jmin2 = abs(is2-lsto2(i))
			jmax2 = is2 + lsto2(i)
			
			do j2 = jmin2, jmax2, 2
				read(12,*) energy(i,j2)
				file_pointer = file_pointer + 1
			enddo
			
			do j2 = jmin2, jmax2, 2
				do jp2 = j2, jmax2, 2
					rnujjp = PC * (energy(i,j2)-energy(i,jp2))
					kmin = (jp2-j2) / 2
					kmax = (jp2+j2) / 2
					do k = kmin, kmax
						do q = -k, k
							if (j2 == jp2 .and. q < 0) then
							else
								do ir = 1, 2
									if (j2 == jp2 .and. q == 0 .and. ir == 2) then
									else
										nrhos = nrhos + 1
									endif
								enddo
							endif
						enddo
					enddo
				enddo
			enddo
		enddo
		
		rewind(12)
		call lb(12,2)


		allocate(ntab(nrhos))
		allocate(j2tab(nrhos))
		allocate(jp2tab(nrhos))
		allocate(ktab(nrhos))
		allocate(qtab(nrhos))
		allocate(irtab(nrhos))
		allocate(rnutab(nrhos))
		
		nrhos = 0
		
		do i = 1, n_terms
			read(12,*) n, lsto2(i)
			jmin2 = abs(is2-lsto2(i))
			jmax2 = is2 + lsto2(i)
			
			do j2 = jmin2, jmax2, 2
				read(12,*) energy(i,j2)
			enddo
			
			do j2 = jmin2, jmax2, 2
				do jp2 = j2, jmax2, 2
					rnujjp = PC * (energy(i,j2)-energy(i,jp2))
					kmin = (jp2-j2) / 2
					kmax = (jp2+j2) / 2
					do k = kmin, kmax
						do q = -k, k
							if (j2 == jp2 .and. q < 0) then
							else
								do ir = 1, 2
									if (j2 == jp2 .and. q == 0 .and. ir == 2) then
									else
										nrhos = nrhos + 1
										ntab(nrhos) = n
										j2tab(nrhos) = j2
										jp2tab(nrhos) = jp2
										ktab(nrhos) = k
										qtab(nrhos) = q
										irtab(nrhos) = ir
										rnutab(nrhos) = rnujjp
										write(13,FMT='(7(1X,I5),6X,E12.5)') nrhos, n, j2, jp2, k, q, ir, rnujjp
									endif
								enddo
							endif
						enddo
					enddo
				enddo
			enddo
		enddo

		close(12)
		
		if (verbose_mode == 1) then
			print *, 'Number of unknowns : ', nrhos
		endif
		
	end subroutine read_model_file
	

!------------------------------------------------------------
! Clean all the allocated memory
!------------------------------------------------------------	
	subroutine clean
		
		deallocate(lsto2)
		deallocate(energy)		
		deallocate(ntab)
		deallocate(ktab)
		deallocate(qtab)
		deallocate(irtab)
		deallocate(j2tab)
		deallocate(jp2tab)
		deallocate(rnutab)		
		
	end subroutine clean

!------------------------------------------------------------
! Number of columns in a line
!------------------------------------------------------------
	function ncolumns(string)
	integer :: ncolumns
	character(len=*) :: string
	character(len=1) :: ch, chprev
	integer :: i, l, blanks
		l = len(trim(string))
		blanks = 1
		chprev = string(1:1)
		do i = 2, l
			ch = string(i:i)			
			if (ch /= ' ' .and. chprev == ' ') then
				blanks = blanks + 1
			endif			
			chprev = ch
		enddo
		ncolumns = blanks - 1
	end function ncolumns


!------------------------------------------------------------
! Read the file with the data for the experiment
!------------------------------------------------------------
	subroutine read_observation(file,in_fixed,in_observation)
	character(len=120) :: file
	character(len=320) :: column
	type(fixed_parameters) :: in_fixed
	type(type_observation) :: in_observation
	integer :: i, j, format_file
	real(kind=8) :: SI_max, SQ_max, SU_max, SV_max, Noise, MinAmpli
		
		open(unit=12,file=file,action='read',status='old')
		
		read(12,*) observation%n		
		
		allocate(observation%wl(observation%n))
		allocate(observation%stokes(0:3,observation%n))
		allocate(observation%sigma(0:3,observation%n))

! Now detect the file type (4 columns for a simple estimation of sigma from the amplitude)
! and 8 columns for explicitly giving the value of sigma for each wavelength)
		read(12,FMT='(A)') column
		format_file = ncolumns(column)
		close(12)

! Re-read the first line of the file
		open(unit=12,file=file,action='read',status='old')
		
		read(12,*) observation%n

		if (format_file == 5) then
				write(*,*) 'Using noise levels to equalize the amplitudes of all Stokes parameters'
			do i = 1, observation%n
				read(12,*) observation%wl(i), (observation%stokes(j,i),j=0,3)
			enddo
		endif

		if (format_file == 9) then
				write(*,*) 'Reading noise levels from file'
			do i = 1, observation%n
				read(12,*) observation%wl(i), (observation%stokes(j,i),j=0,3), (observation%sigma(j,i),j=0,3)
			enddo
		endif
		
		close(12)
		
! Give the same weight to all the Stokes parameters
		if (format_file == 5) then
			
			SI_max = maxval(abs(observation%stokes(0,:)))
			SQ_max = maxval(abs(observation%stokes(1,:)))
			SU_max = maxval(abs(observation%stokes(2,:)))
			SV_max = maxval(abs(observation%stokes(3,:)))
		
			observation%sigma(0,:) = SI_max
			observation%sigma(1,:) = SQ_max
			observation%sigma(2,:) = SU_max
			observation%sigma(3,:) = SV_max
		
			Noise = 1.d-4
			if (SI_max < 1.d-20) SI_max = 1.d0
			if (SQ_max < Noise*20) SQ_max = SI_max / 50.d0
			if (SU_max < Noise*20) SU_max = SI_max / 50.d0
			if (SV_max < Noise*20) SV_max = SI_max / 20.d0
		
			MinAmpli = min(SI_max,SQ_max,SU_max,SV_max)

			observation%sigma = observation%sigma / sqrt(1.d3)
		endif
		
!		observation%sigma(0,:) = Noise * SI_max / MinAmpli
!		observation%sigma(1,:) = Noise * SQ_max / MinAmpli
!		observation%sigma(2,:) = Noise * SU_max / MinAmpli
!		observation%sigma(3,:) = Noise * SV_max / MinAmpli
		
!		observation%sigma = observation%sigma / sum(observation%sigma) * (4*observation%n) * Noise
		
! In the inversion mode, the wavelength axis has to be the same as that introduced in the observation
		if (working_mode == 1) then
			in_fixed%no = observation%n
		endif
				
				
	end subroutine read_observation
	
!------------------------------------------------------------
! Read the file with the parameters to invert
!------------------------------------------------------------
	subroutine read_parameters_to_invert(file,in_params,in_inversion)
	character(len=120) :: file
	type(variable_parameters) :: in_params
	type(type_inversion) :: in_inversion
	integer :: i, j
		
		in_params%inverted = 0
		
! Number of parameters of the inversion		
		in_params%n_total = 10
		
		allocate(in_params%inverted(in_params%n_total))
		
		open(unit=12,file=file,action='read',status='old')
		call lb(12,3)
		call lb(12,2)
		read(12,*) in_inversion%iter_max
		call lb(12,2)
		read(12,*) in_inversion%n_cycles
		
		allocate(in_inversion%cycles(in_params%n_total,in_inversion%n_cycles))
		allocate(in_inversion%stokes_weights(0:3,in_inversion%n_cycles))
		
		do i = 1, in_params%n_total
			call lb(12,2)
			read(12,*) (in_inversion%cycles(i,j),j=1,in_inversion%n_cycles)
		enddo
		
		do i = 0, 3
			call lb(12,2)
			read(12,*) (in_inversion%stokes_weights(i,j),j=1,in_inversion%n_cycles)
		enddo
		
		call lb(12,2)
		allocate(in_inversion%algorithm(in_inversion%n_cycles))
		read(12,*) (in_inversion%algorithm(j),j=1,in_inversion%n_cycles)
				
		close(12)
		
		in_params%inverted = in_inversion%cycles(:,1)
	
		in_params%n_inverted = sum(params%inverted)
				
	end subroutine read_parameters_to_invert
	
!------------------------------------------------------------
! Write a file with the final parameters so that it can be used for initializing again the inversion
!------------------------------------------------------------
	subroutine write_experiment(in_params,in_fixed)
	character(len=120) :: tmp
	type(variable_parameters) :: in_params
	type(fixed_parameters) :: in_fixed
	integer :: i
	
		open(unit=13,file=output_final_parameters,action='write',status='replace')

		write(13,FMT='(A)') '***********************************************'
		write(13,FMT='(A)') '* File defining the specific experiment to solve'
		write(13,FMT='(A)') '***********************************************'
		write(13,*)
		write(13,FMT='(A)') '# Include stimulated emission (0-> no, 1-> yes)'
		write(13,FMT='(I1)') isti
		write(13,*)
		
		write(13,FMT='(A)') '# Include magnetic field (0-> no, 1-> yes)'
		write(13,FMT='(I1)') imag
		write(13,*)

		write(13,FMT='(A)') '# Include depolarization rates (0-> no, 1-> yes)'
		write(13,FMT='(I1)') idep
		write(13,*)

		write(13,FMT='(A)') '# Value of delta if depolarization rates are included (not used if the previous value is 0)'
		write(13,FMT='(E9.3)') in_params%delta_collision
		write(13,*)
		
		write(13,FMT='(A)') '# Include Paschen-Back effect (0-> no, 1-> yes)'
		write(13,FMT='(I1)') use_paschen_back
		write(13,*)
		
		write(13,FMT='(A)') '# Magnetic field strength [G], thetaB [degrees], chiB [degrees]'
		write(13,FMT='(F9.4,3X,F7.2,3X,F7.2)') in_params%bgauss, in_params%thetabd, in_params%chibd
		write(13,*)
		
		write(13,FMT='(A)') '# Apparent height of the He I atoms in arcsec'
		write(13,FMT='(F7.3)') in_params%height
		write(13,*)
		
		write(13,FMT='(A)') '# Optical depth of the slab in the maximum of I (slab) or strength of the line (ME)'
		write(13,FMT='(F8.4)') in_params%dtau
		write(13,*)
		
		write(13,FMT='(A)') '# Source function gradient (only ME)'
		write(13,FMT='(F8.4)') in_params%beta
		write(13,*)
		
		write(13,FMT='(A)') '# Boundary Stokes parameters (I0,Q0,U0,V0)'
		write(13,FMT='(4(E9.3,2X))') (in_fixed%Stokes_incident(i),i=0,3)
		write(13,*)
		
		write(13,FMT='(A)') '# Transition where to compute the emission'
		write(13,FMT='(I1)') in_fixed%nemiss
		write(13,*)
		
		write(13,FMT='(A)') '# Use atomic polarization? (0-> no, 1-> yes)'
		write(13,FMT='(I1)') in_fixed%use_atomic_pol
		write(13,*)
		
		write(13,FMT='(A)') '# Observation angle with respect to the vertical theta,chi,gamma [degrees]'
		write(13,FMT='(3(F7.2,2X))') in_fixed%thetad, in_fixed%chid, in_fixed%gammad
		write(13,*)
		
		write(13,FMT='(A)') '# Wavelength axis: minimum, maximum and number of grid points'
		write(13,*) in_fixed%omin, in_fixed%omax, in_fixed%no
		write(13,FMT='(3(F7.2,2X))')
		
		write(13,FMT='(A)') '# Line wavelength [A], Doppler velocity [km/s] and damping [a]'
		write(13,FMT='(F10.4,3X,F6.3,3X,F7.4)') in_fixed%wl, in_params%vdopp, in_params%damping
		write(13,*)
		
		write(13,FMT='(A)') '# Macroscopic velocity [km/s] (>0 is a redshift)'
		write(13,FMT='(F7.3)') in_params%vmacro
		write(13,*)
		
		write(13,FMT='(A)') '# Include magneto-optical effects in the RT'
		write(13,FMT='(I1)') use_mag_opt_RT
		write(13,*)
		
		write(13,FMT='(A)') '# Include stimulated emission in the RT'
		write(13,FMT='(I1)') use_stim_emission_RT
		write(13,*)
		
		close(13)
	
	end subroutine write_experiment	
	
!------------------------------------------------------------
! Save the inverted profiles
!------------------------------------------------------------
	subroutine write_final_profiles(file,in_observation,in_inversion)
	character(len=120) :: file
	type(type_observation) :: in_observation
	type(type_inversion) :: in_inversion
	integer :: i, j
	
		open(unit=12,file=file,action='write',status='replace')
		
		write(12,*) in_observation%n, in_inversion%chisq
		
		do i = 1, in_observation%n
			write(12,FMT='(5(1X,E15.8))') in_observation%wl(i), &
				(in_inversion%stokes_unperturbed(j,i),j=0,3)
		enddo
		
		close(12)
	end subroutine write_final_profiles	

!------------------------------------------------------------
! Print the parameters depending on the synthesis mode
!------------------------------------------------------------
	subroutine print_parameters(params,text,header)
	type(variable_parameters) :: params
	character(len=20) :: text
	logical :: header
	
! EMISSION
		if (synthesis_mode == 0) then					
			if (params%inverted(6) == 1) then
				if (header) write(*,FMT='(A)') '                        B        thetaB     chiB        v_th       D^(2)    vmacro      a        h'
				write(*,FMT='(A,F9.4,2X,F9.4,2X,F9.4,2X,E8.2,2X,F9.4,2X,F9.4,2X,F9.4)') text, params%bgauss, params%thetabd, &
					params%chibd, 10.d0**params%delta_collision, params%vmacro, params%damping, params%height
			else
				if (header) write(*,FMT='(A)') '                        B        thetaB     chiB        v_th       vmacro       a       h'
				write(*,FMT='(A,F9.4,2X,F9.4,2X,F9.4,2X,F9.4,2X,F9.4,2X,F9.4,2X,F9.4)') text, params%bgauss, params%thetabd, &
					params%chibd, params%vdopp, params%vmacro, params%damping, params%height
			endif				
		endif
! CONSTANT SLAB 
		if (synthesis_mode == 1 .or. synthesis_mode == 3 .or. synthesis_mode == 4 .or. &
			synthesis_mode == 5) then
			if (params%inverted(6) == 1) then
				if (header) write(*,FMT='(A)') '                        B        thetaB     chiB        v_th         tau       D^(2)     vmacro       a       h'
				write(*,FMT='(A,F9.4,2X,F9.4,2X,F9.4,2X,F9.4,2X,F9.4,3X,E9.2,2X,F9.4,2X,F9.4,2X,F9.4)') text, params%bgauss, params%thetabd, &
					params%chibd, params%vdopp, params%dtau, 10.d0**params%delta_collision, params%vmacro, params%damping, params%height
			else
				if (header) write(*,FMT='(A)') '                        B        thetaB     chiB        v_th         tau       vmacro       a        h'
				write(*,FMT='(A,F9.4,2X,F9.4,2X,F9.4,2X,F9.4,2X,F9.4,2X,F9.4,2X,F9.4,2X,F9.4)') text, params%bgauss, params%thetabd, &
					params%chibd, params%vdopp, params%dtau, params%vmacro, params%damping, params%height
			endif				
		endif
! MILNE-EDDINGTON
		if (synthesis_mode == 2) then					
			if (params%inverted(6) == 1) then
				if (header) write(*,FMT='(A)') '                        B        thetaB     chiB        v_th         strength        beta         D^(2)     vmacro       a      h'
				write(*,FMT='(A,F9.4,2X,F9.4,2X,F9.4,2X,F9.4,2X,F9.4,3X,E9.2,2X,F9.4,2X,F9.4,2X,F9.4,2X,F9.4)') text, params%bgauss, &
					params%thetabd, params%chibd, params%vdopp, params%dtau, params%beta, 10.d0**params%delta_collision, params%vmacro, params%damping, params%height
			else
				if (header) write(*,FMT='(A)') '                        B        thetaB     chiB        v_th     strength     beta       vmacro       a       h'
				write(*,FMT='(A,F9.4,2X,F9.4,2X,F9.4,2X,F9.4,2X,F9.4,2X,F9.4,2X,F9.4,2X,F9.4,2X,F9.4)') text, params%bgauss, params%thetabd, &
					params%chibd, params%vdopp, params%dtau, params%beta, params%vmacro, params%damping, params%height
			endif				
		endif
	
		
	end subroutine print_parameters
end module io
