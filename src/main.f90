!=====================================================
! This program calculates the energy distribution of test particles
! to compile and run:
! make fields
! make FC=gfrotran
! serial: OMP_NUM_THREADS=1  ./gc_particles
! paralelo: OMP_NUM_THREADS=8 ./gc_particles
! Update  : 03/06/2026
! Authors: R. Sandez and M. Cécere
!=====================================================

program gc_particles
  use mod_config
  use mod_grid
  use mod_fields_analytic
  use mod_reader
  use mod_derived_fields
  use mod_random
  use mod_particles
  use mod_evolution_guiding_center
  use mod_diagnostics
  use mod_output
  use omp_lib, only: omp_get_wtime, omp_get_max_threads

  implicit none

  type(config_t)    :: cfg
  type(grid_t)      :: grid
  type(particles_t) :: part
  integer           :: istep, ip, nactive, progress_stride
  real(dp)          :: t_start, t_end, t_evolve, time
  character (len=256) :: input_file
  integer :: traj_unit
  

  input_file = 'input.nml'
  if (command_argument_count() >= 1) then
     call get_command_argument(1, input_file)
  end if

  write(*,*) 'Reading configuration from: ', trim(input_file)
  call read_config(trim(input_file), cfg)
  call print_config(cfg)

  select case (trim(cfg%field_source))
  case ('analytic')
     call build_analytic_fields(cfg, grid)
  case ('file')
     call read_fields_from_file(cfg, grid)
  case default
     write(*,*) 'ERROR: unsupported field_source: ', trim(cfg%field_source)
     stop
  end select

  write(*,*) 'Computing derived fields...'
  call compute_derived_fields(grid)

  call init_random_seed(cfg%random_seed)
  write(*,*) 'Initializing particles...'
  call initialize_particles(cfg, grid, part)

  write(*,*) 'Computing initial energies...'
  call compute_energies(cfg, grid, part, initial=.true.)

  write(*,*) 'Evolving particles...'
  progress_stride=max(1,cfg%nsteps/10)

  write(*,*) "OpenMP max threads = ", omp_get_max_threads()

  ! Abrir el archivo de trayectoria si se requiere
  if (cfg%write_trajectories) then
      call open_trajectory_file(cfg, traj_unit)
      time = 0.0_dp
      call write_trajectory(cfg, traj_unit, 0, 0.0_dp, part)
  else
      traj_unit = -1 ! No escribo las trayectorias
  end if

  write(*,*) 'Evolving particles...'
  progress_stride = max (1, cfg%nsteps/10)
  
  t_start = omp_get_wtime()

  do istep=1,cfg%nsteps
     !$omp parallel do default(shared) private(ip) schedule(static)
     do ip=1,part%np
        call advance_guiding_center_rk2(cfg, grid, part, ip)
     end do
     !$omp end parallel do

     if (cfg%write_trajectories) then
         if (mod(istep, max(1, cfg%trajectory_stride)) == 0) then
            time = real(istep, dp) * cfg%dt
            call write_trajectory(cfg, traj_unit, istep, time, part)
         end if
      end if

     if (mod(istep, progress_stride) == 0) then
        nactive=count_active(part)
        write(*,'(A,I8,A,I8,A,I10)') ' step ',istep,' / ',cfg%nsteps,' active=',nactive
     end if
  end do

  if (traj_unit > 0) close(traj_unit)

  t_end = omp_get_wtime()
  t_evolve = t_end - t_start

  write(*,*)
  write(*,*) "Timing:"
  write(*,*) "  Evolution time [s] = ", t_evolve
  write(*,*) "  Time per step [s]  = ", t_evolve / real(cfg%nsteps, dp)
  write(*,*) "  Time per particle-step [s] = ", t_evolve / real(cfg%nsteps * cfg%npart, dp)
  write(*,*)

  write(*,*) 'Computing final energies...'
  call compute_energies(cfg, grid, part, initial=.false.)

  write(*,*) 'Writing outputs...'
  call write_outputs(cfg, part)
  nactive=count_active(part)
  write(*,*) 'Done.'
  write(*,*) 'Active particles: ', nactive, ' / ', part%np
  write(*,*) 'Output prefix: ', trim(cfg%output_dir)//'/'//trim(cfg%output_prefix)
end program gc_particles
