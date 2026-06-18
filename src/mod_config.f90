!=====================================================
! The module contains the configuration values obtained
! from ../input.nml
! Update  : 03/06/2026
! Authors: R. Sandez and M. Cécere
!=====================================================

module mod_config
  use mod_precision
  implicit none

  type :: config_t
     character(len=256) :: field_source = "file"
     character(len=256) :: file_type    = "txt"
     character(len=256) :: field_file   = "fields_uniform.txt"

     integer :: nx = 32, ny = 32, nz = 32
     real(dp) :: xmin = 0.0_dp, xmax = 1.0e5_dp ! 1km en cm
     real(dp) :: ymin = 0.0_dp, ymax = 1.0e5_dp
     real(dp) :: zmin = 0.0_dp, zmax = 1.0e5_dp

     real(dp) :: bx0 = 0.0_dp, by0 = 0.0_dp, bz0 = 1.0_dp ! campo magnético uniforme [gauss]
     real(dp) :: ex0 = 0.1_dp, ey0 = 0.0_dp, ez0 = 0.0_dp ! campo eléctrico uniforme [statvolt/cm]

     integer :: npart = 10000, nsteps = 200
     real(dp) :: dt = 1.0e-3_dp
     real(dp) :: temperature = 1.0_dp 
     real(dp) :: q_particle = 9.109e-28_dp ! carga de e en CGS
     real(dp) :: m_particle = 4.806e-10_dp ! masa de e en CGS
     real(dp) :: kb = 1.381e-16_dp ! constante de Boltzmann en CGS

     logical :: use_exb_drift = .true.
     logical :: use_gradb_drift = .true.
     logical :: use_curvature_drift = .true.
     logical :: use_mirror_force = .true.
     real(dp) :: exb_factor = 2.998e11_dp ! factor de escala para el drift ExB

     character(len=256) :: boundary_condition = "periodic"
     character(len=256) :: output_dir = "output"
     character(len=256) :: output_prefix = "gc_txt"
     integer :: nbins_energy = 80
     integer :: random_seed = 12345
  end type config_t

contains

  subroutine read_config(filename, cfg)
    character(len=*), intent(in) :: filename
    type(config_t), intent(inout) :: cfg

    character(len=256) :: field_source, file_type, field_file
    integer :: nx, ny, nz
    real(dp) :: xmin, xmax, ymin, ymax, zmin, zmax
    real(dp) :: bx0, by0, bz0, ex0, ey0, ez0
    integer :: npart, nsteps
    real(dp) :: dt, temperature, q_particle, m_particle, kb
    logical :: use_exb_drift, use_gradb_drift, use_curvature_drift, use_mirror_force
    real(dp) :: exb_factor
    character(len=256) :: boundary_condition, output_dir, output_prefix
    integer :: nbins_energy, random_seed
    integer :: unit, ios
    logical :: exists

    namelist /run_config/ field_source, file_type, field_file, &
         nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, &
         bx0, by0, bz0, ex0, ey0, ez0, &
         npart, nsteps, dt, temperature, q_particle, m_particle, kb, &
         use_exb_drift, use_gradb_drift, use_curvature_drift, use_mirror_force, &
         exb_factor, boundary_condition, output_dir, output_prefix, &
         nbins_energy, random_seed

    field_source = cfg%field_source
    file_type = cfg%file_type
    field_file = cfg%field_file
    nx = cfg%nx; ny = cfg%ny; nz = cfg%nz
    xmin = cfg%xmin; xmax = cfg%xmax
    ymin = cfg%ymin; ymax = cfg%ymax
    zmin = cfg%zmin; zmax = cfg%zmax
    bx0 = cfg%bx0; by0 = cfg%by0; bz0 = cfg%bz0
    ex0 = cfg%ex0; ey0 = cfg%ey0; ez0 = cfg%ez0
    npart = cfg%npart; nsteps = cfg%nsteps; dt = cfg%dt
    temperature = cfg%temperature
    q_particle = cfg%q_particle
    m_particle = cfg%m_particle
    kb = cfg%kb
    use_exb_drift = cfg%use_exb_drift
    use_gradb_drift = cfg%use_gradb_drift
    use_curvature_drift = cfg%use_curvature_drift
    use_mirror_force = cfg%use_mirror_force
    exb_factor = cfg%exb_factor
    boundary_condition = cfg%boundary_condition
    output_dir = cfg%output_dir
    output_prefix = cfg%output_prefix
    nbins_energy = cfg%nbins_energy
    random_seed = cfg%random_seed

    inquire(file=trim(filename), exist=exists)
    if (.not. exists) then
       write(*,*) "WARNING: input file not found, using defaults: ", trim(filename)
       return
    end if

    open(newunit=unit, file=trim(filename), status="old", action="read", iostat=ios)
    if (ios /= 0) stop "ERROR: could not open input.nml"
    read(unit, nml=run_config, iostat=ios)
    if (ios /= 0) then
       write(*,*) "ERROR while reading namelist. iostat=", ios
       stop
    end if
    close(unit)

    cfg%field_source = adjustl(field_source)
    cfg%file_type = adjustl(file_type)
    cfg%field_file = adjustl(field_file)
    cfg%nx = nx; cfg%ny = ny; cfg%nz = nz
    cfg%xmin = xmin; cfg%xmax = xmax
    cfg%ymin = ymin; cfg%ymax = ymax
    cfg%zmin = zmin; cfg%zmax = zmax
    cfg%bx0 = bx0; cfg%by0 = by0; cfg%bz0 = bz0
    cfg%ex0 = ex0; cfg%ey0 = ey0; cfg%ez0 = ez0
    cfg%npart = npart; cfg%nsteps = nsteps; cfg%dt = dt
    cfg%temperature = temperature
    cfg%q_particle = q_particle
    cfg%m_particle = m_particle
    cfg%kb = kb
    cfg%use_exb_drift = use_exb_drift
    cfg%use_gradb_drift = use_gradb_drift
    cfg%use_curvature_drift = use_curvature_drift
    cfg%use_mirror_force = use_mirror_force
    cfg%exb_factor = exb_factor
    cfg%boundary_condition = adjustl(boundary_condition)
    cfg%output_dir = adjustl(output_dir)
    cfg%output_prefix = adjustl(output_prefix)
    cfg%nbins_energy = nbins_energy
    cfg%random_seed = random_seed
  end subroutine read_config

  subroutine print_config(cfg)
    type(config_t), intent(in) :: cfg
    write(*,*) "Configuration"
    write(*,*) " field_source = ", trim(cfg%field_source)
    write(*,*) " file_type    = ", trim(cfg%file_type)
    write(*,*) " field_file   = ", trim(cfg%field_file)
    write(*,*) " grid         = ", cfg%nx, cfg%ny, cfg%nz
    write(*,*) " npart        = ", cfg%npart
    write(*,*) " nsteps, dt   = ", cfg%nsteps, cfg%dt
    write(*,*) " boundary     = ", trim(cfg%boundary_condition)
  end subroutine print_config
end module mod_config
