!=====================================================
! The subroutine initialize_particles calculates particles position xp, vpar and mu
! Update  : 03/06/2026
! Authors: R. Sandez and M. Cécere
!=====================================================

module mod_boundary
  use mod_precision
  use mod_config
  use mod_grid
  use mod_particles

  implicit none

contains

  subroutine apply_boundary_position(cfg,grid,x,y,z,active)
    type(config_t), intent(in)  :: cfg
    type(grid_t), intent(in)    :: grid
    real(dp), intent(inout)     :: x,y,z
    logical, intent(out)        :: active
    real(dp)                    :: lx,ly,lz

    active=.true.

    select case (trim(cfg%boundary_condition))
    case ('periodic')
       lx=grid%xmax-grid%xmin; ly=grid%ymax-grid%ymin; lz=grid%zmax-grid%zmin
       ! if the particle leaves the box, enter again
       x=grid%xmin+modulo(x-grid%xmin,lx)
       y=grid%ymin+modulo(y-grid%ymin,ly)
       z=grid%zmin+modulo(z-grid%zmin,lz)
       !if (x >= grid%xmax) x=grid%xmin
       !if (y >= grid%ymax) y=grid%ymin
       !if (z >= grid%zmax) z=grid%zmin
    case default
       if (x < grid%xmin .or. x >= grid%xmax .or. y < grid%ymin .or. y >= grid%ymax .or. z < grid%zmin .or. z >= grid%zmax) active=.false.
    end select

  end subroutine apply_boundary_position

  subroutine apply_boundary_particle(cfg,grid,part,ip)
    type(config_t), intent(in)        :: cfg
    type(grid_t), intent(in)          :: grid
    type(particles_t), intent(inout)  :: part
    integer, intent(in)               :: ip
    logical                           :: ok

    call apply_boundary_position(cfg,grid,part%x(ip),part%y(ip),part%z(ip),ok)
    if (.not. ok) part%active(ip)=.false.

  end subroutine apply_boundary_particle
end module mod_boundary
