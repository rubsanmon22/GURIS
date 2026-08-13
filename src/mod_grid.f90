!=====================================================
! The module contains the grid and values on the grid
! x, vecB, vecE
! Bmag, vecb, gradB and kappa
! Update  : 03/06/2026
! Authors: R. Sandez and M. Cécere
!=====================================================

module mod_grid
  use mod_precision
  implicit none

  type :: grid_t
     integer :: nx = 0, ny = 0, nz = 0
     real(dp), allocatable :: x(:), y(:), z(:)
     real(dp), allocatable :: bx(:,:,:), by(:,:,:), bz(:,:,:)
     real(dp), allocatable :: ex(:,:,:), ey(:,:,:), ez(:,:,:)
     real(dp), allocatable :: bmag(:,:,:)
     real(dp), allocatable :: bhx(:,:,:), bhy(:,:,:), bhz(:,:,:)
     real(dp), allocatable :: gradb_x(:,:,:), gradb_y(:,:,:), gradb_z(:,:,:)
     real(dp), allocatable :: kappa_x(:,:,:), kappa_y(:,:,:), kappa_z(:,:,:)
     real(dp) :: dx = 0.0_dp, dy = 0.0_dp, dz = 0.0_dp
     real(dp) :: xmin = 0.0_dp, xmax = 0.0_dp
     real(dp) :: ymin = 0.0_dp, ymax = 0.0_dp
     real(dp) :: zmin = 0.0_dp, zmax = 0.0_dp
  end type grid_t

contains

  subroutine allocate_base_grid(grid, nx, ny, nz)
    type(grid_t), intent(inout) :: grid
    integer, intent(in) :: nx, ny, nz
    grid%nx = nx; grid%ny = ny; grid%nz = nz
    allocate(grid%x(nx), grid%y(ny), grid%z(nz))
    allocate(grid%bx(nx,ny,nz), grid%by(nx,ny,nz), grid%bz(nx,ny,nz))
    allocate(grid%ex(nx,ny,nz), grid%ey(nx,ny,nz), grid%ez(nx,ny,nz))
  end subroutine allocate_base_grid

  subroutine allocate_derived_grid(grid)
    type(grid_t), intent(inout) :: grid
    integer :: nx, ny, nz
    nx = grid%nx; ny = grid%ny; nz = grid%nz
    allocate(grid%bmag(nx,ny,nz))
    allocate(grid%bhx(nx,ny,nz), grid%bhy(nx,ny,nz), grid%bhz(nx,ny,nz))
    allocate(grid%gradb_x(nx,ny,nz), grid%gradb_y(nx,ny,nz), grid%gradb_z(nx,ny,nz))
    allocate(grid%kappa_x(nx,ny,nz), grid%kappa_y(nx,ny,nz), grid%kappa_z(nx,ny,nz))
  end subroutine allocate_derived_grid

  subroutine set_uniform_coordinates(grid, xmin, xmax, ymin, ymax, zmin, zmax)
    type(grid_t), intent(inout) :: grid
    real(dp), intent(in) :: xmin, xmax, ymin, ymax, zmin, zmax
    integer :: i, j, k
    grid%xmin = xmin; grid%xmax = xmax
    grid%ymin = ymin; grid%ymax = ymax
    grid%zmin = zmin; grid%zmax = zmax
    grid%dx = (xmax - xmin) / real(grid%nx - 1, dp)
    grid%dy = (ymax - ymin) / real(grid%ny - 1, dp)
    grid%dz = (zmax - zmin) / real(grid%nz - 1, dp)
    do i=1,grid%nx; grid%x(i) = xmin + real(i-1,dp)*grid%dx; end do
    do j=1,grid%ny; grid%y(j) = ymin + real(j-1,dp)*grid%dy; end do
    do k=1,grid%nz; grid%z(k) = zmin + real(k-1,dp)*grid%dz; end do
  end subroutine set_uniform_coordinates
end module mod_grid
