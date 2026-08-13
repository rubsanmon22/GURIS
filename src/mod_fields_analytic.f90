!=====================================================
! The module contains fields analytic
! Update  : 03/06/2026
! Authors: R. Sandez and M. Cécere
!=====================================================

module mod_fields_analytic
  use mod_precision
  use mod_config
  use mod_grid
  implicit none
contains
  subroutine build_analytic_fields(cfg, grid)
    type(config_t), intent(in) :: cfg
    type(grid_t), intent(out) :: grid
    integer :: i, j, k
    call allocate_base_grid(grid, cfg%nx, cfg%ny, cfg%nz)
    call set_uniform_coordinates(grid, cfg%xmin, cfg%xmax, cfg%ymin, cfg%ymax, cfg%zmin, cfg%zmax)
    do k=1,grid%nz; do j=1,grid%ny; do i=1,grid%nx
       grid%bx(i,j,k)=cfg%bx0; grid%by(i,j,k)=cfg%by0; grid%bz(i,j,k)=cfg%bz0
       grid%ex(i,j,k)=cfg%ex0; grid%ey(i,j,k)=cfg%ey0; grid%ez(i,j,k)=cfg%ez0
    end do; end do; end do
  end subroutine build_analytic_fields
end module mod_fields_analytic
