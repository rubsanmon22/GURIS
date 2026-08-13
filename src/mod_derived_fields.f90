!=====================================================
! This subroutine computes Bmag=bnorm, b=bhat, gradient of B, and vecb . grad vecb
! Update  : 03/06/2026
! Authors: R. Sandez and M. Cécere
!=====================================================

module mod_derived_fields
  use mod_precision
  use mod_grid

  implicit none

contains

  subroutine compute_derived_fields(grid)
    type(grid_t), intent(inout) :: grid
    integer :: i,j,k
    real(dp) :: bnorm, bfloor
    real(dp) :: dbhxdx, dbhxdy, dbhxdz
    real(dp) :: dbhydx, dbhydy, dbhydz
    real(dp) :: dbhzdx, dbhzdy, dbhzdz

    bfloor = 1.0e-30_dp
    call allocate_derived_grid(grid)

    do k=1,grid%nz; do j=1,grid%ny; do i=1,grid%nx
       ! magnetic field magnitude B -> bnorm
       bnorm = sqrt(grid%bx(i,j,k)**2 + grid%by(i,j,k)**2 + grid%bz(i,j,k)**2)
       grid%bmag(i,j,k) = bnorm
       ! unit magnetic field vector -> bh
       if (bnorm > bfloor) then
          grid%bhx(i,j,k)=grid%bx(i,j,k)/bnorm
          grid%bhy(i,j,k)=grid%by(i,j,k)/bnorm
          grid%bhz(i,j,k)=grid%bz(i,j,k)/bnorm
       else
          grid%bhx(i,j,k)=0.0_dp; grid%bhy(i,j,k)=0.0_dp; grid%bhz(i,j,k)=0.0_dp
       end if
    end do; end do; end do

    ! gradient of B
    do k=1,grid%nz; do j=1,grid%ny; do i=1,grid%nx
       grid%gradb_x(i,j,k)=derivative_x(grid%bmag,grid,i,j,k)
       grid%gradb_y(i,j,k)=derivative_y(grid%bmag,grid,i,j,k)
       grid%gradb_z(i,j,k)=derivative_z(grid%bmag,grid,i,j,k)
    end do; end do; end do

    ! kappa = vecb . grad vecb
    do k=1,grid%nz; do j=1,grid%ny; do i=1,grid%nx
       dbhxdx=derivative_x(grid%bhx,grid,i,j,k); dbhxdy=derivative_y(grid%bhx,grid,i,j,k); dbhxdz=derivative_z(grid%bhx,grid,i,j,k)
       dbhydx=derivative_x(grid%bhy,grid,i,j,k); dbhydy=derivative_y(grid%bhy,grid,i,j,k); dbhydz=derivative_z(grid%bhy,grid,i,j,k)
       dbhzdx=derivative_x(grid%bhz,grid,i,j,k); dbhzdy=derivative_y(grid%bhz,grid,i,j,k); dbhzdz=derivative_z(grid%bhz,grid,i,j,k)
       grid%kappa_x(i,j,k)=grid%bhx(i,j,k)*dbhxdx + grid%bhy(i,j,k)*dbhxdy + grid%bhz(i,j,k)*dbhxdz
       grid%kappa_y(i,j,k)=grid%bhx(i,j,k)*dbhydx + grid%bhy(i,j,k)*dbhydy + grid%bhz(i,j,k)*dbhydz
       grid%kappa_z(i,j,k)=grid%bhx(i,j,k)*dbhzdx + grid%bhy(i,j,k)*dbhzdy + grid%bhz(i,j,k)*dbhzdz
    end do; end do; end do
  end subroutine compute_derived_fields

  ! second-order centered finite difference
  function derivative_x(f,grid,i,j,k) result(df)
    real(dp), intent(in) :: f(:,:,:)
    type(grid_t), intent(in) :: grid
    integer, intent(in) :: i,j,k
    real(dp) :: df
    if (i == 1) then
       df=(f(i+1,j,k)-f(i,j,k))/grid%dx
    else if (i == grid%nx) then
       df=(f(i,j,k)-f(i-1,j,k))/grid%dx
    else
       df=(f(i+1,j,k)-f(i-1,j,k))/(2.0_dp*grid%dx)
    end if
  end function derivative_x

  function derivative_y(f,grid,i,j,k) result(df)
    real(dp), intent(in) :: f(:,:,:)
    type(grid_t), intent(in) :: grid
    integer, intent(in) :: i,j,k
    real(dp) :: df
    if (j == 1) then
       df=(f(i,j+1,k)-f(i,j,k))/grid%dy
    else if (j == grid%ny) then
       df=(f(i,j,k)-f(i,j-1,k))/grid%dy
    else
       df=(f(i,j+1,k)-f(i,j-1,k))/(2.0_dp*grid%dy)
    end if
  end function derivative_y

  function derivative_z(f,grid,i,j,k) result(df)
    real(dp), intent(in) :: f(:,:,:)
    type(grid_t), intent(in) :: grid
    integer, intent(in) :: i,j,k
    real(dp) :: df
    if (k == 1) then
       df=(f(i,j,k+1)-f(i,j,k))/grid%dz
    else if (k == grid%nz) then
       df=(f(i,j,k)-f(i,j,k-1))/grid%dz
    else
       df=(f(i,j,k+1)-f(i,j,k-1))/(2.0_dp*grid%dz)
    end if
  end function derivative_z
end module mod_derived_fields
