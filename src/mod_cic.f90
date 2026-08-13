!=====================================================
! The subroutine cic_sample_field calculates the vecE, vecB, Bmag, vecb, grad B, vecb . grab vecb
! in the particle position
! The function cic_scalar uses the cloud in cell aproximation to interpolate a value in the particle position
! Update  : 03/06/2026
! Authors: R. Sandez and M. Cécere
!=====================================================

module mod_cic
  use mod_precision
  use mod_grid
  implicit none
  type :: field_sample_t
     real(dp) :: ex=0.0_dp, ey=0.0_dp, ez=0.0_dp
     real(dp) :: bx=0.0_dp, by=0.0_dp, bz=0.0_dp
     real(dp) :: bmag=0.0_dp
     real(dp) :: bhx=0.0_dp, bhy=0.0_dp, bhz=0.0_dp
     real(dp) :: gradb_x=0.0_dp, gradb_y=0.0_dp, gradb_z=0.0_dp
     real(dp) :: kappa_x=0.0_dp, kappa_y=0.0_dp, kappa_z=0.0_dp
  end type field_sample_t

contains

  subroutine cic_sample_fields(grid,xp,yp,zp,fs,inside)
    type(grid_t), intent(in)          :: grid
    real(dp), intent(in)              :: xp, yp, zp
    type(field_sample_t), intent(out) :: fs
    logical, intent(out)              :: inside
    integer                           :: i,j,k
    real(dp)                          :: fx,fy,fz

    inside = .true.
    if (xp < grid%xmin .or. xp >= grid%xmax .or. yp < grid%ymin .or. yp >= grid%ymax .or. zp < grid%zmin .or. zp >= grid%zmax) then
       inside=.false.; fs=field_sample_t(); return
    end if

    ! finding cell indexes for the particle inside the grid
    i=int((xp-grid%xmin)/grid%dx)+1
    j=int((yp-grid%ymin)/grid%dy)+1
    k=int((zp-grid%zmin)/grid%dz)+1
    if (i < 1 .or. i >= grid%nx .or. j < 1 .or. j >= grid%ny .or. k < 1 .or. k >= grid%nz) then
       inside=.false.; fs=field_sample_t(); return
    end if

    ! normalised length outside the cell
    fx=(xp-grid%x(i))/grid%dx; fy=(yp-grid%y(j))/grid%dy; fz=(zp-grid%z(k))/grid%dz
    ! interpolate vecE
    fs%ex=cic_scalar(grid%ex,i,j,k,fx,fy,fz); fs%ey=cic_scalar(grid%ey,i,j,k,fx,fy,fz); fs%ez=cic_scalar(grid%ez,i,j,k,fx,fy,fz)
    ! interpolate vecB
    fs%bx=cic_scalar(grid%bx,i,j,k,fx,fy,fz); fs%by=cic_scalar(grid%by,i,j,k,fx,fy,fz); fs%bz=cic_scalar(grid%bz,i,j,k,fx,fy,fz)
    ! interpolate Bmag
    fs%bmag=cic_scalar(grid%bmag,i,j,k,fx,fy,fz)
    ! interpolate vecb
    fs%bhx=cic_scalar(grid%bhx,i,j,k,fx,fy,fz); fs%bhy=cic_scalar(grid%bhy,i,j,k,fx,fy,fz); fs%bhz=cic_scalar(grid%bhz,i,j,k,fx,fy,fz)
    ! interpolate grad B
    fs%gradb_x=cic_scalar(grid%gradb_x,i,j,k,fx,fy,fz); fs%gradb_y=cic_scalar(grid%gradb_y,i,j,k,fx,fy,fz); fs%gradb_z=cic_scalar(grid%gradb_z,i,j,k,fx,fy,fz)
    ! interpolate vecb . grab vecb
    fs%kappa_x=cic_scalar(grid%kappa_x,i,j,k,fx,fy,fz); fs%kappa_y=cic_scalar(grid%kappa_y,i,j,k,fx,fy,fz); fs%kappa_z=cic_scalar(grid%kappa_z,i,j,k,fx,fy,fz)
  end subroutine cic_sample_fields

  function cic_scalar(f,i,j,k,wx,wy,wz) result(val)
    real(dp), intent(in)  :: f(:,:,:)
    integer, intent(in)   :: i,j,k
    real(dp), intent(in)  :: wx,wy,wz
    real(dp)              :: val
    val = &
      f(i  ,j  ,k  )*(1.0_dp-wx)*(1.0_dp-wy)*(1.0_dp-wz) + &
      f(i+1,j  ,k  )*wx          *(1.0_dp-wy)*(1.0_dp-wz) + &
      f(i  ,j+1,k  )*(1.0_dp-wx)*wy          *(1.0_dp-wz) + &
      f(i+1,j+1,k  )*wx          *wy          *(1.0_dp-wz) + &
      f(i  ,j  ,k+1)*(1.0_dp-wx)*(1.0_dp-wy)*wz           + &
      f(i+1,j  ,k+1)*wx          *(1.0_dp-wy)*wz           + &
      f(i  ,j+1,k+1)*(1.0_dp-wx)*wy          *wz           + &
      f(i+1,j+1,k+1)*wx          *wy          *wz
  end function cic_scalar

end module mod_cic
