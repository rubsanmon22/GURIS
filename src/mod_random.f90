!=====================================================
! This module contains subroutines and functions to calculate randum values
! Update  : 03/06/2026
! Authors: R. Sandez and M. Cécere
!=====================================================

module mod_random
  use mod_precision
  implicit none
contains
  subroutine init_random_seed(seed)
    integer, intent(in) :: seed
    integer :: n, i
    integer, allocatable :: seed_array(:)
    call random_seed(size=n)
    allocate(seed_array(n))
    do i=1,n
       seed_array(i) = seed + 37*(i-1)
    end do
    call random_seed(put=seed_array)
    deallocate(seed_array)
  end subroutine init_random_seed

  function uniform_random() result(u)
    real(dp) :: u
    call random_number(u)
  end function uniform_random

  function gaussian_random() result(g)
    real(dp) :: g, u1, u2
    real(dp), parameter :: twopi = 6.283185307179586476925286766559_dp
    call random_number(u1); call random_number(u2)
    u1 = max(min(u1,1.0_dp-1.0e-15_dp),1.0e-15_dp)
    g = sqrt(-2.0_dp*log(u1))*cos(twopi*u2)
  end function gaussian_random
end module mod_random
