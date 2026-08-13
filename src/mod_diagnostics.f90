!=====================================================
! The subroutine compute_energies calculates the kinetic energy
! Update  : 03/06/2026
! Authors: R. Sandez and M. Cécere
!=====================================================

module mod_diagnostics
  use mod_precision
  use mod_config
  use mod_grid
  use mod_particles
  use mod_cic
  use mod_evolution_guiding_center !se incorpora para calcular la energia cinetica con las velocidades de los drifts

  implicit none
contains
  
subroutine compute_energies(cfg,grid,part,initial)
    type(config_t), intent(in) :: cfg
    type(grid_t), intent(in) :: grid
    type(particles_t), intent(inout) :: part
    logical, intent(in) :: initial
    integer :: ip
    real(dp) :: ekin
    real(dp) :: dxdt,dydt,dzdt,dvpardt !velocidades que se incorporan a la energia cinetica
    type(field_sample_t) :: fs
    logical :: inside

    do ip=1,part%np
       if (.not. part%active(ip)) then
          ekin=-1.0_dp
       else
          call cic_sample_fields(grid,part%x(ip),part%y(ip),part%z(ip),fs,inside)
          if (.not. inside) then
             ekin=-1.0_dp
          else
            ! llamo las variables necesarias para el cálculo
            call guiding_center_rhs(cfg, part%x(ip), part%y(ip), part%z(ip), part%vpar(ip), part%mu(ip), fs, dxdt, dydt, dzdt, dvpardt)
            ! nuevo valor de ekin
            ekin=0.5_dp*cfg%m_particle*(dxdt**2 + dydt**2 + dzdt**2) + part%mu(ip)*fs%bmag
            ! ekin=0.5_dp*cfg%m_particle*part%vpar(ip)**2 + part%mu(ip)*fs%bmag
          end if
       end if
       if (initial) then
          part%ekin_initial(ip)=ekin
       else
          part%ekin_final(ip)=ekin
       end if
    end do
  end subroutine compute_energies

  function count_active(part) result(nactive)
    type(particles_t), intent(in) :: part
    integer :: nactive, ip
    nactive=0
    do ip=1,part%np
       if (part%active(ip)) nactive=nactive+1
    end do
  end function count_active
end module mod_diagnostics
