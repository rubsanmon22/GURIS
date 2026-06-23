!=====================================================
! The subroutine guiding_center_rhs calculates dx/dt, dy/dt, dz/dt and dvpar/dt
! The subroutine advance_guiding_center_rk2 calculates the evolution of x, y, t and vpar
! Update  : 03/06/2026
! Authors: R. Sandez and M. Cécere
!=====================================================


module mod_evolution_guiding_center
  use mod_precision
  use mod_config
  use mod_grid
  use mod_particles
  use mod_cic
  use mod_boundary
  implicit none
contains
  subroutine guiding_center_rhs(cfg,x,y,z,vpar,mu,fs,dxdt,dydt,dzdt,dvpardt)
    type(config_t), intent(in)        :: cfg
    real(dp), intent(in)              :: x,y,z,vpar,mu
    type(field_sample_t), intent(in)  :: fs
    real(dp), intent(out)             :: dxdt,dydt,dzdt,dvpardt
    real(dp)                          :: b2, inv_b2, epar, vex,vey,vez, vgbx,vgby,vgbz, vcx,vcy,vcz
    real(dp)                          :: bxg,byg,bzg, bxk,byk,bzk, bdgradb, q,m

    ! particle charge and mass
    q=cfg%q_particle; m=cfg%m_particle

    dxdt=0.0_dp; dydt=0.0_dp; dzdt=0.0_dp; dvpardt=0.0_dp
   ! velocidades iniciales de los drifts y aceleración paralela

    if (fs%bmag <= 0.0_dp) return
    ! B square
    b2=fs%bmag*fs%bmag; inv_b2=1.0_dp/b2

    ! parallel vecE
    epar=fs%ex*fs%bhx + fs%ey*fs%bhy + fs%ez*fs%bhz

    ! drift E: ve = vecE x vecB/B^2
    vex=0.0_dp; vey=0.0_dp; vez=0.0_dp
    if (cfg%use_exb_drift) then
       vex=cfg%exb_factor*(fs%ey*fs%bz - fs%ez*fs%by)*inv_b2*(1/m)
       vey=cfg%exb_factor*(fs%ez*fs%bx - fs%ex*fs%bz)*inv_b2*(1/m)
       vez=cfg%exb_factor*(fs%ex*fs%by - fs%ey*fs%bx)*inv_b2*(1/m)
    end if

    ! (bxg, byg, bzg) = vecb x grad B,
    ! drift gradient B: vgb = mu/q B vecb x grad B
    vgbx=0.0_dp; vgby=0.0_dp; vgbz=0.0_dp
    if (cfg%use_gradb_drift) then
       bxg=fs%bhy*fs%gradb_z - fs%bhz*fs%gradb_y
       byg=fs%bhz*fs%gradb_x - fs%bhx*fs%gradb_z
       bzg=fs%bhx*fs%gradb_y - fs%bhy*fs%gradb_x
       vgbx=((mu)/(q*fs%bmag))*bxg; vgby=((mu)/(q*fs%bmag))*byg; vgbz=((mu)/(q*fs%bmag))*bzg
    end if

    ! (bxk, byk, bzk) =  b x kappa: kappa = b . grad b
    ! drift curvature: vc = m*vpar**2/(q*B) vecb x kappa
    vcx=0.0_dp; vcy=0.0_dp; vcz=0.0_dp
    if (cfg%use_curvature_drift) then
       bxk=fs%bhy*fs%kappa_z - fs%bhz*fs%kappa_y
       byk=fs%bhz*fs%kappa_x - fs%bhx*fs%kappa_z
       bzk=fs%bhx*fs%kappa_y - fs%bhy*fs%kappa_x
       vcx=(vpar*vpar/(q*fs%bmag))*bxk
       vcy=(vpar*vpar/(q*fs%bmag))*byk
       vcz=(vpar*vpar/(q*fs%bmag))*bzk
    end if

    ! dr/dt = vpar b + vdrift; vdrift = ve + vgb + vc (electric, gradient of B and curvature)
    ! todo en unidades de velocidad por unidad de masa
    dxdt=vpar*fs%bhx + vex + vgbx + vcx
    dydt=vpar*fs%bhy + vey + vgby + vcy
    dzdt=vpar*fs%bhz + vez + vgbz + vcz
    
    ! dvpar/dt = q/m epar - mu/m bgradb
    if (cfg%use_mirror_force) then
       ! bdgradb = vecb . grad vecb
       bdgradb=fs%bhx*fs%gradb_x + fs%bhy*fs%gradb_y + fs%bhz*fs%gradb_z
       dvpardt=(q/(m*m))*epar - (mu/(m*m))*bdgradb !aceleracion por unidad de masa
    else
       dvpardt=(q/(m*m))*epar !aceleracion por unidad de masa
    end if
  end subroutine guiding_center_rhs

  subroutine advance_guiding_center_rk2(cfg,grid,part,ip)
    type(config_t), intent(in)        :: cfg
    type(grid_t), intent(in)          :: grid
    type(particles_t), intent(inout)  :: part
    integer, intent(in)               :: ip
    type(field_sample_t)              :: fs
    logical                           :: inside, ok
    real(dp)                          :: kx1,ky1,kz1,kv1, kx2,ky2,kz2,kv2, xm,ym,zm,vpm

    if (.not. part%active(ip)) return
    call apply_boundary_particle(cfg,grid,part,ip)

    if (.not. part%active(ip)) return
    call cic_sample_fields(grid,part%x(ip),part%y(ip),part%z(ip),fs,inside)

    if (.not. inside) then; part%active(ip)=.false.; return; end if
    call guiding_center_rhs(cfg,part%x(ip),part%y(ip),part%z(ip),part%vpar(ip),part%mu(ip),fs,kx1,ky1,kz1,kv1)

    ! first step of RK2
    xm=part%x(ip)+0.5_dp*cfg%dt*kx1
    ym=part%y(ip)+0.5_dp*cfg%dt*ky1
    zm=part%z(ip)+0.5_dp*cfg%dt*kz1
    vpm=part%vpar(ip)+0.5_dp*cfg%dt*kv1

    call apply_boundary_position(cfg,grid,xm,ym,zm,ok)

    if (.not. ok) then; part%active(ip)=.false.; return; end if
    call cic_sample_fields(grid,xm,ym,zm,fs,inside)

    if (.not. inside) then; part%active(ip)=.false.; return; end if
    call guiding_center_rhs(cfg,xm,ym,zm,vpm,part%mu(ip),fs,kx2,ky2,kz2,kv2)

    ! last step of RK2
    part%x(ip)=part%x(ip)+cfg%dt*kx2
    part%y(ip)=part%y(ip)+cfg%dt*ky2
    part%z(ip)=part%z(ip)+cfg%dt*kz2
    part%vpar(ip)=part%vpar(ip)+cfg%dt*kv2
    call apply_boundary_particle(cfg,grid,part,ip)
  end subroutine advance_guiding_center_rk2
end module mod_evolution_guiding_center
