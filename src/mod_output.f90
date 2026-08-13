!=====================================================
! The module contains different outputs
! Update  : 03/06/2026
! Authors: R. Sandez and M. Cécere
!=====================================================

module mod_output
  use mod_precision
  use mod_config
  use mod_particles
  implicit none

contains

  subroutine ensure_output_dir(cfg)
    type(config_t), intent(in) :: cfg
    call execute_command_line('mkdir -p '//trim(cfg%output_dir))
  end subroutine ensure_output_dir

  subroutine write_outputs(cfg,part)
    type(config_t), intent(in) :: cfg
    type(particles_t), intent(in) :: part
    call ensure_output_dir(cfg)
    call write_particle_table(cfg,part)
    call write_energy_table(cfg,part)
    call write_energy_histogram(cfg,part)
  end subroutine write_outputs

  subroutine write_particle_table(cfg,part)
    type(config_t), intent(in) :: cfg
    type(particles_t), intent(in) :: part
    character(len=1024) :: filename
    integer :: unit, ip
    filename=trim(cfg%output_dir)//'/'//trim(cfg%output_prefix)//'_particles_final.dat'
    open(newunit=unit,file=trim(filename),status='replace',action='write')
    write(unit,'(A)') '# id x y z vpar mu active'
    do ip=1,part%np
       write(unit,'(I10,1X,5(ES24.15,1X),L1)') ip,part%x(ip),part%y(ip),part%z(ip),part%vpar(ip),part%mu(ip),part%active(ip)
    end do
    close(unit)
  end subroutine write_particle_table

  subroutine write_energy_table(cfg,part)
    type(config_t), intent(in) :: cfg
    type(particles_t), intent(in) :: part
    character(len=1024) :: filename
    integer :: unit, ip
    filename=trim(cfg%output_dir)//'/'//trim(cfg%output_prefix)//'_energy_table.dat'
    open(newunit=unit,file=trim(filename),status='replace',action='write')
    write(unit,'(A)') '# id E_initial E_final active'
    do ip=1,part%np
       write(unit,'(I10,1X,2(ES24.15,1X),L1)') ip,part%ekin_initial(ip),part%ekin_final(ip),part%active(ip)
    end do
    close(unit)
  end subroutine write_energy_table

  subroutine write_energy_histogram(cfg,part)
    type(config_t), intent(in) :: cfg
    type(particles_t), intent(in) :: part
    character(len=1024) :: filename
    integer :: unit, ip, ibin, nbins
    integer, allocatable :: hi(:), hf(:)
    real(dp) :: emin, emax, de, ecenter, ei, ef
    logical :: first
    nbins=cfg%nbins_energy
    allocate(hi(nbins),hf(nbins)); hi=0; hf=0
    first=.true.; emin=0.0_dp; emax=1.0_dp
    do ip=1,part%np
       ei=part%ekin_initial(ip); ef=part%ekin_final(ip)
       if (ei >= 0.0_dp .and. ef >= 0.0_dp) then
          if (first) then
             emin=min(ei,ef); emax=max(ei,ef); first=.false.
          else
             emin=min(emin,min(ei,ef)); emax=max(emax,max(ei,ef))
          end if
       end if
    end do
    if (emax <= emin) then
       if (emin == 0.0_dp) then; emin=-0.5_dp; emax=0.5_dp
       else; emin=emin-0.5_dp*abs(emin); emax=emax+0.5_dp*abs(emax)
       end if
    end if
    de=(emax-emin)/real(nbins,dp)
    do ip=1,part%np
       ei=part%ekin_initial(ip); ef=part%ekin_final(ip)
       if (ei >= 0.0_dp) then
          ibin=int((ei-emin)/de)+1; ibin=min(max(ibin,1),nbins); hi(ibin)=hi(ibin)+1
       end if
       if (ef >= 0.0_dp) then
          ibin=int((ef-emin)/de)+1; ibin=min(max(ibin,1),nbins); hf(ibin)=hf(ibin)+1
       end if
    end do
    filename=trim(cfg%output_dir)//'/'//trim(cfg%output_prefix)//'_energy_histogram.dat'
    open(newunit=unit,file=trim(filename),status='replace',action='write')
    write(unit,'(A)') '# E_center N_initial N_final'
    do ibin=1,nbins
       ecenter=emin+(real(ibin,dp)-0.5_dp)*de
       write(unit,'(ES24.15,1X,I12,1X,I12)') ecenter,hi(ibin),hf(ibin)
    end do
    close(unit)
    deallocate(hi,hf)
  end subroutine write_energy_histogram

  subroutine open_trajectory_file(cfg,unit)
      type(config_t), intent(in) :: cfg
      integer, intent(out) :: unit
      character(1024) :: filename

      call ensure_output_dir(cfg)
      filename=trim(cfg%output_dir)//'/'//trim(cfg%output_prefix)//'_trajectories.dat'
      open(newunit=unit, file=trim(filename), status='replace', action='write')
      write(unit, '(A)') '# id step x y z vpar mu active'
  end subroutine open_trajectory_file

  subroutine write_trajectory(cfg, unit, istep, time, part)
      type(config_t), intent(in) :: cfg
      type(particles_t), intent(in) :: part
      real(dp), intent(in) :: time
      integer, intent(in) :: istep, unit
      integer :: ip, idx

      if (cfg%trajectory_nids >0) then
         do ip=1, cfg%trajectory_nids
            idx = cfg%trajectory_ids(ip)
            if (idx >= 1 .and. idx <= part%np) then
               write(unit, '(I8,1X,ES24.15,1X,I10,1X,3(ES24.15,1X),L1)') &
                  istep, time, idx, part%x(idx), part%y(idx), part%z(idx), part%active(idx)
            end if
         end do
      else
         do idx = 1, min(cfg%trajectory_npart, part%np)
            write(unit, '(I8, 1X, ES24.15,1X,I10,1X,3(ES24.15,1X),L1)') &
               istep, time, idx, part%x(idx), part%y(idx), part%z(idx), part%active(idx)
         end do
      end if
   end subroutine write_trajectory

end module mod_output
