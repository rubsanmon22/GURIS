!=====================================================
! The subroutine selects the type of file to read the fields
! Update  : 03/06/2026
! Authors: R. Sandez and M. Cécere
!=====================================================

module mod_reader
  use mod_config
  use mod_grid
  use mod_reader_txt
  implicit none
contains
  subroutine read_fields_from_file(cfg, grid)
    type(config_t), intent(in) :: cfg
    type(grid_t), intent(out) :: grid
    select case (trim(cfg%file_type))
    case ('txt','dat')
       call read_fields_txt(cfg, grid)
    case default
       write(*,*) "ERROR: unsupported file_type: ", trim(cfg%file_type)
       write(*,*) "This version supports file_type='txt' or 'dat'."
       stop
    end select
  end subroutine read_fields_from_file
end module mod_reader
