!=====================================================
! The subroutine read_fields_txt is an example to read fields from a txt file
! Update  : 03/06/2026
! Authors: R. Sandez and M. Cécere
!=====================================================

module mod_reader_txt
  use mod_precision
  use mod_config
  use mod_grid
  implicit none

contains

  subroutine read_fields_txt(cfg, grid)
    type(config_t), intent(in)  :: cfg
    type(grid_t), intent(out)   :: grid
    integer                     :: unit, ios, n, ntot, i, j, k
    character(len=2048)         :: line
    real(dp)                    :: xx, yy, zz, bx, by, bz, ex, ey, ez
    logical                     :: exists

    ntot = cfg%nx * cfg%ny * cfg%nz
    inquire(file=trim(cfg%field_file), exist=exists)
    if (.not. exists) then
       write(*,*) "ERROR: field file does not exist: ", trim(cfg%field_file)
       stop
    end if
    call allocate_base_grid(grid, cfg%nx, cfg%ny, cfg%nz)
    call set_uniform_coordinates(grid, cfg%xmin, cfg%xmax, cfg%ymin, cfg%ymax, cfg%zmin, cfg%zmax)
    open(newunit=unit, file=trim(cfg%field_file), status="old", action="read", iostat=ios)
    if (ios /= 0) stop "ERROR: could not open TXT field file"
    n = 0
    do
       read(unit, '(A)', iostat=ios) line
       if (ios /= 0) exit
       if (len_trim(line) == 0) cycle
       if (line(1:1) == '#') cycle
       read(line, *, iostat=ios) xx, yy, zz, bx, by, bz, ex, ey, ez
       if (ios /= 0) cycle
       n = n + 1
       if (n > ntot) stop "ERROR: too many data rows in TXT field file"
       i = mod(n-1, grid%nx) + 1
       j = mod((n-1)/grid%nx, grid%ny) + 1
       k = ((n-1)/(grid%nx*grid%ny)) + 1
       grid%bx(i,j,k)=bx; grid%by(i,j,k)=by; grid%bz(i,j,k)=bz
       grid%ex(i,j,k)=ex; grid%ey(i,j,k)=ey; grid%ez(i,j,k)=ez
    end do
    close(unit)
    if (n /= ntot) then
       write(*,*) "ERROR: wrong number of rows in TXT file."
       write(*,*) " Read rows:     ", n
       write(*,*) " Expected rows: ", ntot
       stop
    end if
    write(*,*) "Read TXT field file: ", trim(cfg%field_file)
    write(*,*) "Rows read: ", n
  end subroutine read_fields_txt

end module mod_reader_txt
