program main

!*****************************************************************************80
!
!! MAIN is the main program for LATIN_RANDOM_TEST.
!
!  Discussion:
!
!    LATIN_RANDOM_TEST tests the LATIN_RANDOM library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2014
!
!  Author:
!
!    John Burkardt
!
  use netcdf
  implicit none

  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test


  integer ( kind = 4 ), parameter :: param_dim = 14
  integer ( kind = 4 ), parameter :: sample_dim = 2000
  integer, dimension(2) :: dimids

  integer ( kind = 4 ) j
  real ( kind = 8 ) x(param_dim,sample_dim)
  character(200) :: filename, param_dim_str, sample_dim_str

  integer :: ncid        ! NetCDF file ID
  integer :: varid       ! Variable ID
  integer :: param_dimid, sample_dimid  ! Dimension IDs
  integer :: retval      ! Return code for NetCDF operations

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LATIN_RANDOM_TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the LATIN_RANDOM library.'

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  LATIN_RANDOM chooses a random Latin Square'
  write ( *, '(a)' ) '  cell arrangement, and then returns'
  write ( *, '(a)' ) '  a random point from each cell.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)'  ) '  Spatial dimension =  ', param_dim
  write ( *, '(a,i6)'  ) '  Number of points =   ', sample_dim
  write ( *, '(a,i12)' ) '  Random number SEED = ', seed

  call latin_random ( param_dim, sample_dim, seed, x )

  write(param_dim_str, '(I0)') param_dim
  write(sample_dim_str, '(I0)') sample_dim

  filename = "lhs_out_p"//trim(param_dim_str)//"_n"//trim(sample_dim_str)//".nc"
  retval = nf90_create(trim(filename), NF90_CLOBBER, ncid)
  if (retval /= nf90_noerr) then
    print *, "Error creating file:", nf90_strerror(retval)
    stop
  end if

  retval = nf90_def_dim(ncid, "param_dim", param_dim, param_dimid)
  if (retval /= nf90_noerr) then
    print *, "Error defining dimension param_dim:", nf90_strerror(retval)
    stop
  end if

  retval = nf90_def_dim(ncid, "sample_dim", sample_dim, sample_dimid)
  if (retval /= nf90_noerr) then
    print *, "Error defining dimension sample_dim:", nf90_strerror(retval)
    stop
  end if

  dimids = [param_dimid, sample_dimid]

  retval = nf90_def_var(ncid, "lhs_sample", NF90_FLOAT, dimids, varid)
  if (retval /= nf90_noerr) then
    print *, "Error defining variable:", nf90_strerror(retval)
    stop
  end if

  retval = nf90_enddef(ncid)
  if (retval /= nf90_noerr) then
    print *, "Error ending define mode:", nf90_strerror(retval)
    stop
  end if

  retval = nf90_put_var(ncid, varid, x)
  if (retval /= nf90_noerr) then
    print *, "Error writing data:", nf90_strerror(retval)
    stop
  end if

  retval = nf90_close(ncid)
  if (retval /= nf90_noerr) then
    print *, "Error closing file:", nf90_strerror(retval)
    stop
  end if


!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LATIN_RANDOM_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end

function itoa(i) result(res)
  character(:), allocatable:: res
  integer, intent(in):: i
  character(range(i)+2):: tmp
  write(tmp, '(i0)') i
  res = trim(tmp)
end function
