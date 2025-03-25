module global_fun
use sort_mod

contains

! ------------------------------------------
subroutine concat_mat(mat1, mat2, ndim_opt)
! concatenate mat2 to mat1 to the bottom row by default
! set ndim_opt to 1 to concatenate mat2 to the right
! slow implementation but couldn't find a faster way for fortran at the moment ... 
! there probably exists some fancy way to accelerate this with memory allocation tricks
! but let's worry about it when it becomes the bottleneck

real(8), dimension(:,:), allocatable, intent(inout):: mat1
real(8), dimension(:,:), allocatable:: mat_interm
real(8), dimension(:,:), intent(in):: mat2
real(8), dimension(:), allocatable:: arr
integer:: ncol1, nrow1, ncol2, nrow2, ndim, len_mat1, len_mat2
integer, optional:: ndim_opt

ncol1 = size(mat1, 1)
nrow1 = size(mat1, 2)
ncol2 = size(mat2, 1)
nrow2 = size(mat2, 2)

if (present(ndim_opt)) then
   if (ndim_opt > 2 .or. ndim_opt < 1) then
      print*, '### Error in concat_mat(): ndim_opt out of range'
      print*, 'ndim_opt must be either 1 (append as column) or 2 (append as row)'
      stop
   endif
   ndim = ndim_opt
else
   ndim = 2
endif

! check dimension compatibility
len_mat1 = size(mat1, 3-ndim)
len_mat2 = size(mat2, 3-ndim)
if (len_mat1 /= 0. .and. len_mat1 /= len_mat2) then
   print*, '### Error in concat_mat(): incompatible dimension'
   print*, 'Appended:', len_mat1, 'Appending:', len_mat2
   print*, "Check the dimension of the two matrix you're trying to concatenate"
   stop
endif

if (ndim == 2)  then
   ! concatenate to the bottom
   ! create mat1 if does not exist
   if (ncol1 == 0) ncol1 = ncol2
   allocate(mat_interm(ncol1, nrow1+nrow2))
   mat_interm(1:ncol1, 1:nrow1) = mat1
   mat_interm(:, nrow1+1:nrow1+nrow2) = mat2
else
   ! concatenate to the right
   if (nrow1 == 0) nrow1 = nrow2
   allocate(mat_interm(ncol1+ncol2, nrow1))
   mat_interm(1:ncol1, 1:nrow1) = mat1
   mat_interm(ncol1+1:ncol1+ncol2, :) = mat2
endif

mat1 = mat_interm

deallocate(mat_interm)

end subroutine concat_mat

! ------------------------------------------
subroutine append_mat(mat, arr, arr_len, ndim_opt)
! similar to concat_mat but the only append a 1D array instead of a 2D matrix
integer, intent(in) :: arr_len
real(8), dimension(:,:), allocatable, intent(inout):: mat
real(8), dimension(arr_len), intent(in):: arr
real(8), dimension(:,:), allocatable:: mat1d
integer, optional, intent(in):: ndim_opt
integer:: ndim



if (present(ndim_opt)) then
   if (ndim_opt > 2 .or. ndim_opt < 1) then
      print*, '### Error in append_mat(): ndim_opt out of range'
      print*, 'ndim_opt must be either 1 (append as column) or 2 (append as row)'
      stop
   endif
   ndim = ndim_opt
else
   ndim = 2
endif

if (ndim == 2) then
   ! append at the bottom
   allocate(mat1d(arr_len, 1))
else
   allocate(mat1d(1, arr_len))
endif

mat1d = reshape(arr, shape(mat1d))
call concat_mat(mat, mat1d, ndim)
deallocate(mat1d)

end subroutine append_mat

! ------------------------------------------
subroutine sort_col_prio(mat, ncol, nrow, prio_col_opt, last_col_opt, tol_sigfig_opt)
! sort rows of matrix with descending column priority (or as specified by prio_col_opt)
! e.g.,
! [ 1 5 2 
!   2 3 5
!   1 3 4 ]
! becomes: 
! [ 1 3 4
!   1 5 2
!   2 3 5 ]
! prio_col_opt has to be allocatable when being passed in to the subroutine
! one can also pass in last_col_opt (instead of prio_col_opt) to have sorted columns 
! from 1 to last_col_opt
! pass in tol_sigfig_opt to specifies a sigfig below which order doesn't matter
! (good for real number comparison)

real(8), dimension(ncol,nrow), intent(inout):: mat
integer, intent(in):: ncol, nrow
integer, intent(in), optional:: last_col_opt, tol_sigfig_opt
integer, intent(in), allocatable, dimension(:), optional:: prio_col_opt
integer, dimension(:), allocatable:: idxs, sort_order
integer:: i, icol, last_col, ncol_sort, tol_sigfig

allocate(idxs(nrow))

if (present(last_col_opt)) then
   last_col = last_col_opt
else
   last_col = ncol
endif

if (present(tol_sigfig_opt)) tol_sigfig = tol_sigfig_opt

if (present(prio_col_opt)) then
   ! print*, size(prio_col_opt)
   ! debug before making changes
   if (maxval(prio_col_opt) > last_col) then
      print*, maxval(prio_col_opt), last_col
      print*, '### ERROR in sort_col_prio()'
      print*, 'Entry for column priority contains column number greater than the'
      print*, 'last (right most) column interested'
      stop
   elseif (has_duplicate(prio_col_opt, size(prio_col_opt))) then
      print*, '### ERROR in sort_col_prio()'
      print*, 'Duplicate column number found in prio_col_opt, please fix!'
      stop
   else
      ncol_sort = size(prio_col_opt)
      allocate(sort_order(ncol_sort))
      do icol = 1, ncol_sort
         sort_order(icol) = prio_col_opt(ncol_sort-icol+1)
      enddo
   endif
else
   ncol_sort = last_col
   allocate(sort_order(ncol_sort))
   do icol = 1, last_col
      sort_order(icol) = last_col-icol+1
   enddo
endif

do i = 1, ncol_sort
   icol = sort_order(i)
   print*, 'sorting column...', icol
   if (present(tol_sigfig_opt)) then
      call getSortArrayIdx(nrow, mat(icol, :), idxs, tol_sigfig_opt=tol_sigfig)
   else
      call indexArrayReal(nrow, mat(icol, :), idxs)
   endif
   mat = mat(:, idxs)
   ! print*, 'within sorting, icol=', icol
   ! print*, mat(icol,1:6)
   ! print*, idxs
enddo


deallocate(idxs)
deallocate(sort_order)

end subroutine sort_col_prio

! ------------------------------------------
logical function has_duplicate(arr, nval)
! does not have to be sorted

integer:: idx, nval
integer, dimension(nval):: arr
integer, allocatable, dimension(:):: arr_dupcheck

do idx = 1, nval
   if (size(arr_dupcheck) > 0) then
      if (any(arr_dupcheck == arr(idx))) then
         has_duplicate = .true.
         return
      endif
   else
      if (size(arr_dupcheck) > 0) then
         arr_dupcheck = [arr_dupcheck, arr(idx)]
      else
         arr_dupcheck = arr(idx)
      endif
   endif   
enddo
has_duplicate = .false.

end function

! ------------------------------------------
subroutine dedup_sorted(mat, last_col_opt, ndim_opt, tol_sigfig_opt)
! checks duplicated row (default) or column of a sorted matrix
! two rows that have all elements sharing the same 3 sig figs will be considered
! duplicated by default, unless otherwise specified by `tol_sigfig_opt`

! WARNING: MAKE SURE the input matrix is sorted (column priority does not matter)
! otherwise the output matrix won't be dedup'd WITHOUT WARNING

real(8), dimension(:,:), allocatable, intent(inout):: mat
integer, intent(in), optional:: tol_sigfig_opt, ndim_opt, last_col_opt
integer, dimension(:), allocatable:: idx_incl
integer:: ndim, irow, icol, tol_sigfig, last_col
character(len=100):: str_val1, str_val2, sigfig_fmt

if (present(tol_sigfig_opt)) then
   tol_sigfig = tol_sigfig_opt
else
   tol_sigfig = 4
endif

if (present(last_col_opt)) then
   last_col = last_col_opt
else
   last_col = size(mat, 1)
endif

if (present(ndim_opt)) then
   if (ndim_opt > 2 .or. ndim_opt < 1) then
      print*, '### Error in concat_mat(): ndim_opt out of range'
      print*, 'ndim_opt must be either 1 (append as column) or 2 (append as row)'
      stop
   endif
   ndim = ndim_opt
else
   ndim = 2
endif

ncol = size(mat, 1)
nrow = size(mat, 2)

if (ndim == 2) then
   sigfig_fmt = '('//itoa(ncol)//'e'//itoa(tol_sigfig+8)//'.'//itoa(tol_sigfig)//')'
   do irow = 1, nrow-1
      write(str_val1, trim(sigfig_fmt)) mat(1:last_col,irow)
      write(str_val2, trim(sigfig_fmt)) mat(1:last_col,irow+1)
      if (str_val1 .ne. str_val2) then
         idx_incl = [idx_incl, irow]
      endif
   enddo
   idx_incl = [idx_incl, nrow]
   mat = mat(:,idx_incl)
else
   sigfig_fmt = '('//itoa(nrow)//'e'//itoa(tol_sigfig+8)//'.'//itoa(tol_sigfig)//')'
   do icol = 1, ncol-1
      write(str_val1, trim(sigfig_fmt)) mat(icol, 1:last_col)
      write(str_val2, trim(sigfig_fmt)) mat(icol+1, 1:last_col)
      if (str_val1 .ne. str_val2) then
         idx_incl = [idx_incl, icol]
      endif
   enddo
   idx_incl = [idx_incl, ncol]
   mat = mat(idx_incl, :)
endif

deallocate(idx_incl)
end subroutine dedup_sorted

! ------------------------------------------
integer function find_closest_arr_in_mat(mat, arr, ncol, nrow, tol_sigfig, prio_col_opt) result (idx_found)
! find the row idx with vector closest to `arr` in a matrix `mat`
! tolin: fractional tolerance when looking for closest value
! prio_col_opt: specifies which column to find first (from the most important to the least)
!     defaults to descending order of importance if not specified.
real(8), dimension(ncol, nrow):: mat 
real(8), dimension(ncol):: arr
real(8), allocatable, dimension(:):: mat_slice
integer, allocatable, dimension(:), optional:: prio_col_opt
integer, allocatable, dimension(:):: prio_col
integer:: idx, i, icol, idx_s, idx_e, tol_sigfig, idx_prev

idx_found = 1
idx_s = 1
idx_e = nrow
idx_prev = 1

if (present(prio_col_opt)) then
   ! debug before making changes
   if (maxval(prio_col_opt) > ncol) then
      print*, '### ERROR in find_closest_arr_in_mat()'
      print*, 'Entry for column priority contains column number greater than the'
      print*, 'last (right most) column interested'
      stop
   elseif (has_duplicate(prio_col_opt, size(prio_col_opt))) then
      print*, '### ERROR in find_closest_arr_in_mat()'
      print*, 'Duplicate column number found in prio_col_opt, please fix!'
      stop
   else
      allocate(prio_col(ncol))
      prio_col = prio_col_opt
   endif
else
   allocate(prio_col(ncol))
   do icol = 1, ncol
      prio_col(icol) = icol
   enddo
endif

! this should narrow down the possible index to one given that the matrix is sorted and 
! and dedup'd with the same tol_sigfig
do i = 1, ncol
   idx_prev = idx_s
   icol = prio_col(i)
   mat_slice = mat(icol, idx_s:idx_e)
   call find_closest_val_in_arr(mat_slice, arr(icol), idx_s, idx_e, &
                                'sigfig', tol_sigfig_opt=tol_sigfig)
   ! exit as long as we've narrowed down to one row
   idx_s = idx_s + idx_prev - 1
   idx_e = idx_e + idx_prev - 1
   if (idx_s == idx_e) exit
enddo

idx_found = idx_s

deallocate(prio_col)

end function

! ------------------------------------------
Subroutine read_sctab
use micro_prm
use parameters, only: h_shape, num_h_moments
use namelists, only: KiD_outdir

character(len=100) :: cnuval
logical :: file_exist

if (num_h_moments(1) > 3 .and. num_h_moments(2) > 3) then
   if (num_h_moments(1) .ne. num_h_moments(2)) then
      print*, '### FATAL ###'
      print*, 'the number of moments predicted has to be the same for single category AMP'
      print*, 'change the values of `num_h_moments` in the namelist to skip this warning'
      stop
   endif
   if (h_shape(1) .ne. h_shape(2)) then
      print*, '### FATAL ###'
      print*, 'the input shape parameter has to be the same for cloud and rain'
      print*, 'change the values of `h_shape` in the namelist to skip this warning'
      print*, 'or change the code to incorporate this scenario'
      stop
   endif
endif

if (num_h_moments(1) == 4) then
   prio_c = column_priority_4m
   pmomsc=(/3,0,imomc1,imomc2,imomr1,imomr2/)
   folder_sclut='./src/input_data/singcat_lutables/'

   write(cnuval, '(I0)'), int(h_shape(1))
   if (imomr1 >= 10) then
      write(momstr4, '(A, I1, I1, I2, I2)') 'M',imomc1, imomc2, imomr1, imomr2
   elseif (imomc2 >= 10 .and. imomr2 >= 10) then
      write(momstr4, '(A, I1, I2, I1, I2)') 'M',imomc1, imomc2, imomr1, imomr2
   elseif (imomr2 >= 10) then
      write(momstr4, '(A, I1, I1, I1, I2)') 'M',imomc1, imomc2, imomr1, imomr2
   else
      write(momstr4, '(A, I1, I1, I1, I1)') 'M',imomc1, imomc2, imomr1, imomr2
   endif

   sc4m_lufile='nu'//trim(cnuval)//'_'//trim(momstr4)//'.txt'
   sc4m_lu_abspath = trim(folder_sclut)//trim(sc4m_lufile)
   inquire(file = sc4m_lu_abspath, exist = file_exist, size = fsize)

   if (.not. file_exist .or. fsize == 0) then
      tab4m_nrow = 0
      tab4m_ncol = 0
      print*, 'No existing table for the selected moments'
      print*, 'Creating a new one...'
   else
      ! find the dimension of the lutable and allocate array
      CALL execute_command_line('wc -l < '//trim(sc4m_lu_abspath)//' > '//&
         trim(KiD_outdir)//'rowcount.txt')
      OPEN(unit = 20, file=trim(KiD_outdir)//'rowcount.txt') 
      READ(20, *) tab4m_nrow
      close(20, status='delete')
      CALL execute_command_line("awk '{print NF}' "//trim(sc4m_lu_abspath)//&
         " | sort -nu | tail -n 1 > "//trim(KiD_outdir)//"colcount.txt")
      OPEN(unit = 20, file=trim(KiD_outdir)//'colcount.txt') 
      READ(20, *) tab4m_ncol
      close(20, status='delete')
      allocate(sc4m_tab(tab4m_ncol, tab4m_nrow))
      open(20, file = sc4m_lu_abspath)
      read(20, *) sc4m_tab
      close(20)

      sc4m_tab_moms = sc4m_tab(1:4, :)
      sc4m_tab_dns = sc4m_tab(5:6, :)
      sc4m_tab_giveup = sc4m_tab(7, :)
      print*, 'Existing lookup table found...'
   endif
endif

End Subroutine read_sctab 

! ------------------------------------------
Subroutine  write_sctab
use micro_prm

if (npm == 4) then
   print*, 'Updating lookup table'
   ! call append_mat(sc4m_tab_dns, sc4m_tab_giveup, *array_length_variable*, ndim_opt = 1)
   call concat_mat(sc4m_tab_moms, sc4m_tab_dns, ndim_opt = 1)
   call concat_mat(sc4m_tab, sc4m_tab_moms, ndim_opt = 2)
   call sort_col_prio(sc4m_tab, size(sc4m_tab,1), size(sc4m_tab,2), &
      prio_col_opt=prio_c, tol_sigfig_opt=p_sigfig)
   call dedup_sorted(sc4m_tab, last_col_opt=4, tol_sigfig_opt=p_sigfig)
   if (size(sc4m_tab,2) > tab4m_nrow) then
      open(20, file = sc4m_lu_abspath)
      sigfig_fmt = '(7e'//itoa(p_sigfig+8)//'.'//itoa(p_sigfig)//')'
      write(20, trim(sigfig_fmt)) sc4m_tab
      close(20)
      print*, 'Lookup table updated'
   endif
   if (allocated(sc4m_tab)) deallocate(sc4m_tab)
endif

End Subroutine write_sctab 

! ------------------------------------------
Function  get_meandiam(m3, m0)
real(8) :: m3, m0
real(8) :: get_meandiam

get_meandiam = (m3/m0)**(1/3.)*1.e6

End Function get_meandiam 

! ------------------------------------------
! mapRange2R: {{{
Subroutine mapRange2R(n, x, rangeMin, rangeMax, direc)
! map a ratio [0, 1] to the entire real number space if direc == 'f'
! opposite if direc == 'b'
use physconst, only: pi
implicit none

integer, intent(in) :: n
real(8), dimension(n), intent(inout) :: x
character(1), intent(in) :: direc
real(8), intent(in) :: rangeMin, rangeMax
real(8) :: a, b, d
integer :: i

a = LOG10(rangeMin)
b = LOG10(rangeMax)
d = b - a
if (direc .eq. 'f') then
   do i=1,n
      if (x(i)<rangeMin) x(i) = rangeMin*1.00001
      if (x(i)>rangeMin) x(i) = rangeMax*0.99999
      ! if (x(i)<0 .or. x(i)>1) then
      !    print*, '### error in mapRange2R'
      !    print*, 'ratio must be between 0 and 1'
      !    print*, 
      !    stop
      ! endif
      x(i) = dtan(((LOG10(x(i)) - a)/d - 0.5)*pi)
   enddo
elseif (direc .eq. 'b') then
   do i=1,n
      x(i) = 10.**((datan(x(i))/pi+0.5)*d+a)
   enddo
else
   stop "the 2nd argument of mapRange2R has to be either 'f' or 'b'."
endif

End Subroutine mapRange2R 
!}}}

! ------------------------------------------

! mapRatio2R : {{{
Subroutine  mapRatio2R(n, x, direc)
integer, intent(in) :: n
real(8), dimension(n), intent(inout) :: x
integer :: i
character(1) :: direc

if (direc .eq. 'f') then
   do i=1,n
      if (x(i)<0 .or. x(i)>1) then
         print*, '### error in mapRange2R'
         print*, 'ratio must be between 0 and 1'
         stop
      endif
      x(i) = dtan((x(i)-0.5)*pi)
   enddo
elseif (direc .eq. 'b') then
   do i=1,n
      x(i) = datan(x(i))/pi + 0.5
   enddo
else
   stop "the 2nd argument of mapRange2R has to be either 'f' or 'b'."
endif

End Subroutine mapRatio2R 

!}}}

! ------------------------------------------
! linspace: {{{
! adopted from: https://gist.github.com/ivan-pi/f4b4741d7ed54ceff787c85d6ba22a5a

function linspace(start,end,num,endpoint,step) result(samples)
   ! PARAMETERS
   real(8), intent(in) :: start 
   !! The starting value of the sequence.
   real(8), intent(in) :: end
   !! The end value of the sequence, unless `endpoint` is set to `.false.`. 
   !! In that case, the sequence consists of all but the last of `num + 1` 
   !! evenly spaced samples, so that `end` is excluded. Note that the 
   !! step size changes when `endpoint` is `.false.`.
   integer, intent(in), optional :: num
   !! Number of samples to generate. Default value is 50.
   logical, intent(in), optional :: endpoint
   !! If `.true.`, `end` is the last sample. Otherwise, it is not included. Default is `.true.`.
   real(8), intent(out), optional :: step
   !! If present, `step` is the size of spacing between samples.

   ! RETURNS
   real(8), allocatable :: samples(:)
   !! There are `num` equally spaced samples in the closed interval `[start, stop]` or 
   !! the half-open interval `[start, stop)` (depending on whether `endpoint` is `.true.` or `.false.`).

   integer :: num_, i
   logical :: endpoint_
   real(8) :: step_

   num_ = 50
   if (present(num)) num_ = num

   endpoint_ = .true.
   if (present(endpoint)) endpoint_ = endpoint

   ! find step size
   if (endpoint_) then
      step_ = (end - start)/real(num_-1,8)
   else
      step_ = (end - start)/real(num_,8)
   end if

   if (present(step)) step = step_

   allocate(samples(num_))
   do i = 1, num_
      samples(i) = start + (i-1)*step_
   end do
end function linspace
! }}}

! ------------------------------------------
! Incomplete gamma function                                                                 
! from Numerical Recipes in Fortran 77: The Art of                                          
! Scientific Computing                                                                      

      function gammq(a,x)

      double precision a,gammq,x

! USES gcf,gser                                                                             
! Returns the incomplete gamma function Q(a,x) = 1-P(a,x)                                   

      double precision gammcf,gamser,gln
!     if (x.lt.0..or.a.le.0) pause 'bad argument in gammq'                                  
      ! if (x.lt.0..or.a.le.0) print*, 'bad argument in gammq'
      if (x.lt.a+1.) then
         call gser(gamser,a,x,gln)
         gammq=1.-gamser
      else
         call gcf(gammcf,a,x,gln)
         gammq=gammcf
      end if
      return
      end

! ------------------------------------------

      subroutine gser(gamser,a,x,gln)
      integer itmax
      double precision a,gamser,gln,x,eps
      parameter(itmax=100,eps=3.e-7)
      integer n
      double precision ap,del,sum,gamma
      gln = log(gamma(a))
      if (x.le.0.) then
!        if (x.lt.0.) pause 'x < 0 in gser'                                                 
         ! if (x.lt.0.) print*, 'x < 0 in gser'
         gamser = 0.
         return
      end if
      ap=a
      sum=1./a
      del=sum
      do n=1,itmax
         ap=ap+1.
         del=del*x/ap
         sum=sum+del
         if (abs(del).lt.abs(sum)*eps) goto 1
      end do
!     pause 'a too large, itmax too small in gser'                                          
      ! print*, 'a too large, itmax too small in gser'
 1    gamser=sum*exp(-x+a*log(x)-gln)
      return
      end

! ------------------------------------------

      subroutine gcf(gammcf,a,x,gln)
      integer itmax
      double precision a,gammcf,gln,x,eps,fpmin
      parameter(itmax=100,eps=3.e-7,fpmin=1.e-30)
      integer i
      double precision an,b,c,d,del,h,gamma
      gln=log(gamma(a))
      b=x+1.-a
      c=1./fpmin
      d=1./b
      h=d
      do i=1,itmax
         an=-i*(i-a)
         b=b+2.
         d=an*d+b
         if(abs(d).lt.fpmin) d=fpmin
         c=b+an/c
         if(abs(c).lt.fpmin) c=fpmin
         d=1./d
         del=d*c
         h = h*del
         if(abs(del-1.).lt.eps)goto 1
      end do
!     pause 'a too large, itmax too small in gcf'                                           
      ! print*, 'a too large, itmax too small in gcf'
 1    gammcf=exp(-x+a*log(x)-gln)*h
      return
    end subroutine gcf

! ------------------------------------------

! Diagnose any moment for a bin from M0 and M3 in that bin.
elemental function diagnose_moment(x, m0, m3) result(mx)
  real(8), intent(in) :: x, m0, m3
  real(8) :: mx
  real(8) :: m3o0

  if (m0<=1d-14) then
    mx = 0.
  else
    m3o0 = m3 / m0
    mx = m0 * m3o0**(x/3.)
  endif

end function diagnose_moment

! ----------------------------------------------------------

subroutine random_stduniform(u)
   implicit none
   double precision,intent(out) :: u
   double precision :: r
   call random_number(r)
   u = 1 - r
end subroutine random_stduniform

! Randomly sample from a normal distribution
subroutine random_stdnormal(x)
   implicit none
   double precision,intent(out) :: x
   double precision,parameter :: pi=3.14159265
   double precision :: u1,u2
   call random_stduniform(u1)
   call random_stduniform(u2)
   x = sqrt(-2*log(u1))*cos(2*pi*u2)
end subroutine random_stdnormal

! ----------------------------------------------------------
 subroutine read_netcdf(variable, filename, varname)
     use netcdf
     implicit none

     ! Variable declarations
     integer :: ncid         ! NetCDF file ID
     integer :: varid        ! Variable ID
     integer :: retval       ! Return value for error handling
     integer :: dimids(2)    ! Dimension IDs for the variable
     integer :: var_dims(2)  ! Dimensions of the variable
     integer :: nx, ny, ndims! Dimensions sizes
     real(8), dimension(:,:), allocatable :: variable   ! Array to hold the variable data

     character(len=200) :: filename
     character(len=100) :: varname

     ! Open the NetCDF file (read-only mode)
     retval = nf90_open(filename, nf90_nowrite, ncid)
     if (retval /= nf90_noerr) then
       print *, 'Error opening NetCDF file'
       stop
     end if

     ! Get dimension lengths (assuming a 2D variable here)
     dimids = (/1,2/)
     retval = nf90_inquire_dimension(ncid, dimids(1), len=nx)
     if (retval /= nf90_noerr) then
       print *, 'Error getting dimension length for dim 1:', nf90_strerror(retval)
       retval = nf90_close(ncid)
       stop
     endif

     retval = nf90_inquire_dimension(ncid, dimids(2), len=ny)
     if (retval /= nf90_noerr) then
       print *, 'Error getting dimension length for dim 2:', nf90_strerror(retval)
       retval = nf90_close(ncid)
       stop
     endif

     retval = nf90_inq_varid(ncid, varname, varid)
     if (retval /= nf90_noerr) then
       print *, 'Error getting variable ID'
       retval = nf90_close(ncid)
       stop
     end if


     ! Allocate the array to hold the data
     allocate(variable(nx, ny))

     ! Read the variable data
     retval = nf90_get_var(ncid, varid, variable)
     if (retval /= nf90_noerr) then
       print *, 'Error reading variable data'
       retval = nf90_close(ncid)
       stop
     end if

     ! Close the NetCDF file
     retval = nf90_close(ncid)
     if (retval /= nf90_noerr) then
       print *, 'Error closing NetCDF file'
       stop
     end if

     ! Output the data (or process as needed)
     ! print *, 'Data read successfully:', variable

 end subroutine read_netcdf

! ----------------------------------------------------------

subroutine load_latinhc(n_perturbed_param, n_ppe)
use parameters, only: lsample, aero_N_init
use switches, only: wctrl
use namelists, only: irealz, Na_max, Na_min, w_max, w_min
  
integer, intent(in) :: n_perturbed_param, n_ppe
character(len=200) :: filename
character(len=100) :: varname
character(len=6) :: n_ppe_str, n_perturbed_param_str

write(n_perturbed_param_str,'(I0)') n_perturbed_param+2
write(n_ppe_str,'(I0)') n_ppe

filename = '/home/arthurhu/KiD_AMP/KiD_1mode_gam/lhs_module/lhs_out_p'// &
  trim(n_perturbed_param_str)//'_n'//trim(n_ppe_str)//'.nc'
print*, 'loading LHS:', filename
varname = 'lhs_sample'
call read_netcdf(lsample, filename, varname)

! draw a aero_N_init between Na_min and Na_max
aero_N_init(1) = lsample(1, irealz) * (Na_max - Na_min) + Na_min
! draw a wctrl(1) between w_min and w_max
wctrl(1) = lsample(2, irealz) * (w_max - w_min) + w_min

end subroutine load_latinhc
! ----------------------------------------------------------
end module global_fun
