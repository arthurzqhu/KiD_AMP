module sort_mod
   use, intrinsic :: iso_fortran_env, only: IK=>int32, RK=>real64
   implicit none
contains
   ! modified from https://stackoverflow.com/questions/54860293/ordering-function-in-fortran
   Subroutine  getSortArrayIdx(n, array, idx, tol_sigfig_opt)
      integer(IK), intent(in)           :: n
      real(RK)   , intent(in)           :: array(n)
      integer    , intent(in), optional :: tol_sigfig_opt
      integer(IK), intent(out)          :: idx(n)
      integer                           :: tol_sigfig, iidx, idx_s, idx_e
      real(RK)   , allocatable          :: array_alloc(:)
      integer(IK), allocatable          :: idx_array(:), idx_idx(:)

      if (.not. present(tol_sigfig_opt)) then
         tol_sigfig = 0
      else
         tol_sigfig = tol_sigfig_opt
      endif

      call indexArrayReal(n, array, idx)

      iidx = 1
      array_alloc = array(idx)
      do while (iidx < n)
         call find_closest_val_in_arr(array_alloc, array_alloc(iidx), &
            idx_s, idx_e, cmp_method='sigfig', tol_sigfig_opt=tol_sigfig)

         if (idx_e > idx_s) then
            idx_array = idx(idx_s:idx_e)
            allocate(idx_idx(idx_e - idx_s + 1))
            call indexArrayReal(idx_e - idx_s + 1, dble(idx_array), idx_idx)
            idx_array = idx_array(idx_idx)
            idx(idx_s:idx_e) = idx_array
            deallocate(idx_idx)
         endif

         iidx = idx_e + 1
      end do


   End Subroutine getSortArrayIdx 

   ! ========================================
   subroutine indexArrayReal(n, Array, Index)
      implicit none
      integer(IK), intent(in)           :: n
      real(RK)   , intent(in)           :: Array(n)
      integer(IK), intent(out)          :: Index(n)
      integer(IK), parameter            :: nn=15, nstack=50
      integer(IK)                       :: k, i, j, indext, jstack, l, r
      integer(IK)                       :: istack(nstack)
      real(RK)                          :: a

      do j = 1, n
         Index(j) = j
      end do

      jstack = 0
      l = 1
      r = n

      do
         if (r-l < nn) then
            do j = l+1, r
               indext = Index(j)
               a = Array(indext)
               do i = j-1, l, -1
                  if (Array(Index(i)) <= a) exit
                  Index(i+1) = Index(i)
               end do
               Index(i+1) = indext
            end do
            if (jstack == 0) return
            r = istack(jstack)
            l = istack(jstack-1)
            jstack = jstack-2
         else
            k = (l+r)/2
            call swap(Index(k), Index(l+1))
            call exchangeIndex(Index(l), Index(r))
            call exchangeIndex(Index(l+1), Index(r))
            call exchangeIndex(Index(l), Index(l+1))
            i = l+1
            j = r
            indext = Index(l+1)
            a = Array(indext)
            do
               do
                  i = i+1
                  if (Array(Index(i)) >= a) exit
               end do
               do
                  j = j-1
                  if (Array(Index(j)) <= a) exit
               end do
               if (j < i) exit
               call swap(Index(i), Index(j))
            end do
            Index(l+1) = Index(j)
            Index(j) = indext
            jstack = jstack+2
            if (jstack > nstack) then
               write(*, *) 'NSTACK too small in indexArrayReal()'   ! xxx
               error stop
            end if
            if (r-i+1 >= j-l) then
               istack(jstack) = r
               istack(jstack-1) = i
               r = j-1
            else
               istack(jstack) = j-1
               istack(jstack-1) = l
               l = i
            end if
         end if
      end do
   contains
      subroutine exchangeIndex(i, j)
         integer(IK), intent(inout) :: i, j
         integer(IK)                :: swp
         if (Array(j) < Array(i)) then
            swp = i
            i = j
            j = swp
         end if
      end subroutine exchangeIndex
      pure elemental subroutine swap(a, b)
         implicit none
         integer(IK), intent(inout) :: a, b
         integer(IK) :: dum
         dum = a
         a = b
         b = dum
      end subroutine swap
   end subroutine indexArrayReal

   ! ------------------------------------------
   logical function are_close(val1, val2, cmp_method, tol_rat_opt, tol_sigfig_opt)
      ! compares val1 and var2 and returns true if they are close enough in value
      ! cmp_method = 'sigfig' or 'ratio'
      ! val1 here is treated as the standard/benchmark if cmp_method == 'ratio'
      implicit none

      real(8), intent(in):: val1, val2
      real(8), intent(in), optional:: tol_rat_opt
      integer, intent(in), optional:: tol_sigfig_opt
      real(8):: tol_rat
      integer:: tol_sigfig
      character(len=*), intent(in):: cmp_method
      character(len=50):: sigfig_fmt, str_val1, str_val2

      if (present(tol_sigfig_opt)) tol_sigfig = tol_sigfig_opt
      if (present(tol_rat_opt)) tol_rat = tol_rat_opt

      if (trim(cmp_method) .eq. 'ratio') then
         if (.not. present(tol_rat_opt)) then
            print*, '### Error in are_close()'
            print*, 'Must pass in optional argument tol_rat_opt = [] as real'
            stop
         endif

         if (abs(log(val2/val1)) <= log(1. + tol_rat)) then
            are_close = .true.
         else
            are_close = .false.
         endif

      elseif (trim(cmp_method) .eq. 'sigfig') then
         if (.not. present(tol_sigfig_opt)) then
            print*, '### Error in are_close()'
            print*, 'Must pass in optional argument tol_sigfig_opt = [] as integer'
            stop
         endif

         sigfig_fmt = '('//'e'//itoa(tol_sigfig+8)//'.'//itoa(tol_sigfig)//')'
         write(str_val1, trim(sigfig_fmt)) val1
         write(str_val2, trim(sigfig_fmt)) val2
         if (str_val1 .eq. str_val2) then
            are_close = .true.
         else
            are_close = .false.
         endif
      else
         print*, 'cmp_method = ', trim(cmp_method)
         print*, '### Error in are_close()'
         print*, 'the cmp_method argument must be either "ratio" or "sigfig"'
         stop
      endif

   end function

   ! ------------------------------------------
   function itoa(i) result(res)
      character(:), allocatable:: res
      integer, intent(in):: i
      character(range(i)+2):: tmp
      write(tmp, '(i0)') i
      res = trim(tmp)
   end function

   ! ------------------------------------------
   subroutine find_closest_val_in_arr(arr, val, idx_s, idx_e, cmp_method, tol_rat_opt, tol_sigfig_opt)

   ! if cmp_method == 'ratio', find the two indices in a sorted array 
   ! that bound the target value +/- tolerance with binary search.
   ! if cmp_method == 'sigfig', find the two indices that bound all 
   ! values with the same sig fig.
   ! if val is outside the tolerance range, or when tol_rat == 0, 
   ! idx_s = idx_e = the index for the closest value

   ! WARNING: input array `arr` must be sorted

   real(8), allocatable, dimension(:), intent(in):: arr
   real(8), intent(in):: val
   character(*), intent(in):: cmp_method
   real(8), intent(in), optional:: tol_rat_opt
   integer, intent(in), optional:: tol_sigfig_opt
   integer, intent(out):: idx_s, idx_e
   integer:: idxL, idxR, idxM, idx_T 
   integer:: tol_sigfig
   logical:: l_close
   real(8):: tol_rat
   real(8):: minv, maxv  ! min and max value within the tolerated range 
                         ! so that arr(idx_s-1) < minv <= arr(idx_s) and 
                         ! arr(idx_e) <= maxv < arr(idx_e+1)
   character(len=20):: valstr_sigfig, sigfig_fmt

   if (trim(cmp_method) .eq. 'ratio') then
      if (.not. present(tol_rat_opt)) then
         print*, '### Error in find_closest_val_in_arr()'
         print*, 'Must pass in optional argument tol_rat_opt=[] as real'
         stop
      endif
      tol_rat = tol_rat_opt

      minv = val*(1-tol_rat)
      maxv = val*(1+tol_rat)
   elseif (trim(cmp_method) .eq. 'sigfig') then
      if (.not. present(tol_sigfig_opt)) then
         print*, '### Error in find_closest_val_in_arr()'
         print*, 'Must pass in optional argument tol_sigfig_opt=[] as integer'
         stop
      endif
      tol_sigfig = tol_sigfig_opt

      sigfig_fmt = '(e'//itoa(tol_sigfig+8)//'.'//itoa(tol_sigfig)//')'
      write(valstr_sigfig, sigfig_fmt) val
   else
      print*, '### Error in find_closest_val_in_arr()'
      print*, "cmp_method argument must be either 'ratio' or 'sigfig'"
      stop
   endif

   idxL = 1
   idxR = size(arr)
   idxM = floor((idxL+idxR)/2.)


   if (trim(cmp_method) .eq. 'sigfig') then

      ! otherwise, first find the closest index
      do while (idxR > idxL+1)  ! while idxL and idxR are more than 1 apart
         if (arr(idxM) < val) then
            idxL = idxM
            idxM = ceiling((idxL+idxR)/2.)  ! ceiling here and floor below to avoid skipping indices
         else
            idxR = idxM
            idxM = floor((idxL+idxR)/2.)
         endif
         idx_T = idxR
      enddo

      ! then make sure idx_T is closest to the target value
      if (abs(arr(idxL)-val) >= abs(arr(idxR)-val)) then
         idx_T = idxR
      else
         idx_T = idxL
      endif

      ! finally find the first and the last index that have the same sigfig as
      ! the target value
      idxL = 1; idxR = idx_T; idxM = floor((idxL+idxR)/2.)
      do while (idxR > idxL+1)
         ! close enough to the target value, or as close as the best we can find in the array
         l_close = are_close(arr(idxM), val, 'sigfig', tol_sigfig_opt=tol_sigfig) & 
                        .or. are_close(arr(idxM), arr(idx_T), 'sigfig', tol_sigfig_opt=tol_sigfig)
         if (.not. l_close) then
            idxL = idxM
            idxM = ceiling((idxL+idxR)/2.)  ! ceiling here and floor below to avoid skipping indices
         else
            idxR = idxM
            idxM = floor((idxL+idxR)/2.)
         endif
         idx_s = idxR
      end do   
      ! check idxL so that it's not left out, again, close enough or as close as we got
      if (are_close(arr(idxL), val, 'sigfig', tol_sigfig_opt=tol_sigfig) .or. &
            are_close(arr(idxL), arr(idx_T), 'sigfig', tol_sigfig_opt=tol_sigfig)) then
         idx_s = idxL
      else ! in case val is exactly in between arr(idxL) and arr(idxR), which might skip the do loop above
         idx_s = idxR
      endif

      idxL = idx_T; idxR = size(arr); idxM = floor((idxL+idxR)/2.)
      do while (idxR > idxL+1)
         l_close = are_close(arr(idxM), val, 'sigfig', tol_sigfig_opt=tol_sigfig) & 
                        .or. are_close(arr(idxM), arr(idx_T), 'sigfig', tol_sigfig_opt=tol_sigfig)
         if (l_close) then
            idxL = idxM
            idxM = ceiling((idxL+idxR)/2.)  ! ceiling here and floor below to avoid skipping indices
         else
            idxR = idxM
            idxM = floor((idxL+idxR)/2.)
         endif
         idx_e = idxL
      enddo
      ! check idxR so that it's not left out
      if (are_close(arr(idxR), val, 'sigfig', tol_sigfig_opt=tol_sigfig) .or. &
            are_close(arr(idxR), arr(idx_T), 'sigfig', tol_sigfig_opt=tol_sigfig)) then
         idx_e = idxR
      else
         idx_e = idxL
      endif

      return

   endif ! comparing sigfig


   if (trim(cmp_method) .eq. 'ratio') then
      if (arr(idxL) >= minv) then
         idx_s = idxL
      else
         do while (idxR > idxL+1)  ! while idxL and idxR are more than 1 apart
            if (arr(idxM) < minv) then
               idxL = idxM
               idxM = ceiling((idxL+idxR)/2.)  ! ceiling here and floor below to avoid skipping indices
            else
               idxR = idxM
               idxM = floor((idxL+idxR)/2.)
            endif
            idx_s = idxR
         enddo
      endif


      if (tol_rat > 0) then
         idxL = idx_s
         idxR = size(arr)
         idxM = floor((idxL+idxR)/2.)

         if (idxL == idxR) then
            idx_e = idxR
            return
         endif

         do while (idxR > idxL+1)
            if (arr(idxM) <= maxv) then
               idxL = idxM
               idxM = ceiling((idxL+idxR)/2.)  ! ceiling here and floor below to avoid skipping indices
            else
               idxR = idxM
               idxM = floor((idxL+idxR)/2.)
            endif
            idx_e = idxL
         enddo 

      elseif (tol_rat == 0.) then
         ! make final comparison of idxL and idxR see which one is closer to minv
         if (abs(arr(idxL)-val) > abs(arr(idxR)-val)) then
            idx_s = idxR
            idx_e = idxR
         else
            idx_s = idxL
            idx_e = idxL
         endif
      endif
   endif

   end subroutine find_closest_val_in_arr
end module sort_mod
