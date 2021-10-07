module misc_fun
implicit none

contains
   subroutine diag_nu(nu_diag, momvalx, momvaly, momx, momy)
      implicit none
      integer, intent(in) :: momx,momy
      real(8), intent(in) :: momvalx, momvaly
      real(8), intent(out) :: nu_diag

      nu_diag=11.8*(1000* (momvaly/momvalx)**(1/real(momy-momx)) -0.7)**2 + 2.
      !nu_diag=max(minnu,min(maxnu,(0.0005714e-6*cx(k,1)+0.2714)**(-2.)))
   end subroutine diag_nu


end module misc_fun
