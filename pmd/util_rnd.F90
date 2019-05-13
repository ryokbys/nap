function urnd(dseed0)
!
!  Uniform random number generator
!      
  implicit none
  real(8),intent(in),optional:: dseed0
  real(8):: urnd
  real(8),save:: dseed= 12345d0
  real(8),save:: d2p31m,d2p31
  data d2p31m/2147483647d0/
  data d2p31 /2147483648d0/

  if( present(dseed0) ) dseed = dseed0

  dseed=dmod(16807d0*dseed,d2p31m)
  urnd=dseed/d2p31
  return
end function urnd
!=======================================================================
subroutine sub_urnd(dseed,urnd)
!     
!     Subroutine version of uniform random number generator.
!     The input, DSEED, must be kept untouched outside this routine.
!     
  implicit none
  real(8),intent(inout):: dseed
  real(8),intent(out):: urnd

  real(8),save:: d2p31m,d2p31
  data d2p31m/2147483647d0/
  data d2p31 /2147483648d0/

  dseed=dmod(16807d0*dseed,d2p31m)
  urnd=dseed/d2p31
  return
end subroutine sub_urnd
!=======================================================================
function box_muller()
!
!  Generate Gaussian distribution from two uniform random number.
!  Only one of two dependent random numbers is returned.
!
  implicit none
  real(8):: box_muller
  real(8),parameter:: pi= 3.14159265358979d0
  real(8):: r1,r2
  interface
    function urnd(dseed0)
      real(8),intent(in),optional:: dseed0
      real(8):: urnd
    end function urnd
  end interface

  r1= urnd()
  r2= urnd()
  box_muller= sqrt(-2d0*dlog(1d0-r1)) *cos(2d0*pi*r2)
  return
end function box_muller
