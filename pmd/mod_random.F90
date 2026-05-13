module random
  use mod_precision
  implicit none
  save
! rseed must be real(8): the LCG requires exact integer arithmetic up to ~2^31;
! in single precision 2147483647 and 2147483648 both round to 2^31, corrupting
! the generator and causing urnd() to return 1.0, which makes log(1-r1) = -Inf.
  real(8):: rseed  = 12345.0d0
  real(rp),parameter:: pi= 3.14159265358979_rp

contains
!=======================================================================
  subroutine set_seed(seed)
    real(rp),intent(in):: seed
    rseed = real(seed, 8)
  end subroutine set_seed
!=======================================================================
  function get_seed()
    real(rp):: get_seed
    get_seed = real(rseed, rp)
    return
  end function get_seed
!=======================================================================
  function urnd()
!
!  Uniform random number generator (Lehmer LCG, Park-Miller).
!  Internal arithmetic is always real(8): in single precision,
!  2147483647 and 2147483648 both round to 2^31, breaking the generator.
!
    real(rp):: urnd
    real(8),parameter:: d2p31m = 2147483647.0d0
    real(8),parameter:: d2p31  = 2147483648.0d0

    rseed = mod(16807.0d0*rseed, d2p31m)
    urnd  = real(rseed/d2p31, rp)
! In single precision, values near 1.0 round up to 1.0; clamp to avoid
! log(0) in callers such as box_muller.
    if (urnd >= 1.0_rp) urnd = nearest(1.0_rp, -1.0_rp)
    return
  end function urnd
!=======================================================================
  function polarbm()
!
!  Polar form of Box-Muller method.
!  See: http://www.design.caltech.edu/erik/Misc/Gaussian.html
!  Generate Gaussian distribution from two uniform random numbers.
!  Only one of two dependent random numbers is returned.
!
    real(rp):: polarbm
    real(rp):: x1,x2,w,y

    w = 2.0_rp
    do while( w.ge.1.0_rp )
      x1 = 2.0_rp*urnd()-1.0_rp
      x2 = 2.0_rp*urnd()-1.0_rp
      w = x1*x1 +x2*x2
    enddo
    y = sqrt( (-2.0_rp*log(w))/w)
    polarbm = x1*y
    return
  end function polarbm
!=======================================================================
  function box_muller()
!
!  Generate Gaussian distribution from two uniform random number.
!  Only one of two dependent random numbers is returned.
!
    real(rp):: box_muller
    real(rp):: r1,r2

    r1= urnd()
    r2= urnd()
    box_muller= sqrt(-2.0_rp*log(1.0_rp-r1)) *cos(2.0_rp*pi*r2)
    return
  end function box_muller
  
end module random
