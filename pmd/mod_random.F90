module random
  use mod_precision
  implicit none
  save
  real(rp):: rseed  = 12345.0_rp
  real(rp),parameter:: pi= 3.14159265358979_rp

contains
!=======================================================================
  subroutine set_seed(seed)
    real(rp),intent(in):: seed
    rseed = seed
  end subroutine set_seed
!=======================================================================
  function get_seed()
    real(rp):: get_seed
    get_seed = rseed
    return
  end function get_seed
!=======================================================================
  function urnd()
!
!  Uniform random number generator
!      
    real(rp):: urnd
    real(rp),save:: d2p31m,d2p31
    data d2p31m/2147483647.0_rp/
    data d2p31 /2147483648.0_rp/

    rseed=mod(16807.0_rp*rseed,d2p31m)
    urnd=rseed/d2p31
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
