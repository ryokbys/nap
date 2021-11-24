module random
  implicit none
  save
  real(8):: rseed  = 12345d0
  real(8),parameter:: pi= 3.14159265358979d0

contains
!=======================================================================
  subroutine set_seed(seed)
    real(8),intent(in):: seed
    rseed = seed
  end subroutine set_seed
!=======================================================================
  function get_seed()
    real(8):: get_seed
    get_seed = rseed
    return
  end function get_seed
!=======================================================================
  function urnd()
!
!  Uniform random number generator
!      
    real(8):: urnd
    real(8),save:: d2p31m,d2p31
    data d2p31m/2147483647d0/
    data d2p31 /2147483648d0/

    rseed=dmod(16807d0*rseed,d2p31m)
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
    real(8):: polarbm
    real(8):: x1,x2,w,y

    w = 2d0
    do while( w.ge.1d0 )
      x1 = 2d0*urnd()-1d0
      x2 = 2d0*urnd()-1d0
      w = x1*x1 +x2*x2
    enddo
    y = sqrt( (-2d0*log(w))/w)
    polarbm = x1*y
    return
  end function polarbm
!=======================================================================
  function box_muller()
!
!  Generate Gaussian distribution from two uniform random number.
!  Only one of two dependent random numbers is returned.
!
    real(8):: box_muller
    real(8):: r1,r2

    r1= urnd()
    r2= urnd()
    box_muller= sqrt(-2d0*log(1d0-r1)) *cos(2d0*pi*r2)
    return
  end function box_muller
  
end module random
