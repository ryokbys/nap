      function urnd()
c
c  Uniform random number generator
c      
      implicit none 
      real(8):: urnd
      real(8),save:: dseed= 12345d0
      real(8),save:: d2p31m,d2p31
      data d2p31m/2147483647d0/
      data d2p31 /2147483648d0/
      
      dseed=dmod(16807d0*dseed,d2p31m)
      urnd=dseed/d2p31
      return
      end function urnd
c=======================================================================
      function box_muller()
c
c  Generate Gaussian distribution from two uniform random number.
c  Only one of two dependent random numbers is returned.
c
      implicit none
      real(8):: box_muller
      real(8),parameter:: pi= 3.14159265358979d0
      real(8):: r1,r2
      real(8),external:: urnd

      r1= urnd()
      r2= urnd()
      box_muller= sqrt(-2d0*dlog(r1)) *cos(2d0*pi*r2)
      return
      end function box_muller
