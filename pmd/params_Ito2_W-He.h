c.....Mass
      real(rp),parameter:: am_W   = real(183.84d0, rp)
      real(rp),parameter:: am_He  =   real(4.0026d0, rp)

c.....energy scaling
      real(rp),parameter:: sfac  = real(1d0, rp) !/2

c.....length scaling for hybrid QMCL calculation
      real(rp),parameter:: slen  = real(1d0, rp)!/3.204d0 *3.172d0
      real(rp),parameter:: aa2bs = slen

c.....2-body
c.....prefactor in eV*A, see the note on 2013-06-26
      real(rp),parameter:: p_fac = real(14.3997065885d0, rp) *aa2bs
c.....Z
      real(rp),parameter:: p_Z(1:2)= (/ real(74d0, rp), real(2d0, rp) /)
      real(rp),parameter:: p_alpha(1:2,1:2)=
     &     (/ real(3.942657d0, rp) /aa2bs,
     &        real(4.610618d0, rp) /aa2bs,
     &        real(4.610618d0, rp) /aa2bs,
     &        real(4.217128d0, rp) /aa2bs /)
      real(rp),parameter:: p_beta(1:2,1:2)=
     &     (/  real(0d0, rp) /aa2bs**2,
     &        -real(0.985010d0, rp) /aa2bs**2,
     &        -real(0.985010d0, rp) /aa2bs**2,
     &         real(7.056687d0, rp) /aa2bs**2 /)
      real(rp),parameter:: p_gamma(1:2,1:2)=
     &     (/  real(0d0, rp) /aa2bs**3,
     &         real(1.129080d0, rp) /aa2bs**3,
     &         real(1.129080d0, rp) /aa2bs**3,
     &        -real(3.096151d0, rp) /aa2bs**3 /)


c.....Many-body
      real(rp),parameter:: p_rs  = 3.270 *aa2bs
      real(rp),parameter:: p_rc  = 4.070 *aa2bs
c      real(rp),parameter:: p_B   = 2106.421590d0
      real(rp),parameter:: p_B   = real(45.8957687592d0, rp)
      real(rp),parameter:: p_c   = real(1.896978d0, rp) /aa2bs
