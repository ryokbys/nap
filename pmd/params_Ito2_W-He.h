c.....Mass
      real(rp),parameter:: am_W   = 183.84_rp
      real(rp),parameter:: am_He  =   4.0026_rp

c.....energy scaling
      real(rp),parameter:: sfac  = 1.0_rp !/2

c.....length scaling for hybrid QMCL calculation
      real(rp),parameter:: slen  = 1.0_rp!/3.204d0 *3.172d0
      real(rp),parameter:: aa2bs = slen

c.....2-body
c.....prefactor in eV*A, see the note on 2013-06-26
      real(rp),parameter:: p_fac = 14.3997065885_rp *aa2bs
c.....Z
      real(rp),parameter:: p_Z(1:2)= (/ 74.0_rp, 2.0_rp /)
      real(rp),parameter:: p_alpha(1:2,1:2)=
     &     (/ 3.942657_rp /aa2bs,
     &        4.610618_rp /aa2bs,
     &        4.610618_rp /aa2bs,
     &        4.217128_rp /aa2bs /)
      real(rp),parameter:: p_beta(1:2,1:2)=
     &     (/  0.0_rp /aa2bs**2,
     &        -0.985010_rp /aa2bs**2,
     &        -0.985010_rp /aa2bs**2,
     &         7.056687_rp /aa2bs**2 /)
      real(rp),parameter:: p_gamma(1:2,1:2)=
     &     (/  0.0_rp /aa2bs**3,
     &         1.129080_rp /aa2bs**3,
     &         1.129080_rp /aa2bs**3,
     &        -3.096151_rp /aa2bs**3 /)


c.....Many-body
      real(rp),parameter:: p_rs  = 3.270 *aa2bs
      real(rp),parameter:: p_rc  = 4.070 *aa2bs
c      real(rp),parameter:: p_B   = 2106.421590_rp
      real(rp),parameter:: p_B   = 45.8957687592_rp
      real(rp),parameter:: p_c   = 1.896978_rp /aa2bs
