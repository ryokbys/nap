c.....Mass
      real(8),parameter:: am_W   = 183.84d0 *aump/aume
      real(8),parameter:: am_He  =   4.0026d0 *aump/aume

c.....energy scaling
      real(8),parameter:: sfac  = 1d0 !/2

c.....length scaling for hybrid QMCL calculation
      real(8),parameter:: slen  = 1d0!/3.204d0 *3.172d0
      real(8),parameter:: aa2bs = aa2bohr*slen

c.....2-body
c.....prefactor in eV*A, see the note on 2013-06-26
      real(8),parameter:: p_fac = 14.3997065885d0 *ev2hrt *aa2bs
c.....Z
      real(8),parameter:: p_Z(1:2)= (/ 74d0, 2d0 /)
      real(8),parameter:: p_alpha(1:2,1:2)=
     &     (/ 3.942657d0 /aa2bs,
     &        4.610618d0 /aa2bs,
     &        4.610618d0 /aa2bs,
     &        4.217128d0 /aa2bs /)
      real(8),parameter:: p_beta(1:2,1:2)=
     &     (/  0d0 /aa2bs**2,
     &        -0.985010d0 /aa2bs**2,
     &        -0.985010d0 /aa2bs**2,
     &         7.056687d0 /aa2bs**2 /)
      real(8),parameter:: p_gamma(1:2,1:2)=
     &     (/  0d0 /aa2bs**3,
     &         1.129080d0 /aa2bs**3,
     &         1.129080d0 /aa2bs**3,
     &        -3.096151d0 /aa2bs**3 /)


c.....Many-body
      real(8),parameter:: p_rs  = 3.270 *aa2bs
      real(8),parameter:: p_rc  = 4.070 *aa2bs
c      real(8),parameter:: p_B   = 2106.421590d0
      real(8),parameter:: p_B   = 45.8957687592d0 *ev2hrt
      real(8),parameter:: p_c   = 1.896978d0 /aa2bs
