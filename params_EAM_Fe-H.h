c=======================================================================
c  EAM potential of Ackland model for Fe-Fe (iron)
c   See Philos. Mag. 83(35), 3977-3994 (2003)
c  Note that raw values are in AA and eV.
c-----------------------------------------------------------------------
c  This will be extended to Fe-H system
c  See Ramasubramaniam et al. PRB 79, 174101 (2009)
c=======================================================================
      real(8),parameter:: am_fe = 55.847d0 *aump/aume
c      real(8),parameter:: alcfe = 2.8665d-10 /bohr
      real(8),parameter:: alcfe = 2.8553d-10 /bohr
      real(8),parameter:: z_fe  = 26d0
c.....prefactor for rs
      real(8),parameter:: a_rs  = 0.88534d0
c.....prefactor for varphi r < r1
      real(8),parameter:: z2q2  = 9.7342365892908E+03 /ehrt *aa2bohr
c.....Range in varphi
      real(8),parameter:: r1  = 1.0d-10 /bohr
      real(8),parameter:: r2  = 2.0d-10 /bohr
c.....Bs for r1 < r <= r2
      real(8),parameter:: b0  = 6.4265260576348d0
      real(8),parameter:: b1  = 1.7900488524286d0 /aa2bohr
      real(8),parameter:: b2  =-4.5108316729807d0 /aa2bohr**2
      real(8),parameter:: b3  = 1.0866199373306d0 /aa2bohr**3
c.....varphi for r2 < r
      real(8),parameter:: a_vphi(1:15)=
     &     (/  0.0d0
     &     , -24.028204854115d0   /ehrt /aa2bohr**3 ! 2
     &     ,  11.300691696477d0   /ehrt /aa2bohr**3 ! 3
     &     ,   5.3144495820462d0  /ehrt /aa2bohr**3 ! 4
     &     ,  -4.6659532856049d0  /ehrt /aa2bohr**3 ! 5
     &     ,   5.9637758529194d0  /ehrt /aa2bohr**3
     &     ,  -1.7710262006061d0  /ehrt /aa2bohr**3
     &     ,   0.85913830768731d0 /ehrt /aa2bohr**3
     &     ,  -2.1845362968261d0  /ehrt /aa2bohr**3
     &     ,   2.6424377007466d0  /ehrt /aa2bohr**3 ! 10
     &     ,  -1.0358345370208d0  /ehrt /aa2bohr**3
     &     ,   0.33548264951582d0 /ehrt /aa2bohr**3
     &     ,  -4.6448582149334d-2 /ehrt /aa2bohr**3
     &     ,  -7.0294963048689d-3 /ehrt /aa2bohr**3 ! 14
     &     ,   0.0d0 /)
      real(8),parameter:: r_vphi(1:15)=
     &     (/  0.0d0
     &     ,   2.2d0 *aa2bohr ! 2
     &     ,   2.3d0 *aa2bohr
     &     ,   2.4d0 *aa2bohr
     &     ,   2.5d0 *aa2bohr ! 5
     &     ,   2.6d0 *aa2bohr
     &     ,   2.7d0 *aa2bohr
     &     ,   2.8d0 *aa2bohr
     &     ,   3.0d0 *aa2bohr
     &     ,   3.3d0 *aa2bohr ! 10
     &     ,   3.7d0 *aa2bohr
     &     ,   4.2d0 *aa2bohr
     &     ,   4.7d0 *aa2bohr
     &     ,   5.3d0 *aa2bohr ! 14
     &     ,   0.0d0 /)
c.....psi for calculation of rho
      real(8),parameter:: a_psi(1:3)=
     &     (/ 11.686859407970d0     /aa2bohr**3
     &     ,  -0.014710740098830d0  /aa2bohr**3
     &     ,   0.47193527075943d0   /aa2bohr**3  /)
      real(8),parameter:: r_psi(1:3)=
     &     (/  2.4d0 *aa2bohr
     &     ,   3.2d0 *aa2bohr
     &     ,   4.2d0 *aa2bohr /)
c.....coeff for the embedded func
      real(8),parameter:: a_emb  = -3.5387096579929d-4
c.....cutoffs
      real(8),parameter:: rc_rho = 4.2d0 *aa2bohr
      real(8),parameter:: rc_vphi= 5.3d0 *aa2bohr
