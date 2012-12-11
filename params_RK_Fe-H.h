c=======================================================================
c  EAM potential of Ackland model for Fe-Fe (iron)
c   See Philos. Mag. 83(35), 3977-3994 (2003)
c  Note that raw values are in AA and eV.
c-----------------------------------------------------------------------
c  Extended to Fe-H system
c    See Ramasubramaniam et al. PRB 79, 174101 (2009), (Potential A)
c-----------------------------------------------------------------------
c  Parameters were modified to fit to lattice constant obtained by VASP
c  for the hybrid QM-CL calculation.
c=======================================================================
      real(8),parameter:: am_fe = 55.847d0 *aump/aume
      real(8),parameter:: am_h  =  1.008d0 *aump/aume
c      real(8),parameter:: alcfe = 2.8553d-10 /bohr ! Ackland original
      real(8),parameter:: alcfe = 2.835d-10 /bohr
c.....r correction factor for RK modification
      real(8),parameter:: rfac  = 2.835d0 /2.8553d0
      real(8),parameter:: rfaci = 1d0 /rfac
      real(8),parameter:: z_fe  = 26d0
      real(8),parameter:: z_h   = 1d0
c.....prefactor for rs
      real(8),parameter:: a_rs  = 0.88534d0
c.....prefactor for varphi r < r1
      real(8),parameter:: qe    = 3.794701096299472d0 /ehrt *aa2bohr
      real(8),parameter:: z2q2  = 9.7342365892908E+03 /ehrt *aa2bohr
      real(8),parameter:: z2q2_feh= z_fe*z_h*qe*qe /ehrt *aa2bohr
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

c-----------------------------------------------------------------------
c  For Fe-H potential parameters
c-----------------------------------------------------------------------
c.....Range in phi_FeH
      real(8),parameter:: r1_feh =  0.6d-10 /bohr
      real(8),parameter:: r2_feh =  1.2d-10 /bohr
c.....Bs for FeH
      real(8),parameter:: b0_feh = 853.4769964964161d0
      real(8),parameter:: b1_feh =-4206.406420131467d0 /aa2bohr
      real(8),parameter:: b2_feh = 8686.215689507188d0 /aa2bohr**2
      real(8),parameter:: b3_feh =-9137.341019760202d0 /aa2bohr**3
      real(8),parameter:: b4_feh = 4807.823405345844d0 /aa2bohr**4
      real(8),parameter:: b5_feh =-1002.904058496011d0 /aa2bohr**5
c.....a_phi_feh
      real(8),parameter:: a_phi_feh(1:7)=
     &     (/ 14.0786236789212005d0  /ehrt /aa2bohr**3 !1
     &     ,  -4.4526835887173704d0  /ehrt /aa2bohr**3
     &     ,   5.5025121262565992d0  /ehrt /aa2bohr**3
     &     ,  -1.0687489808214079d0  /ehrt /aa2bohr**3
     &     ,  -0.3461498208163201d0  /ehrt /aa2bohr**3
     &     ,  -0.0064991947759021d0  /ehrt /aa2bohr**3
     &     ,  -0.0357435602984102d0  /ehrt /aa2bohr**3 !7
     &     /)
c.....r_phi_feh
      real(8),parameter:: r_phi_feh(1:7)=
     &     (/ 1.6d0 *aa2bohr !1
     &     ,  1.7d0 *aa2bohr
     &     ,  1.8d0 *aa2bohr
     &     ,  2.0d0 *aa2bohr
     &     ,  2.5d0 *aa2bohr
     &     ,  3.2d0 *aa2bohr
     &     ,  4.2d0 *aa2bohr !7
     &     /)
c.....Cut off for phi_HH
      real(8),parameter:: rc_phi_hh  = 2.3d-10 /bohr
      real(8),parameter:: r0_hh      = 0.74d-10 /bohr
c.....S(r) for phi_HH
      real(8),parameter:: r_tanh_hh  = 0.9d0 *aa2bohr
      real(8),parameter:: a_tanh_hh  = 25d0 /aa2bohr
c.....Coeff for phi_HH
      real(8),parameter:: c1_phi_hh  = 0.0d0 /ehrt
      real(8),parameter:: c2_phi_hh  = 0.0d0 /ehrt
      real(8),parameter:: eb_hh      = 2.37d0 /ehrt
      real(8),parameter:: almbd_hh   = 0.4899d0
c.....a^F parameters
      real(8),parameter:: a_f(1:6)=
     &     (/ -0.0581256120818134d0 /ehrt
     &     ,   0.0022854552833736d0 /ehrt
     &     ,  -0.0000314202805805d0 /ehrt
     &     ,   0.0000013764132084d0 /ehrt
     &     ,  -0.0000000253707731d0 /ehrt
     &     ,   0.0000000001483685d0 /ehrt
     &     /)
c.....Rho_FeH
      real(8),parameter:: a_rho_feh(1:6)=
     &     (/ 10.0073629216300581d0 /aa2bohr**3
     &     ,  32.4861983261490295d0 /aa2bohr**3
     &     ,  -0.9494226032063788d0 /aa2bohr**3
     &     ,  11.6659812262450338d0 /aa2bohr**3
     &     ,  -0.0147080251458273d0 /aa2bohr**3
     &     ,   0.4943383753319843d0 /aa2bohr**3
     &     /)
      real(8),parameter:: r_rho_feh(1:6)=
     &     (/ 1.6d0 *aa2bohr
     &     ,  1.8d0 *aa2bohr
     &     ,  2.0d0 *aa2bohr
     &     ,  2.4d0 *aa2bohr
     &     ,  3.2d0 *aa2bohr
     &     ,  4.2d0 *aa2bohr
     &     /)
c.....Rho_HFe
      real(8),parameter:: a_rho_hfe(1:5)=
     &     (/ 11.1667357634216433d0 /aa2bohr**3
     &     ,  -3.0351307365078730d0 /aa2bohr**3
     &     ,   3.6096144794370653d0 /aa2bohr**3
     &     ,   0.0212509034775648d0 /aa2bohr**3
     &     ,   0.030391493994625d0  /aa2bohr**3
     &     /)
      real(8),parameter:: r_rho_hfe(1:5)=
     &     (/ 1.5d0 *aa2bohr
     &     ,  2.0d0 *aa2bohr
     &     ,  2.5d0 *aa2bohr
     &     ,  3.0d0 *aa2bohr
     &     ,  4.2d0 *aa2bohr
     &     /)
c.....Rho_HH
      real(8),parameter:: c_rho_hh = 1800d0 /aa2bohr**2
