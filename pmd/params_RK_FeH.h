!=======================================================================
!  EAM potential of Ackland model for Fe-Fe (iron)
!   See Philos. Mag. 83(35), 3977-3994 (2003)
!  Note that raw values are in Ang and eV.
!-----------------------------------------------------------------------
!  Extended to Fe-H system
!    See Ramasubramaniam et al. PRB 79, 174101 (2009), (Potential A)
!-----------------------------------------------------------------------
!  Parameters were modified to fit to lattice constant obtained by VASP
!  for the hybrid QM-CL calculation.
!  Cutoff radius should be 5.26232Ang. (originally 5.3 Ang.)
!=======================================================================
      real(rp),parameter:: am_fe = 55.847_rp
      real(rp),parameter:: am_h  =  1.008_rp
!      real(rp),parameter:: alcfe = 2.8553d0 ! Ackland original
      real(rp),parameter:: alcfe = 2.835_rp
!.....length scaling factor for hybrid QMCL calc
      real(rp),parameter:: slen  = 1.0_rp /2.8553_rp *2.835_rp
      real(rp),parameter:: rfac  = 2.835_rp /2.8553_rp
      real(rp),parameter:: rfaci = 1.0_rp /rfac
      real(rp),parameter:: z_fe  = 26.0_rp
      real(rp),parameter:: z_h   = 1.0_rp
!.....prefactor for rs
      real(rp),parameter:: a_rs  = 0.88534_rp
!.....prefactor for varphi r < r1
      real(rp),parameter:: qe    = 3.794701096299472_rp *slen
      real(rp),parameter:: z2q2  = 9.7342365892908E+03 *slen
      real(rp),parameter:: z2q2_feh= z_fe*z_h*qe*qe *slen
!.....Range in varphi
      real(rp),parameter:: r1  = 1.0_rp *slen
      real(rp),parameter:: r2  = 2.0_rp *slen
!.....Bs for r1 < r <= r2
      real(rp),parameter:: b0  = 6.4265260576348_rp
      real(rp),parameter:: b1  = 1.7900488524286_rp /slen
      real(rp),parameter:: b2  =-4.5108316729807_rp /slen**2
      real(rp),parameter:: b3  = 1.0866199373306_rp /slen**3
!.....varphi for r2 < r
      real(rp),parameter:: a_vphi(1:15)= &
     &     (/  0.0_rp &
     &     , -24.028204854115_rp   /slen**3 &! 2
     &     ,  11.300691696477_rp   /slen**3 &! 3
     &     ,   5.3144495820462_rp  /slen**3 &! 4
     &     ,  -4.6659532856049_rp  /slen**3 &! 5
     &     ,   5.9637758529194_rp  /slen**3 &
     &     ,  -1.7710262006061_rp  /slen**3 &
     &     ,   0.85913830768731_rp /slen**3 &
     &     ,  -2.1845362968261_rp  /slen**3 &
     &     ,   2.6424377007466_rp  /slen**3 &! 10
     &     ,  -1.0358345370208_rp  /slen**3 &
     &     ,   0.33548264951582_rp /slen**3 &
     &     ,  -4.6448582149334e-2_rp /slen**3 &
     &     ,  -7.0294963048689e-3_rp /slen**3 &! 14
     &     ,   0.0_rp /)
      real(rp),parameter:: r_vphi(1:15)= &
     &     (/  0.0_rp &
     &     ,   2.2_rp *slen & ! 2
     &     ,   2.3_rp *slen & 
     &     ,   2.4_rp *slen & 
     &     ,   2.5_rp *slen & ! 5
     &     ,   2.6_rp *slen & 
     &     ,   2.7_rp *slen & 
     &     ,   2.8_rp *slen & 
     &     ,   3.0_rp *slen & 
     &     ,   3.3_rp *slen & ! 10
     &     ,   3.7_rp *slen & 
     &     ,   4.2_rp *slen & 
     &     ,   4.7_rp *slen & 
     &     ,   5.3_rp *slen & ! 14
     &     ,   0.0_rp /)
!.....psi for calculation of rho
      real(rp),parameter:: a_psi(1:3)= &
     &     (/ 11.686859407970_rp     /slen**3 &
     &     ,  -0.014710740098830_rp  /slen**3 &
     &     ,   0.47193527075943_rp   /slen**3  /)
      real(rp),parameter:: r_psi(1:3)= &
     &     (/  2.4_rp *slen &
     &     ,   3.2_rp *slen &
     &     ,   4.2_rp *slen /)
!.....coeff for the embedded func
      real(rp),parameter:: a_emb  = -3.5387096579929e-4_rp
!.....cutoffs
      real(rp),parameter:: rc_rho = 4.2_rp *slen
      real(rp),parameter:: rc_vphi= 5.3_rp *slen

!-----------------------------------------------------------------------
!  For Fe-H potential parameters
!-----------------------------------------------------------------------
!.....Range in phi_FeH
      real(rp),parameter:: r1_feh =  0.6_rp
      real(rp),parameter:: r2_feh =  1.2_rp
!.....Bs for FeH
      real(rp),parameter:: b0_feh = 1242.1709168218642_rp
      real(rp),parameter:: b1_feh =-6013.566711223783_rp 
      real(rp),parameter:: b2_feh = 12339.540893927151_rp
      real(rp),parameter:: b3_feh =-12959.66163724488_rp 
      real(rp),parameter:: b4_feh = 6817.850021676971_rp 
      real(rp),parameter:: b5_feh =-1422.1723964897117_rp
!      real(rp),parameter:: b0_feh = 853.4769964964161d0
!      real(rp),parameter:: b1_feh =-4206.406420131467d0
!      real(rp),parameter:: b2_feh = 8686.215689507188d0
!      real(rp),parameter:: b3_feh =-9137.341019760202d0
!      real(rp),parameter:: b4_feh = 4807.823405345844d0
!      real(rp),parameter:: b5_feh =-1002.904058496011d0
!.....a_phi_feh
      real(rp),parameter:: a_phi_feh(1:7)= &
     &     (/ 14.0786236789212005_rp & !1
     &     ,  -4.4526835887173704_rp & 
     &     ,   5.5025121262565992_rp & 
     &     ,  -1.0687489808214079_rp & 
     &     ,  -0.3461498208163201_rp & 
     &     ,  -0.0064991947759021_rp & 
     &     ,  -0.0357435602984102_rp & !7
     &     /)
!.....r_phi_feh
      real(rp),parameter:: r_phi_feh(1:7)= &
     &     (/ 1.6_rp & !1
     &     ,  1.7_rp & 
     &     ,  1.8_rp & 
     &     ,  2.0_rp & 
     &     ,  2.5_rp & 
     &     ,  3.2_rp & 
     &     ,  4.2_rp & !7
     &     /)
!.....Cut off for phi_HH
      real(rp),parameter:: rc_phi_hh  = 2.3_rp
      real(rp),parameter:: r0_hh      = 0.74_rp
!.....S(r) for phi_HH
      real(rp),parameter:: r_tanh_hh  = 0.9_rp
      real(rp),parameter:: a_tanh_hh  = 25.0_rp
!.....Coeff for phi_HH
      real(rp),parameter:: c1_phi_hh  = 0.0_rp
      real(rp),parameter:: c2_phi_hh  = 0.0_rp
      real(rp),parameter:: eb_hh      = 2.37_rp
      real(rp),parameter:: almbd_hh   = 0.4899_rp
!.....a^F parameters
      real(rp),parameter:: a_f(1:6)= &
     &     (/ -0.0581256120818134_rp &
     &     ,   0.0022854552833736_rp &
     &     ,  -0.0000314202805805_rp &
     &     ,   0.0000013764132084_rp &
     &     ,  -0.0000000253707731_rp &
     &     ,   0.0000000001483685_rp &
     &     /)
!.....Rho_FeH
      real(rp),parameter:: a_rho_feh(1:6)= &
     &     (/ 10.0073629216300581_rp & 
     &     ,  32.4861983261490295_rp & 
     &     ,  -0.9494226032063788_rp & 
     &     ,  11.6659812262450338_rp & 
     &     ,  -0.0147080251458273_rp & 
     &     ,   0.4943383753319843_rp & 
     &     /)
      real(rp),parameter:: r_rho_feh(1:6)= &
     &     (/ 1.6_rp &
     &     ,  1.8_rp &
     &     ,  2.0_rp &
     &     ,  2.4_rp &
     &     ,  3.2_rp &
     &     ,  4.2_rp &
     &     /)
!.....Rho_HFe
      real(rp),parameter:: a_rho_hfe(1:5)= &
     &     (/ 11.1667357634216433_rp &
     &     ,  -3.0351307365078730_rp &
     &     ,   3.6096144794370653_rp &
     &     ,   0.0212509034775648_rp &
     &     ,   0.030391493994625_rp  &
     &     /)
      real(rp),parameter:: r_rho_hfe(1:5)= &
     &     (/ 1.5_rp &
     &     ,  2.0_rp &
     &     ,  2.5_rp &
     &     ,  3.0_rp &
     &     ,  4.2_rp &
     &     /)
!.....Rho_HH
      real(rp),parameter:: c_rho_hh = 1800.0_rp

!.....correction for H-H potential
!.....see Song and Curtin, nature materials 3479
!.....potential A
      real(rp),parameter:: k_hh_corr  = 1.5_rp
      real(rp),parameter:: r1_hh_corr  = 3.0_rp
      real(rp),parameter:: r0_hh_corr  = 0.9_rp
      real(rp),parameter:: lmbd_hh_corr= 1.0_rp
      real(rp),parameter:: b0_hh_corr  = 1.44_rp
      real(rp),parameter:: c0_hh_corr  = 0.19_rp
!.....potential B
!!$      real(rp),parameter:: k_hh_corr  = 1.5d0
!!$      real(rp),parameter:: r1_hh_corr  = 2.5d0
!!$      real(rp),parameter:: r0_hh_corr  = 0.9d0
!!$      real(rp),parameter:: lmbd_hh_corr= 1d0
!!$      real(rp),parameter:: b0_hh_corr  = 3.24d0
!!$      real(rp),parameter:: c0_hh_corr  = 0.239d0
