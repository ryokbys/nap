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
      real(rp),parameter:: am_fe = real(55.847d0, rp)
      real(rp),parameter:: am_h  =  real(1.008d0, rp)
!      real(rp),parameter:: alcfe = 2.8553d0 ! Ackland original
      real(rp),parameter:: alcfe = real(2.835d0, rp)
!.....length scaling factor for hybrid QMCL calc
      real(rp),parameter:: slen  = real(1d0, rp) /real(2.8553d0, rp) *real(2.835d0, rp)
      real(rp),parameter:: rfac  = real(2.835d0, rp) /real(2.8553d0, rp)
      real(rp),parameter:: rfaci = real(1d0, rp) /rfac
      real(rp),parameter:: z_fe  = real(26d0, rp)
      real(rp),parameter:: z_h   = real(1d0, rp)
!.....prefactor for rs
      real(rp),parameter:: a_rs  = real(0.88534d0, rp)
!.....prefactor for varphi r < r1
      real(rp),parameter:: qe    = real(3.794701096299472d0, rp) *slen
      real(rp),parameter:: z2q2  = 9.7342365892908E+03 *slen
      real(rp),parameter:: z2q2_feh= z_fe*z_h*qe*qe *slen
!.....Range in varphi
      real(rp),parameter:: r1  = real(1.0d0, rp) *slen
      real(rp),parameter:: r2  = real(2.0d0, rp) *slen
!.....Bs for r1 < r <= r2
      real(rp),parameter:: b0  = real(6.4265260576348d0, rp)
      real(rp),parameter:: b1  = real(1.7900488524286d0, rp) /slen
      real(rp),parameter:: b2  =-real(4.5108316729807d0, rp) /slen**2
      real(rp),parameter:: b3  = real(1.0866199373306d0, rp) /slen**3
!.....varphi for r2 < r
      real(rp),parameter:: a_vphi(1:15)= &
     &     (/  real(0.0d0, rp) &
     &     , -real(24.028204854115d0, rp)   /slen**3 &! 2
     &     ,  real(11.300691696477d0, rp)   /slen**3 &! 3
     &     ,   real(5.3144495820462d0, rp)  /slen**3 &! 4
     &     ,  -real(4.6659532856049d0, rp)  /slen**3 &! 5
     &     ,   real(5.9637758529194d0, rp)  /slen**3 &
     &     ,  -real(1.7710262006061d0, rp)  /slen**3 &
     &     ,   real(0.85913830768731d0, rp) /slen**3 &
     &     ,  -real(2.1845362968261d0, rp)  /slen**3 &
     &     ,   real(2.6424377007466d0, rp)  /slen**3 &! 10
     &     ,  -real(1.0358345370208d0, rp)  /slen**3 &
     &     ,   real(0.33548264951582d0, rp) /slen**3 &
     &     ,  -real(4.6448582149334d-2, rp) /slen**3 &
     &     ,  -real(7.0294963048689d-3, rp) /slen**3 &! 14
     &     ,   real(0.0d0, rp) /)
      real(rp),parameter:: r_vphi(1:15)= &
     &     (/  real(0.0d0, rp) &
     &     ,   real(2.2d0, rp) *slen & ! 2
     &     ,   real(2.3d0, rp) *slen & 
     &     ,   real(2.4d0, rp) *slen & 
     &     ,   real(2.5d0, rp) *slen & ! 5
     &     ,   real(2.6d0, rp) *slen & 
     &     ,   real(2.7d0, rp) *slen & 
     &     ,   real(2.8d0, rp) *slen & 
     &     ,   real(3.0d0, rp) *slen & 
     &     ,   real(3.3d0, rp) *slen & ! 10
     &     ,   real(3.7d0, rp) *slen & 
     &     ,   real(4.2d0, rp) *slen & 
     &     ,   real(4.7d0, rp) *slen & 
     &     ,   real(5.3d0, rp) *slen & ! 14
     &     ,   real(0.0d0, rp) /)
!.....psi for calculation of rho
      real(rp),parameter:: a_psi(1:3)= &
     &     (/ real(11.686859407970d0, rp)     /slen**3 &
     &     ,  -real(0.014710740098830d0, rp)  /slen**3 &
     &     ,   real(0.47193527075943d0, rp)   /slen**3  /)
      real(rp),parameter:: r_psi(1:3)= &
     &     (/  real(2.4d0, rp) *slen &
     &     ,   real(3.2d0, rp) *slen &
     &     ,   real(4.2d0, rp) *slen /)
!.....coeff for the embedded func
      real(rp),parameter:: a_emb  = -real(3.5387096579929d-4, rp)
!.....cutoffs
      real(rp),parameter:: rc_rho = real(4.2d0, rp) *slen
      real(rp),parameter:: rc_vphi= real(5.3d0, rp) *slen

!-----------------------------------------------------------------------
!  For Fe-H potential parameters
!-----------------------------------------------------------------------
!.....Range in phi_FeH
      real(rp),parameter:: r1_feh =  real(0.6d0, rp)
      real(rp),parameter:: r2_feh =  real(1.2d0, rp)
!.....Bs for FeH
      real(rp),parameter:: b0_feh = real(1242.1709168218642d0, rp)
      real(rp),parameter:: b1_feh =-real(6013.566711223783d0, rp) 
      real(rp),parameter:: b2_feh = real(12339.540893927151d0, rp)
      real(rp),parameter:: b3_feh =-real(12959.66163724488d0, rp) 
      real(rp),parameter:: b4_feh = real(6817.850021676971d0, rp) 
      real(rp),parameter:: b5_feh =-real(1422.1723964897117d0, rp)
!      real(rp),parameter:: b0_feh = 853.4769964964161d0
!      real(rp),parameter:: b1_feh =-4206.406420131467d0
!      real(rp),parameter:: b2_feh = 8686.215689507188d0
!      real(rp),parameter:: b3_feh =-9137.341019760202d0
!      real(rp),parameter:: b4_feh = 4807.823405345844d0
!      real(rp),parameter:: b5_feh =-1002.904058496011d0
!.....a_phi_feh
      real(rp),parameter:: a_phi_feh(1:7)= &
     &     (/ real(14.0786236789212005d0, rp) & !1
     &     ,  -real(4.4526835887173704d0, rp) & 
     &     ,   real(5.5025121262565992d0, rp) & 
     &     ,  -real(1.0687489808214079d0, rp) & 
     &     ,  -real(0.3461498208163201d0, rp) & 
     &     ,  -real(0.0064991947759021d0, rp) & 
     &     ,  -real(0.0357435602984102d0, rp) & !7
     &     /)
!.....r_phi_feh
      real(rp),parameter:: r_phi_feh(1:7)= &
     &     (/ real(1.6d0, rp) & !1
     &     ,  real(1.7d0, rp) & 
     &     ,  real(1.8d0, rp) & 
     &     ,  real(2.0d0, rp) & 
     &     ,  real(2.5d0, rp) & 
     &     ,  real(3.2d0, rp) & 
     &     ,  real(4.2d0, rp) & !7
     &     /)
!.....Cut off for phi_HH
      real(rp),parameter:: rc_phi_hh  = real(2.3d0, rp)
      real(rp),parameter:: r0_hh      = real(0.74d0, rp)
!.....S(r) for phi_HH
      real(rp),parameter:: r_tanh_hh  = real(0.9d0, rp)
      real(rp),parameter:: a_tanh_hh  = real(25d0, rp)
!.....Coeff for phi_HH
      real(rp),parameter:: c1_phi_hh  = real(0.0d0, rp)
      real(rp),parameter:: c2_phi_hh  = real(0.0d0, rp)
      real(rp),parameter:: eb_hh      = real(2.37d0, rp)
      real(rp),parameter:: almbd_hh   = real(0.4899d0, rp)
!.....a^F parameters
      real(rp),parameter:: a_f(1:6)= &
     &     (/ -real(0.0581256120818134d0, rp) &
     &     ,   real(0.0022854552833736d0, rp) &
     &     ,  -real(0.0000314202805805d0, rp) &
     &     ,   real(0.0000013764132084d0, rp) &
     &     ,  -real(0.0000000253707731d0, rp) &
     &     ,   real(0.0000000001483685d0, rp) &
     &     /)
!.....Rho_FeH
      real(rp),parameter:: a_rho_feh(1:6)= &
     &     (/ real(10.0073629216300581d0, rp) & 
     &     ,  real(32.4861983261490295d0, rp) & 
     &     ,  -real(0.9494226032063788d0, rp) & 
     &     ,  real(11.6659812262450338d0, rp) & 
     &     ,  -real(0.0147080251458273d0, rp) & 
     &     ,   real(0.4943383753319843d0, rp) & 
     &     /)
      real(rp),parameter:: r_rho_feh(1:6)= &
     &     (/ real(1.6d0, rp) &
     &     ,  real(1.8d0, rp) &
     &     ,  real(2.0d0, rp) &
     &     ,  real(2.4d0, rp) &
     &     ,  real(3.2d0, rp) &
     &     ,  real(4.2d0, rp) &
     &     /)
!.....Rho_HFe
      real(rp),parameter:: a_rho_hfe(1:5)= &
     &     (/ real(11.1667357634216433d0, rp) &
     &     ,  -real(3.0351307365078730d0, rp) &
     &     ,   real(3.6096144794370653d0, rp) &
     &     ,   real(0.0212509034775648d0, rp) &
     &     ,   real(0.030391493994625d0, rp)  &
     &     /)
      real(rp),parameter:: r_rho_hfe(1:5)= &
     &     (/ real(1.5d0, rp) &
     &     ,  real(2.0d0, rp) &
     &     ,  real(2.5d0, rp) &
     &     ,  real(3.0d0, rp) &
     &     ,  real(4.2d0, rp) &
     &     /)
!.....Rho_HH
      real(rp),parameter:: c_rho_hh = real(1800d0, rp)

!.....correction for H-H potential
!.....see Song and Curtin, nature materials 3479
!.....potential A
      real(rp),parameter:: k_hh_corr  = real(1.5d0, rp)
      real(rp),parameter:: r1_hh_corr  = real(3.0d0, rp)
      real(rp),parameter:: r0_hh_corr  = real(0.9d0, rp)
      real(rp),parameter:: lmbd_hh_corr= real(1d0, rp)
      real(rp),parameter:: b0_hh_corr  = real(1.44d0, rp)
      real(rp),parameter:: c0_hh_corr  = real(0.19d0, rp)
!.....potential B
!!$      real(rp),parameter:: k_hh_corr  = 1.5d0
!!$      real(rp),parameter:: r1_hh_corr  = 2.5d0
!!$      real(rp),parameter:: r0_hh_corr  = 0.9d0
!!$      real(rp),parameter:: lmbd_hh_corr= 1d0
!!$      real(rp),parameter:: b0_hh_corr  = 3.24d0
!!$      real(rp),parameter:: c0_hh_corr  = 0.239d0
