!-----mass
      real(rp),parameter:: am_W  = 183.84
      real(rp),parameter:: am_He =   4.0026
!
!  cut off should be larger than R+D (=4.13158A+0.93018A =5.06176A)
!
!.....length scaling for hybrid QMCL calculation
      real(rp),parameter:: slen= 1.0_rp /3.165_rp *3.172_rp
!.....W-W and W-He parameters
      real(rp),parameter:: p_D0(1:2,1:2)= reshape([ &
           2.87454_rp, &
           0.377472_rp, &
           0.377472_rp, &
           0.0_rp],[2,2])
      real(rp),parameter:: p_r0(1:2,1:2)=reshape([  &
           2.38631_rp *slen, &
           3.40623_rp *slen, &
           3.40623_rp *slen, &
           0.0_rp *slen],[2,2])
      real(rp),parameter:: p_beta(1:2,1:2)=reshape([  &
           1.33682_rp  /slen, &
           0.602677_rp /slen, &
           0.602677_rp /slen, &
           0.0_rp /slen],[2,2])
      real(rp),parameter:: p_S(1:2,1:2)=reshape([ &
           1.25348_rp, &
           4.399296_rp, &
           4.399296_rp, &
           0.0_rp],[2,2])
      real(rp),parameter:: p_gamma(1:2,1:2)=reshape([ &
           8.3879e-4_rp, &
           2.623239e-2_rp, &
           2.623239e-2_rp, &
           0.0_rp],[2,2])
      real(rp),parameter:: p_c(1:2,1:2)=reshape([ &
           0.850284_rp, &
           0.371185_rp, &
           0.371185_rp, &
           0.0_rp],[2,2])
      real(rp),parameter:: p_d(1:2,1:2)=reshape([ &
           0.144317_rp, &
           0.252377_rp, &
           0.252377_rp, &
           0.0_rp],[2,2])
      real(rp),parameter:: p_h(1:2,1:2)=reshape([ &
           -0.36846_rp, &
           0.741345_rp, &
           0.741345_rp, &
           0.0_rp],[2,2])
      real(rp),parameter:: p_R1(1:2,1:2)=reshape([ &
           4.131580_rp *slen, &
           3.00_rp *slen, &
           3.00_rp *slen, &
           0.0_rp *slen],[2,2])
      real(rp),parameter:: p_D1(1:2,1:2)=reshape([ &
           0.930180_rp *slen, &
           0.35_rp *slen, &
           0.35_rp *slen, &
           0.0_rp *slen],[2,2])
!.....He-He parameters
      real(rp),parameter:: p_HeHe_rc    = p_R1(1,1) +p_D1(1,1)
      real(rp),parameter:: p_HeHe_rm    = 2.970_rp *slen
      real(rp),parameter:: p_HeHe_A     = 1.9211529e+5_rp
      real(rp),parameter:: p_HeHe_alpha = 10.73520708_rp
      real(rp),parameter:: p_HeHe_c6    = 1.34920045_rp
      real(rp),parameter:: p_HeHe_c8    = 0.41365922_rp
      real(rp),parameter:: p_HeHe_c10   = 0.17078164_rp
      real(rp),parameter:: p_HeHe_beta  = -1.89296514_rp
      real(rp),parameter:: p_HeHe_D     = 1.4132_rp
      real(rp),parameter:: p_HeHe_eps   = 10.94_rp *fkb

      
