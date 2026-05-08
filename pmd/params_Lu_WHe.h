!-----mass
      real(rp),parameter:: am_W  = 183.84
      real(rp),parameter:: am_He =   4.0026
!
!  cut off should be larger than R+D (=4.13158A+0.93018A =5.06176A)
!
!.....length scaling for hybrid QMCL calculation
      real(rp),parameter:: slen= real(1d0, rp) /real(3.165d0, rp) *real(3.172d0, rp)
!.....W-W and W-He parameters
      real(rp),parameter:: p_D0(1:2,1:2)= reshape([ &
           real(2.87454d0, rp), &
           real(0.377472d0, rp), &
           real(0.377472d0, rp), &
           real(0d0, rp)],[2,2])
      real(rp),parameter:: p_r0(1:2,1:2)=reshape([  &
           real(2.38631d0, rp) *slen, &
           real(3.40623d0, rp) *slen, &
           real(3.40623d0, rp) *slen, &
           real(0d0, rp) *slen],[2,2])
      real(rp),parameter:: p_beta(1:2,1:2)=reshape([  &
           real(1.33682d0, rp)  /slen, &
           real(0.602677d0, rp) /slen, &
           real(0.602677d0, rp) /slen, &
           real(0d0, rp) /slen],[2,2])
      real(rp),parameter:: p_S(1:2,1:2)=reshape([ &
           real(1.25348d0, rp), &
           real(4.399296d0, rp), &
           real(4.399296d0, rp), &
           real(0d0, rp)],[2,2])
      real(rp),parameter:: p_gamma(1:2,1:2)=reshape([ &
           real(8.3879d-4, rp), &
           real(2.623239d-2, rp), &
           real(2.623239d-2, rp), &
           real(0d0, rp)],[2,2])
      real(rp),parameter:: p_c(1:2,1:2)=reshape([ &
           real(0.850284d0, rp), &
           real(0.371185d0, rp), &
           real(0.371185d0, rp), &
           real(0d0, rp)],[2,2])
      real(rp),parameter:: p_d(1:2,1:2)=reshape([ &
           real(0.144317d0, rp), &
           real(0.252377d0, rp), &
           real(0.252377d0, rp), &
           real(0d0, rp)],[2,2])
      real(rp),parameter:: p_h(1:2,1:2)=reshape([ &
           -real(0.36846d0, rp), &
           real(0.741345d0, rp), &
           real(0.741345d0, rp), &
           real(0d0, rp)],[2,2])
      real(rp),parameter:: p_R1(1:2,1:2)=reshape([ &
           real(4.131580d0, rp) *slen, &
           real(3.00d0, rp) *slen, &
           real(3.00d0, rp) *slen, &
           real(0d0, rp) *slen],[2,2])
      real(rp),parameter:: p_D1(1:2,1:2)=reshape([ &
           real(0.930180d0, rp) *slen, &
           real(0.35d0, rp) *slen, &
           real(0.35d0, rp) *slen, &
           real(0d0, rp) *slen],[2,2])
!.....He-He parameters
      real(rp),parameter:: p_HeHe_rc    = p_R1(1,1) +p_D1(1,1)
      real(rp),parameter:: p_HeHe_rm    = real(2.970d0, rp) *slen
      real(rp),parameter:: p_HeHe_A     = real(1.9211529d+5, rp)
      real(rp),parameter:: p_HeHe_alpha = real(10.73520708d0, rp)
      real(rp),parameter:: p_HeHe_c6    = real(1.34920045d0, rp)
      real(rp),parameter:: p_HeHe_c8    = real(0.41365922d0, rp)
      real(rp),parameter:: p_HeHe_c10   = real(0.17078164d0, rp)
      real(rp),parameter:: p_HeHe_beta  = -real(1.89296514d0, rp)
      real(rp),parameter:: p_HeHe_D     = real(1.4132d0, rp)
      real(rp),parameter:: p_HeHe_eps   = real(10.94d0, rp) *fkb

      
