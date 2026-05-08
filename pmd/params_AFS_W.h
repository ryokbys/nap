!.....Mass
      real(rp),parameter:: am_W   = real(183.84d0, rp)

!.....energy scaling
      real(rp),parameter:: sfac  = real(1d0, rp) !/2

!.....length scaling for hybrid QMCL calculation
      real(rp),parameter:: slen  = real(1d0, rp) 

!.....2-body
!.....W-W
      real(rp),parameter:: p_WW_c    =  real(3.25d0, rp) *slen
      real(rp),parameter:: p_WW_c0   = real(47.1346499d0, rp) 
      real(rp),parameter:: p_WW_c1   =-real(33.7665655d0, rp)  /slen
      real(rp),parameter:: p_WW_c2   =  real(6.2541999d0, rp)  /slen**2
      real(rp),parameter:: p_WW_B    = real(90.3d0, rp)  /slen**3
      real(rp),parameter:: p_WW_alpha=  real(1.2d0, rp) /slen
      real(rp),parameter:: p_WW_b0   =  real(2.7411d0, rp) *slen

!.....Many-body
!.....W
      real(rp),parameter:: p_W_A    = real(1.896373d0, rp) 
      real(rp),parameter:: p_W_d    = real(4.400224d0, rp) *slen
      real(rp),parameter:: p_W_beta = real(0.0d0, rp)
