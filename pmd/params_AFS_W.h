!.....Mass
      real(rp),parameter:: am_W   = 183.84_rp

!.....energy scaling
      real(rp),parameter:: sfac  = 1.0_rp !/2

!.....length scaling for hybrid QMCL calculation
      real(rp),parameter:: slen  = 1.0_rp 

!.....2-body
!.....W-W
      real(rp),parameter:: p_WW_c    =  3.25_rp *slen
      real(rp),parameter:: p_WW_c0   = 47.1346499_rp 
      real(rp),parameter:: p_WW_c1   =-33.7665655_rp  /slen
      real(rp),parameter:: p_WW_c2   =  6.2541999_rp  /slen**2
      real(rp),parameter:: p_WW_B    = 90.3_rp  /slen**3
      real(rp),parameter:: p_WW_alpha=  1.2_rp /slen
      real(rp),parameter:: p_WW_b0   =  2.7411_rp *slen

!.....Many-body
!.....W
      real(rp),parameter:: p_W_A    = 1.896373_rp 
      real(rp),parameter:: p_W_d    = 4.400224_rp *slen
      real(rp),parameter:: p_W_beta = 0.0_rp
