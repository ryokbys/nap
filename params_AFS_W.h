c.....Mass
      real(8),parameter:: am_W   = 183.84d0 *aump/aume

c.....energy scaling
      real(8),parameter:: sfac  = 1d0 !/2

c.....length scaling for hybrid QMCL calculation
      real(8),parameter:: slen  = 1d0 
c/3.204d0 *3.172d0
      real(8),parameter:: aa2bs = aa2bohr*slen

c.....2-body
c.....W-W
      real(8),parameter:: p_WW_c    =  3.25d0 *aa2bs
      real(8),parameter:: p_WW_c0   = 47.1346499d0 *ev2hrt
      real(8),parameter:: p_WW_c1   =-33.7665655d0 *ev2hrt /aa2bs
      real(8),parameter:: p_WW_c2   =  6.2541999d0 *ev2hrt /aa2bs**2
      real(8),parameter:: p_WW_B    = 90.3d0 *ev2hrt /aa2bs**3
      real(8),parameter:: p_WW_alpha=  1.2d0 /aa2bs
      real(8),parameter:: p_WW_b0   =  2.7411d0 *aa2bs

c.....Many-body
c.....W
      real(8),parameter:: p_W_A    = 1.896373d0 *ev2hrt
      real(8),parameter:: p_W_d    = 4.400224d0 *aa2bs
      real(8),parameter:: p_W_beta = 0.0d0
