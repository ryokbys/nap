c.....Mass
      real(8),parameter:: am_W   = 183.84d0 *aump/aume
      real(8),parameter:: am_He  =   4.0026d0 *aump/aume

c.....energy scaling
      real(8),parameter:: sfac  = 1d0 !/2

c.....2-body
c.....W-W
      real(8),parameter:: p_WW_c    = 3.25d0 *aa2bohr
      real(8),parameter:: p_WW_c0   = 47.042512d0 *ev2hrt /aa2bohr**2
      real(8),parameter:: p_WW_c1   = -33.987390d0 *ev2hrt /aa2bohr**3
      real(8),parameter:: p_WW_c2   = 6.319348d0 *ev2hrt /aa2bohr**4
      real(8),parameter:: p_WW_B    = 90.300865d0 *ev2hrt /aa2bohr**3
      real(8),parameter:: p_WW_alpha= 0.997985d0 /aa2bohr
      real(8),parameter:: p_WW_b0   = 2.951811d0 *aa2bohr
c.....W-He
      real(8),parameter:: p_WHe_rc = 4.0d0 *aa2bohr
      real(8),parameter:: p_WHe_rs = 2.2d0 *aa2bohr
      real(8),parameter:: p_WHe_cs = p_WHe_rc -p_WHe_rs
      real(8),parameter:: p_WHe_z1 = 1488.31d0 *ev2hrt *aa2bohr
      real(8),parameter:: p_WHe_a1 = 4.129d0 /aa2bohr
      real(8),parameter:: p_WHe_b0 = -0.0033d0
c.....He-He
      real(8),parameter:: p_HeHe_rc = 2.4d0 *aa2bohr
      real(8),parameter:: p_HeHe_rs = 1.1d0 *aa2bohr
      real(8),parameter:: p_HeHe_cs = p_HeHe_rc -p_HeHe_rs
      real(8),parameter:: p_HeHe_z1 = 87.8831d0 *ev2hrt *aa2bohr
      real(8),parameter:: p_HeHe_a1 = 3.01517d0 /aa2bohr
      real(8),parameter:: p_HeHe_b0 = 0.0d0

c.....Many-body
c.....W
      real(8),parameter:: p_W_w    = 1.786367d0 *ev2hrt /aa2bohr
      real(8),parameter:: p_W_d    = 4.40d0 *aa2bohr
      real(8),parameter:: p_W_beta = -0.23433d0
c.....He
      real(8),parameter:: p_He_A = -0.154896d0 /ev2hrt
      real(8),parameter:: p_He_w = 0.812281d0 *ev2hrt /aa2bohr
      real(8),parameter:: p_He_d = p_W_d
      real(8),parameter:: p_He_beta = p_W_beta

