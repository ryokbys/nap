c.....Mass
      real(8),parameter:: am_W   = 183.84d0 *aump/aume
      real(8),parameter:: am_He  =   4.0026d0 *aump/aume
c.....W-W, 2-body
      real(8),parameter:: p_WW_c  = 3.25d0 *aa2bohr
      real(8),parameter:: p_WW_c0 = 47.062906d0 *ev2hrt
      real(8),parameter:: p_WW_c1 = -33.967128d0 *ev2hrt /aa2bohr
      real(8),parameter:: p_WW_c2 = 6.2542d0 *ev2hrt /aa2bohr /aa2bohr
      real(8),parameter:: p_WW_B  = 90.300452d0 *ev2hrt
      real(8),parameter:: p_WW_alpha= 1.093995d0 /aa2bohr
      real(8),parameter:: p_WW_b0   = 2.982916d0 *aa2bohr
c.....W-W, N-body
      real(8),parameter:: p_WW_A  = 1.773037d0 *ev2hrt /aa2bohr
      real(8),parameter:: p_WW_d  = 4.400224d0 *aa2bohr
      real(8),parameter:: p_WW_beta = -0.083332d0

c.....W-He, 2-body
      real(8),parameter:: p_WHe_rc = 4.0d0 *aa2bohr
      real(8),parameter:: p_WHe_z0 = 7.40929d0 *ev2hrt
      real(8),parameter:: p_WHe_z1 = 106.783d0 *ev2hrt *aa2bohr
      real(8),parameter:: p_WHe_z2 = 314.141d0 *ev2hrt *aa2bohr*aa2bohr
      real(8),parameter:: p_WHe_z3 = 16.491d0 *ev2hrt /aa2bohr
      real(8),parameter:: p_WHe_e0 = 12.5676d0 /aa2bohr
      real(8),parameter:: p_WHe_e1 = 3.75022d0 /aa2bohr
      real(8),parameter:: p_WHe_e2 = 2.78994d0 /aa2bohr
      real(8),parameter:: p_WHe_e3 = 5.43841d0 /aa2bohr
c.....W-He, N-body
      real(8),parameter:: p_WHe_d = 3.5d0 *aa2bohr
      real(8),parameter:: p_WHe_A = 0.0d0 *ev2hrt

c.....He-He, 2-body
      real(8),parameter:: p_HeHe_rc = 4.0d0 *aa2bohr
      real(8),parameter:: p_HeHe_z0 = 155.491d0 *ev2hrt
      real(8),parameter:: p_HeHe_z1 = 21.6251d0 *ev2hrt *aa2bohr
      real(8),parameter:: p_HeHe_z2 = 32.4872d0 *ev2hrt *aa2bohr*aa2bohr
      real(8),parameter:: p_HeHe_z3 = -90.4287d0 *ev2hrt /aa2bohr
      real(8),parameter:: p_HeHe_e0 = 3.17058d0 /aa2bohr
      real(8),parameter:: p_HeHe_e1 = 2.86124d0 /aa2bohr
      real(8),parameter:: p_HeHe_e2 = 12.9308d0 /aa2bohr
      real(8),parameter:: p_HeHe_e3 = 3.20055d0 /aa2bohr
