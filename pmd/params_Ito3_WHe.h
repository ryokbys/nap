!.....Mass
      real(rp),parameter:: am_W   = 183.84_rp
      real(rp),parameter:: am_He  =   4.0026_rp

!.....energy scaling
      real(rp),parameter:: sfac  = 1.0_rp !/2

!.....lattice constant of this potential
      real(rp),parameter:: alcfe = 3.202_rp 

!.....2-body
!.....prefactor in eV*A, see the note on 2013-06-26
      real(rp),parameter:: p_fac = 14.3997065885_rp
!.....Z
      real(rp),parameter:: p_Z(1:2)= (/ 74.0_rp, 2.0_rp /)
      real(rp),parameter:: p_alpha(1:2,1:2)= &
     &     reshape([ 4.081955_rp, &
     &        3.170542_rp, &
     &        3.170542_rp, &
     &        3.887407_rp ],[2,2])
      real(rp),parameter:: p_beta(1:2,1:2)= &
     &     reshape([  0.0_rp, &
     &        -0.984959_rp, &
     &        -0.984959_rp, &
     &         5.302379_rp ],[2,2])
      real(rp),parameter:: p_gamma(1:2,1:2)= &
     &     reshape([  0.0_rp, &
     &         0.436908_rp, &
     &         0.436908_rp, &
     &        -2.619318_rp ],[2,2])


!.....Many-body
      real(rp),parameter:: p_rs(1:2,1:2)= &
     &     reshape([ 3.235582_rp, &
     &        3.5_rp, &
     &        3.5_rp, &
     &        3.8_rp ],[2,2])
      real(rp),parameter:: p_rl(1:2,1:2)= &
     &     reshape([ 4.274613_rp, &
     &        4.0_rp, &
     &        4.0_rp, &
     &        4.3_rp ],[2,2]) ! cutoff radius of this potential
      real(rp),parameter:: p_rsp  = 3.151734
      real(rp),parameter:: p_rlp  = 4.069256
      real(rp),parameter:: p_B   = 192.0050_rp
      real(rp),parameter:: p_c   = 1.418559_rp
      real(rp),parameter:: p_d   = 24.05130_rp
!.....cutoff radius of this potential
      real(rp),parameter:: rc_pot = p_rlp
