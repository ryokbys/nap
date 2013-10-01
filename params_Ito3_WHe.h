!.....Mass
      real(8),parameter:: am_W   = 183.84d0 *aump/aume
      real(8),parameter:: am_He  =   4.0026d0 *aump/aume

!.....energy scaling
      real(8),parameter:: sfac  = 1d0 !/2

!.....lattice constant of this potential
      real(8),parameter:: alcfe = 3.202d0 *aa2bohr

!.....2-body
!.....prefactor in eV*A, see the note on 2013-06-26
      real(8),parameter:: p_fac = 14.3997065885d0 *ev2hrt *aa2bohr
!.....Z
      real(8),parameter:: p_Z(1:2)= (/ 74d0, 2d0 /)
      real(8),parameter:: p_alpha(1:2,1:2)= &
     &     (/ 4.081955d0 /aa2bohr, &
     &        3.170542d0 /aa2bohr, &
     &        3.170542d0 /aa2bohr, &
     &        3.887407d0 /aa2bohr /)
      real(8),parameter:: p_beta(1:2,1:2)= &
     &     (/  0d0 /aa2bohr**2, &
     &        -0.984959d0 /aa2bohr**2, &
     &        -0.984959d0 /aa2bohr**2, &
     &         5.302379d0 /aa2bohr**2 /)
      real(8),parameter:: p_gamma(1:2,1:2)= &
     &     (/  0.0d0 /aa2bohr**3, &
     &         0.436908d0 /aa2bohr**3, &
     &         0.436908d0 /aa2bohr**3, &
     &        -2.619318d0 /aa2bohr**3 /)


!.....Many-body
      real(8),parameter:: p_rs(1:2,1:2)= &
     &     (/ 3.235582d0 *aa2bohr, &
     &        3.5d0 *aa2bohr, &
     &        3.5d0 *aa2bohr, &
     &        3.8d0 *aa2bohr /)
      real(8),parameter:: p_rl(1:2,1:2)= &
     &     (/ 4.274613d0 *aa2bohr, &
     &        4.0d0 *aa2bohr, &
     &        4.0d0 *aa2bohr, &
     &        4.3d0 *aa2bohr /) ! cutoff radius of this potential
      real(8),parameter:: p_rsp  = 3.151734 *aa2bohr
      real(8),parameter:: p_rlp  = 4.069256 *aa2bohr
      real(8),parameter:: p_B   = 192.0050d0 *ev2hrt**2/aa2bohr
      real(8),parameter:: p_c   = 1.418559d0 /aa2bohr
      real(8),parameter:: p_d   = 24.05130d0 *ev2hrt**2
!.....cutoff radius of this potential
      real(8),parameter:: rc_pot = p_rlp
