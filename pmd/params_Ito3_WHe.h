!.....Mass
      real(rp),parameter:: am_W   = real(183.84d0, rp)
      real(rp),parameter:: am_He  =   real(4.0026d0, rp)

!.....energy scaling
      real(rp),parameter:: sfac  = real(1d0, rp) !/2

!.....lattice constant of this potential
      real(rp),parameter:: alcfe = real(3.202d0, rp) 

!.....2-body
!.....prefactor in eV*A, see the note on 2013-06-26
      real(rp),parameter:: p_fac = real(14.3997065885d0, rp)
!.....Z
      real(rp),parameter:: p_Z(1:2)= (/ real(74d0, rp), real(2d0, rp) /)
      real(rp),parameter:: p_alpha(1:2,1:2)= &
     &     reshape([ real(4.081955d0, rp), &
     &        real(3.170542d0, rp), &
     &        real(3.170542d0, rp), &
     &        real(3.887407d0, rp) ],[2,2])
      real(rp),parameter:: p_beta(1:2,1:2)= &
     &     reshape([  real(0d0, rp), &
     &        -real(0.984959d0, rp), &
     &        -real(0.984959d0, rp), &
     &         real(5.302379d0, rp) ],[2,2])
      real(rp),parameter:: p_gamma(1:2,1:2)= &
     &     reshape([  real(0.0d0, rp), &
     &         real(0.436908d0, rp), &
     &         real(0.436908d0, rp), &
     &        -real(2.619318d0, rp) ],[2,2])


!.....Many-body
      real(rp),parameter:: p_rs(1:2,1:2)= &
     &     reshape([ real(3.235582d0, rp), &
     &        real(3.5d0, rp), &
     &        real(3.5d0, rp), &
     &        real(3.8d0, rp) ],[2,2])
      real(rp),parameter:: p_rl(1:2,1:2)= &
     &     reshape([ real(4.274613d0, rp), &
     &        real(4.0d0, rp), &
     &        real(4.0d0, rp), &
     &        real(4.3d0, rp) ],[2,2]) ! cutoff radius of this potential
      real(rp),parameter:: p_rsp  = 3.151734
      real(rp),parameter:: p_rlp  = 4.069256
      real(rp),parameter:: p_B   = real(192.0050d0, rp)
      real(rp),parameter:: p_c   = real(1.418559d0, rp)
      real(rp),parameter:: p_d   = real(24.05130d0, rp)
!.....cutoff radius of this potential
      real(rp),parameter:: rc_pot = p_rlp
