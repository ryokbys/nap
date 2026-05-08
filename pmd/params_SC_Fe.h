!-----EAM Sutton-Chen model for Fe (iron)
      real(rp),parameter:: am_fe = real(55.847d0, rp)
!      real(rp),parameter:: alcfe = 2.8665d-10 /ang
      real(rp),parameter:: alcfe = real(2.9782935d-10, rp) /ang
!-----Sutton-Chen model parameters for Fe (in eV,Ang)
      real(rp),parameter:: sc_eps= real(0.017306d0, rp)
      real(rp),parameter:: sc_a  = real(3.471392d-10, rp) /ang
      real(rp),parameter:: sc_n  = real(8.137381d0, rp)
      real(rp),parameter:: sc_m  = real(4.7877d0, rp)
      real(rp),parameter:: sc_c  = real(24.9390d0, rp)
!      real(rp),parameter:: sc_rc = 6d-10 /ang
