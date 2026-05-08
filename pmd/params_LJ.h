!-----LJ parameters for Argon
      real(rp),parameter:: epslj = real(120d0, rp) *fkb
      real(rp),parameter:: sgmlj = real(3.41d-10, rp) /ang
      real(rp),parameter:: am_ar = real(39.948d0, rp)
      real(rp),parameter:: alcar = real(2d0, rp)**(real(1d0, rp)/6)*sgmlj &
          *real(1.41421356d0, rp)*real(0.996d0, rp)
