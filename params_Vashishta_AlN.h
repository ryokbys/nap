c-----Al,N mass (to be multiplied by amu)
      real(8),parameter:: am_Al = 26.9815d0 *amu
      real(8),parameter:: am_N  = 14.0067d0 *amu
c-----J*m/C**2 ---> Ht*Bohr/e**2
      real(8),parameter:: cfct= 8.9876d+9 *j2ehrt /bohr
     &     /(6.2414d+18)**2
c      real(8),parameter:: cfct= 1d0
c-----Vashishta potential parameters for AlN
      real(8),parameter:: v_z(1:2)=
     &     (/ 1.0366d0,
     &       -1.0366d0 /)
      real(8),parameter:: v_alp(1:2)= (/ 0.0d0 *aa2bohr**3,
     &                                   3.0d0 *aa2bohr**3 /)
      real(8),parameter:: v_n(1:2,1:2)= (/ 7d0, 9d0, 9d0, 7d0 /)
      real(8),parameter:: v_h(1:2,1:2)=
     &     (/ 8.7162d-17 *j2ehrt *aa2bohr**v_n(1,1),
     &        2.6824d-18 *j2ehrt *aa2bohr**v_n(2,1),
     &        2.6824d-18 *j2ehrt *aa2bohr**v_n(1,2),
c     &        2.6824d-17 *j2ehrt *aa2bohr**v_n(2,1),
c     &        2.6824d-17 *j2ehrt *aa2bohr**v_n(1,2),
     &        2.1790d-17 *j2ehrt *aa2bohr**v_n(2,2) /)
      real(8),parameter:: v_w(1:2,1:2)=
     &     (/ 0d0,
     &        9.7901d-18 *j2ehrt *aa2bohr**6,
     &        9.7901d-18 *j2ehrt *aa2bohr**6,
     &        0d0 /)
      real(8),parameter:: v_r1s(1:2,1:2)=
     &     (/ 5d0 *aa2bohr,
     &        5d0 *aa2bohr,
     &        5d0 *aa2bohr,
     &        5d0 *aa2bohr /)
      real(8),parameter:: v_r4s(1:2,1:2)=
     &     (/ 3.75d0 *aa2bohr,
     &        3.75d0 *aa2bohr,
     &        3.75d0 *aa2bohr,
     &        3.75d0 *aa2bohr /)
      real(8),parameter:: v_b  = 40d-19 *j2ehrt
      real(8),parameter:: v_c  = 30d0
      real(8),parameter:: v_xi = 1d0 *aa2bohr
      real(8),parameter:: v_r0 = 2.6d0 *aa2bohr
      real(8),parameter:: v_tht= 109.47122d0 /180d0 *3.14159265358979d0
c-----num of division for 2-body term table
      integer,parameter:: nd2b = 2048
