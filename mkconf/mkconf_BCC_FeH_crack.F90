program mkconf_BCC_FeH_crack
!-----------------------------------------------------------------------
!  Create a BCC crack model with Fe and H atoms.
!  Crack runs towards x-direction and tensile stress in on z-direction,
!  and periodic in y-direction.
!-----------------------------------------------------------------------
!  OUTPUT:
!    - pmd00000
!    - akr0000
!-----------------------------------------------------------------------
  implicit real*8(a-h,o-z),integer(i-n)
  include '../pmd/params_unit.h'
  include '../pmd/params_RK_FeH.h'
!      include '../pmd/params_Ramas_FeH.h'
!-----max # of atoms
  integer,parameter::nmax=100000
!-----# of unit cells
!      integer,parameter:: nuc(1:3)= (/ 1, 1, 1 /)
  integer,parameter:: nuc(1:3)= (/ 300,1,100 /)
!-----vacuum width in unit of cell
!      integer,parameter:: nvac(1:3)= (/ 5, 5, 5 /)
  integer,parameter:: nvac(1:3)= (/ 4, 10, 10 /)
  real(8):: ua(3,10)
  real(8):: tag(nmax),ra(3,nmax),va(3,nmax),eki(nmax),epi(nmax) &
       ,h(3,3,0:1),strs(3,3,nmax)

  small=1d-7

!.....Lattice constant of Fe, see Philos. Mag. 83 (2003) 3977
!      cunit= 2.835d0
  cunit= 2.8553d0

!-----simulation box size
  h(1:3,1:3,0:1)= 0d0
  h(1,1,0)= cunit*(nuc(1)+nvac(1))
  h(2,2,0)= cunit*(nuc(2)+nvac(2))
  h(3,3,0)= cunit*(nuc(3)+nvac(3))

!-----unit cell, BCC
  ua(1:3,1)= (/ 0.0d0, 0.0d0, 0.0d0 /)
  ua(1:3,2)= (/ 0.5d0, 0.5d0, 0.5d0 /)

!.....oval parameters
  a= 1d0/9
  b= a/10
  inc=0
  x0= (0.5d0*nvac(1)) /(nuc(1)+nvac(1))
  y0= (0.5d0*nvac(2)) /(nuc(2)+nvac(2))
  z0= (0.5d0*nvac(3)) /(nuc(3)+nvac(3))
  do ix=0,nuc(1)-1
    do iy=0,nuc(2)-1
      do iz=0,nuc(3)-1
!!$!.....remove center layer of z of 1/3 of x as a crack tip
!!$        if( iz.eq. nuc(3)/2 .and. ix.lt.nuc(1)/3 ) cycle
        do m=1,2
          x=(ua(1,m)+dble(ix))/(nuc(1)+nvac(1)) +small +x0
          y=(ua(2,m)+dble(iy))/(nuc(2)+nvac(2)) +small +y0
          z=(ua(3,m)+dble(iz))/(nuc(3)+nvac(3)) +small +z0
!.....remove an eval region at 1/3 of x as a crack tip
          if( x.lt.a .and. abs(z-0.5d0).lt.b*sqrt(1d0-(x/a)**2) ) cycle
          inc=inc+1
          if(inc.gt.nmax)then
            write(*,*)'Error inc>nmax',inc,nmax
            stop
          endif
          ra(1,inc)= x
          ra(2,inc)= y
          ra(3,inc)= z
!              ra(1:3,inc)= ra(1:3,inc)
!     &             +dble(nvac(1:3))/(nuc(1:3)+nvac(1:3))/2
          is= 1
          ifmv= 1
          if( ix.eq.0 .or. ix.eq.nuc(1)-1 ) ifmv= 3
          if( iz.eq.0 .or. iz.eq.nuc(3)-1 ) ifmv= 2
          tag(inc)= 1d0*is +0.1d0*ifmv +1d-14*inc
        enddo
      enddo
    enddo
  enddo

  write(6,'(a,i10)') " num of atoms=",inc

  va(1:3,1:inc)= 0d0

  call write_pmd0_ascii(15,'pmd00000','replace',inc,tag &
       ,ra,va,h,cunit,eki,epi,strs)
!!$      call write_pmd0_bin(15,'pmd00000','replace',inc,tag,ra,va,h &
!!$           ,eki,epi,strs)

!-----output 'akr000' for Akira visualization
  call write_akr(15,'akr0000',inc,h,cunit,tag,ra,va)


end program mkconf_BCC_FeH_crack
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make mkconf_BCC_FeH_crack"
!     End:
