      program mkconf_BCC_FeH
!-----------------------------------------------------------------------
!  Create a BCC crystal with Fe and H atoms
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
     integer,parameter:: nuc(1:3)= (/ 5,5,5 /)
!-----vacuum width in unit of cell
!      integer,parameter:: nvac(1:3)= (/ 5, 5, 5 /)
      integer,parameter:: nvac(1:3)= (/ 0, 0, 0 /)
      real(8):: ua(3,10)
      real(8):: tag(nmax),ra(3,nmax),va(3,nmax),eki(nmax),epi(nmax) &
           ,h(3,3,0:1),strs(3,3,nmax)

      small=1d-7

!.....Lattice constant of Fe, see Philos. Mag. 83 (2003) 3977
      cunit= 2.835d0
!      cunit= 2.8553d0

!-----simulation box size
      h(1:3,1:3,0:1)= 0d0
      h(1,1,0)= cunit*(nuc(1)+nvac(1))
      h(2,2,0)= cunit*(nuc(2)+nvac(2))
      h(3,3,0)= cunit*(nuc(3)+nvac(3))

!-----unit cell, BCC
      ua(1:3,1)= (/ 0.0d0, 0.0d0, 0.0d0 /)
      ua(1:3,2)= (/ 0.5d0, 0.5d0, 0.5d0 /)
      
      inc=0 
      do ix=0,nuc(1)-1
        do iy=0,nuc(2)-1
          do iz=0,nuc(3)-1
            do m=1,2
              x=(ua(1,m)+dble(ix))/(nuc(1)+nvac(1)) +small
              y=(ua(2,m)+dble(iy))/(nuc(2)+nvac(2)) +small
              z=(ua(3,m)+dble(iz))/(nuc(3)+nvac(3)) +small
!              if( .not. (x.gt.0.2d0 .and. x.lt.0.4d0.and.
!     &             y.gt.0.2d0 .and. y.lt.0.4d0) ) cycle
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
              tag(inc)= 1d0*is +0.1d0*ifmv +1d-14*inc
            enddo
          enddo
        enddo
      enddo

!!$!.....Add H atoms
!!$      ix= nuc(1)/2
!!$      iy= nuc(2)/2
!!$      iz= nuc(3)/2
!!$      inc=inc+1
!!$      is= 2
!!$      ifmv= 1
!!$      tag(inc)= 1d0*is +0.1d0*ifmv +1d-14*inc
!!$!.....O-site
!!$      ra(1,inc)= (0.5d0 +dble(ix))/(nuc(1)+nvac(1)) +small
!!$      ra(2,inc)= (0.5d0 +dble(iy))/(nuc(2)+nvac(2)) +small
!!$      ra(3,inc)= (0.0d0 +dble(iz))/(nuc(3)+nvac(3)) +small
!!$!.....T-site
!!$      ra(1,inc)= (0.5d0 +dble(ix))/(nuc(1)+nvac(1)) +small
!!$      ra(2,inc)= (0.25d0 +dble(iy))/(nuc(2)+nvac(2)) +small
!!$      ra(3,inc)= (0.0d0 +dble(iz))/(nuc(3)+nvac(3)) +small

      write(6,'(a,i10)') " num of atoms=",inc
!      write(6,'(a,i10)') " id of inc=",nint(mod(tag(inc)*1d14,1d13))

      va(1:3,1:inc)= 0d0

      call write_pmd0_ascii(15,'pmd00000','replace',inc,tag,ra,va,h &
           ,eki,epi,strs)
!!$      call write_pmd0_bin(15,'pmd00000','replace',inc,tag,ra,va,h &
!!$           ,eki,epi,strs)
      
!-----output 'akr000' for Akira visualization
      open(15,file='akr0000',form='formatted',status='replace')
      write(15,'(i10,3i5)') inc, 3, 0, 0
      write(15,'(3es11.3)') ((h(ia,ib,0),ia=1,3),ib=1,3)
      do i=1,inc
        write(15,'(i3,6es11.3)') int(tag(i)),ra(1:3,i),va(1:3,i)
      enddo
      close(15)
      
      end program mkconf_BCC_FeH
!=======================================================================
      subroutine myrnd(rnd,dseed)
      real*8 rnd,dseed
      real*8 d2p31m,d2p31
      save d2p31m,d2p31
      data d2p31m/2147483647d0/
      data d2p31 /2147483648d0/
      
      dseed=dmod(16807d0*dseed,d2p31m)
      rnd=dseed/d2p31
      return
      end subroutine myrnd
!=======================================================================
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make 10mkconf"
!     End:
