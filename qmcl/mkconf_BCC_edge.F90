      program mkconf_BCC
!-----------------------------------------------------------------------
!  Create a BCC crystal for edge dislocation calculation.
!-----------------------------------------------------------------------
!  Output
!  ------
!    * pmd000-000
!    * akr000
!-----------------------------------------------------------------------
      implicit real*8(a-h,o-z),integer(i-n)
      include './params_unit.h'
!-----max # of atoms
      integer,parameter::nmax=100000
      integer,parameter:: nuamax= 100
!-----# of unit cells
!.....for system1
      integer,parameter:: nuc(1:3)= (/ 25, 20, 1 /) ! default
!      integer,parameter:: nuc(1:3)= (/ 101, 52, 1 /) ! default
!!$!.....for system2
!!$      integer,parameter:: nuc(1:3)= (/ 51, 62, 1 /)
!-----vacuum width in unit of cell
!      integer,parameter:: nvac(1:3)= (/ 5, 5, 5 /)
      integer,parameter:: nvac(1:3)= (/ 0, 2, 0 /)
      integer:: ix,iy,iz,inc,nh
      real(8):: ua(3,nuamax),uh(3),hunit
      real(8):: tag(nmax),ra(3,nmax),va(3,nmax),eki(nmax),epi(nmax) &
           ,h(3,3,0:1),strs(3,3,nmax),h0(3,3),s(3),ymax,ymin,dseed,rnd
!.....H density in atomic ppm
      real(8),parameter:: hppm    = 200d0
!-----alcfe
      real(8),parameter:: alcfe   = 2.835d0
!.....small value
      real(8),parameter:: small   = 1d-8

!.....num of cells in along x must be odd
      if( mod(nuc(1),2).ne.1 ) then
        write(6,'(a)') ' [Error] mod(nuc(1),2)ne.1 !!!'
        stop
      endif
!.....check if num of cells along y is even number
      if( mod(nuc(2),2).ne.0 .or. mod(nvac(2),2).ne.0 ) then
        write(6,'(a)') ' [Error] mod(nuc(2),2).ne.0 !!!'
        write(6,'(a)') ' [Error] mod(nvac(2),2).ne.0 !!!'
        stop
      endif

!.....lattice constant of Fe, see Philos. Mag. 83 (2003) 3977
      alc= alcfe
      daa = sqrt(3d0)/2 *alcfe
      write(6,'(a,es12.4,a)') " Lattice constant =",alc," [Bohr]"
      write(6,'(a,es12.4,a)') " Fe-Fe bond length=",daa," [Bohr]"
      hunit= alc

      call set_system_1(h0,nuamax,ua,hunit,uh)
!      call set_system_2(h0,nuamax,ua,hunit,uh)

!.....simulation box size
      h(1:3,1:3,0:1)= 0d0
      h(1,1,0)= h0(1,1)*(nuc(1)+nvac(1))
      h(2,2,0)= h0(2,2)*(nuc(2)+nvac(2))
      h(3,3,0)= h0(3,3)*(nuc(3)+nvac(3))

!.....width of to-be-removed half plane
      daah= daa/h(1,1,0) +small
      write(6,'(a,es12.4)') ' Burgers vector in reduce unit='&
           ,daah -small

      nh= 0
      inc=0 
      ymin=1d0
      ymax=0d0
      do ix=0,nuc(1)-1
        do iy=0,nuc(2)-1
          do iz=0,nuc(3)-1
            do m=1,12
              s(1)= (ua(1,m)+dble(ix))/(nuc(1)+nvac(1)) +small
              s(2)= (ua(2,m)+dble(iy))/(nuc(2)+nvac(2)) +small
              s(3)= (ua(3,m)+dble(iz))/(nuc(3)+nvac(3)) +small
              s(1:3)= s(1:3) &
                   +dble(nvac(1:3))/(nuc(1:3)+nvac(1:3))/2
!.....remove an atomic plane along yz-plane of bottom half
              if( iy .lt. nuc(2)/2 ) then
                if(  s(1).ge.0.5d0-daah/2 .and. &
                     s(1).lt.0.5d0+daah/2 ) then
                  cycle
!.....dislocation at left
                elseif( s(1).ge.0.5d0-daah*3/2 .and. &
                     s(1).lt.0.5d0-daah/2 ) then
                  s(1)= s(1) +daah
!!$!.....dislocation at right
!!$                elseif( s(1).ge.0.5d0+daah/2 .and. &
!!$                     s(1).lt.0.5d0+daah*3/2 ) then
!!$                  s(1)= s(1) -daah
                endif
              endif
!!$!.....shift linearly to bury removed space
!!$              if( iy .lt. nuc(2)/2 ) then
!!$                if( s(1).lt.0.5d0 ) then
!!$                  s(1)=s(1) +(daah*0.4)/2*s(1)/0.5d0
!!$                else
!!$                  s(1)=s(1) -(daah*0.6)/2*(1d0-s(1))/0.5d0
!!$                endif
!!$              endif
              inc=inc+1
              if(inc.gt.nmax)then
                write(*,*)' [Error] inc>nmax, inc,nmax=',inc,nmax
                stop
              endif
              is= 1
              ifmv= 1
              ra(1:3,inc)= s(1:3)
              tag(inc)= 1d0*is +0.1d0*ifmv +1d-14*inc
!.....search top and bottom in y-direction
              ymin=min(ymin,s(2))
              ymax=max(ymax,s(2))
            enddo
!!$!.....add H at the dislocation core
!!$            if( iy.eq.nuc(2)/2 .and. ix.eq.nuc(1)/2 ) then
!!$              nh= nh +1
!!$              inc=inc+1
!!$              is= 2
!!$              ifmv= 1
!!$              s(1)= (uh(1)+dble(ix))/(nuc(1)+nvac(1)) +small
!!$              s(2)= (uh(2)+dble(iy))/(nuc(2)+nvac(2)) +small
!!$              s(3)= (uh(3)+dble(iz))/(nuc(3)+nvac(3)) +small
!!$              s(1:3)= s(1:3) &
!!$                   +dble(nvac(1:3))/(nuc(1:3)+nvac(1:3))/2
!!$              ra(1:3,inc)= s(1:3)
!!$              tag(inc)= 1d0*is +0.1d0*ifmv +1d-14*inc
!!$            endif
          enddo
        enddo
      enddo
!.....Set top and bottom atoms to ifmv=2
      do i=1,inc
        if( ra(2,i).lt.ymin+daah ) then
          is= int(tag(i))
          ifmv=2
          tag(i)= 1d0*is +0.1d0*ifmv +1d-14*i
        elseif( ra(2,i).gt.ymax-daah ) then
          is= int(tag(i))
          ifmv=2
          tag(i)= 1d0*is +0.1d0*ifmv +1d-14*i
        endif
      enddo

      write(6,'(a,i10)') " num of atom=",inc
      write(6,'(a,i10)') " num of H   =",nh

      va(1:3,1:inc)= 0d0

      call write_pmd_ascii(15,'pmd000-000',inc,h,hunit,tag,ra,va &
           ,eki,epi,strs)
      
!-----output 'akr000' for Akira visualization
      open(15,file='akr000',form='formatted',status='replace')
      write(15,'(es15.7)') hunit
      write(15,'(3es11.3)') ((h(ia,ib,0)/hunit,ia=1,3),ib=1,3)
      write(15,'(i10,3i5)') inc, 3, 0, 0
      do i=1,inc
        write(15,'(i3,6es11.3)') int(tag(i)),ra(1:3,i),va(1:3,i)
      enddo
      close(15)
      
      end program mkconf_BCC
!=======================================================================
      subroutine set_system_1(h0,nuamax,ua,alc,uh)
!    a1=[-1,-1,-1]
!    a2=[-1,-1, 2]
!    a3=[-1, 1, 0]
        implicit none
        integer,intent(in):: nuamax
        real(8),intent(out):: h0(3,3),ua(3,nuamax),uh(3)
        real(8),intent(in):: alc

!.....unit vectors, h0= (a1,a2,a3) where a1,a2,a3 are column vectors
        h0(1:3,1:3)= 0d0
        h0(1,1)=  1.73205080756888d0 *alc
        h0(2,2)=  2.44948974278318d0 *alc
        h0(3,3)=  1.41421356237310d0 *alc
!.....atom positions in the unit cell
        ua(1:3,1) =(/ 1.666666d-01, 3.333333d-01, 0.000d+00 /)
        ua(1:3,2) =(/ 0.000000d+00, 0.000000d+00, 0.000d+00 /)
        ua(1:3,3) =(/ 5.000000d-01, 0.000000d+00, 0.000d+00 /)
        ua(1:3,4) =(/ 1.666666d-01, 8.333333d-01, 5.000d-01 /)
        ua(1:3,5) =(/ 0.000000d+00, 5.000000d-01, 5.000d-01 /)
        ua(1:3,6) =(/ 5.000000d-01, 5.000000d-01, 5.000d-01 /)
        ua(1:3,7) =(/ 3.333333d-01, 1.666666d-01, 5.000d-01 /)
        ua(1:3,8) =(/ 8.333333d-01, 1.666666d-01, 5.000d-01 /)
        ua(1:3,9) =(/ 3.333333d-01, 6.666666d-01, 0.000d+00 /)
        ua(1:3,10)=(/ 8.333333d-01, 6.666666d-01, 0.000d+00 /)
        ua(1:3,11)=(/ 6.666666d-01, 3.333333d-01, 0.000d+00 /)
        ua(1:3,12)=(/ 6.666666d-01, 8.333333d-01, 5.000d-01 /)

!.....hydrogen position, T-site
        uh(1:3)= (/ 7d0/12, 1d0/24, 3d0/8 /)


      end subroutine set_system_1
!=======================================================================
      subroutine set_system_2(h0,nuamax,ua,alc,uh)
!    a1=[-1,-1,-1]
!    a2=[-1, 1, 0]
!    a3=[-1,-1, 2]
        implicit none
        integer,intent(in):: nuamax
        real(8),intent(out):: h0(3,3),ua(3,nuamax),uh(3)
        real(8),intent(in):: alc

!.....unit vectors, h0= (a1,a2,a3) where a1,a2,a3 are column vectors
        h0(1:3,1:3)= 0d0
        h0(1,1)=  1.73205080756888d0 *alc
        h0(2,2)=  1.41421356237310d0 *alc
        h0(3,3)=  2.44948974278318d0 *alc
!.....atom positions in the unit cell
        ua(1:3,1) =(/ 1.666666d-01, 0.000d+00, 3.333333d-01 /)
        ua(1:3,2) =(/ 0.000000d+00, 0.000d+00, 0.000000d+00 /)
        ua(1:3,3) =(/ 5.000000d-01, 0.000d+00, 0.000000d+00 /)
        ua(1:3,4) =(/ 1.666666d-01, 5.000d-01, 8.333333d-01 /)
        ua(1:3,5) =(/ 0.000000d+00, 5.000d-01, 5.000000d-01 /)
        ua(1:3,6) =(/ 5.000000d-01, 5.000d-01, 5.000000d-01 /)
        ua(1:3,7) =(/ 3.333333d-01, 5.000d-01, 1.666666d-01 /)
        ua(1:3,8) =(/ 8.333333d-01, 5.000d-01, 1.666666d-01 /)
        ua(1:3,9) =(/ 3.333333d-01, 0.000d+00, 6.666666d-01 /)
        ua(1:3,10)=(/ 8.333333d-01, 0.000d+00, 6.666666d-01 /)
        ua(1:3,11)=(/ 6.666666d-01, 0.000d+00, 3.333333d-01 /)
        ua(1:3,12)=(/ 6.666666d-01, 5.000d-01, 8.333333d-01 /)

!.....hydrogen position, T-site
        uh(1:3)= (/ 7d0/12, 1d0/24, 3d0/8 /)

      end subroutine set_system_2
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
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make 10mkconf"
!     End:
