      program mkconf
c-----------------------------------------------------------------------
c     Make atom configuration graphene sheet
c-----------------------------------------------------------------------
      implicit none
      real(8),parameter:: pi = 3.14159265358979d0
      real(8),parameter:: pi2 = 2d0*pi
      real(8),parameter:: umass   = 1822.9d0
      real(8),parameter:: tempau= 3.15785d5
c-----initial temperature
      real(8),parameter:: tini = 300d0 /tempau
c-----mass of C
      real(8),parameter:: am   = 12.0107d0 *umass

      integer:: i,j,k,l,m,n,ii,i00,i0,iz,iy,ix,itemp,inc,ia,ib
      integer:: natm,np(3),ncell(3),nsize1
      real(8):: alc,h(3,3,0:1),dcc,radius,x0,ua(3,4),x1(1:3)
      integer,allocatable:: ifmv(:),is(:)
      real(8),allocatable,dimension(:,:):: ra,va
      real(8):: dseed,rnd1,rnd2,sumvx,sumvy,sumvz,facv

c-----decide ncell: num. of unit cells in each direction
      ncell(1)= 10
      ncell(2)= 10
      ncell(3)= 2

c-----set system size
      np(1:3)= ncell(1:3)
c      np(3)= 1
      natm= np(1)*np(2)*np(3)*4
c-----equilibrium length of C-C in graphene
c      dcc= 1.46d0 *0.975615d0
      dcc= 1.65d0 *0.975615d0
      alc= dcc*sqrt(3d0)
      h(1:3,1:3,0:1)= 0d0
      h(1,1,0)= np(1) *alc
      h(2,2,0)= np(2) *sqrt(3d0)*alc
      h(3,3,0)= h(1,1,0)
      
      write(6,'(a,3i4)') " size:",np(1:3)
      write(6,'(a,i10)') " natm=",natm
      write(6,'(a,es12.4)') " alc=",alc
      write(6,'(a)') " h(1:3,1:3,0):"
      write(6,'(3es12.4)') h(1:3,1:3,0)

      allocate(is(natm),ifmv(natm),ra(3,natm),va(3,natm))

c-----make graphene sheet
      ua(1:3,1)= (/ 0.0d0,       0.0d0,   0.5d0 /)
      ua(1:3,2)= (/ 0.5d0,       1d0/6,   0.5d0 /)
      ua(1:3,3)= (/ 0.5d0,       0.5d0,   0.5d0 /)
      ua(1:3,4)= (/ 0.0d0, 0.5d0+1d0/6,   0.5d0 /)

c-----extend system
      inc=0
      do i=0,np(1)-1
        do j=0,np(2)-1
          do k=0,np(3)-1
            do m=1,4
              inc=inc+1
              ra(1,inc)= (ua(1,m) + i)/np(1)
              ra(2,inc)= (ua(2,m) + j)/np(2)
              ra(3,inc)= (ua(3,m) + k)/np(3)
            enddo
          enddo
        enddo
      enddo
      if(inc.ne.natm) stop "error: inc.ne.natm!!"

c-----check pos
      open(90,file="pos_graphene.dat")
      write(90,'(3es15.7)') (ra(1:3,i),i=1,natm)
      close(90)

c-----set is,ifmv
      is(1:natm)= 1
      ifmv(1:natm)= 1

c-----set initial velocity as zero
      va(1:3,natm)= 0d0

c-----velocities in Maxwell-Boltzmann distribution
      facv= dsqrt(2d0*tini/am)
      dseed=12345
      do i=1,natm
        do l=1,3
          call myrnd(rnd1,dseed)
          call myrnd(rnd2,dseed)
          va(l,i)=facv*dsqrt(-dlog(rnd1))*dcos(2d0*pi*rnd2)
        enddo
c-------collision test
        if(ra(3,i).gt.0.5d0) then
          va(3,i)= -facv/2
        else
          va(3,i)= facv/2
        endif
      enddo
c-----set center of mass motion to zero
      sumvx=0d0
      sumvy=0d0
      sumvz=0d0
      do i=1,natm
        sumvx=sumvx+va(1,i)
        sumvy=sumvy+va(2,i)
        sumvz=sumvz+va(3,i)
      enddo
      do i=1,natm
        va(1,i)=va(1,i)-sumvx/dble(natm)
        va(2,i)=va(2,i)-sumvy/dble(natm)
        va(3,i)=va(3,i)-sumvz/dble(natm)
      enddo
c-----scale velocities to reduced unit
      do i=1,natm
        va(1,i)=va(1,i) /h(1,1,0)
        va(2,i)=va(2,i) /h(2,2,0)
        va(3,i)=va(3,i) /h(3,3,0)
      enddo
      
c-----output to a file
      open(10,file='ini000',form='formatted',status='replace')
c-----num. of atoms
      write(10,'(i10)') natm
c-----h matrix
      write(10,'(3es22.14)') (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
c-----positions, velocities
      do i=1,natm
        write(10,'(2i4,10es22.14)') is(i),ifmv(i),ra(1:3,i),va(1:3,i)
      enddo
      close(10)

c-----output 'kvs000' for KVS visualization
      open(15,file='kvs000',form='formatted',status='replace')
      write(15,'(i10)') natm
      write(15,'(3es11.3)') (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
      do i=1,inc
        write(15,'(i3,6es11.3)') is(i),ra(1:3,i),va(1:3,i)
      enddo
      close(15)

      deallocate(is,ifmv,ra,va)
      end program mkconf
c=======================================================================
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
c-----for emacs---------------------------------------------------------
c     Local Variables:
c     compile-command: "make 10mkconf_graphene"
c     End:
