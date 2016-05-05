      program mkconf
c-----------------------------------------------------------------------
c     Make atom configuration graphene sheet
c-----------------------------------------------------------------------
      implicit none
      integer:: i,j,k,l,m,n,ii,i00,i0,iz,iy,ix,itemp,inc,ia,ib
      integer:: natm,np(3),ncell(3),nsize1,ifmv,is
      real(8):: alc,h(3,3,0:1),dcc,radius,x0,ua(3,4),x1(1:3)
      real(8),allocatable,dimension(:,:):: ra,va
      real(8),allocatable:: strs(:,:,:),tag(:),eki(:),epi(:)
      real(8):: dseed,rnd1,rnd2,sumvx,sumvy,sumvz,facv

c-----decide ncell: num. of unit cells in each direction
      ncell(1)= 5
      ncell(2)= 5
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
      h(3,3,0)= ncell(3) *3.7d0 *5
      
      write(6,'(a,3i4)') " size:",np(1:3)
      write(6,'(a,i10)') " natm=",natm
      write(6,'(a,es12.4)') " alc=",alc
      write(6,'(a)') " h(1:3,1:3,0):"
      write(6,'(3es12.4)') h(1:3,1:3,0)

      allocate(ra(3,natm),va(3,natm),tag(natm),strs(3,3,natm)
     &     ,eki(natm),epi(natm))

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
              is= 1
              ifmv= 1
              tag(inc)= 1d0*is +0.1d0*ifmv +1d-14*inc
            enddo
          enddo
        enddo
      enddo
      if(inc.ne.natm) stop "error: inc.ne.natm!!"

c-----check pos
      open(90,file="pos_graphene.dat")
      write(90,'(3es15.7)') (ra(1:3,i),i=1,natm)
      close(90)

c-----set initial velocity as zero
      va(1:3,natm)= 0d0

      call write_pmd0_ascii(15,'pmd0000','replace',inc,tag
     &     ,ra,va,h,alc,eki,epi,strs)

c-----output for Akira
      open(15,file='akr0000',form='formatted',status='replace')
      write(15,'(es15.7)') alc
      write(15,'(3es11.3)') (((h(ia,ib,l)/alc,ia=1,3),ib=1,3),l=0,1)
      write(15,'(i10)') natm
      do i=1,inc
        write(15,'(i3,6es11.3)') int(tag(i)),ra(1:3,i),va(1:3,i)
      enddo
      close(15)

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
c     compile-command: "make 10mkconf"
c     End:
