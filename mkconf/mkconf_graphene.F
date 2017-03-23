      program mkconf
c-----------------------------------------------------------------------
c     Make atom configuration graphene sheet
c-----------------------------------------------------------------------
      implicit none
      integer:: i,j,k,l,m,n,ii,i00,i0,iz,iy,ix,itemp,inc,ia,ib
      integer:: natm,np(3),ncell(3),nsize1,ifmv,is,nvac(3)
      real(8):: alc,h(3,3,0:1),dcc,radius,x0,ua(3,4),x1(1:3)
      real(8),allocatable,dimension(:,:):: ra,va
      real(8),allocatable:: strs(:,:,:),tag(:),eki(:),epi(:)

c-----decide ncell: num. of unit cells in each direction
      ncell(1)= 40
      ncell(2)= 1
      ncell(3)= 20

c.....vacuum region
      nvac(1) = 2
      nvac(2) = 0
      nvac(3) = 10
      
c-----set system size
      np(1:3)= ncell(1:3) +nvac(1:3)
c      np(3)= 1
      natm= ncell(1)*ncell(2)*ncell(3)*4
c-----equilibrium length of C-C in graphene
      dcc= 1.46d0 *0.975615d0
c      dcc= 1.65d0 *0.975615d0
      alc= dcc*sqrt(3d0)
      h(1:3,1:3,0:1)= 0d0
      h(1,1,0)= np(1) *alc
      h(2,2,0)= np(2) *3.7d0 *5
      h(3,3,0)= np(3) *sqrt(3d0)*alc
      
      write(6,'(a,3i4)') " size:",ncell(1:3)
      write(6,'(a,i10)') " natm=",natm
      write(6,'(a,es12.4)') " alc=",alc
      write(6,'(a)') " h(1:3,1:3,0):"
      write(6,'(3es12.4)') h(1:3,1:3,0)

      allocate(ra(3,natm),va(3,natm),tag(natm),strs(3,3,natm)
     &     ,eki(natm),epi(natm))

c-----make graphene sheet
      ua(1:3,1)= (/ 0.0d0,   0.5d0,       0.0d0 /)
      ua(1:3,2)= (/ 0.5d0,   0.5d0,       1d0/6 /)
      ua(1:3,3)= (/ 0.5d0,   0.5d0,       0.5d0 /)
      ua(1:3,4)= (/ 0.0d0,   0.5d0, 0.5d0+1d0/6 /)

c-----extend system
      inc=0
      do i=0,ncell(1)-1
        do j=0,ncell(2)-1
          do k=0,ncell(3)-1
            do m=1,4
              inc=inc+1
              ra(1,inc)= (ua(1,m) + i)/np(1) +dble(nvac(1))/2/np(1)
              ra(2,inc)= (ua(2,m) + j)/np(2) +dble(nvac(2))/2/np(2)
              ra(3,inc)= (ua(3,m) + k)/np(3) +dble(nvac(3))/2/np(3)
              is= 1
              ifmv= 1
              tag(inc)= 1d0*is +0.1d0*ifmv +1d-14*inc
            enddo
          enddo
        enddo
      enddo
      if(inc.ne.natm) stop "error: inc.ne.natm!!"

c-----set initial velocity as zero
      va(1:3,natm)= 0d0

      call write_pmd0_ascii(15,'pmd0000','replace',inc,tag
     &     ,ra,va,h,alc,eki,epi,strs)

      end program mkconf
c-----for emacs---------------------------------------------------------
c     Local Variables:
c     compile-command: "make mkconf_graphene"
c     End:
