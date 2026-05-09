      program cna2disl
c-----------------------------------------------------------------------
c  Create dislocation line profile from CNA data
c-----------------------------------------------------------------------
c  USAGE:
c    $ ./cna2disl < cna000 > disl000
c-----------------------------------------------------------------------
      use mod_precision
      implicit none
      real(rp),parameter:: bohr = 0.5291772e-10_rp
      real(rp),parameter:: rcz  = 2.8553e-10_rp /bohr

      integer:: i,iz,natm,naux,nz,nc
      real(rp):: h(3,3),dz,xav,yav,xref,yref,xi,yi,zi
      real(rp),allocatable:: ra(:,:),aux(:,:)
      integer,allocatable:: itype(:),lshd(:),lscl(:)

c-----read cna file from STDIN
      read(5,*) natm,naux
      read(5,*) h(1,1:3)
      read(5,*) h(2,1:3)
      read(5,*) h(3,1:3)
      allocate(ra(3,natm),aux(naux,natm),itype(natm))
      do i=1,natm
        read(5,*) itype(i),ra(1:3,i),aux(1:naux,i)
      enddo

      nz= int(h(3,3)/rcz)+1
      dz= 1.0_rp/nz
      allocate(lshd(nz),lscl(natm))

c.....Make lshd and lscl
      lshd(1:nz)= 0
      lscl(1:natm)= 0
      do i=1,natm
c.....Only dislocation core atoms
        if( ra(2,i).gt.0.75_rp .or. ra(2,i).lt.0.25_rp )  cycle
c.....Assign PBC
        if( ra(1,i).ge.1.0_rp ) ra(1,i)=ra(1,i) -1.0_rp
        if( ra(1,i).lt.0.0_rp ) ra(1,i)=ra(1,i) +1.0_rp
        if( ra(2,i).ge.1.0_rp ) ra(2,i)=ra(2,i) -1.0_rp
        if( ra(2,i).lt.0.0_rp ) ra(2,i)=ra(2,i) +1.0_rp
        if( ra(3,i).ge.1.0_rp ) ra(3,i)=ra(3,i) -1.0_rp
        if( ra(3,i).lt.0.0_rp ) ra(3,i)=ra(3,i) +1.0_rp
c.....Store atoms to the lists
        iz= ra(3,i)/dz +1
        lscl(i)= lshd(iz)
        lshd(iz)= i
      enddo

c.....Search and list up dislocation core atoms
      do iz=1,nz
        nc = 0
        xav= 0.0_rp
        yav= 0.0_rp
        i= lshd(iz)
        if( iz.eq.1 ) then
          xref= ra(1,i)
          yref= ra(2,i)
        endif
        do
          if(i.eq.0) exit
          xi= ra(1,i)
          yi= ra(2,i)

          xav= xav +xi -anint(xi-xref)
          yav= yav +yi -anint(yi-yref)
          nc= nc +1

          i= lscl(i)
        enddo
        if( nc.eq.0 ) stop ' [Error] nc.eq.0 !!!'
        xav= xav/nc
        yav= yav/nc
        zi = dz*(iz-1) +dz/2
        write(6,'(4es15.7)') zi,zi*h(3,3)*0.5291772_rp,xav,yav
      enddo

      end program cna2disl
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make cna2disl"
c     End:
