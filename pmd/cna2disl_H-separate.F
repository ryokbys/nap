      program cna2disl
c-----------------------------------------------------------------------
c.....Create dislocation line profile from CNA data
c-----------------------------------------------------------------------
c.....USAGE:
c.....$ ./cna2disl < cna000 > disl000
c-----------------------------------------------------------------------
      implicit none
      real(8),parameter:: bohr = 0.5291772d-10
      real(8),parameter:: rcz  = 2.8553d-10 /bohr
c.....real(8),parameter:: rcx  = sqrt(3d0)/2 *rcz
      real(8),parameter:: rcx  = rcz
      integer,parameter:: nnmax= 100
      integer,parameter:: ngmax= 5

      integer:: i,j,k,iz,natm,naux,nz,nc,igh,ig,ng
      real(8):: h(3,3),dz,xref,yref,xi,yi,zi,dx,nx,x
      real(8),allocatable:: ra(:,:),aux(:,:),xav(:),yav(:)
      integer,allocatable:: itype(:),lshd(:),lscl(:),lspr(:,:),igrp(:)

c-----read cna file from STDIN
      read(5,*) natm,naux
      read(5,*) h(1,1:3)
      read(5,*) h(2,1:3)
      read(5,*) h(3,1:3)
      allocate(ra(3,natm),aux(naux,natm),itype(natm))
      do i=1,natm
        read(5,*) itype(i),ra(1:3,i),aux(1:naux,i)
      enddo

      nx= int(h(1,1)/rcx)+1
      dx= 1d0/nx
      nz= int(h(3,3)/rcz)+1
      dz= 1d0/nz
      allocate(lshd(nz),lscl(natm),lspr(0:nnmax,natm),igrp(natm)
     &     ,xav(ngmax),yav(ngmax))

c.....Make lshd and lscl
      lshd(1:nz)= 0
      lscl(1:natm)= 0
      do i=1,natm
c.....Only dislocation core atoms
        if( ra(2,i).gt.0.75d0 .or. ra(2,i).lt.0.25d0 )  cycle
c.....Assign PBC
        if( ra(1,i).ge.1d0 ) ra(1,i)=ra(1,i) -1d0
        if( ra(1,i).lt.0d0 ) ra(1,i)=ra(1,i) +1d0
        if( ra(2,i).ge.1d0 ) ra(2,i)=ra(2,i) -1d0
        if( ra(2,i).lt.0d0 ) ra(2,i)=ra(2,i) +1d0
        if( ra(3,i).ge.1d0 ) ra(3,i)=ra(3,i) -1d0
        if( ra(3,i).lt.0d0 ) ra(3,i)=ra(3,i) +1d0
c.....Store atoms to the lists
        iz= ra(3,i)/dz +1
        lscl(i)= lshd(iz)
        lshd(iz)= i
      enddo

c.....Detect disl core clusters
      lspr(0:nnmax,1:natm)= 0
      do i=1,natm
c.....Only dislocation core atoms
        if( ra(2,i).gt.0.75d0 .or. ra(2,i).lt.0.25d0 )  cycle
        iz= ra(3,i)/dz +1
        j= lshd(iz)
        do
          if(j.eq.0) exit
          x= h(1,1)*abs(ra(1,j)-ra(1,i) -anint(ra(1,j)-ra(1,i)))
          if( x.lt.rcx ) then
            lspr(0,i)=lspr(0,i) +1
            lspr(lspr(0,i),i)= j
          endif
          j=lscl(j)
        enddo
      enddo

      igrp(1:natm)= 0
      do iz=1,nz
        ig=0
        i= lshd(iz)
        do
          if(i.eq.0) exit
          call make_bonded_cluster(i,lspr,igrp,ig,natm,nnmax)
          i= lscl(i)
        enddo

c.....Check which cluster includes H atoms
        igh= 0
        i=lshd(iz)
        do
          if(i.eq.0) exit
          if(itype(i).eq.2) then ! H atom
            igh= igrp(i)
            exit
          endif
          i=lscl(i)
        enddo

        ng=ig
        if(ng.gt.ngmax) stop ' [Error] ng.gt.ngmax !!!'
        xav(1:ngmax)= 0d0
        yav(1:ngmax)= 0d0
        do ig=1,ng
          nc = 0
          xav(ig)= 0d0
          yav(ig)= 0d0
          i= lshd(iz)
          do
            if(i.eq.0) exit
            if(igrp(i).eq.ig) exit
            i=lscl(i)
          enddo
          xref= ra(1,i)
          yref= ra(2,i)
          do
            if(i.eq.0) exit
            if(igrp(i).ne.ig) then
              i=lscl(i)
              cycle
            endif
            xi= ra(1,i)
            yi= ra(2,i)

            xav(ig)= xav(ig) +xi -anint(xi-xref)
            yav(ig)= yav(ig) +yi -anint(yi-yref)
            nc= nc +1

            i= lscl(i)
          enddo
          if( nc.eq.0 ) stop ' [Error] nc.eq.0 !!!'
          xav(ig)= xav(ig)/nc
          yav(ig)= yav(ig)/nc
        enddo
        zi = dz*(iz-1) +dz/2
        if(ng.gt.1) then
          if(igh.eq.1) then
            write(6,'(10es15.7)') zi,zi*h(3,3)*0.5291772d0,xav(1),yav(1)
     &           ,xav(2),yav(2)
          else
            write(6,'(10es15.7)') zi,zi*h(3,3)*0.5291772d0,xav(2),yav(2)
     &           ,xav(1),yav(1)
          endif
        elseif(ng.eq.1) then
          write(6,'(10es15.7)') zi,zi*h(3,3)*0.5291772d0,xav(1),yav(1)
     &         ,xav(1),yav(1)
        endif
      enddo

      end program cna2disl
c=======================================================================
      recursive subroutine make_bonded_cluster(i,lspr,igrp,ig
     &     ,natm,nnmax)
      implicit none
      integer,intent(in):: i,nnmax,natm,lspr(0:nnmax,natm)
      integer,intent(inout):: igrp(natm),ig
      integer:: j,k

      if(igrp(i).eq.0) then
        ig=ig+1
        igrp(i)=ig
      endif

      do k=1,lspr(0,i)
        j=lspr(k,i)
        if(igrp(j).eq.0) then
          igrp(j)= ig
          call make_bonded_cluster(j,lspr,igrp,ig,natm,nnmax)
        endif
      enddo
      
      end  subroutine make_bonded_cluster
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make cna2disl_H"
c     End:
