      program pmd2cna
c-----------------------------------------------------------------------
c Reduce atoms from mts file using Common Neighbor Analysis.
c   - FCC crystal only?
c   - change cutoff radius, rc, to apply this to your system
c  
c-----------------------------------------------------------------------
c
c  << Input files >>
c    pmd.in: simulation config file
c    pmd###-???: MD-coordination files,
c                ### is the node number and ??? is the frame number
c
c  << Output files >>
c    cna???: cna files
c
c-----------------------------------------------------------------------
      implicit none
      include "mpif.h"
c-----max. num. of atoms
      integer,parameter:: namax = 1000000
c-----max. num. of species
      integer,parameter:: nismax= 10
c-----max. num. of boundary-particles
      integer,parameter:: nbmax = 500000
c-----max. num. of neighbors
      integer,parameter:: nnmax = 30
c-----max. num. of buffer data
      integer,parameter:: maxbuf= 100000
c-----cutoff radius
      real(8),parameter:: rc1   = 8.d-10  /0.5291772d-10
      real(8),parameter:: rc    = 3.5d-10 /0.5291772d-10

      integer:: i,j,k,l,m,n,ia,ib
      integer:: nstp,nouterg,noutpmd,iout_count,istp,nerg,ifpmd,npmd
      integer:: natm,nb,ntot,nis,itmp,ipmd
      real(8):: dt
      real(8):: tcpu,tcpu1,tcpu2,tcom
      real(8):: tmp,dmp,treq,ttol,epot,ekin,epot0
c-----parallel-related variables
      integer:: nx,ny,nz,nxyz
      integer:: nn(6),myparity(3),lsrc(6),lsb(0:nbmax,6)
     &     ,myid_md,nodes_md,mpi_md_world,myx,myy,myz,ierr,idis
      real(8):: sv(3,6),anxi,anyi,anzi,sorg(3)
c-----simulation box
      real(8):: hunit,h(3,3,0:1),hi(3,3),volume,sgm(3,3)
c-----positions, velocities, and accelerations
      integer:: itype(namax)
      real(8):: ra(3,namax),va(3,namax),op(10,namax),tag(namax)
c-----buffer data
      integer:: ibuf(maxbuf)
      real(8):: rbuf(12,maxbuf)
c-----cna data
      integer:: ncna,ncnatot
      integer:: icna(namax),lscna(namax)
      integer,allocatable:: lcna(:),nrecvs(:),idisps(:)
      real(8):: xcna(12,namax)
c-----species, index of FMV
      integer:: lspr(0:nnmax,namax)
c-----output file names
      character:: cnode*3,cnum*3
      
c-----initialize the MPI environment
      call mpi_init(ierr)
c-----total number of MD-nodes
      call mpi_comm_size(MPI_COMM_WORLD, nodes_md, ierr)
c-----my rank in MD-nodes
      call mpi_comm_rank(MPI_COMM_WORLD, myid_md, ierr)

c-----read input parameters, this should be after mpi_{init,size,rank}
c  nx,ny,nz: num of divisions in each direction
c  dt:     time step in sec.
c  nstp:   num. of steps
c  rc:     cutoff radius in Angstrom
c  dmp:    damping factor
c  treq:   required temperature in K
c  ttol:   temperature tolerence in percent
c  nerg:   num. of energy outputs
c  ifpmd:  where or not to write pmd files
c  npmd:   num. of pmd outputs
      open(10,file="in.control",status="old")
c      read(10,*) nx,ny,nz
      nx=1
      ny=1
      nz=1
      read(10,*) itmp,tmp       ! dt, nstp
      read(10,*) nerg,npmd
      close(10)
c-----parallel configuration
      nxyz= nx*ny*nz
      anxi= 1d0/nx
      anyi= 1d0/ny
      anzi= 1d0/nz
c-----error trap
      if(nodes_md.ne.nxyz) then
        write(6,'(a)') " error: nodes_md .ne. nxyz!!"
        call mpi_finalize(ierr)
        stop
      endif
      mpi_md_world= MPI_COMM_WORLD
      
c-----allocation
      allocate(lcna(nxyz),nrecvs(0:nxyz-1),idisps(0:nxyz-1))

c-----vector node indices: range [0:nx-1]
      myx=myid_md/(ny*nz)
      myy=mod(myid_md/nz,ny)
      myz=mod(myid_md,nz)
c-----reduced node origin
      sorg(1)= anxi*myx
      sorg(2)= anyi*myy
      sorg(3)= anzi*myz
c-----ntset
      call ntset(myx,myy,myz,nx,ny,nz,nn,sv,myparity,anxi,anyi,anzi)
      call get_cnum(cnode,myid_md)

c-----loop over num. of pmd outputs
      do ipmd=0,npmd
        if(myid_md.eq.0) write(6,'(a,i5)') " ipmd=",ipmd
c-------each node opens one MD-configuration file
        call get_cnum(cnum,ipmd)
        open(20,file="pmd"//cnode//"-"//cnum,status="old")
c        open(20,file="fin"//cnode,status="old")
c-------natm: num. of particles in this node
        read(20,*) hunit
        read(20,*) (((h(ia,ib,l)*hunit,ia=1,3),ib=1,3),l=0,1)
        read(20,*) natm
        nis= 0
        do i=1,natm
c---------species, positions, velocities, index of motion
c          read(20,*) itype(i),itmp,ra(1:3,i),va(1:3,i)
          read(20,*) tag(i),ra(1:3,i),va(1:3,i),op(1:5,i)
          ra(1:3,i)= ra(1:3,i) -sorg(1:3)
        enddo
        close(20)

c=======make neighbor list
c-------set HI and SGM
        call boxmat(h,hi,volume,sgm)
c-------copy RA of boundary atoms
        call bacopy(tcom,sgm,volume,lsb,nbmax,ra,namax,natm,nb
     &       ,anxi,anyi,anzi,nn,itype,rc1,myid_md,myparity,lsrc,sv
     &       ,mpi_md_world)
c-------make pair list, LSPR
        call mk_lspr(namax,natm,nb,nnmax,itype,ra,rc,h,hi
     &       ,anxi,anyi,anzi,lspr,myid_md)
c-------common neighbor analysis
        call cna(namax,natm,nb,nnmax,lspr,ra,va,tag,ncna,lscna)
cc-------reduce num. of atoms
c        call z_reduce(namax,natm,ra,ncna,lscna)
c-------gather num. of cna atoms
        call mpi_reduce(ncna,ncnatot,1,MPI_INTEGER,MPI_SUM
     &       ,0,mpi_md_world,ierr)
        if(myid_md.eq.0) write(6,'(a,i8)') " ncnatot=",ncnatot
c-------error trap
        if(ncnatot.gt.namax) then
          write(*,*) " !!!ncnatot.gt.namax!!!"
          call mpi_finalize(ierr)
          stop 
        endif
        call mpi_allgather(ncna,1,MPI_INTEGER
     &       ,lcna,1,MPI_INTEGER,mpi_md_world,ierr)
        idisps(0:nodes_md-1)=0
        do i=1,nodes_md-1
          idisps(i)= idisps(i-1) +lcna(i)
        enddo
        do i=0,nodes_md-1
          nrecvs(i)= lcna(i+1)
        enddo
c-------gather
        do i=1,ncna
          itmp= lscna(i)
          rbuf(1,i)   = tag(itmp)
          rbuf(2:4,i) = ra(1:3,itmp) +sorg(1:3)
          rbuf(5:7,i) = va(1:3,itmp)
          rbuf(8:12,i)= op(1:5,itmp)
        enddo
c        if(myid_md.eq.0)write(6,'(a)') " 1st gatherv..."
c        call mpi_gatherv(ibuf,ncna,MPI_INTEGER
c     &       ,icna,nrecvs,idisps,MPI_INTEGER
c     &       ,0,mpi_md_world,ierr)
        do i=0,nodes_md-1
          idisps(i)= idisps(i) *12
          nrecvs(i)= nrecvs(i) *12
        enddo
c        if(myid_md.eq.0)write(6,'(a)') " 2nd gatherv..."
        if(myid_md.eq.0)write(6,'(a)') " gatherv..."
        call mpi_gatherv(rbuf,12*ncna,MPI_DOUBLE_PRECISION
     &       ,xcna,nrecvs,idisps,MPI_DOUBLE_PRECISION
     &       ,0,mpi_md_world,ierr)

c-------output as cna file
        if(myid_md.eq.0) then
          write(6,'(a)') " writing file cna"//cnum//"..."
          open(21,file="cna"//cnum,status="replace")
          write(21,'(es15.7)') hunit
          write(21,'(3es12.4)') ((h(ia,ib,0)/hunit,ia=1,3),ib=1,3)
          write(21,'(i10,3i5)') ncnatot, 8,0,0
          do i=1,ncnatot
            write(21,'(i4,20es12.4)') int(xcna(1,i)),xcna(2:12,i)
          enddo
          close(21)
        endif
      enddo

      call mpi_finalize(ierr)

      end program pmd2cna
c=======================================================================
      subroutine setup(nismax,am,acon,fack,dt)
      implicit none
      include "params_unit.h"
      include "params_eam_al.h"
c      include "params_lj_ar.h"
      integer,intent(in):: nismax
      real(8),intent(in):: dt
      real(8),intent(out):: am(nismax),acon(nismax),fack(nismax)

      integer:: i

c-----atomic masses
      am(1:nismax)= am_al
c      am(1:nismax)= am_ar

c-----prefactors for normalized accelerration
      do i=1,nismax
        acon(i)= 0.5d0*dt**2/am(i)
      enddo

c-----prefactors for kinetic energy, FACK
      do i=1,nismax
        fack(i)= 0.5d0*am(i)/dt**2
      enddo

      end subroutine setup
c=======================================================================
      subroutine boxmat(h,hi,volume,sgm)
c-----------------------------------------------------------------------
c  setup matrices of MD-box
c    H:   MD-box matrix
c    HI:  inverse MD-box matrix
c    SGM: cofactor matrix
c-----------------------------------------------------------------------
      implicit none
      real(8),intent(in):: h(3,3,0:1)
      real(8),intent(out):: volume,sgm(3,3),hi(3,3)

      real(8):: hit(3,3)
      integer:: i,j,im,ip,jm,jp

c-----cofactor matrix, SGM
      do j=1,3
        jm=mod(j+1,3)+1
        jp=mod(j,  3)+1
        do i=1,3
          im=mod(i+1,3)+1
          ip=mod(i,  3)+1
          sgm(i,j)=h(ip,jp,0)*h(im,jm,0)-h(im,jp,0)*h(ip,jm,0)
        enddo
      enddo
c-----MD-box volume
      volume=h(1,1,0)*sgm(1,1)+h(2,1,0)*sgm(2,1)+h(3,1,0)*sgm(3,1)
      do j=1,3
        do i=1,3
          hit(i,j)= sgm(i,j)/volume
        enddo
      enddo
c-----transpose
      do j=1,3
        do i=1,3
          hi(i,j)= hit(j,i)
        enddo
      enddo

      return
      end
c=======================================================================
      subroutine ntset(myx,myy,myz,nx,ny,nz,nn,sv,myparity,anxi,anyi
     &     ,anzi)
c-----------------------------------------------------------------------
c  Preparation for network related properties
c-----------------------------------------------------------------------
      implicit none
      integer,intent(in):: myx,myy,myz,nx,ny,nz
      real(8),intent(in):: anxi,anyi,anzi
      integer,intent(out):: nn(6),myparity(3)
      real(8),intent(out):: sv(3,6)
      integer:: iv(3,6),ku,k1x,k1y,k1z

      iv(1:3,1)= (/ -1, 0, 0 /)
      iv(1:3,2)= (/  1, 0, 0 /)
      iv(1:3,3)= (/  0,-1, 0 /)
      iv(1:3,4)= (/  0, 1, 0 /)
      iv(1:3,5)= (/  0, 0,-1 /)
      iv(1:3,6)= (/  0, 0, 1 /)

      do ku=1,6
        k1x=mod(myx+iv(1,ku)+nx,nx)
        k1y=mod(myy+iv(2,ku)+ny,ny)
        k1z=mod(myz+iv(3,ku)+nz,nz)
c-------scalar neighbor ID, nn
        nn(ku)=k1x*(ny*nz)+k1y*nz+k1z
c-------shift vector for exchnaging ra
        sv(1,ku)= anxi*iv(1,ku)
        sv(2,ku)= anyi*iv(2,ku)
        sv(3,ku)= anzi*iv(3,ku)
      enddo
      
c-----Set up the node parity table
      if (nx.eq.1) then
        myparity(1)=2
      else
        myparity(1)=mod(myx,2)
      endif

      if (ny.eq.1) then
        myparity(2)=2
      else
        myparity(2)=mod(myy,2)
      endif

      if (nz.eq.1) then
        myparity(3)=2
      else
        myparity(3)=mod(myz,2)
      endif

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine bacopy(tcom,sgm,volume,lsb,nbmax,ra,namax,natm,nb
     &     ,anxi,anyi,anzi,nn,is,rc,myid_md,myparity,lsrc,sv
     &     ,mpi_md_world)
c-----------------------------------------------------------------------
c     Exchanges boundary-atom data among neighbor nodes: is and ra
c-----------------------------------------------------------------------
      implicit none
      include 'mpif.h'
      integer,intent(in):: namax,nbmax,nn(6)
     &     ,myid_md,myparity(3),mpi_md_world
      integer,intent(inout):: natm,nb,lsb(0:nbmax,6),is(namax),lsrc(6)
      real(8),intent(in):: sv(3,6),sgm(3,3),volume,anxi,anyi,anzi,rc
      real(8),intent(inout):: ra(3,namax),tcom

c      integer:: status(MPI_STATUS_SIZE)
      integer:: i,j,kd,kdd,kul,kuh,ku,ierr
      integer:: nav,maxna,maxb,inode,nsd,nrc,nbnew
      real(8):: tcom1,tcom2,vala,valb,valc,rca,rcb,rcc,xi(3)
      logical,external:: bbd
      integer,allocatable:: ibuf(:),ibufr(:)
      real(8),allocatable:: dbuf(:,:),dbufr(:,:)

      allocate(ibuf(nbmax),ibufr(nbmax),dbuf(3,nbmax),dbufr(3,nbmax))

c-----reset the num of "received" boundary atoms
      nbnew=0

c-----calculate the cut-off lengths
      vala=dsqrt(sgm(1,1)**2+sgm(2,1)**2+sgm(3,1)**2)/volume
      valb=dsqrt(sgm(1,2)**2+sgm(2,2)**2+sgm(3,2)**2)/volume
      valc=dsqrt(sgm(1,3)**2+sgm(2,3)**2+sgm(3,3)**2)/volume
      rca=rc*vala
      rcb=rc*valb
      rcc=rc*valc
c      write(6,'(a,3es12.4)') "rca,rcb,rcc=",rca,rcb,rcc
c      write(6,'(i5,3es12.4)') myid_md,rca,rcb,rcc

c-----loop over x, y, & z directions
      do 100 kd=1,3

c-------No. of to-be-copied atoms, LSB(0,)
        do kdd= -1,0
          lsb(0,2*kd+kdd)=0
        enddo
c-------Scan all the residents & copies
        do i=1,natm+nbnew
          xi(1:3)= ra(1:3,i)
c---------For low & high directions
          kul=2*kd-1
          kuh=2*kd
c---------Lower neighbor
          if (bbd(xi(1),xi(2),xi(3),rca,rcb,rcc
     &         ,kul,anxi,anyi,anzi)) then
            lsb(0,kul)=lsb(0,kul)+1
            lsb(lsb(0,kul),kul)=i
          endif
c---------Higher neighbor
          if(bbd(xi(1),xi(2),xi(3),rca,rcb,rcc
     &           ,kuh,anxi,anyi,anzi)) then
            lsb(0,kuh)=lsb(0,kuh)+1
            lsb(lsb(0,kuh),kuh)=i
          endif
        enddo
c-------Error trap
        do kdd= -1,0
          ku=2*kd+kdd
c          write(*,*) " myid_md,ku,lsb(0)=",myid_md,ku,lsb(0,ku)
          call MPI_ALLREDUCE(lsb(0,ku),maxb,1,MPI_INTEGER,MPI_MAX,
     &         mpi_md_world,ierr)
          if (maxb.gt.nbmax) then
            if (myid_md.eq.0) then
              write(*,*)'Buffer or list overflowed at bacopy'
              write(*,*)'LSB(0) NBMAX =',maxb,nbmax
              write(*,*)'KU=',ku
            endif
            call MPI_FINALIZE(ierr)
            stop
          endif
        enddo

c-------To calculate the communication time
        tcom1=MPI_WTIME()

        do 200 kdd= -1,0
          ku=2*kd+kdd
          inode=nn(ku)
          nsd=lsb(0,ku)
          ibuf(1)=nsd
          call mespasi(inode,myparity(kd),ibuf,ibufr,1,1,10
     &         ,mpi_md_world)
          nrc=ibufr(1)
c---------Store the # of received boundary atoms in LSRC
          lsrc(ku)=nrc

c---------Exchange IS
          do i=1,nsd
            ibuf(i)=is(lsb(i,ku))
          enddo
          call mespasi(inode,myparity(kd),ibuf,ibufr,nsd,nrc,20
     &         ,mpi_md_world)
          do i=1,nrc
            is(natm+nbnew+i)=ibufr(i)
          enddo

c---------Exchange ra
          do i=1,nsd
            j= lsb(i,ku)
            dbuf(1:3,i)= ra(1:3,j) -sv(1:3,ku)
          enddo
          call mespasd(inode,myparity(kd),dbuf,dbufr,nsd*3,nrc*3,21
     &         ,mpi_md_world)
          do i=1,nrc
            ra(1:3,natm+nbnew+i)= dbufr(1:3,i)
          enddo

          call MPI_BARRIER(mpi_md_world,ierr)
c---------increase the # of received boundary atoms
          nbnew=nbnew+nrc
c          write(6,'(a,2i8)') " nbnew,nrc=",nbnew,nrc
200     continue

c-------Add the communication time to COMT
        tcom2=MPI_WTIME()
        tcom=tcom+tcom2-tcom1

c-------Error trap
        nav=natm+nbnew
        call MPI_ALLREDUCE(nav,maxna,1,MPI_INTEGER,MPI_MAX
     &       ,mpi_md_world,ierr)
        if (maxna.gt.namax) then
          if (myid_md.eq.0) then
            write(*,*)'NAMAX overflowed at bacopy'
            write(*,*)'N+NB NAMAX = ',maxna,namax
          endif
          call MPI_FINALIZE(ierr)
          stop
        endif

100   continue

c-----num. of received boundary atoms
      nb=nbnew

      deallocate(ibuf,ibufr,dbuf,dbufr)
      end subroutine bacopy
c=======================================================================
      function bbd(xv,yv,zv,rcav,rcbv,rccv,ku,anxi,anyi,anzi)
c-----------------------------------------------------------------------
c  BBD = .true. if the coordinates are in the boundary to neighbor ku
c-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      logical:: bbd
      
      if (ku.eq.1) then
        bbd = xv.lt.rcav
      else if (ku.eq.2) then
        bbd = anxi-rcav.lt.xv
      else if (ku.eq.3) then
        bbd = yv.lt.rcbv
      else if (ku.eq.4) then
        bbd = anyi-rcbv.lt.yv
      else if (ku.eq.5) then
        bbd = zv.lt.rccv
      else if (ku.eq.6) then
        bbd = anzi-rccv.lt.zv
      else
        write(*,*)'BBD call is out of range'
      endif
      return
      end function bbd
c=======================================================================
      function bmv(xv,yv,zv,ku,anxi,anyi,anzi)
c-----------------------------------------------------------------------
c  BMV = .true. if the coordinates should belong to neighbor ku
c------------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      logical bmv

      if (ku.eq.1) then
        bmv = xv.lt.0d0
      else if (ku.eq.2) then
        bmv = anxi.lt.xv
      else if (ku.eq.3) then
        bmv = yv.lt.0d0
      else if (ku.eq.4) then
        bmv = anyi.lt.yv
      else if (ku.eq.5) then
        bmv = zv.lt.0d0
      else if (ku.eq.6) then
        bmv = anzi.lt.zv
      else
        write(*,*)'BMV call is out of range'
      endif
      return
      end function bmv
c=======================================================================
      subroutine mk_lspr(namax,natm,nb,nnmax,is,ra,rc,h,hi
     &     ,anxi,anyi,anzi,lspr,myid_md)
      implicit none
      integer,intent(in):: namax,natm,nb,nnmax,is(namax),myid_md
      real(8),intent(in):: ra(3,namax),rc,anxi,anyi,anzi
     &     ,hi(3,3),h(3,3,0:1)
      integer,intent(out):: lspr(0:nnmax,namax)

      integer:: i,j,k,l,m,n
      integer:: lcx,lcy,lcz,lcx2,lcy2,lcz2,lcyz2,lcxyz2,mx,my,mz
     &     ,kux,kuy,kuz,m1x,m1y,m1z,m1,ic,jc,ierr
      real(8):: rcx,rcy,rcz,rcxi,rcyi,rczi,xi(3),xij(3),rij(3),rij2
      integer,allocatable:: lscl(:),lshd(:)

c-----reset pair list, LSPR
      lspr(0,1:natm+nb)= 0

c-----make a linked cell list, LSCL
      lcx=anxi/dsqrt(hi(1,1)**2+hi(1,2)**2+hi(1,3)**2)/rc
      lcy=anyi/dsqrt(hi(2,1)**2+hi(2,2)**2+hi(2,3)**2)/rc
      lcz=anzi/dsqrt(hi(3,1)**2+hi(3,2)**2+hi(3,3)**2)/rc
      lcx2= lcx +2
      lcy2= lcy +2
      lcz2= lcz +2
      lcyz2=lcy2*lcz2
      lcxyz2=lcx2*lcyz2
      rcx= anxi/lcx
      rcy= anyi/lcy
      rcz= anzi/lcz
      rcxi=1d0/rcx
      rcyi=1d0/rcy
      rczi=1d0/rcz
c-----allocate LSCL & LSHD after obtaining lcxyz2
      allocate(lscl(natm+nb),lshd(lcxyz2))

c-----reset headers
      lshd(1:lcxyz2)= 0

c-----construct a linked-cell list, LSCL, & a header list, LSHD
      do i=1,natm+nb
c-------assign a vector cell index
        mx=(ra(1,i)+rcx)*rcxi
        my=(ra(2,i)+rcy)*rcyi
        mz=(ra(3,i)+rcz)*rczi
c-------classify residents in inner cells even if they are not
        if(i.le.natm) then
          mx= min(max(mx,1),lcx)
          my= min(max(my,1),lcy)
          mz= min(max(mz,1),lcz)
c-------copied atoms are either in inner or surface cells
        else
          mx= min(max(mx,0),lcx+1)
          my= min(max(my,0),lcy+1)
          mz= min(max(mz,0),lcz+1)
        endif
        m= mx*lcyz2 +my*lcz2 +mz +1
        lscl(i)= lshd(m)
c-------the last one goes to the header
        lshd(m)= i
      enddo

c-----make a pair list, LSPR
c-----Scan resident cells
c      do 100 mz=1,lcz
c      do 100 my=1,lcy
c      do 100 mx=1,lcx
c-----Scan all cells, not only resident but also surface cells
      do 100 mz=0,lcz+1
      do 100 my=0,lcy+1
      do 100 mx=0,lcx+1
        m= mx*lcyz2 +my*lcz2 +mz +1
        if(m.gt.lcxyz2) then
          write(6,'(a)') " !!!m.gt.lcxyz2!!!"
          call mpi_finalize(ierr)
          stop
        endif
        if (lshd(m).eq.0) goto 100
        do 102 kuz= -1,1
        do 102 kuy= -1,1
        do 102 kux= -1,1
          m1x= mx +kux
          m1y= my +kuy
          m1z= mz +kuz
          if(m1x.lt.0 .or. m1x.gt.lcx+1
     &         .or. m1y.lt.0 .or. m1y.gt.lcy+1
     &         .or. m1z.lt.0 .or. m1z.gt.lcz+1) goto 102
          m1=m1x*lcyz2 +m1y*lcz2 +m1z +1
          if(m1.gt.lcxyz2) then
            write(6,'(a)') " !!!m1.gt.lcxyz2!!!"
            call mpi_finalize(ierr)
          stop
        endif
          if (lshd(m1).eq.0) goto 102

          i=lshd(m)
1         continue
c          if (natm.lt.i) goto 4

          ic=is(i)
          xi(1:3)= ra(1:3,i)

          j=lshd(m1)

2         continue
          if (j.eq.i) goto 3
          jc=is(j)
          xij(1:3)= ra(1:3,j) -xi(1:3)
          rij(1)= h(1,1,0)*xij(1) +h(1,2,0)*xij(2) +h(1,3,0)*xij(3)
          rij(2)= h(2,1,0)*xij(1) +h(2,2,0)*xij(2) +h(2,3,0)*xij(3)
          rij(3)= h(3,1,0)*xij(1) +h(3,2,0)*xij(2) +h(3,3,0)*xij(3)
          rij2= rij(1)**2 +rij(2)**2 +rij(3)**2

c          if (rij2.lt.rcij2(ic,jc)) then
          if(rij2.lt.rc**2) then
            lspr(0,i)=lspr(0,i)+1
            lspr(lspr(0,i),i)=j
          endif

c---------Continue until j= 0
3         j=lscl(j)
          if (j.gt.0) goto 2

c---------Continue until i= 0
4         i=lscl(i)
          if (i.gt.0) goto 1

102     continue
100   continue

      deallocate(lscl,lshd)
      end subroutine mk_lspr
c=======================================================================
      subroutine mespasi(inode,parity,ibufs,ibufr,nsd,nrc,tag,
     &     mpi_md_world)
c-----------------------------------------------------------------------
c     Integer message passing.  
c-----------------------------------------------------------------------
      include 'mpif.h'
      integer,intent(in):: inode,parity,nsd,nrc,tag
      integer,intent(in):: ibufs(nsd)
      integer,intent(out):: ibufr(nrc)
c-----locals
      integer status(MPI_STATUS_SIZE),ierr

c-----Even: send & recv
      if (parity.eq.0) then
        call MPI_SEND(ibufs,nsd,MPI_INTEGER,inode,tag,
     &       mpi_md_world,ierr)
        call MPI_RECV(ibufr,nrc,MPI_INTEGER,MPI_ANY_SOURCE,tag,
     &       mpi_md_world,status,ierr)
c-----Odd: recv & send
      else if (parity.eq.1) then
        call MPI_RECV(ibufr,nrc,MPI_INTEGER,MPI_ANY_SOURCE,tag,
     &       mpi_md_world,status,ierr)
        call MPI_SEND(ibufs,nsd,MPI_INTEGER,inode,tag,
     &       mpi_md_world,ierr)
c-----Exchange information with myself
      else
        do i=1,nrc
          ibufr(i)=ibufs(i)
        enddo
      endif
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mespasd(inode,parity,bufs,bufr,nsd,nrc,tag,
     & mpi_md_world)
c-----------------------------------------------------------------------
c     Real*8 message passing.
c-----------------------------------------------------------------------
      include 'mpif.h'
      integer,intent(in):: inode,parity,nsd,nrc,tag
      real(8),intent(in):: bufs(nsd)
      real(8),intent(out):: bufr(nrc)
      integer:: status(MPI_STATUS_SIZE),ierr

c-----Even: send & recv
      if (parity.eq.0) then
        call MPI_SEND(bufs,nsd,MPI_DOUBLE_PRECISION,inode,tag,
     &       mpi_md_world,ierr) 
        call MPI_RECV(bufr,nrc,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,tag,
     &       mpi_md_world,status,ierr) 
c-----Odd: recv & send
      else if (parity.eq.1) then
        call MPI_RECV(bufr,nrc,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,tag,
     &       mpi_md_world,status,ierr) 
        call MPI_SEND(bufs,nsd,MPI_DOUBLE_PRECISION,inode,tag,
     &       mpi_md_world,ierr)
c-----Exchange information with myself
      else
        do i=1,nrc
          bufr(i)=bufs(i)
        enddo
      endif
      return
      end
c=======================================================================
      subroutine get_cnum(chr,num)
      implicit none
      integer,intent(in):: num
      character,intent(inout):: chr*3

      chr="000"
      if(num.lt.10) then
        write(chr(3:3),'(i1)') num
      elseif(num.lt.100) then
        write(chr(2:3),'(i2)') num
      else
        write(chr(1:3),'(i3)') num
      endif
      
      end subroutine get_cnum
c=======================================================================
      subroutine cna(namax,natm0,nbndr,nnmax,lspr,ra0,va0,tag,ncna
     &     ,lscna)
c-----------------------------------------------------------------------
c Reduce natm0 to natm by using Common Neighbor Analysis
c-----------------------------------------------------------------------
      implicit none
      integer,intent(in):: namax,natm0,nbndr,nnmax,lspr(0:nnmax,namax)
      real(8),intent(in):: ra0(3,namax),va0(3,namax)
      integer,intent(inout):: ncna,lscna(namax)
      real(8),intent(inout):: tag(namax)

      integer:: i,j,l,m,n,ii,iii,ni,jj,nj,il,jl,n1,n2,iil,nn1,im,iim
     &     ,ib1,ib2,iib1,iib2,n421,n422
      integer,allocatable:: icommon(:),ibond(:,:),nb(:),idc(:,:,:)
      integer,parameter:: lmax= 12
      integer,parameter:: mmax= 12

cc-----for serial code
c      allocate(icommon(lmax),ibond(2,mmax),nb(mmax)
c     &     ,idc(3,nnmax,natm0))
c-----for parallel code
      allocate(icommon(lmax),ibond(2,mmax),nb(mmax)
     &     ,idc(3,nnmax,natm0+nbndr))

c-----init three indices
      idc(1:3,1:nnmax,1:natm0+nbndr)= 0

c-----for each atom-i, store three indices (LMN)
      do i=1,natm0+nbndr

c-------for each 1st n.n.
        do ii=1,lspr(0,i)
          j=lspr(ii,i)
c---------j>i only
          if(j.le.i) cycle

c---------count num of common neighbors: L
          l= 0
          icommon(1:lmax)= 0
          do iii=1,lspr(0,i)
            ni=lspr(iii,i)
            if(ni.eq.j) cycle
            do jj=1,lspr(0,j)
              nj=lspr(jj,j)
              if(nj.eq.ni) then
                l=l+1
                if(l.gt.lmax) stop " l.gt.lmax!!!"
                icommon(l)= ni
                exit
              endif
            enddo
c---------end of counting L
          enddo
          
c---------count num of bonds between common neighbors: M
          m= 0
          ibond(1:2,1:mmax)= 0
c---------for each common neighbor-n1
          do il=1,l
            n1=icommon(il)
c-----------for common neighbor-n2 which must be larger than n1
            do jl=1,l
              n2=icommon(jl)
              if(n2.le.n1) cycle
c-------------scan 1st n.n. of n1
              do iil=1,lspr(0,n1)
                nn1=lspr(iil,n1)
                if(nn1.eq.n2) then
                  m=m+1
                  if(m.gt.mmax) stop " m.gt.mmax!!"
                  ibond(1:2,m)= (/ n1,n2 /)
                  exit
                endif
              enddo
            enddo
          enddo

c---------count max num of continuous bonds: N
          nb(1:mmax)= 1
c---------this does not distinguish star and chain connections
          do im=1,m-1
            ib1= ibond(1,im)
            ib2= ibond(2,im)
            do iim=im+1,m
              iib1= ibond(1,iim)
              iib2= ibond(2,iim)
c-------------if two bonds are connected, up nb
              if(iib1.eq.ib1 .or. iib2.eq.ib1
     &             .or. iib1.eq.ib2 .or. iib2.eq.ib2) then
                nb(im)=nb(im) +1
                nb(iim)=nb(iim) +1
              endif
            enddo
          enddo
c---------maximum nb
          n= 0
          do im=1,m
            n= max(nb(im),n)
          enddo

c---------store (LMN) to i
          idc(1:3,ii,i)= (/ l,m,n /)
c---------store (LMN) to j, too
          do jj=1,lspr(0,j)
            if(lspr(jj,j).eq.i) then
              idc(1:3,jj,j)= (/ l,m,n /)
              exit
            endif
          enddo
c-------end of 1st n.n. of i: j
        enddo
c-----end of atom-i
      enddo

c-----reduce atoms
      ncna= 0
      do i=1,natm0
        n421= 0
        n422= 0
        do ii=1,lspr(0,i)
          l=idc(1,ii,i)
          m=idc(2,ii,i)
          n=idc(3,ii,i)
          if(l.eq.4 .and. m.eq.2 .and. n.eq.1 ) n421=n421 +1
          if(l.eq.4 .and. m.eq.2 .and. n.eq.2 ) n422=n422 +1
        enddo
c-------if perfect FCC structure, reduce the atom-i, so skip
        if(n421.eq.12 .and. n422.eq.0) cycle
c        if(mod(i,10).ne.0) cycle
c-------if not FCC
        ncna=ncna +1
        lscna(ncna)= i
c-------if HCP structure, itype(i)= 2
        if(n421.eq.6 .and. n422.eq.6) then
          tag(i)= (tag(i)-int(tag(i))) +2d0
c          itype(i)= 2
c          itype(ncna)= 2
c-------otherwise, itype(i)= 1
        else
          tag(i)= (tag(i)-int(tag(i))) +1d0
c          itype(i)= 1
c          itype(ncna)= 1
        endif
cc-------copy to ra,va
c        ra(1:3,ncna)= ra0(1:3,i)
c        va(1:3,ncna)= va0(1:3,i)
      enddo

      deallocate(icommon,ibond,nb,idc)
      end subroutine cna
c=======================================================================
      subroutine z_reduce(namax,natm,ra,ncna,lscna)
c-----------------------------------------------------------------------
c  Reduce num. of atoms by eliminating according to z-coordinate
c-----------------------------------------------------------------------
      implicit none
      integer,intent(in):: namax,natm
      real(8),intent(in):: ra(3,namax)
      integer,intent(inout):: ncna,lscna(namax)

      integer:: i,n

      n=0
      do i=1,natm
        if(ra(1,i).gt.0.05d0 .and. ra(1,i).lt.0.2d0) then
          n=n+1
          lscna(n)= i
        endif
      enddo
      ncna= n

      end subroutine z_reduce
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make pmd2cna"
c     End:
