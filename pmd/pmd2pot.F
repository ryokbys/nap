      program pmd2pot
c-----------------------------------------------------------------------
c  Reduce atoms from pmd file by evaluating potential energy per atom
c    - for any crystal structures
c-----------------------------------------------------------------------
c
c  << Input files >>
c    pmd.in: simulation config file
c    pmd###-???: MD-coordination files,
c                ### is the node number and ??? is the frame number
c
c  << Output files >>
c    pot???: pot files
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
c-----criteria or threshold for potential energy per atom
      real(8),parameter:: epic  = -6.8d10

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
      real(8):: ra(3,namax),va(3,namax),op(2,namax)
c-----buffer data
      integer:: ibuf(maxbuf)
      real(8):: rbuf(8,maxbuf)
c-----pot data
      integer:: npot,npottot
      integer:: ipot(namax),lspot(namax)
      integer,allocatable:: lpot(:),nrecvs(:),idisps(:)
      real(8):: xpot(8,namax)
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
      open(10,file="pmd.in",status="old")
      read(10,*) nx,ny,nz
      read(10,*) tmp,itmp
      read(10,*) tmp,tmp
      read(10,*) tmp,tmp
      read(10,*) itmp
      read(10,*) ifpmd,npmd
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
      allocate(lpot(nxyz),nrecvs(0:nxyz-1),idisps(0:nxyz-1))

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
          read(20,*) itype(i),ra(1:3,i),va(1:3,i),op(1:2,i)
          ra(1:3,i)= ra(1:3,i) -sorg(1:3)
        enddo
        close(20)

c=======extract atom whose epi > epic
        npot= 0
        do i=1,natm
          if(op(2,i).gt.epic) then
            npot= npot +1
            lspot(npot)= i
          endif
        enddo

c-------gather num. of pot atoms
        call mpi_reduce(npot,npottot,1,MPI_INTEGER,MPI_SUM
     &       ,0,mpi_md_world,ierr)
        if(myid_md.eq.0) write(6,'(a,i8)') " npottot=",npottot
c-------error trap
        if(npottot.gt.namax) then
          write(*,*) " !!!npottot.gt.namax!!!"
          call mpi_finalize(ierr)
          stop 
        endif
        call mpi_allgather(npot,1,MPI_INTEGER
     &       ,lpot,1,MPI_INTEGER,mpi_md_world,ierr)
        idisps(0:nodes_md-1)=0
        do i=1,nodes_md-1
          idisps(i)= idisps(i-1) +lpot(i)
        enddo
        do i=0,nodes_md-1
          nrecvs(i)= lpot(i+1)
        enddo
c-------gather
        do i=1,npot
          itmp= lspot(i)
          rbuf(1:3,i)= ra(1:3,itmp) +sorg(1:3)
          rbuf(4:6,i)= va(1:3,itmp)
          rbuf(7:8,i)= op(1:2,itmp)
          ibuf(i)= itype(itmp)
        enddo
        if(myid_md.eq.0)write(6,'(a)') " 1st gatherv..."
        call mpi_gatherv(ibuf,npot,MPI_INTEGER
     &       ,ipot,nrecvs,idisps,MPI_INTEGER
     &       ,0,mpi_md_world,ierr)
        do i=0,nodes_md-1
          idisps(i)= idisps(i) *8
          nrecvs(i)= nrecvs(i) *8
        enddo
        if(myid_md.eq.0)write(6,'(a)') " 2nd gatherv..."
        call mpi_gatherv(rbuf,8*npot,MPI_DOUBLE_PRECISION
     &       ,xpot,nrecvs,idisps,MPI_DOUBLE_PRECISION
     &       ,0,mpi_md_world,ierr)

c-------output as pot file
        if(myid_md.eq.0) then
          write(6,'(a)') " writing file pot"//cnum//"..."
          open(21,file="pot"//cnum,status="replace")
          write(21,'(es15.7)') hunit
          write(21,'(3es12.4)') (((h(ia,ib,l)/hunit,ia=1,3)
     &         ,ib=1,3),l=0,1)
          write(21,'(i10)') npottot
          do i=1,npottot
            write(21,'(i4,8es13.4e3)') ipot(i),xpot(1:8,i)
          enddo
          close(21)
        endif
      enddo

      call mpi_finalize(ierr)

      end program pmd2pot
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
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make pmd2pot"
c     End:
