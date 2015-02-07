      subroutine copy_rho_ba(tcom,namax,natm,nb,nbmax,lsb
     &     ,lsrc,myparity,nn,sv,mpi_md_world,rho)
c-----------------------------------------------------------------------
c     Exchanges boundary-atom data among neighbor nodes
c-----------------------------------------------------------------------
      implicit none
      include "mpif.h"
      integer:: status(MPI_STATUS_SIZE)
c-----in
      integer,intent(in):: namax,natm,nb,nbmax,mpi_md_world
      integer,intent(in):: lsb(0:nbmax,6),lsrc(6),myparity(3),nn(6)
      real(8),intent(in):: sv(3,6)
c-----out
      real(8),intent(inout):: rho(natm+nb),tcom

c-----locals
      integer:: i,j,k,l,m,n,kd,kdd,ku,inode,nsd,nsd3,nrc,nrc3,nbnew,ierr
      real(8):: tcom1,tcom2
      logical,save:: l1st=.true.
      real(8),allocatable,save:: dbuf(:),dbufr(:)

      if( l1st ) then
        allocate(dbuf(nbmax),dbufr(nbmax))
        l1st=.false.
      endif

      nbnew= 0

c-----loop over z, y, & x directions
      do kd=1,3
        tcom1= mpi_wtime()
        do kdd=-1,0
          ku= 2*kd +kdd
          inode= nn(ku)
c---------num. of to-be-sent particles
          nsd= lsb(0,ku)
c---------num. of to-be-recieved particles
          nrc= lsrc(ku)

c---------exchange x
          do i=1,nsd
            j=lsb(i,ku)
            dbuf(i)= rho(j)
          enddo
          call mespasd(inode,myparity(kd),dbuf,dbufr,nsd,nrc,21
     &         ,mpi_md_world)
          do i=1,nrc
            rho(natm+nbnew+i)= dbufr(i)
          enddo

c---------mpi barrier
          call mpi_barrier(mpi_md_world,ierr)
c---------accumulate num. of boundary particles
c          write(6,'(a,2i8)') "nbnew,nrc=",nbnew,nrc
          nbnew=nbnew +nrc
        enddo
        tcom2= mpi_wtime()
        tcom= tcom +tcom2-tcom1
      enddo

      if(nbnew.ne.nb) then
        write(6,'(a,2i8)') "nbnew,(natm+nb)=",nbnew,natm+nb
        stop "error: nbnew.ne.(natm+nb)!!"
      endif
      
      end subroutine copy_rho_ba
c=======================================================================
      subroutine copy_strs_ba(tcom,namax,natm,nb,nbmax,lsb
     &     ,lsrc,myparity,nn,sv,mpi_md_world,strs)
c-----------------------------------------------------------------------
c  Exchanges boundary-atom data among neighbor nodes
c-----------------------------------------------------------------------
      implicit none
      include "mpif.h"
      integer:: status(MPI_STATUS_SIZE)
c-----in
      integer,intent(in):: namax,natm,nb,nbmax,mpi_md_world
      integer,intent(in):: lsb(0:nbmax,6),lsrc(6),myparity(3),nn(6)
      real(8),intent(in):: sv(3,6)
c-----out
      real(8),intent(inout):: strs(9,natm+nb),tcom

c-----locals
      integer:: i,j,k,l,m,n,kd,kdd,ku,inode,nsd,nrc,nbnew,ierr
      real(8):: tcom1,tcom2
      
      logical,save:: l1st=.true.
      real(8),save,allocatable:: dbuf(:,:),dbufr(:,:)

      if( l1st ) then
        allocate(dbuf(9,nbmax),dbufr(9,nbmax))
        l1st=.false.
      endif

      nbnew= 0

c-----loop over z, y, & x directions
      do kd=1,3
        tcom1= mpi_wtime()
        do kdd=-1,0
          ku= 2*kd +kdd
          inode= nn(ku)
c---------num. of to-be-sent particles
          nsd= lsb(0,ku)
c---------num. of to-be-recieved particles
          nrc= lsrc(ku)

c---------exchange strs
          do i=1,nsd
            j=lsb(i,ku)
            dbuf(1:9,i)= strs(1:9,j)
          enddo
          call mespasd(inode,myparity(kd),dbuf,dbufr,9*nsd,9*nrc,21
     &         ,mpi_md_world)
          do i=1,nrc
            strs(1:9,natm+nbnew+i)= dbufr(1:9,i)
          enddo

c---------mpi barrier
          call mpi_barrier(mpi_md_world,ierr)
c---------accumulate num. of boundary particles
c          write(6,'(a,2i8)') "nbnew,nrc=",nbnew,nrc
          nbnew=nbnew +nrc
        enddo
        tcom2= mpi_wtime()
        tcom= tcom +tcom2-tcom1
      enddo

      if(nbnew.ne.nb) then
        write(6,'(a,2i8)') "nbnew,(natm+nb)=",nbnew,natm+nb
        stop "error: nbnew.ne.(natm+nb)!!"
      endif
      
      end subroutine copy_strs_ba
c=======================================================================
      subroutine copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb
     &     ,lsrc,myparity,nn,mpi_md_world,x,ndim)
c-----------------------------------------------------------------------
c     Send-back & receive reaction on cached-atoms
c-----------------------------------------------------------------------
      implicit none
      include "mpif.h"
      integer,intent(in):: namax,natm,nbmax,nb,mpi_md_world,ndim
      integer,intent(in):: lsb(0:nbmax,6),lsrc(6),myparity(3),nn(6)
      real(8),intent(inout):: x(ndim,namax),tcom

      integer:: status(MPI_STATUS_SIZE)
      integer:: i,j,k,l,m,n,kd,kdd,ku,kuc,ibkwd,nsd,nsd3,nrc,nrc3,nsdbk
     &     ,ierr,natmx
      real(8):: tcom1,tcom2
      real(8),save,allocatable:: dbuf(:,:),dbufr(:,:)
      logical,save:: l1st=.true.
      integer,save:: mdim

      if( l1st ) then
        mdim= ndim
        allocate(dbuf(mdim,nbmax),dbufr(mdim,nbmax))
        l1st=.false.
      endif

      if( ndim.gt.mdim ) then
        deallocate(dbuf,dbufr)
        mdim= ndim
        allocate(dbuf(mdim,nbmax),dbufr(mdim,nbmax))
      endif

c-----natmx
      natmx= natm +nb

c-----num. of sent-back reactions
      nsdbk= 0

c-----send-back reactions in the reverse orer, z, y, & x
      do kd=3,1,-1

c-------To calculate the communication time
        tcom1=MPI_WTIME()

c-------higher & lower directions
        do kdd=0,-1,-1
          ku= 2*kd +kdd
          if(mod(ku,2).eq.0) then
            kuc= ku-1
          else
            kuc= ku+1
          endif
          ibkwd= nn(kuc)
c---------num. of to-be-sent particles
          nsd= lsrc(ku)
c          nsd3= ndim*nsd
          nsd3= mdim*nsd
c---------num. of to-be-recieved particles
          nrc= lsb(0,ku)
c          nrc3= ndim*nrc
          nrc3= mdim*nrc

c---------to-be-sent-back particles
          do i=1,nsd
            dbuf(1:ndim,i)= x(1:ndim,natmx-nsdbk-nsd+i)
          enddo
          call mespasd(ibkwd,myparity(kd),dbuf,dbufr,nsd3,nrc3,500
     &         ,mpi_md_world)
          do k=1,nrc
            i=lsb(k,ku)
            x(1:ndim,i)= x(1:ndim,i) +dbufr(1:ndim,k)
          enddo

c---------mpi barrier
          call mpi_barrier(mpi_md_world,ierr)
c---------accumulate num. of already sent-back-particles
          nsdbk=nsdbk +nsd
        enddo

c-------Add the communication time to COMT
        tcom2=MPI_WTIME()
        tcom=tcom+tcom2-tcom1

      enddo

c-----check
      if(nsdbk.ne.nb) then
        write(6,'(a,2i8)') "nsdbk,nb=",nsdbk,nb
        stop "error: nsdbk.ne.nb!!"
      endif

c      deallocate(dbuf,dbufr)
      end subroutine copy_dba_bk
c=======================================================================
      function hvsd(x)
c
c  Heaviside's stepwise function
c
      implicit none
      real(8),intent(in):: x
      real(8):: hvsd

      hvsd= 0d0
      if( x.ge.0 ) then
        hvsd= 1d0
        return
      endif
      return 

      end function hvsd
c=======================================================================
