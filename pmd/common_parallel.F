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
