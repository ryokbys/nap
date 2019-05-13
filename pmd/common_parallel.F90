subroutine mespasi(inode,parity,ibufs,ibufr,nsd,nrc,tag, &
     mpi_md_world)
!-----------------------------------------------------------------------
!     Integer message passing.  
!-----------------------------------------------------------------------
  include 'mpif.h'
  integer,intent(in):: inode,parity,nsd,nrc,tag
  integer,intent(in):: ibufs(nsd)
  integer,intent(out):: ibufr(nrc)
!-----locals
  integer status(MPI_STATUS_SIZE),ierr

!-----Even: send & recv
  if (parity.eq.0) then
    call mpi_send(ibufs,nsd,MPI_INTEGER,inode,tag, &
         mpi_md_world,ierr)
    call mpi_recv(ibufr,nrc,MPI_INTEGER,MPI_ANY_SOURCE,tag, &
         mpi_md_world,status,ierr)
!-----Odd: recv & send
  else if (parity.eq.1) then
    call mpi_recv(ibufr,nrc,MPI_INTEGER,MPI_ANY_SOURCE,tag, &
         mpi_md_world,status,ierr)
    call mpi_send(ibufs,nsd,MPI_INTEGER,inode,tag, &
         mpi_md_world,ierr)
!-----Exchange information with myself
  else
    do i=1,nrc
      ibufr(i)=ibufs(i)
    enddo
  endif
  return
end subroutine mespasi
!=======================================================================
subroutine mespasd(inode,parity,bufs,bufr,nsd,nrc,tag, &
     mpi_md_world)
!-----------------------------------------------------------------------
!     Real*8 message passing.
!-----------------------------------------------------------------------
  include 'mpif.h'
  integer,intent(in):: inode,parity,nsd,nrc,tag
  real(8),intent(in):: bufs(nsd)
  real(8),intent(out):: bufr(nrc)
  integer:: status(MPI_STATUS_SIZE),ierr

!-----Even: send & recv
  if (parity.eq.0) then
    call mpi_send(bufs,nsd,MPI_DOUBLE_PRECISION,inode,tag, &
         mpi_md_world,ierr) 
    call mpi_recv(bufr,nrc,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,tag, &
         mpi_md_world,status,ierr) 
!-----Odd: recv & send
  else if (parity.eq.1) then
    call mpi_recv(bufr,nrc,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,tag, &
         mpi_md_world,status,ierr) 
    call mpi_send(bufs,nsd,MPI_DOUBLE_PRECISION,inode,tag, &
         mpi_md_world,ierr)
!-----Exchange information with myself
  else
    do i=1,nrc
      bufr(i)=bufs(i)
    enddo
  endif
  return
end subroutine mespasd
!=======================================================================
