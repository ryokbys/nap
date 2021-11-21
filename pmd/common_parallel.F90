subroutine mespasi(inode,ipar,ibufs,ibufr,nsd,nrc,tag,mpi_world)
!-----------------------------------------------------------------------
!     Integer message passing.  
!-----------------------------------------------------------------------
  implicit none 
  include 'mpif.h'
  integer,intent(in):: inode,ipar,nsd,nrc,tag,mpi_world
  integer,intent(in):: ibufs(nsd)
  integer,intent(out):: ibufr(nrc)
!-----locals
  integer status(MPI_STATUS_SIZE),ierr,i

!-----Even: send & recv
  if (ipar.eq.0) then
    call mpi_send(ibufs,nsd,MPI_INTEGER,inode,tag,mpi_world,ierr)
    call mpi_recv(ibufr,nrc,MPI_INTEGER,MPI_ANY_SOURCE,tag,mpi_world,status,ierr)
!-----Odd: recv & send
  else if (ipar.eq.1) then
    call mpi_recv(ibufr,nrc,MPI_INTEGER,MPI_ANY_SOURCE,tag,mpi_world,status,ierr)
    call mpi_send(ibufs,nsd,MPI_INTEGER,inode,tag,mpi_world,ierr)
!-----Exchange information with myself
  else
    do i=1,nrc
      ibufr(i)=ibufs(i)
    enddo
  endif
  return
end subroutine mespasi
!=======================================================================
subroutine mespasd(inode,ipar,bufs,bufr,nsd,nrc,tag,mpi_world)
!-----------------------------------------------------------------------
! Real(8) message passing.
!-----------------------------------------------------------------------
  implicit none
  include 'mpif.h'
  integer,intent(in):: inode,ipar,nsd,nrc,tag,mpi_world
  real(8),intent(in):: bufs(nsd)
  real(8),intent(out):: bufr(nrc)
  integer:: status(MPI_STATUS_SIZE),ierr,i

!-----Even: send & recv
  if (ipar.eq.0) then
    call mpi_send(bufs,nsd,mpi_real8,inode,tag,mpi_world,ierr) 
    call mpi_recv(bufr,nrc,mpi_real8,MPI_ANY_SOURCE,tag, &
         mpi_world,status,ierr) 
!-----Odd: recv & send
  else if (ipar.eq.1) then
    call mpi_recv(bufr,nrc,mpi_real8,MPI_ANY_SOURCE,tag, &
         mpi_world,status,ierr) 
    call mpi_send(bufs,nsd,mpi_real8,inode,tag,mpi_world,ierr)
!-----Exchange information with myself
  else
    do i=1,nrc
      bufr(i)=bufs(i)
    enddo
  endif
  return
end subroutine mespasd
!=======================================================================
