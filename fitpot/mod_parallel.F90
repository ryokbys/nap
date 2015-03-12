  module parallel
    include 'mpif.h'
    integer,save:: myid,mpi_world,ierr
    integer,save:: nnode,isid0,isid1
    integer,allocatable,save:: nspn(:),ispn(:)
    integer,save:: mynsmpl
  end module parallel
