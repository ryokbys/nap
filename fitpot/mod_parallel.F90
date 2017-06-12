  module parallel
    include 'mpif.h'
    integer,save:: myid,mpi_world,ierr
    integer,save:: nnode
    integer,save:: isid0,isid1
    integer,allocatable,save:: nspn(:),ispn(:)
    integer,save:: mynsmpl,myntrn,myntst,maxmynsmpl,maxmyntrn

    integer:: mpi_comm_pmd, myid_pmd, nnode_pmd
  end module parallel
