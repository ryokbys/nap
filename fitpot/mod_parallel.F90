module parallel
  save
  include 'mpif.h'
  integer:: myid,mpi_world,ierr
  integer:: nnode
  integer:: isid0,isid1
  integer,allocatable:: nspn(:),ispn(:)
  integer:: mynsmpl,myntrn,myntst,maxmynsmpl,maxmyntrn

  integer:: mpi_comm_pmd, myid_pmd, nnode_pmd
end module parallel
