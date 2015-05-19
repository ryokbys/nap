  module parallel
    include 'mpif.h'
    integer,save:: myid,mpi_world,ierr
    integer,save:: nnode
    integer,save:: isid0_trn,isid1_trn
    integer,save:: isid0_tst,isid1_tst
    integer,allocatable,save:: nspn_trn(:),ispn_trn(:)
    integer,allocatable,save:: nspn_tst(:),ispn_tst(:)
    integer,save:: mynsmpl_trn,maxmynsmpl_trn
    integer,save:: mynsmpl_tst,maxmynsmpl_tst
  end module parallel
