module pmdmpi
!-----------------------------------------------------------------------
!                     Last modified: <2016-11-02 18:32:36 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
! Module that includes variables and parameters used for parallel
! computation with mpi for spatial decomposition MD simulation.
!-----------------------------------------------------------------------
  implicit none
  save

  integer:: myid_md,nodes_md,mpi_md_world
  integer:: nx = 1
  integer:: ny = 1
  integer:: nz = 1

end module pmdmpi
