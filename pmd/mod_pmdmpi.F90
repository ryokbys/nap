module pmdmpi
!-----------------------------------------------------------------------
!                     Last modified: <2017-10-19 11:01:44 Ryo KOBAYASHI>
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

contains
  function get_factor(n)
    implicit none
    integer,intent(in):: n
    integer:: imax,i
    integer:: get_factor

    if( n.eq.1 ) then
      get_factor = 1
      return
    endif

    imax = int(sqrt(dble(n)))
    if( mod(n,2).eq.0 ) then
      get_factor = 2
      return
    else if( mod(n,3).eq.0 ) then
      get_factor = 3
      return
    else
      i = 6
      do while( i.le.imax )
        if( mod(n,i-1).eq.0 ) then
          get_factor = i-1
          return
        else if( mod(n,i+1).eq.0 ) then
          get_factor = i+1
          return
        endif
        i = i + 6
      enddo
    endif
  end function get_factor
!=======================================================================
  subroutine assign_num_nodes(al1,al2,al3)
    implicit none
    real(8),intent(in):: al1,al2,al3

    integer:: nfac,f,n
    integer:: maxfac = 100
    integer,allocatable:: factors(:)

    if( .not.allocated(factors)) allocate(factors(maxfac))

    n = nodes_md
    nfac = 1
    factors(nfac) = 1
    do while( .true. )
      f = get_factor(n)
      nfac = nfac + 1
      factors(nfac) = f
      if( f.eq.n ) exit
      n = n/f
    enddo

    do while( nfac.lt.3 )
      nfac = nfac + 1
      factors(nfac) = 1
    enddo

!.....TODO: assign nodes taking the ratio al1:al2:al3 into account...
    
  end subroutine assign_num_nodes
end module pmdmpi
