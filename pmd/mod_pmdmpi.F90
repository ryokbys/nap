module pmdmpi
!-----------------------------------------------------------------------
!                     Last modified: <2021-02-27 10:12:32 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
! Module that includes variables and parameters used for parallel
! computation with mpi for spatial decomposition MD simulation.
!-----------------------------------------------------------------------
  implicit none
  save

contains
  function get_factor(n) result(ifac)
    implicit none
    integer,intent(in):: n
    integer:: imax,i
    integer:: ifac

    ifac = 1
    if( n.eq.1 ) return

    imax = int(sqrt(dble(n)))
    if( mod(n,2).eq.0 ) then
      ifac = 2
      return
    else if( mod(n,3).eq.0 ) then
      ifac = 3
      return
    else
      i = 6
      do while( i.le.imax )
        if( mod(n,i-1).eq.0 ) then
          ifac = i-1
          return
        else if( mod(n,i+1).eq.0 ) then
          ifac = i+1
          return
        endif
        i = i + 6
      enddo
    endif
  end function get_factor
!=======================================================================
  subroutine assign_num_nodes(al1,al2,al3)
    use pmdvars,only: nodes_md
    real(8),intent(in):: al1,al2,al3

    integer:: nfac,f,n
    integer:: maxfac = 100
    integer,allocatable,save:: factors(:)

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
!=======================================================================
  subroutine nid2xyz(id,ix,iy,iz)
!
!  Convert continuous node-id and cell position (ix,iy,iz)
!  Note: id,ix,iy,iz start from 0
!
    use pmdvars,only: nx,ny,nz
    integer,intent(in):: id
    integer,intent(out):: ix,iy,iz

    ix = id/(ny*nz)
    iy = mod(id/nz,ny)
    iz = mod(id,nz)
    return
  end subroutine nid2xyz
!=======================================================================
  subroutine xyz2nid(ix,iy,iz,id)
!
!     Convert cell position (kx,ky,kz) to continuous node-id
!
    use pmdvars,only: nx,ny,nz
    integer,intent(in):: ix,iy,iz
    integer,intent(out):: id

    id = ix*(ny*nz) +iy*nz +iz
    return
  end subroutine xyz2nid
!=======================================================================
end module pmdmpi
