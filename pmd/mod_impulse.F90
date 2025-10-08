module impulse
!
!  Module for evaluation of probability density.
!
  use memory,only: accum_mem
  use pmdvars,only: nspmax
  use util,only: itotOf
  implicit none
  include "mpif.h"
  include "./const.h"
  save

  character(len=128):: paramsdir = '.'
  character(len=128):: cfout_impls = 'out.impulse'
  integer,parameter:: io_impls = 69

  logical:: l_impls = .false.  ! Toggle for impulse analysis
  logical:: l_macro_impls = .false.  ! For checking macro ON/OFF
  integer:: itot_impls = 0  ! Atom index that must be given in in.pmd
  integer:: ia_impls  ! Atom index in local ordering
  real(8):: tau_impls(3)  ! Hopping direction
  real(8):: ftau(nspmax), ftaul(nspmax)  ! sum_j F_{ji}(t).tau
  real(8):: ptau, ptaul  ! p_i(t).tau

contains
!=======================================================================
  subroutine init_impulse(myid)
    use pmdvars,only: specorder,nsp
    integer,intent(in):: myid
    real(8):: dtau
    integer:: is

    if( .not. l_macro_impls ) then
      if( myid == 0 ) print *,'ERROR@init_impulse: ' &
           //'IMPULSE macro not set.\n' &
           //'You need to add -DIMPULSE in pmd makefile as,\n' &
           //'  CPPFLAGS= -DIMPULSE '
      stop
    endif

    if( itot_impls < 1 ) then
      if( myid == 0 ) print *,'ERROR@init_impulse: itot_impls<1,'&
           //' which should not happen!'
      stop
    endif

!.....normalize tau vector
    dtau = dsqrt(tau_impls(1)**2 +tau_impls(2)**2 +tau_impls(3)**2)
    tau_impls(:) = tau_impls(:)/dtau

!.....Open file
    if( myid == 0 ) then
      open(io_impls, file=trim(cfout_impls), status='replace')
      write(io_impls,'(a,9a)') "# istp, simtime, ptau,",&
           (specorder(is), is=1,nsp)
    endif
    return
  end subroutine init_impulse
!=======================================================================
  subroutine bcast_impulse(myid, mpi_world)
    integer,intent(in):: myid, mpi_world
    
    integer:: ierr

    call mpi_bcast(l_impls, 1, mpi_logical, 0, mpi_world, ierr)
    call mpi_bcast(itot_impls, 1, mpi_integer, 0, mpi_world, ierr)
    call mpi_bcast(tau_impls, 3, mpi_real8, 0, mpi_world, ierr)

  end subroutine bcast_impulse
!=======================================================================
  subroutine set_ia_impls(natm,tag)
    integer,intent(in):: natm
    real(8),intent(in):: tag(natm)

    integer:: ia
    
    ia_impls = 0
    do ia=1,natm
      if( itotOf(tag(ia)) == itot_impls ) then
        ia_impls = ia
        exit
      endif
    enddo

  end subroutine set_ia_impls
!=======================================================================
  subroutine comp_ptau(natm,tag,va,h)
    use pmdvars,only: am,fa2v,nsp
    integer,intent(in):: natm
    real(8),intent(in):: tag(natm),va(3,natm),h(3,3)

    integer:: is, ixyz
    real(8):: vi(3)

    if( ia_impls > 0 ) then  ! only the node contains itot_impls does this
      is = int(tag(ia_impls))
      ptaul = 0d0
      do ixyz=1,3
        ptaul = ptaul +va(ixyz,ia_impls)*tau_impls(ixyz)
      enddo
      ptaul = ptaul*am(is)
    endif
    return
  end subroutine comp_ptau
!=======================================================================
  subroutine write_impulse(istp,simtime,myid,mpi_world,iprint)
!
!  Write out impulse information into the file.
!  This will be called after velocity update is done.
!
    use pmdvars,only: fa2v,am,dt,nsp
    integer,intent(in):: istp, myid, mpi_world, iprint
    real(8),intent(in):: simtime
    integer:: ierr, is

!.....Scale forces
    do is=1,nsp
      ftaul(is) = ftaul(is) *fa2v(is)*dt*am(is)
    enddo
    
!.....Gather information from child nodes
    ptau= 0d0
    call mpi_reduce(ptaul,ptau,1,mpi_real8,mpi_sum, &
         0,mpi_world,ierr)
    ftau(:) = 0d0
    call mpi_reduce(ftaul,ftau,nsp,mpi_real8,mpi_sum, &
         0,mpi_world,ierr)

!.....Write to the file only at node-0
    if( myid == 0 ) then
      write(io_impls,'(i10, f10.1, 11es12.4)') istp, simtime, ptau, &
           (ftau(is), is=1,nsp)
    endif
    
  end subroutine write_impulse
!=======================================================================
  subroutine final_impulse(myid)
    integer,intent(in):: myid
    if( myid == 0 ) close(io_impls)
  end subroutine final_impulse
end module impulse
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd lib"
!     End:
