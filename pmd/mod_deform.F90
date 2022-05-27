module deform
!-----------------------------------------------------------------------
!  Module for applying deformation to the simulation cell.
!                     Last-modified: <2022-05-27 22:37:00 KOBAYASHI Ryo>
!-----------------------------------------------------------------------
  implicit none
  include "./const.h"
  save

  character(len=20):: cdeform= 'none'
  real(8):: trlx_deform = -1d0  ! default: minus, THIS MUST BE SPECIFIED BY THE USER
  real(8):: dhmat(1:3,1:3) = 0d0  ! deviation of h-matrix in given trlx_deform
  real(8):: ddhmat(1:3,1:3) = 0d0  ! deviation of h-matrix per fs
  
contains
!=======================================================================
  subroutine init_deform(hmat,myid,iprint)
    integer,intent(in):: myid,iprint
    real(8),intent(in):: hmat(3,3)

    integer:: i,j

    if( trlx_deform.lt.0d0 ) then
      if( myid.eq.0 ) print *,'ERROR: trlx_deform must be specified'&
           //' if you use deformation.'
      stop 1
    endif

    ddhmat(1:3,1:3) = dhmat(1:3,1:3) /trlx_deform

    if( myid.eq.0 .and. iprint.ge.ipl_basic ) then
      print '(a)',' Deformation'
      print '(a)','   Initial h-matrix:'
      print '(3x,3es12.3)',hmat(1,1:3)
      print '(3x,3es12.3)',hmat(2,1:3)
      print '(3x,3es12.3)',hmat(3,1:3)
      print '(a)','   Final h-matrix after deformation:'
      print '(3x,3es12.3)',hmat(1,1:3) + dhmat(1,1:3)
      print '(3x,3es12.3)',hmat(2,1:3) + dhmat(2,1:3)
      print '(3x,3es12.3)',hmat(3,1:3) + dhmat(3,1:3)
      print '(a)','   Deviation of h-matrix per fs:'
      print '(3x,3es12.3)',ddhmat(1,1:3)
      print '(3x,3es12.3)',ddhmat(2,1:3)
      print '(3x,3es12.3)',ddhmat(3,1:3)
    endif
    
    return
  end subroutine init_deform
!=======================================================================
  subroutine apply_deform(hmat,dt,simtime)
!
!  Apply deformation to the cell
!
    real(8),intent(inout):: hmat(3,3)
    real(8),intent(in):: dt,simtime

    if( simtime.gt.trlx_deform ) return
    hmat(1:3,1:3) = hmat(1:3,1:3) +ddhmat(1:3,1:3)*dt
    return
  end subroutine apply_deform
end module deform
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
