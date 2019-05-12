module deform
!-----------------------------------------------------------------------
!  Module for applying deformation to the simulation cell.
!                     Last-modified: <2019-05-12 22:16:08 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
  implicit none
  save

!.....Increment of hmat per step
  real(8):: dhmat(3,3)
  
contains
!=======================================================================
  subroutine init_deform(nstp,hmat,dhratio,myid,iprint)
    integer,intent(in):: nstp,myid,iprint
    real(8),intent(in):: hmat(3,3),dhratio(3,3)

    integer:: i,j
    real(8):: hmatfin(3,3)

    if( myid.eq.0 .and. iprint.gt.0 ) then
      print '(a)',' Deformation'
      print '(a)','   Ratio of h-matrix deformation:'
      print '(3x,3es12.3)',dhratio(1,1:3)
      print '(3x,3es12.3)',dhratio(2,1:3)
      print '(3x,3es12.3)',dhratio(3,1:3)
    endif
    
!.....Determine deformation increament from the initial and final size
    if( nstp.eq.0 ) then
      dhmat(:,:) = 0d0
    else
      do i=1,3
        do j=1,3
          hmatfin(j,i) = hmat(j,i)*dhratio(j,i)
          dhmat(j,i) = (hmatfin(j,i) -hmat(j,i)) /nstp
        enddo
      enddo
    endif

    if( myid.eq.0 .and. iprint.gt.0 ) then
      print '(a)','   Increment of h-matrix deformation:'
      print '(3x,3es12.3)',dhmat(1,1:3)
      print '(3x,3es12.3)',dhmat(2,1:3)
      print '(3x,3es12.3)',dhmat(3,1:3)
    endif
    return
  end subroutine init_deform
!=======================================================================
  subroutine apply_deform(hmat)
!
!  Apply deformation to the cell
!
    real(8),intent(inout):: hmat(3,3)

    hmat(1:3,1:3) = hmat(1:3,1:3) +dhmat(1:3,1:3)
    return
  end subroutine apply_deform
end module deform
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
