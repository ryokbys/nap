module virtual_wall
!-----------------------------------------------------------------------
!                     Last modified: <2025-05-16 15:26:42 KOBAYASHI Ryo>
!-----------------------------------------------------------------------
! Module for virtual wall in MD simulation
!-----------------------------------------------------------------------
  use pmdvars,only: nvwall, ivwall, spos_vwall, iside_vwall, iovwall, &
       frc_vwall, dt,h,hi,fa2v,sorg, myid_md
  use util,only: itotOf
  implicit none
  private
  save

  public:: correct_pos_vwall, write_frc_vwall,bcast_vwall

  character(len=13),parameter:: cfrcname = 'out.frc_vwall'
       
contains
!=======================================================================
  subroutine bcast_vwall(mpicomm)
    include 'mpif.h'
    integer,intent(in):: mpicomm
    integer:: ierr

    call mpi_bcast(nvwall,1,mpi_integer,0,mpicomm,ierr)
    if( .not. allocated(ivwall) ) then
      allocate(ivwall(nvwall),spos_vwall(nvwall),iside_vwall(nvwall),&
           frc_vwall(nvwall))
    endif
    call mpi_bcast(ivwall,nvwall,mpi_integer,0,mpicomm,ierr)
    call mpi_bcast(spos_vwall,nvwall,mpi_real8,0,mpicomm,ierr)
    call mpi_bcast(iside_vwall,nvwall,mpi_integer,0,mpicomm,ierr)

  end subroutine bcast_vwall
!=======================================================================
  subroutine correct_pos_vwall(natm,tag,ra,va)
!
!  Correct positions using velocities with the following corrections
!  for atoms that cross the virtual wall:
!    - shift position to reflected place
!    - invert velocity
!    - add momentum change to force on wall
!
    implicit none
    integer,intent(in):: natm
    real(8),intent(in):: tag(natm)
    real(8),intent(inout):: ra(3,natm),va(3,natm)

    integer:: ivw,i123,iside,is,ia
    real(8):: spw,vt(3),rt

    frc_vwall(:) = 0d0
    do ivw=1,nvwall
      i123 = ivwall(ivw)
      spw = spos_vwall(ivw)
      iside = iside_vwall(ivw)
      do ia=1,natm
        vt(1:3) = (hi(1:3,1)*va(1,ia) &
             +hi(1:3,2)*va(2,ia) +hi(1:3,3)*va(3,ia))*dt
        rt = ra(i123,ia) +sorg(i123)
!.....Position is already updated before this routine is called.
!!$        ra(1:3,ia) = ra(1:3,ia) +vt(1:3)
!!$        if( myid_md.eq.1 ) print *,'itot,r(t-dt),r(t)=',itotOf(tag(ia)), &
!!$             iside*(rt-vt(i123)-spw), iside*(rt-spw)
        if( iside*(rt-vt(i123)-spw) > 0d0 .and. &
             iside*(rt-spw) < 0d0 )  then
          ra(i123,ia) = spw -(rt -spw) -sorg(i123)
          vt(i123) = -vt(i123)
          va(1:3,ia) = (h(1:3,1,0)*vt(1) +h(1:3,2,0)*vt(2) &
               +h(1:3,3,0)*vt(3)) /dt
!.....Note: this vt is in scaled by h-mat, thus frc_vwall should be scaled later.
          is = int(tag(ia))
          frc_vwall(ivw) = frc_vwall(ivw) +2d0*vt(i123)/dt /(fa2v(is)*2d0)
        endif
      enddo
    enddo

  end subroutine correct_pos_vwall
!=======================================================================
  subroutine write_frc_vwall(istp)
!
!  Write out force on virtual walls.
!
    integer,intent(in):: istp
    
    integer:: num,ivw,i123,iside
    real(8):: frci(3),frc
    logical:: lopen = .false.
    
    inquire(file=trim(cfrcname), number=num, opened=lopen)
    if( .not.lopen ) then
      open(iovwall, file=trim(cfrcname), status='replace')
      write(iovwall,'(a)',advance='no') '# istp, frcs on'
      do ivw=1,nvwall
        write(iovwall,'(a,i0)',advance='no') ' wall-',ivw
      enddo
      write(iovwall,'(a)') ''
    endif

    write(iovwall,'(i8)',advance='no') istp
    do ivw=1,nvwall
      i123 = ivwall(ivw)
      iside= iside_vwall(ivw)
      frci(1:3) = h(1:3,i123,0)*frc_vwall(ivw)
      frc = sqrt(frci(1)*frci(1) +frci(2)*frci(2) +frci(3)*frci(3))
      write(iovwall,'(e15.4)',advance='no') frc
    enddo
    write(iovwall,'(a)') ''
    
    
  end subroutine write_frc_vwall
end module virtual_wall
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd lib"
!     End:

