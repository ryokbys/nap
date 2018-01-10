module force
!-----------------------------------------------------------------------
!                     Last-modified: <2018-01-10 12:41:42 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
  implicit none
  save
  integer:: num_forces = -1
  character(len=128),allocatable:: force_list(:)

contains
!=======================================================================
  function use_force(force_name)
    implicit none
    character(len=*),intent(in):: force_name
    logical:: use_force
    integer:: i

    use_force = .false.
    do i=1,num_forces
      if( trim(force_name).eq.trim(force_list(i)) ) then
        use_force = .true.
        return
      endif
    enddo
    return
  end function use_force
!=======================================================================
  subroutine write_forces(myid)
    implicit none
    integer,intent(in):: myid
    integer:: i

    if( myid.eq.0 ) then
      write(6,'(/,a)') ' Use the following force-fields:'
      do i=1,num_forces
        write(6,'(a)') '   '//trim(force_list(i))
      enddo
    endif
  end subroutine write_forces
end module force
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
