module force
!-----------------------------------------------------------------------
!                     Last-modified: <2017-12-13 20:51:44 Ryo KOBAYASHI>
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

end module force
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
