module version
!-----------------------------------------------------------------------
!                     Last-modified: <2025-05-03 22:37:58 KOBAYASHI Ryo>
!-----------------------------------------------------------------------
! A module for version/revision.
!-----------------------------------------------------------------------
  implicit none
  private
  save

  public:: write_revision, write_authors, write_runtime
  
  character(len=128),parameter:: crevision = 'rev260208'
  character(len=128),parameter:: cgithash = GIT_HASH
  character(len=128),parameter:: cauthors(1) = &
       (/ 'Ryo KOBAYASHI <kobayashi.ryo@nitech.ac.jp>' /)
  
contains
!=======================================================================
  subroutine write_revision()
!!$    write(6,'(a)') '   Revision: '//trim(crevision)
    write(6,'(a)') '   Git commit-ID: '//trim(cgithash)
    return
  end subroutine write_revision
!=======================================================================
  subroutine write_authors()
    integer:: i

    write(6,'(a)') '   Contributors: '
    do i=1,size(cauthors)
      write(6,'(a)') '     - '//trim(cauthors(i))
    enddo
    return
  end subroutine write_authors
!=======================================================================
  subroutine write_runtime()
    use,intrinsic:: iso_fortran_env, only: compiler_version
    
    write(6,'(/a)') '   Compiler: '//trim(compiler_version())
    write(6,'(a)',advance='no') '   Machine: '
    call execute_command_line("uname -snmr")
    return
  end subroutine write_runtime
end module version
