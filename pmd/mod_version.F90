module version
!-----------------------------------------------------------------------
!                     Last-modified: <2025-05-03 22:37:58 KOBAYASHI Ryo>
!-----------------------------------------------------------------------
! A module for version/revision.
!-----------------------------------------------------------------------
  implicit none
  private
  save

  public:: write_revision, write_authors
  
  character(len=128),parameter:: crevision = 'rev250503'

  character(len=128),parameter:: cauthors(1) = &
       (/ 'Ryo KOBAYASHI <kobayashi.ryo@nitech.ac.jp>' /)
  
contains
  subroutine write_revision()
    write(6,'(a)') '   Revision: '//trim(crevision)
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
end module version
