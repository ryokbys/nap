module version
!-----------------------------------------------------------------------
!                     Last-modified: <2022-01-18 23:08:08 KOBAYASHI Ryo>
!-----------------------------------------------------------------------
! A module for version/revision.
!-----------------------------------------------------------------------
  implicit none
  private
  save

  public:: write_version, write_authors
  
  character(len=128),parameter:: cversion = 'rev220118'

  character(len=128),parameter:: cauthors(1) = &
       (/ 'Ryo KOBAYASHI <kobayashi.ryo@nitech.ac.jp>' /)
  
contains
  subroutine write_version()
    write(6,'(a)') '   Revision: '//trim(cversion)
    return
  end subroutine write_version
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
