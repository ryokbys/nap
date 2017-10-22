module version
!-----------------------------------------------------------------------
!                     Last-modified: <2017-10-21 11:08:58 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
! A module for version/revision.
!-----------------------------------------------------------------------
  implicit none
  character(len=128),parameter:: cversion = 'rev171021'

  character(len=128),parameter:: cauthors(1) = &
       (/ 'Ryo KOBAYASHI <kobayashi.ryo@nitech.ac.jp>' /)
  
contains
  subroutine write_version()
    write(6,'(a)') ' Revision: '//trim(cversion)
    return
  end subroutine write_version
!=======================================================================
  subroutine write_authors()
    integer:: i

    write(6,'(a)') ' Contributors: '
    do i=1,size(cauthors)
      write(6,'(a)') '   - '//trim(cauthors(i))
    enddo
    return
  end subroutine write_authors
end module version
