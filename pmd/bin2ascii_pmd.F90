program bin2ascii_pmd
!-----------------------------------------------------------------------
!                     Last-modified: <2021-02-08 10:21:16 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
! Convert binary pmd format file to ascii one.
!-----------------------------------------------------------------------
! Usage:
!   $ /path/to/b2a pmd_bin_## pmd_ascii_##
!
! Output:
!   - pmd_ascii_##
!-----------------------------------------------------------------------
  use pmdio
  implicit none

  integer:: nargc
  character(len=128):: ciname, coname

  nargc= command_argument_count()
  if( nargc.ne.2 ) then
    stop 'Usage: $ /path/to/b2a pmd_bin_## pmd_ascii_##'
  endif
  call getarg(1,ciname)
  call getarg(2,coname)

  call read_pmdtot_bin(10,trim(ciname))
  call write_pmdtot_ascii(20,trim(coname))

  return
end program bin2ascii_pmd
