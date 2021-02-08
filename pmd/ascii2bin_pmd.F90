program ascii2bin_pmd
!-----------------------------------------------------------------------
!                     Last-modified: <2021-02-08 10:26:29 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
! Convert ascii pmd format file to binary one.
!-----------------------------------------------------------------------
! Usage:
!   $ /path/to/a2b pmd_ascii_## pmd_bin_##
!
! Output:
!   - pmd_bin_##
!-----------------------------------------------------------------------
  use pmdio
  implicit none

  integer:: nargc
  character(len=128):: ciname, coname

  nargc= command_argument_count()
  if( nargc.ne.2 ) then
    stop 'Usage: $ /path/to/b2a pmd_ascii_## pmd_bin_##'
  endif
  call getarg(1,ciname)
  call getarg(2,coname)

  call read_pmdtot_ascii(20,trim(ciname))
  call write_pmdtot_bin(10,trim(coname))

  return
end program ascii2bin_pmd
