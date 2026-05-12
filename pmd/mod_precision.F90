module mod_precision
!-----------------------------------------------------------------------
!                     Last modified: <2026-05-12 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
! Module that defines the kind parameter rp (real precision) used
! throughout pmd. Compile with -D__SINGLE__ to use single precision.
! rpfmt  -- es format string for screen output (no exponent-width spec)
! rpfmt3 -- es format string for file I/O (3-digit exponent, e3)
!-----------------------------------------------------------------------
  implicit none
  save
#ifdef __SINGLE__
  integer,parameter:: rp = 4
  character(len=*),parameter:: rpfmt  = 'es16.8'
  character(len=*),parameter:: rpfmt3 = 'es18.8e3'
#else
  integer,parameter:: rp = 8
  character(len=*),parameter:: rpfmt  = 'es22.14'
  character(len=*),parameter:: rpfmt3 = 'es23.14e3'
#endif

end module mod_precision
