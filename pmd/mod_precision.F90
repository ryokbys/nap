module mod_precision
!-----------------------------------------------------------------------
!                     Last modified: <2026-05-08 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
! Module that defines the kind parameter rp (real precision) used
! throughout pmd. Compile with -D__SINGLE__ to use single precision.
!-----------------------------------------------------------------------
  implicit none
  save
#ifdef __SINGLE__
  integer,parameter:: rp = 4
#else
  integer,parameter:: rp = 8
#endif

end module mod_precision
