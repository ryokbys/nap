      subroutine getarg(arg,argt)
!
!     get a command line argument using gfortran
!
      implicit none
!
! input
!
!     number of argument to get
      integer arg
!
! output
!
!     string
      character*(*) argt
!
! executable
!
      call get_command_argument(arg,argt)

      return
      end
