module conditions
  use pmdvars,only: nspmax
  implicit none
  save
  private

!.....Condition information is loaded from "in.vars.conditions"
!  The file format is like the following:
!-----------------------------------------------------------------------  
!  #  var-ID, operator, operand
!     30    ! num of conditions
!     v2    < v1
!     v3    <= v2
!     ...
!     v10   == v9
!-----------------------------------------------------------------------
!  - "v10" -- v + digits means 10-th variable
!  - Conditions are written as "operator" + "operand",
!    and the operand can be either only "v##" or a value.
!  - NOTE: "##" in RHS of operator must be smaller than that in LHS.
!  - Available operators:
!    1) "<"
!    2) ">"
!    3) "<="
!    4) ">="
!    5) "=="
!  
  character(len=128):: cfname = 'in.vars.conditions'

  integer:: ionum = 81
  integer:: nconds  ! number of conditions

contains
!=======================================================================
  subroutine read_conditions(myid, mpi_world, iprint)
    integer,intent(in):: myid, mpi_world, iprint

    integer:: istat, ierr

    istat = 0
    if( myid.eq.0 ) then
      open(ionum,file=trim(cffname),status='old')
!.....1st detect number of conditions
      do while(.true.)
        read(ionum, *, end=99) c1
        if( c1(1:1).eq.'!' .or. c1(1:1).eq.'#') cycle
        if(  c1(1:1).eq.'1' .or. &
             c1(1:1).eq.'2' .or. &
             c1(1:1).eq.'3' .or. &
             c1(1:1).eq.'4' .or. &
             c1(1:1).eq.'5' .or. &
             c1(1:1).eq.'6' .or. &
             c1(1:1).eq.'7' .or. &
             c1(1:1).eq.'8' .or. &
             c1(1:1).eq.'9' ) then
          read(c1,*) nconds
10        exit
        endif
      enddo
      rewind(ionum)
!.....Then, read conditions
      inc = 0
      do while(.true.)
        read(ionum, *, end=99) c1
        if( .not. c1(1:1).eq.'v' ) cycle
        inc = inc + 1
        if( inc .gt. nconds ) exit
        backspace(ionum)
        read(ionum, *, end=99) c1, c2, c3
        if( trim(c2).ne.'<' .and. trim(c2).ne.'>' .and. &
             trim(c2).ne.'<=' .and. trim(c2).ne.'>=' &
             .and. trim(c2).ne.'==' ) then
          print *,' ERROR: operator must be either <, >, <=, >=, or ==,' &
               //' but '//trim(c2)//' is given.'
          istat = 1
          goto 99
        endif
        
      enddo
99    close(ionum)
    endif
    call mpi_bcast(istat, 1, mpi_integer, 0, mpi_world, ierr)

    if( istat.gt.0 ) stop 1
    
  end subroutine read_conditions
end module conditions
!-----------------------------------------------------------------------
! Local Variables:
! compile-command: "make fitpot"
! End:
