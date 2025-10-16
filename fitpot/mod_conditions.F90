module conditions
!=======================================================================
!  Condition information is loaded from "in.vars.conditions"
!  The file format is like the following:
!-----------------------------------------------------------------------  
!  #  var-ID-LHS, operator,  var-ID-RHS
!     30    ! num of conditions
!     2    <  1
!     3    <= 2
!     ...
!     10   == 9
!-----------------------------------------------------------------------
!  - var-ID-LHS/RHS: var-IDs of LHS and RHS of the operator
!  - Conditions are written as "operator" + "operand",
!    and the operand is the var-ID-RHS.
!  - NOTE: var-ID-RHS of operator must be smaller than var-ID-LHS.
!  - Available operators (with operator-ID):
!    1) "<"
!    2) ">"
!    3) "<="
!    4) ">="
!    5) "=="
!=======================================================================
  implicit none
  include 'mpif.h'
  save
  private
  public:: lconds, read_conds, calc_fpenal_conds, calc_gpenal_conds
  
  character(len=128):: cfname = 'in.vars.conditions'

  integer:: ionum = 81
  logical:: lconds = .false.
  integer:: nconds = 0  ! number of conditions
  real(8):: wgt = 1d0   ! weight for the condition
  integer,allocatable:: ivls(:), ivrs(:), iops(:)

contains
!=======================================================================
  subroutine read_conds(myid, mpi_world, iprint)
    integer,intent(in):: myid, mpi_world, iprint

    integer:: istat, ierr, inc
    character(len=128):: c1, c2, c3

    istat = 0
    if( myid.eq.0 ) then
      open(ionum,file=trim(cfname),status='old')
!.....1st detect number of conditions
      do while(.true.)
        read(ionum, *, end=99) c1
        if( c1(1:1).eq.'!' .or. c1(1:1).eq.'#') cycle
        backspace(ionum)
        read(ionum, *) nconds, wgt
!.....Exit this loop after reading the 1st non-comment line
!     that should be "num and weight of conditions".
10      exit
      enddo
      if( nconds.eq.0 ) then
        print *,' ERROR: nconds must not be 0.'
        istat = 1
        goto 99
      endif
      rewind(ionum)
!.....Allocate some
      allocate( ivls(nconds), ivrs(nconds), iops(nconds) )
!.....Then, read conditions
      inc = -1
      do while(.true.)
        read(ionum, *, end=99) c1
        if( c1(1:1).eq.'!' .or. c1(1:1).eq.'#') cycle
        inc = inc + 1
        if( inc.lt.1 ) cycle
        if( inc.gt.nconds ) goto 99
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
        read(c1,*) ivls(inc)
        read(c3,*) ivrs(inc)
        if( ivrs(inc).ge.ivls(inc) ) then
          print *,' ERROR: var-ID-RHS must be smaller than var-ID-LHS' &
               //' in in.vars.conditions'
          istat = 1
          goto 99
        endif
        if( trim(c2).eq.'<' ) iops(inc) = 1
        if( trim(c2).eq.'>' ) iops(inc) = 2
        if( trim(c2).eq.'<=' ) iops(inc) = 3
        if( trim(c2).eq.'>=' ) iops(inc) = 4
        if( trim(c2).eq.'==' ) iops(inc) = 5
      enddo
99    close(ionum)
    endif
    call mpi_bcast(istat, 1, mpi_integer, 0, mpi_world, ierr)
    if( istat.gt.0 ) stop 1
    call mpi_bcast(nconds, 1, mpi_integer, 0, mpi_world, ierr)
    if( .not. allocated(ivls)) then
      allocate( ivls(nconds), ivrs(nconds), iops(nconds) )
    endif
    call mpi_bcast(ivls, nconds, mpi_integer, 0, mpi_world, ierr)
    call mpi_bcast(ivrs, nconds, mpi_integer, 0, mpi_world, ierr)
    call mpi_bcast(iops, nconds, mpi_integer, 0, mpi_world, ierr)
    return
  end subroutine read_conds
!=======================================================================
  subroutine calc_fpenal_conds(ndim,x,fp)
!
!  Calculate f of penalty w.r.t. conditions.
!  This is supposed to be called in a serial process.
!
    integer,intent(in):: ndim
    real(8),intent(in):: x(ndim)
    real(8),intent(inout):: fp

    integer:: ic, iop
    real(8):: vl, vr

!.....fp is initialized before this is called.    
!    fp = 0d0
    do ic=1,nconds
      vl = x(ivls(ic))
      vr = x(ivrs(ic))
      iop = iops(ic)
      if( iop.eq.1 ) then  ! '<'
        if( vl.ge.vr ) fp = fp + wgt*(vl-vr)**2
      else if( iop.eq.2 ) then  ! '>'
        if( vl.le.vr ) fp = fp + wgt*(vr-vl)**2
      else if( iop.eq.3 ) then  ! '<='
        if( vl.gt.vr ) fp = fp + wgt*(vl-vr)**2
      else if( iop.eq.4 ) then  ! '>='
        if( vl.lt.vr ) fp = fp + wgt*(vr-vl)**2
      else if( iop.eq.5 ) then  ! '=='
        fp = fp + wgt*(vl-vr)**2
      endif
    enddo
    
    return
  end subroutine calc_fpenal_conds
!=======================================================================
  subroutine calc_gpenal_conds(ndim,x,gp)
!
!  Calculate g of penalty w.r.t. conditions
!  This is supposed to be called in a serial process.
!
    integer,intent(in):: ndim
    real(8),intent(in):: x(ndim)
    real(8),intent(inout):: gp(ndim)

    integer:: ic, iop, ivl, ivr
    real(8):: vl, vr

!.....gp is initialized before this is called.
!    gp(:) = 0d0
    do ic=1,nconds
      ivl = ivls(ic)
      ivr = ivrs(ic)
      iop = iops(ic)
      vl = x(ivl)
      vr = x(ivr)
      if( iop.eq.1 ) then  ! '<'
        if( vl.ge.vr ) then
          gp(ivl) = gp(ivl) + 2d0*wgt*(vl-vr)
          gp(ivr) = gp(ivr) + 2d0*wgt*(vr-vl)
        endif
      else if( iop.eq.2 ) then  ! '>'
        if( vl.le.vr ) then
          gp(ivl) = gp(ivl) + 2d0*wgt*(vr-vl)
          gp(ivr) = gp(ivr) + 2d0*wgt*(vl-vr)
        endif
      else if( iop.eq.3 ) then  ! '<='
        if( vl.gt.vr ) then
          gp(ivl) = gp(ivl) + 2d0*wgt*(vl-vr)
          gp(ivr) = gp(ivr) + 2d0*wgt*(vr-vl)
        endif
      else if( iop.eq.4 ) then  ! '>='
        if( vl.lt.vr ) then
          gp(ivl) = gp(ivl) + 2d0*wgt*(vr-vl)
          gp(ivr) = gp(ivr) + 2d0*wgt*(vl-vr)
        endif
      else if( iop.eq.5 ) then  ! '=='
        gp(ivl) = gp(ivl) + 2d0*wgt*(vl-vr)
        gp(ivr) = gp(ivr) + 2d0*wgt*(vr-vl)
      endif
    enddo
    
    return
  end subroutine calc_gpenal_conds
end module conditions
!-----------------------------------------------------------------------
! Local Variables:
! compile-command: "make fitpot"
! End:
