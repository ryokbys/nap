module util
!-----------------------------------------------------------------------
!                     Last modified: <2019-05-17 13:45:59 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
  implicit none
  save

contains
!=======================================================================
  function num_data(str,delim)
!
!  Return number of entries within str separated by delim.
!
    implicit none
    character(len=*),intent(in):: str
    character(len=1),intent(in):: delim
    integer:: num_data

    integer:: i

    i=1
    num_data = 0
    do
      if( i.gt.len(str) ) exit
      if( str(i:i).ne.delim ) then
        num_data = num_data + 1
        do
          i = i + 1
          if( i.gt.len(str) ) exit
          if( str(i:i).eq.delim ) exit
        end do
      end if
      i = i + 1
    end do
    return
  end function num_data
!=======================================================================
  function is_numeric(str)
!
!  Check if the given string is numeric or not.
!
    character(len=*),intent(in):: str
    logical:: is_numeric
    integer:: i
    i = verify(trim(str),'0123456789')
    is_numeric = i==0
    return
  end function is_numeric
!=======================================================================
  subroutine time_stamp(prefix)
    implicit none
    character(len=*),intent(in):: prefix
    character(len=10):: c1,c2,c3
    integer:: date(8)
    character(len=128):: cdate,ctime

    call date_and_time(c1,c2,c3,date)

    write(ctime,'(i0.2,a,i0.2,a,i0.2)') date(5),':',date(6) &
         ,':',date(7)
    write(cdate,'(i4,a,i0.2,a,i0.2)') date(1),'-',date(2),'-',date(3)
    write(6,'(a,1x,a,1x,a)') prefix, 'at '//trim(ctime), &
         'on '//trim(cdate)
    return
  end subroutine time_stamp
!=======================================================================
  function ispOf(tag)
    implicit none
    real(8),intent(in):: tag
    integer:: ispOf
    ispOf= int(tag)
    return
  end function ispOf
!=======================================================================
  function ifmvOf(tag)
    implicit none
    real(8),intent(in):: tag
    integer:: ifmvOf
    ifmvOf= int(mod(tag*10,10d0))
    return
  end function ifmvOf
!=======================================================================
  function itotOf(tag)
    implicit none
    real(8),intent(in):: tag
    integer:: itotOf
    real(8):: tmp

    tmp= tag -ispOf(tag) -ifmvOf(tag)*1d-1
    itotOf= nint(tmp*1d+14)
    return
  end function itotOf
!=======================================================================
  subroutine replaceTag(ctx,ival,tag)
    implicit none
    character(len=*),intent(in):: ctx
    integer,intent(in):: ival
    real(8),intent(inout):: tag

    integer:: ifmv

    if( trim(ctx) .eq. 'isp' ) then
! not implemented
    else if( trim(ctx) .eq. 'ifmv' ) then
      if( ival.lt.0 .or. ival.gt.9 ) then
        stop 'Error @replaceTag: ival.lt.0 .or. ival.gt.9'
      endif
      ifmv = ifmvOf(tag)
      tag = tag -ifmv*0.1 +ival*0.1
    else if( trim(ctx) .eq. 'itot' ) then
! not implemented
    endif

  end subroutine replaceTag
!=======================================================================
end module util
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:

