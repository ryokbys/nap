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
