function ispOf(tag)
  use mod_precision
  implicit none
  real(rp),intent(in):: tag
  integer:: ispOf
  ispOf= int(tag)
  return
end function ispOf
!=======================================================================
function ifmvOf(tag)
  use mod_precision
  implicit none
  real(rp),intent(in):: tag
  integer:: ifmvOf
  ifmvOf= int(mod(tag*10,10.0_rp))
  return
end function ifmvOf
!=======================================================================
function itotOf(tag)
  use mod_precision
  implicit none
  real(rp),intent(in):: tag
  integer:: itotOf
  real(rp):: tmp

  tmp= tag -ispOf(tag) -ifmvOf(tag)*1e-1_rp
  itotOf= nint(tmp*1e+14_rp)
  return
end function itotOf
!=======================================================================
subroutine replaceTag(ctx,ival,tag)
  use mod_precision
  implicit none
  character(len=*),intent(in):: ctx
  integer,intent(in):: ival
  real(rp),intent(inout):: tag

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
