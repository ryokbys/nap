function ndat_in_line(ionum,delim) result(ndat)
  implicit none
  integer,intent(in):: ionum
  character(len=1),intent(in):: delim
  integer:: ndat
  integer,external:: num_data
  character:: ctmp*128
  
  read(ionum,'(a)') ctmp
  ndat = num_data(trim(ctmp),delim)
  backspace(ionum)
  return
  
end function ndat_in_line
!=======================================================================
function num_data(str,delim)
  implicit none
  character(len=*),intent(in):: str
  character(len=1),intent(in):: delim
  integer:: num_data

  integer:: i
!    print *, 'len(str), str = ',len(str),str
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
