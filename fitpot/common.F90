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
