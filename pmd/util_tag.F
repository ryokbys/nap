      function ispOf(tag)
      implicit none
      real(8),intent(in):: tag
      integer:: ispOf
      ispOf= int(tag)
      return
      end function ispOf
c=======================================================================
      function ifmvOf(tag)
      implicit none
      real(8),intent(in):: tag
      integer:: ifmvOf
      ifmvOf= int(mod(tag*10,10d0))
      return
      end function ifmvOf
c=======================================================================
      function itotOf(tag)
      implicit none
      real(8),intent(in):: tag
      integer:: itotOf
      real(8):: tmp
      integer,external:: ispOf,ifmvOf
      tmp= tag -ispOf(tag) -ifmvOf(tag)*1d-1
      itotOf= nint(tmp*1d+14)
      return
      end function itotOf
