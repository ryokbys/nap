      subroutine heapsort(n,nmax,tag,ndim,arr)
c-----------------------------------------------------------------------
c  Heap sort
c    - See Numerical Recipes in Fortran, Chap.8
c    - sort array TAG and ARR according to total id in TAG
c      by ascending order
c    - assuming tag includes only total id information
c-----------------------------------------------------------------------
      implicit none
      integer,intent(in):: n,nmax,ndim
      real(8),intent(inout):: tag(nmax),arr(ndim,nmax)

      integer:: i,ir,j,l,irra
      real(8):: rtag,rarr(ndim)

c-----check size
      if( n.lt.2 ) return

c  The index l will be decremented to 1 during the "hiring" phase.
c  Once it reaches 1, the index ir will be decremented to 1 during
c    the "retirement-and-promotion" phase.

      l= n/2+1
      ir=n
 10   continue
      if( l.gt.1) then          ! still in hiring phase
        l=l-1
        rtag= tag(l)
        rarr(1:ndim)= arr(1:ndim,l)
      else                      ! retirement and promotion phase
        rtag=tag(ir)
        rarr(1:ndim)= arr(1:ndim,ir)
c
        tag(ir)= tag(1)
        arr(1:ndim,ir)= arr(1:ndim,1)
c
        ir=ir-1
        if(ir.eq.1)then
          tag(1)= rtag
          arr(1:ndim,1)= rarr(1:ndim)
          return
        endif
      endif
      i=l
      j=l+l
 20   if( j.le.ir ) then
        if( j.lt.ir ) then
          if( tag(j).lt.tag(j+1) ) j=j+1
        endif
        if( rtag.lt.tag(j) ) then
          tag(i)=tag(j)
          arr(1:ndim,i)= arr(1:ndim,j)
          i=j
          j=j+j
        else
          j=ir+1
        endif
        goto 20
      endif
      tag(i)=rtag
      arr(1:ndim,i)= rarr(1:ndim)
      goto 10

      end subroutine heapsort
c=======================================================================
      subroutine heapsort_i(n,nmax,tag,iarr)
c-----------------------------------------------------------------------
c  Heap sort
c    - See Numerical Recipes in Fortran, Chap.8
c    - sort array TAG and IARR according to total id in TAG
c      by ascending order
c    - assuming tag includes only total id information
c-----------------------------------------------------------------------
      implicit none
      integer,intent(in):: n,nmax
      integer,intent(inout):: iarr(nmax)
      real(8),intent(inout):: tag(nmax)

      integer:: i,ir,j,l,irra,irarr
      real(8):: rtag

c-----check size
      if( n.lt.2 ) return

c  The index l will be decremented to 1 during the "hiring" phase.
c  Once it reaches 1, the index ir will be decremented to 1 during
c    the "retirement-and-promotion" phase.

      l= n/2+1
      ir=n
 10   continue
      if( l.gt.1) then          ! still in hiring phase
        l=l-1
        rtag= tag(l)
        irarr=iarr(l)
c        rarr(1:ndim)= arr(1:ndim,l)
      else                      ! retirement and promotion phase
        rtag=tag(ir)
        irarr= iarr(ir)
c
        tag(ir)= tag(1)
        iarr(ir)= iarr(1)
c        arr(1:ndim,ir)= arr(1:ndim,1)
c
        ir=ir-1
        if(ir.eq.1)then
          tag(1)= rtag
          iarr(1)= irarr
c          arr(1:ndim,1)= rarr(1:ndim)
          return
        endif
      endif
      i=l
      j=l+l
 20   if( j.le.ir ) then
        if( j.lt.ir ) then
          if( tag(j).lt.tag(j+1) ) j=j+1
        endif
        if( rtag.lt.tag(j) ) then
          tag(i)=tag(j)
          iarr(i)= iarr(j)
c          arr(1:ndim,i)= arr(1:ndim,j)
          i=j
          j=j+j
        else
          j=ir+1
        endif
        goto 20
      endif
      tag(i)=rtag
      iarr(i)= irarr
c      arr(1:ndim,i)= rarr(1:ndim)
      goto 10

      end subroutine heapsort_i
