subroutine heapsort(n,nmax,tag,ndim,arr)
!-----------------------------------------------------------------------
!  Heap sort
!    - See Numerical Recipes in Fortran, Chap.8
!    - sort array TAG and ARR according to total id in TAG
!      by ascending order
!    - assuming tag includes only total id information
!-----------------------------------------------------------------------
  implicit none
  integer,intent(in):: n,nmax,ndim
  real(8),intent(inout):: tag(nmax),arr(ndim,nmax)

  integer:: i,ir,j,l,irra
  real(8):: rtag,rarr(ndim)

!-----check size
  if( n.lt.2 ) return

!  The index l will be decremented to 1 during the "hiring" phase.
!  Once it reaches 1, the index ir will be decremented to 1 during
!    the "retirement-and-promotion" phase.

  l= n/2+1
  ir=n
10 continue
  if( l.gt.1) then          ! still in hiring phase
    l=l-1
    rtag= tag(l)
    rarr(1:ndim)= arr(1:ndim,l)
  else                      ! retirement and promotion phase
    rtag=tag(ir)
    rarr(1:ndim)= arr(1:ndim,ir)
!
    tag(ir)= tag(1)
    arr(1:ndim,ir)= arr(1:ndim,1)
!
    ir=ir-1
    if(ir.eq.1)then
      tag(1)= rtag
      arr(1:ndim,1)= rarr(1:ndim)
      return
    endif
  endif
  i=l
  j=l+l
20 if( j.le.ir ) then
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
!=======================================================================
subroutine heapsort_itag(n,nmax,itag,ndim,arr)
!-----------------------------------------------------------------------
!  Heap sort
!    - See Numerical Recipes in Fortran, Chap.8
!    - sort array TAG and ARR according to total id in TAG
!      by ascending order
!    - assuming tag includes only total id information
!-----------------------------------------------------------------------
  implicit none
  integer,intent(in):: n,nmax,ndim
  integer,intent(inout):: itag(nmax)
  real(8),intent(inout):: arr(ndim,nmax)

  integer:: i,ir,j,l,irra,jtag
  real(8):: rarr(ndim)

!-----check size
  if( n.lt.2 ) return

!  The index l will be decremented to 1 during the "hiring" phase.
!  Once it reaches 1, the index ir will be decremented to 1 during
!    the "retirement-and-promotion" phase.

  l= n/2+1
  ir=n
10 continue
  if( l.gt.1) then          ! still in hiring phase
    l=l-1
    jtag= itag(l)
    rarr(1:ndim)= arr(1:ndim,l)
  else                      ! retirement and promotion phase
    jtag=itag(ir)
    rarr(1:ndim)= arr(1:ndim,ir)
!
    itag(ir)= itag(1)
    arr(1:ndim,ir)= arr(1:ndim,1)
!
    ir=ir-1
    if(ir.eq.1)then
      itag(1)= jtag
      arr(1:ndim,1)= rarr(1:ndim)
      return
    endif
  endif
  i=l
  j=l+l
20 if( j.le.ir ) then
    if( j.lt.ir ) then
      if( itag(j).lt.itag(j+1) ) j=j+1
    endif
    if( jtag.lt.itag(j) ) then
      itag(i)=itag(j)
      arr(1:ndim,i)= arr(1:ndim,j)
      i=j
      j=j+j
    else
      j=ir+1
    endif
    goto 20
  endif
  itag(i)=jtag
  arr(1:ndim,i)= rarr(1:ndim)
  goto 10

end subroutine heapsort_itag
!=======================================================================
subroutine arg_heapsort_itag(n,nmax,itag,idxarr)
!-----------------------------------------------------------------------
!  Heap sort
!    - See Numerical Recipes in Fortran, Chap.8
!    - sort array TAG and ARR according to total id in TAG
!      by ascending order
!    - assuming tag includes only total id information
!-----------------------------------------------------------------------
  implicit none
  integer,intent(in):: n,nmax
  integer,intent(inout):: itag(nmax)
  integer,intent(out):: idxarr(nmax)

  integer:: i,ir,j,l,irra,jtag,jdx
!!$  real(8):: rarr(ndim)

!-----check size
  if( n.lt.2 ) return

  do i=1,nmax
    idxarr(i) = i
  enddo

!  The index l will be decremented to 1 during the "hiring" phase.
!  Once it reaches 1, the index ir will be decremented to 1 during
!    the "retirement-and-promotion" phase.

  l= n/2+1
  ir=n
10 continue
  if( l.gt.1) then          ! still in hiring phase
    l=l-1
    jtag= itag(l)
    jdx = idxarr(l)
  else                      ! retirement and promotion phase
    jtag= itag(ir)
    jdx = idxarr(ir)
!
    itag(ir)= itag(1)
    idxarr(ir)= idxarr(1)
!
    ir=ir-1
    if(ir.eq.1)then
      itag(1)= jtag
      idxarr(1)= jdx
      return
    endif
  endif
  i=l
  j=l+l
20 if( j.le.ir ) then
    if( j.lt.ir ) then
      if( itag(j).lt.itag(j+1) ) j=j+1
    endif
    if( jtag.lt.itag(j) ) then
      itag(i)=itag(j)
      idxarr(i)= idxarr(j)
      i=j
      j=j+j
    else
      j=ir+1
    endif
    goto 20
  endif
  itag(i)= jtag
  idxarr(i)= jdx
  goto 10

end subroutine arg_heapsort_itag
!=======================================================================
subroutine heapsort_i(n,nmax,tag,iarr)
!-----------------------------------------------------------------------
!  Heap sort
!    - See Numerical Recipes in Fortran, Chap.8
!    - sort array TAG and IARR according to total id in TAG
!      by ascending order
!    - assuming tag includes only total id information
!-----------------------------------------------------------------------
  implicit none
  integer,intent(in):: n,nmax
  integer,intent(inout):: iarr(nmax)
  real(8),intent(inout):: tag(nmax)

  integer:: i,ir,j,l,irra,irarr
  real(8):: rtag

!-----check size
  if( n.lt.2 ) return

!  The index l will be decremented to 1 during the "hiring" phase.
!  Once it reaches 1, the index ir will be decremented to 1 during
!    the "retirement-and-promotion" phase.

  l= n/2+1
  ir=n
10 continue
  if( l.gt.1) then          ! still in hiring phase
    l=l-1
    rtag= tag(l)
    irarr=iarr(l)
!        rarr(1:ndim)= arr(1:ndim,l)
  else                      ! retirement and promotion phase
    rtag=tag(ir)
    irarr= iarr(ir)
!
    tag(ir)= tag(1)
    iarr(ir)= iarr(1)
!        arr(1:ndim,ir)= arr(1:ndim,1)
!
    ir=ir-1
    if(ir.eq.1)then
      tag(1)= rtag
      iarr(1)= irarr
!          arr(1:ndim,1)= rarr(1:ndim)
      return
    endif
  endif
  i=l
  j=l+l
20 if( j.le.ir ) then
    if( j.lt.ir ) then
      if( tag(j).lt.tag(j+1) ) j=j+1
    endif
    if( rtag.lt.tag(j) ) then
      tag(i)=tag(j)
      iarr(i)= iarr(j)
!          arr(1:ndim,i)= arr(1:ndim,j)
      i=j
      j=j+j
    else
      j=ir+1
    endif
    goto 20
  endif
  tag(i)=rtag
  iarr(i)= irarr
!      arr(1:ndim,i)= rarr(1:ndim)
  goto 10

end subroutine heapsort_i
!=======================================================================
recursive subroutine qsort_itag(nmax,il,ir,itag,ndim,arr)
!
!  Sort arr using itag by quick sort algorithm.
!
  implicit none
  integer,intent(in):: nmax,ndim,il,ir
  integer,intent(inout):: itag(nmax)
  real(8),intent(inout):: arr(ndim,nmax)

  integer:: ip,ipiv,i,j

  if( ir-il.lt.1 ) return
  ip = int((il+ir)/2)
  ipiv = itag(ip)
  call swap_itag(nmax,ip,ir,itag,ndim,arr)
  i = il
  do j=il,ir-1
    if( itag(j).lt.ipiv ) then
      call swap_itag(nmax,i,j,itag,ndim,arr)
      i = i + 1
    endif
  enddo
  call swap_itag(nmax,i,ir,itag,ndim,arr)
  call qsort_itag(nmax,il,i,itag,ndim,arr)
  call qsort_itag(nmax,i+1,ir,itag,ndim,arr)

end subroutine qsort_itag
!=======================================================================
recursive subroutine arg_qsort_itag(nmax,il,ir,itag,idxarr)
!
!  Sort itag array by quick sort algorithm and return sorted index array.
!
  implicit none
  integer,intent(in):: nmax,il,ir
  integer,intent(inout):: itag(nmax)
  integer,intent(out):: idxarr(nmax)

  integer:: ip,ipiv,i,j

  do i=1,nmax
    idxarr(i) = i
  enddo
  if( ir-il.lt.1 ) return
  ip = int((il+ir)/2)
  ipiv = itag(ip)
  call swap_itagidx(nmax,ip,ir,itag,idxarr)
  i = il
  do j=il,ir-1
    if( itag(j).lt.ipiv ) then
      call swap_itagidx(nmax,i,j,itag,idxarr)
      i = i + 1
    endif
  enddo
  call swap_itagidx(nmax,i,ir,itag,idxarr)
  call arg_qsort_itag(nmax,il,i,itag,idxarr)
  call arg_qsort_itag(nmax,i+1,ir,itag,idxarr)

end subroutine arg_qsort_itag
!=======================================================================
subroutine swap_itag(nmax,i,j,itag,ndim,arr)
  implicit none 
  integer,intent(in):: nmax,i,j,ndim
  integer,intent(inout):: itag(nmax)
  real(8),intent(inout):: arr(ndim,nmax)

  integer:: itmp
  real(8):: tmp(ndim)

  itmp = itag(i)
  itag(i) = itag(j)
  itag(j) = itmp
  
  tmp(1:ndim) = arr(1:ndim,i)
  arr(1:ndim,i) = arr(1:ndim,j)
  arr(1:ndim,j) = tmp(1:ndim)
  return
end subroutine swap_itag
!=======================================================================
subroutine swap_itagidx(nmax,i,j,itag,idx)
  implicit none 
  integer,intent(in):: nmax,i,j
  integer,intent(inout):: itag(nmax),idx(nmax)

  integer:: itmp

  itmp = itag(i)
  itag(i) = itag(j)
  itag(j) = itmp

  itmp = idx(i)
  idx(i) = idx(j)
  idx(j) = itmp

  return
end subroutine swap_itagidx
