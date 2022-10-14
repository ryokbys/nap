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
subroutine arg_heapsort_iarr(n,nmax,iarr,idxarr)
!-----------------------------------------------------------------------
!  Heap sort
!    - See Numerical Recipes in Fortran, Chap.8
!    - sort indices according to integer array IARR
!      by ascending order
!-----------------------------------------------------------------------
  implicit none
  integer,intent(in):: n,nmax
  integer,intent(inout):: iarr(nmax)
  integer,intent(out):: idxarr(nmax)

  integer:: i,ir,j,l,irra,jarr,jdx
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
    jarr= iarr(l)
    jdx = idxarr(l)
  else                      ! retirement and promotion phase
    jarr= iarr(ir)
    jdx = idxarr(ir)
!
    iarr(ir)= iarr(1)
    idxarr(ir)= idxarr(1)
!
    ir=ir-1
    if(ir.eq.1)then
      iarr(1)= jarr
      idxarr(1)= jdx
      return
    endif
  endif
  i=l
  j=l+l
20 if( j.le.ir ) then
    if( j.lt.ir ) then
      if( iarr(j).lt.iarr(j+1) ) j=j+1
    endif
    if( jarr.lt.iarr(j) ) then
      iarr(i)=iarr(j)
      idxarr(i)= idxarr(j)
      i=j
      j=j+j
    else
      j=ir+1
    endif
    goto 20
  endif
  iarr(i)= jarr
  idxarr(i)= jdx
  goto 10

end subroutine arg_heapsort_iarr
!=======================================================================
subroutine arg_heapsort_arr(n,nmax,arr,idxarr)
!-----------------------------------------------------------------------
!  Heap sort
!    - See Numerical Recipes in Fortran, Chap.8
!    - sort indices according to REAL8 array ARR
!      by ascending order
!-----------------------------------------------------------------------
  implicit none
  integer,intent(in):: n,nmax
  real(8),intent(inout):: arr(nmax)
  integer,intent(out):: idxarr(nmax)

  integer:: i,ir,j,l,irra,jdx
  real(8):: tmp
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
    tmp= arr(l)
    jdx = idxarr(l)
  else                      ! retirement and promotion phase
    tmp= arr(ir)
    jdx = idxarr(ir)
!
    arr(ir)= arr(1)
    idxarr(ir)= idxarr(1)
!
    ir=ir-1
    if(ir.eq.1)then
      arr(1)= tmp
      idxarr(1)= jdx
      return
    endif
  endif
  i=l
  j=l+l
20 if( j.le.ir ) then
    if( j.lt.ir ) then
      if( arr(j).lt.arr(j+1) ) j=j+1
    endif
    if( tmp.lt.arr(j) ) then
      arr(i)=arr(j)
      idxarr(i)= idxarr(j)
      i=j
      j=j+j
    else
      j=ir+1
    endif
    goto 20
  endif
  arr(i)= tmp
  idxarr(i)= jdx
  goto 10

end subroutine arg_heapsort_arr
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
recursive subroutine qsort_iarr(nmax,il,ir,iarr,ndim,arr)
!
!  Sort arr using iarr by quick sort algorithm.
!
  implicit none
  integer,intent(in):: nmax,ndim,il,ir
  integer,intent(inout):: iarr(nmax)
  real(8),intent(inout):: arr(ndim,nmax)

  integer:: ip,ipiv,i,j

  if( ir-il.lt.1 ) return
  ip = int((il+ir)/2)
  ipiv = iarr(ip)
  call swap_iarr(nmax,ip,ir,iarr,ndim,arr)
  i = il
  do j=il,ir-1
    if( iarr(j).lt.ipiv ) then
      call swap_iarr(nmax,i,j,iarr,ndim,arr)
      i = i + 1
    endif
  enddo
  call swap_iarr(nmax,i,ir,iarr,ndim,arr)
  call qsort_iarr(nmax,il,i,iarr,ndim,arr)
  call qsort_iarr(nmax,i+1,ir,iarr,ndim,arr)

end subroutine qsort_iarr
!=======================================================================
recursive subroutine arg_qsort_iarr(nmax,il,ir,iarr,idxarr)
!
!  Sort iarr array by quick sort algorithm and return sorted index array.
!
  implicit none
  integer,intent(in):: nmax,il,ir
  integer,intent(inout):: iarr(nmax)
  integer,intent(out):: idxarr(nmax)

  integer:: ip,ipiv,i,j

  do i=1,nmax
    idxarr(i) = i
  enddo
  if( ir-il.lt.1 ) return
  ip = int((il+ir)/2)
  ipiv = iarr(ip)
  call swap_iarridx(nmax,ip,ir,iarr,idxarr)
  i = il
  do j=il,ir-1
    if( iarr(j).lt.ipiv ) then
      call swap_iarridx(nmax,i,j,iarr,idxarr)
      i = i + 1
    endif
  enddo
  call swap_iarridx(nmax,i,ir,iarr,idxarr)
  call arg_qsort_iarr(nmax,il,i,iarr,idxarr)
  call arg_qsort_iarr(nmax,i+1,ir,iarr,idxarr)

end subroutine arg_qsort_iarr
!=======================================================================
recursive subroutine arg_qsort_arr(nmax,il,ir,arr,idxarr)
!
!  Sort arr array by quick sort algorithm and return sorted index array.
!
  implicit none
  integer,intent(in):: nmax,il,ir
  real(8),intent(inout):: arr(nmax)
  integer,intent(out):: idxarr(nmax)

  integer:: ip,i,j
  real(8):: aip

  do i=1,nmax
    idxarr(i) = i
  enddo
  if( ir-il.lt.1 ) return
  ip = int((il+ir)/2)
  aip = arr(ip)
  call swap_arridx(nmax,ip,ir,arr,idxarr)
  i = il
  do j=il,ir-1
    if( arr(j).lt.aip ) then
      call swap_arridx(nmax,i,j,arr,idxarr)
      i = i + 1
    endif
  enddo
  call swap_arridx(nmax,i,ir,arr,idxarr)
  call arg_qsort_arr(nmax,il,i,arr,idxarr)
  call arg_qsort_arr(nmax,i+1,ir,arr,idxarr)

end subroutine arg_qsort_arr
!=======================================================================
subroutine swap_iarr(nmax,i,j,iarr,ndim,arr)
  implicit none 
  integer,intent(in):: nmax,i,j,ndim
  integer,intent(inout):: iarr(nmax)
  real(8),intent(inout):: arr(ndim,nmax)

  integer:: itmp
  real(8):: tmp(ndim)

  itmp = iarr(i)
  iarr(i) = iarr(j)
  iarr(j) = itmp
  
  tmp(1:ndim) = arr(1:ndim,i)
  arr(1:ndim,i) = arr(1:ndim,j)
  arr(1:ndim,j) = tmp(1:ndim)
  return
end subroutine swap_iarr
!=======================================================================
subroutine swap_iarridx(nmax,i,j,iarr,idx)
  implicit none 
  integer,intent(in):: nmax,i,j
  integer,intent(inout):: iarr(nmax),idx(nmax)

  integer:: itmp

  itmp = iarr(i)
  iarr(i) = iarr(j)
  iarr(j) = itmp

  itmp = idx(i)
  idx(i) = idx(j)
  idx(j) = itmp

  return
end subroutine swap_iarridx
!=======================================================================
subroutine swap_arridx(nmax,i,j,arr,idx)
  implicit none 
  integer,intent(in):: nmax,i,j
  real(8),intent(inout):: arr(nmax)
  integer,intent(inout):: idx(nmax)

  integer:: itmp
  real(8):: tmp

  tmp = arr(i)
  arr(i) = arr(j)
  arr(j) = tmp

  itmp = idx(i)
  idx(i) = idx(j)
  idx(j) = itmp

  return
end subroutine swap_arridx
