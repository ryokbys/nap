module util
!-----------------------------------------------------------------------
!                     Last modified: <2025-05-16 13:43:26 KOBAYASHI Ryo>
!-----------------------------------------------------------------------
!  Utility functions/subroutines used in nap.
!-----------------------------------------------------------------------
  implicit none
  save

  character(*),private,parameter:: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
  character(*),private,parameter:: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
contains
!=======================================================================
  function num_data(str,delim) result(num)
!
!  Return number of entries within str separated by delim.
!
    implicit none
    character(len=*),intent(in):: str
    character(len=1),intent(in):: delim
    integer:: num

    integer:: i

    i=1
    num = 0
    do
      if( i.gt.len(str) ) exit
      if( str(i:i).ne.delim ) then
        num = num + 1
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
!    use,intrinsic:: ieee_arithmetic
    character(len=*),intent(in):: str
    logical:: is_numeric
    real(8):: x
    integer:: e

    read(str,*,iostat=e) x
!    is_numeric = ( e == 0 .and. .not.isnan(x) )
!    is_numeric = ( e == 0 .and. .not.ieee_is_nan(x) )
    is_numeric = ( e == 0 .and. (x*0d0.eq.0d0) )
    return
  end function is_numeric
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
  function ithOf(tag,ith)
    implicit none
    real(8),intent(in):: tag
    integer,intent(in):: ith
    integer:: ithOf
    ithOf= int(mod(tag*10d0**ith,10d0))
    return
  end function ithOf
!=======================================================================
  function itotOf(tag)
!
!  Get total-ID of an atom of the given TAG.
!  Assuming that 1st 5 floating-number digits are for ifmv and igroup,
!  and the last 9 floating-number digits are for the total-ID.
!  For example, in the following tag,
!     1.20000000001234
!            ^^^^^^^^^ <-- these digits are for the total-ID,
!  and the other digits are reserved for species, ifmv, and igroups.
!
    implicit none
    real(8),intent(in):: tag
    integer:: itotOf
    real(8):: tmp

!!$    tmp= tag -ispOf(tag) -ifmvOf(tag)*1d-1
!!$    itotOf= nint(tmp*1d+14)
    tmp = tag*1d+5
    tmp = tmp -int(tmp)
    itotOf = nint(tmp*1d+9)
    return
  end function itotOf
!=======================================================================
  function igvarOf(tag,gid)
!
!  Get the value of a given group-ID (GID) from 2nd-5th digits of 
!  floating-number, e.g.,
!     1.21234000001234
!        ^^^^   <------ these are groups
!  And, 1st-4th groups correspond to 2nd-5th digits, respectively.
!  Thus, the max number of groups in this implementation is 4.
!
    real(8),intent(in):: tag
    integer,intent(in):: gid
    integer:: igvarOf

    if( gid.lt.1 .or. gid.gt.4 ) then
      print *,'ERROR: gid must be 1<=(gid)<=4, but gid=',gid
      stop
    endif
    igvarOf = ithOf(tag,gid+1)
    return
  end function igvarOf
!=======================================================================
  function iauxof(cauxname)
    use pmdvars,only: cauxarr
    character(len=*),intent(in):: cauxname
    integer:: iauxof
    integer:: i
    
    iauxof = 0
    do i=1,size(cauxarr)
      if( trim(cauxname).eq.trim(cauxarr(i)) ) then
        iauxof = i
        return
      endif
    enddo
    if( iauxof.eq.0 ) then
      print *,'ERROR @iauxof: No such auxname, '//trim(cauxname)
      stop 
    endif
    return
  end function iauxof
!=======================================================================
  function csp2isp(csp) result(isp)
!
!  Convert csp to isp.
!  If not found, return -1.
!
    use pmdvars,only: specorder,nspmax
    character(len=*),intent(in):: csp
    integer:: isp

    integer:: i

    isp = -1
    do i=1,nspmax
      if( trim(csp).eq.trim(specorder(i)) ) then
        isp = i
        return
      endif
    enddo
    return
  end function csp2isp
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
  subroutine replace_igvar(tagi,gid,igvar)
!
!  Replace the current group variable of the given GID with the given IGVAR.
!
    real(8),intent(inout):: tagi
    integer,intent(in):: gid,igvar
    integer:: igvar0
    
    igvar0 = ithOf(tagi,gid+1)
    tagi = tagi +(igvar -igvar0) *10d0**(-(gid+1))
    return
  end subroutine replace_igvar
!=======================================================================
  subroutine cell_info(h)
    use vector,only: dot
    implicit none
    real(8),intent(in):: h(3,3)

    integer:: i,j,jm,jp,im,ip
    real(8):: a,b,c,alpha,beta,gamma,sgm(3,3),vol
    real(8),parameter:: pi = 3.14159265358979d0

    write(6,*) ''
    write(6,'(a)') " Lattice vectors:"
    write(6,'(a,"[ ",3f12.3," ]")') '   a = ',h(1:3,1)
    write(6,'(a,"[ ",3f12.3," ]")') '   b = ',h(1:3,2)
    write(6,'(a,"[ ",3f12.3," ]")') '   c = ',h(1:3,3)

    a = dsqrt(dot(h(1:3,1),h(1:3,1)))
    b = dsqrt(dot(h(1:3,2),h(1:3,2)))
    c = dsqrt(dot(h(1:3,3),h(1:3,3)))
    alpha = acos(dot(h(1:3,2),h(1:3,3))/b/c) /pi *180d0
    beta  = acos(dot(h(1:3,1),h(1:3,3))/a/c) /pi *180d0
    gamma = acos(dot(h(1:3,1),h(1:3,2))/a/b) /pi *180d0

    write(6,'(a)') ' Lattice parameters:'
    write(6,'(a,f10.3,a,f7.2,a)') '   |a| = ',a,' Ang.,  alpha = ' &
         ,alpha,' deg.'
    write(6,'(a,f10.3,a,f7.2,a)') '   |b| = ',b,' Ang.,  beta  = ' &
         ,beta,' deg.'
    write(6,'(a,f10.3,a,f7.2,a)') '   |c| = ',c,' Ang.,  gamma = ' &
         ,gamma,' deg.'

    vol = get_vol(h)
    write(6,'(a,f0.2,a)') ' Cell volume = ',vol,' Ang^3'

  end subroutine cell_info
!=======================================================================
  subroutine spcs_info(ntot,tagtot)
    use pmdvars,only: nspmax,specorder
    integer,intent(in):: ntot
    real(8),intent(in):: tagtot(ntot)
    integer:: nsps(nspmax),ia,is

    nsps(:) = 0
    do ia=1,ntot
      is = int(tagtot(ia))
      nsps(is) = nsps(is) +1
    enddo
    write(6,'(a)') ' Number of each species in the initial configuration'
    do is=1,nspmax
      if( trim(specorder(is)).eq.'x' ) cycle
      write(6,'(a,i0)') '   '//specorder(is)//':  ',nsps(is)
    enddo

  end subroutine spcs_info
!=======================================================================
  function get_vol(h) result(vol)
    real(8),intent(in):: h(3,3)

    integer:: i,j,jm,jp,im,ip
    real(8):: sgm(3,3)
    real(8):: vol

!.....cofactor matrix, SGM
    do j=1,3
      jm=mod(j+1,3)+1
      jp=mod(j,  3)+1
      do i=1,3
        im=mod(i+1,3)+1
        ip=mod(i,  3)+1
        sgm(i,j)=h(ip,jp)*h(im,jm)-h(im,jp)*h(ip,jm)
      enddo
    enddo
!.....MD-box volume
    vol= h(1,1)*sgm(1,1) +h(2,1)*sgm(2,1) +h(3,1)*sgm(3,1)
    return
  end function get_vol
!=======================================================================
  subroutine make_cdumpauxarr()
!
!  Builld cdumpauxarr if ifpmd==2 (dump output)
!
    use pmdvars,only: ndumpaux,cdumpaux,cdumpauxarr
    integer:: i,ivx,ivy,ivz
    character(len=6):: ctmp

    ndumpaux = num_data(trim(cdumpaux),' ')
    if( allocated(cdumpauxarr) ) deallocate(cdumpauxarr)
    allocate(cdumpauxarr(ndumpaux))
    read(cdumpaux,*) (cdumpauxarr(i),i=1,ndumpaux)
!.....If cdumpauxarr contains vx,vy,vz, bring them to the beginning of the array
    ivz = idumpauxof('vz')
    if( ivz.gt.0 ) then
      ctmp = cdumpauxarr(ivz)
      do i=ivz-1,1,-1
        cdumpauxarr(i+1) = cdumpauxarr(i)
      enddo
      cdumpauxarr(1) = ctmp
    endif
    ivy = idumpauxof('vy')
    if( ivy.gt.0 ) then
      ctmp = cdumpauxarr(ivy)
      do i=ivy-1,1,-1
        cdumpauxarr(i+1) = cdumpauxarr(i)
      enddo
      cdumpauxarr(1) = ctmp
    endif
    ivx = idumpauxof('vx')
    if( ivx.gt.0 ) then
      ctmp = cdumpauxarr(ivx)
      do i=ivx-1,1,-1
        cdumpauxarr(i+1) = cdumpauxarr(i)
      enddo
      cdumpauxarr(1) = ctmp
    endif
    return
  end subroutine make_cdumpauxarr
!=======================================================================
  function idumpauxof(cauxname) result(idumpaux)
    use pmdvars,only: ndumpaux,cdumpauxarr
    character(len=*),intent(in):: cauxname
    integer:: idumpaux
    integer:: i
    
    idumpaux = -1
    do i=1,ndumpaux
      if( trim(cdumpauxarr(i)).eq.trim(cauxname) ) then
        idumpaux = i
        return
      endif
    enddo
    return
  end function idumpauxof
!=======================================================================
  subroutine calc_nfmv(ntot,tagtot)
    use pmdvars,only: myid_md,mpi_md_world,iprint,nfmv
    include 'mpif.h'
    include './const.h'
    integer,intent(in):: ntot
    real(8),intent(in):: tagtot(ntot)

    integer:: ia,ierr,igrp
!!$    integer,external:: ifmvOf

    if( myid_md.eq.0 ) then
      igrp = 1  ! ifmv is assigned to group #1
      nfmv = 0
      do ia=1,ntot
!!$        nfmv = max(nfmv,ifmvOf(tagtot(ia)))
        nfmv = max(nfmv,ithOf(tagtot(ia),igrp))
      enddo
      if( iprint.ge.ipl_basic ) then
        print *,''
        print '(a,i0)',' Number of ifmvs = ',nfmv
      endif
    endif
!!$    call mpi_bcast(nfmv,1,mpi_integer,0,mpi_md_world,ierr)
    
  end subroutine calc_nfmv
!=======================================================================
  function to_lower(instr) result(outstr)
!-----------------------------------------------------------------------
!   Ref. pg 80, "Upgrading to Fortran 90", by Cooper Redwine,
!   1995 Springer-Verlag, New York.
!-----------------------------------------------------------------------
    character(len=*),intent(in):: instr
    character(len=len(instr)):: outstr
    integer:: i,n

    outstr = instr
    do i=1,len(instr)
!.....Find location of letter in upper case constant string
      n = index(UPPER_CASE,outstr(i:i))
!.....If current substring is a upper case letter, make it lower case
      if( n.ne.0 ) outstr(i:i) = LOWER_CASE(n:n)
    enddo
    return
  end function to_lower
!=======================================================================
  function to_upper(instr) result(outstr)
    character(len=*),intent(in):: instr
    character(len=len(instr)):: outstr
    integer:: i,n

    outstr = instr
    do i=1,len(instr)
!.....Find location of letter in lower case constant string
      n = index(LOWER_CASE,outstr(i:i))
!.....If current substring is a lower case letter, make it upper case
      if( n.ne.0 ) outstr(i:i) = UPPER_CASE(n:n)
    enddo
    return
  end function to_upper
!=======================================================================
  subroutine resize_iarr(iarr, new_size, default_value)
!
!  Resize 1D int array.
!
    integer,allocatable:: iarr(:)
    integer,intent(in):: new_size, default_value
    integer,allocatable:: itemp(:)

    allocate(itemp(new_size))

    if( size(iarr) > 0 ) then
      itemp(:min(size(iarr), new_size)) = iarr(:min(size(iarr),new_size))
    endif
    if( new_size > size(iarr) ) then
      itemp(size(iarr)+1:new_size) = default_value
    endif

    call move_alloc(itemp, iarr)
  end subroutine resize_iarr
!=======================================================================
  subroutine resize_darr(arr, new_size, default_value)
!
!  Resize 1D int array.
!
    real(8),allocatable:: arr(:)
    integer,intent(in):: new_size
    real(8),allocatable:: temp(:)
    real(8):: default_value

    allocate(temp(new_size))

    if( size(arr) > 0 ) then
      temp(:min(size(arr), new_size)) = arr(:min(size(arr),new_size))
    endif
    if( new_size > size(arr) ) then
      temp(size(arr)+1:new_size) = default_value
    endif

    call move_alloc(temp, arr)
  end subroutine resize_darr
!=======================================================================
  function basename(path) result(name)
!
!  Extract basename from full path (by removing strings before "/").
!
    character(len=*), intent(in) :: path
    character(len=:), allocatable :: name
    integer :: p
    character(len=:), allocatable :: tmp

    tmp = trim(path)

    if (len(tmp) > 1 .and. tmp(len(tmp):len(tmp)) == "/") then
      tmp = tmp(:len(tmp)-1)
    end if

    p = index(tmp, "/", back=.true.)

    if (p > 0) then
      name = tmp(p+1:)
    else
      name = tmp
    end if
  end function basename

end module util
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd lib"
!     End:

