module pmdio
!-----------------------------------------------------------------------
!                     Last modified: <2021-11-24 16:05:50 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
  implicit none
  save

contains
!=======================================================================
  function get_ntot_ascii(ionum,cfname) result(ntot)
    integer,intent(in):: ionum
    character(len=*),intent(in):: cfname
    integer:: ntot

    integer:: ia,ib,l,i
    real(8):: hunit,h(3,3,0:1)
    character(len=128):: ctmp 

    open(ionum,file=trim(cfname),status='old')
    do while(.true.)   ! skip comment lines
      read(ionum,'(a)') ctmp
      if( .not. (ctmp(1:1).eq.'!' .or. ctmp(1:1).eq.'#') ) then
        backspace(ionum)
        exit
      endif
    enddo
    read(ionum,*) hunit
    read(ionum,*) (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
    read(ionum,*) ntot
    close(ionum)
    return

  end function get_ntot_ascii
!=======================================================================
  function get_ntot_bin(ionum,cfname) result(ntot)
    integer,intent(in):: ionum
    character(len=*),intent(in):: cfname
    integer:: ntot

    integer:: ia,ib,l,i,msp,naux
    real(8):: h(3,3,0:1),hunit
    character:: ctmp*3

    open(ionum,file=trim(cfname),form='unformatted',status='old')
    read(ionum) msp
    read(ionum) (ctmp,i=1,msp)
    read(ionum) hunit
    read(ionum) (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
    read(ionum) ntot, naux
    close(ionum)
    return

  end function get_ntot_bin
!=======================================================================
  subroutine read_pmdtot_ascii(ionum,cfname,ntot,hunit,h,tagtot, &
       rtot,vtot)
    integer,intent(in):: ionum,ntot
    character(len=*),intent(in):: cfname
    real(8),intent(out):: hunit,h(3,3,0:1)
    real(8),intent(out):: tagtot(ntot),rtot(3,ntot),vtot(3,ntot)

    integer:: ia,ib,l,i,itmp
    character(len=128):: ctmp 

    open(ionum,file=trim(cfname),status='old')
!.....Comment lines at the top of pmdini could contain options.
    do while(.true.)
      read(ionum,'(a)') ctmp
      if( ctmp(1:1).eq.'!' .or. ctmp(1:1).eq.'#' ) then
        call parse_option(ctmp)
      else
        backspace(ionum)
        exit
      endif
    enddo
    read(ionum,*) hunit
    read(ionum,*) (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
    h(1:3,1:3,0:1)= h(1:3,1:3,0:1)*hunit
    read(ionum,*) itmp
    if( itmp.ne.ntot ) then
      print *,' ERROR: itmp.ne.ntot'
      stop
    endif
    do i=1,ntot
      read(ionum,*) tagtot(i),rtot(1:3,i),vtot(1:3,i)
    enddo
    close(ionum)

  end subroutine read_pmdtot_ascii
!=======================================================================
  subroutine write_pmdtot_ascii(ionum,cfname,ntot,hunit,h,tagtot, &
       rtot,vtot)
    use pmdvars,only: has_specorder,specorder
    include './params_unit.h'
    integer,intent(in):: ionum,ntot
    character(len=*),intent(in) :: cfname
    real(8),intent(in):: hunit,h(3,3,0:1)
    real(8),intent(in):: tagtot(ntot),rtot(3,ntot),vtot(3,ntot)

    integer:: ia,ib,l,i,msp

    open(ionum,file=cfname,status='replace')
    if( has_specorder ) then
      msp = 0
      do i=1,ntot
        msp = max(msp,int(tagtot(i)))
      enddo
      write(ionum,'(a)') '!'
      write(ionum,'(a,9(2x,a))') '!  specorder: ',(trim(specorder(i)),i=1,msp)
      write(ionum,'(a)') '!'
    endif
    write(ionum,'(es23.14e3)') hunit
    write(ionum,'(3es23.14e3)') (((h(ia,ib,l)/hunit,ia=1,3) &
         ,ib=1,3),l=0,1)
    write(ionum,'(i10)') ntot
    do i=1,ntot
      write(ionum,'(7es23.14e3,11es13.4e3)') tagtot(i) &
           ,rtot(1:3,i),vtot(1:3,i)    ! dt
    enddo
    close(ionum)
    return
  end subroutine write_pmdtot_ascii
!=======================================================================
  subroutine read_pmdtot_bin(ionum,cfname,ntot,hunit,h,tagtot, &
       rtot,vtot)
    use pmdvars,only: specorder
    integer,intent(in):: ionum,ntot
    character(len=*),intent(in):: cfname
    real(8),intent(out):: hunit,h(3,3,0:1)
    real(8),intent(out):: tagtot(ntot),rtot(3,ntot),vtot(3,ntot)

    integer:: ia,ib,l,i,msp,itmp

    open(ionum,file=trim(cfname),form='unformatted',status='old')
!-----natm: num. of particles in this node
    read(ionum) msp
    read(ionum) (specorder(i),i=1,msp)
    read(ionum) hunit
    read(ionum) (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
    h(1:3,1:3,0:1)= h(1:3,1:3,0:1)*hunit
    read(ionum) itmp
    if( itmp.ne.ntot ) then
      print *,' ERROR: itmp.ne.ntot !'
      stop
    endif
    read(ionum) tagtot(1:ntot)
    read(ionum) rtot(1:3,1:ntot)
    read(ionum) vtot(1:3,1:ntot)
    close(ionum)

  end subroutine read_pmdtot_bin
!=======================================================================
  subroutine write_pmdtot_bin(ionum,cfname,ntot,hunit,h,tagtot, &
       rtot,vtot)
    use pmdvars,only: specorder
    include './params_unit.h'
    integer,intent(in):: ionum,ntot
    character(len=*),intent(in) :: cfname
    real(8),intent(in):: hunit,h(3,3,0:1)
    real(8),intent(in):: tagtot(ntot),rtot(3,ntot),vtot(3,ntot)

    integer:: ia,ib,l,i,msp

    open(ionum,file=cfname,form='unformatted',status='replace')
    msp = 0
    do ia=1,ntot
      msp = max(msp,int(tagtot(ia)))
    enddo
    write(ionum) msp
    write(ionum) (specorder(i),i=1,msp)
    write(ionum) hunit
    write(ionum) (((h(ia,ib,l)/hunit,ia=1,3),ib=1,3),l=0,1)
    write(ionum) ntot
    write(ionum) tagtot(1:ntot)
    write(ionum) rtot(1:3,1:ntot)
    write(ionum) vtot(1:3,1:ntot)
    close(ionum)

  end subroutine write_pmdtot_bin
!=======================================================================
  subroutine write_dump(ionum,cfname,ntot,hunit,h,tagtot,rtot,vtot, &
       stot,ekitot,epitot,naux,auxtot)
!
!     Write atomic configuration in LAMMPS-dump format file.
!
    use pmdvars,only: ndumpaux,cdumpauxarr,specorder,has_specorder
    use util,only: itotOf,iauxof
    use time,only: accum_time
    implicit none
    include "mpif.h"
    integer,intent(in):: ionum,ntot,naux
    character(len=*),intent(in) :: cfname
    real(8),intent(in):: hunit,h(3,3,0:1)
    real(8),intent(in):: tagtot(ntot),rtot(3,ntot),vtot(3,ntot), &
         stot(3,3,ntot),ekitot(3,3,ntot),epitot(ntot),auxtot(naux,ntot)

    integer:: i,j,k,l,is,idlmp
    real(8):: xi(3),ri(3),xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz, &
         xlo_bound,xhi_bound,ylo_bound,yhi_bound, &
         zlo_bound,zhi_bound,st(3,3)
!!$    integer,external:: itotOf
!!$    real(8),allocatable,save:: rlmp(:,:),vlmp(:,:)
    integer,save:: ndlmp
    real(8),allocatable,save:: dlmp(:,:)
    character(len=3),save:: cndlmp
    character(len=3):: csp
    character(len=6):: caux
    real(8),parameter:: tiny = 1d-14
    logical,save:: l1st = .true.

    if( l1st ) then
      ndlmp = 3 +ndumpaux
!.....If vx,vy,vz are included in dumpauxarr, take that into account
!!$      if( idumpauxof('vx').gt.0 ) ndlmp = ndlmp -1
!!$      if( idumpauxof('vy').gt.0 ) ndlmp = ndlmp -1
!!$      if( idumpauxof('vz').gt.0 ) ndlmp = ndlmp -1
      write(cndlmp,'(i0)') ndlmp
      allocate(dlmp(ndlmp,ntot))
      l1st = .false.
    endif

    if( .not.allocated(dlmp) ) then
      allocate(dlmp(ndlmp,ntot))
    else if( size(dlmp).ne.ndlmp*ntot ) then
      deallocate(dlmp)
      allocate(dlmp(ndlmp,ntot))
    endif

    call pmd2lammps(h,ntot,rtot,dlmp(1:3,:),vtot,dlmp(4:6,:) &
         ,xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz)
    xlo_bound = xlo +min(0d0, xy, xz, xy+xz)
    xhi_bound = xhi +max(0d0, xy, xz, xy+xz)
    ylo_bound = ylo +min(0d0, yz)
    yhi_bound = yhi +max(0d0, yz)
    zlo_bound = zlo
    zhi_bound = zhi
!.....Make a data array to be written out
    idlmp = 3  ! skip x,y,z pos data
    do j=1,ndumpaux
      caux = cdumpauxarr(j)
      idlmp = idlmp +1
      if( trim(caux).eq.'ekin'.or.trim(caux).eq.'eki' ) then
        dlmp(idlmp,:) = ekitot(1,1,:) +ekitot(2,2,:) +ekitot(3,3,:)
      else if( trim(caux).eq.'epot'.or.trim(caux).eq.'epi' ) then
        dlmp(idlmp,:) = epitot(:)
      else if( trim(caux).eq.'vx' .or. trim(caux).eq.'vy' .or. trim(caux).eq.'vz' ) then
!.....velocity data are already stored at dlmp(4:6,:) and they must be always at 4:6
        cycle
      else if( trim(caux).eq.'sxx' ) then
        dlmp(idlmp,:) = stot(1,1,:)
      else if( trim(caux).eq.'syy' ) then
        dlmp(idlmp,:) = stot(2,2,:)
      else if( trim(caux).eq.'szz' ) then
        dlmp(idlmp,:) = stot(3,3,:)
      else if( trim(caux).eq.'syz' .or. trim(caux).eq.'szy' ) then
        dlmp(idlmp,:) = stot(2,3,:)
      else if( trim(caux).eq.'sxz' .or. trim(caux).eq.'szx' ) then
        dlmp(idlmp,:) = stot(1,3,:)
      else if( trim(caux).eq.'sxy' .or. trim(caux).eq.'syx' ) then
        dlmp(idlmp,:) = stot(1,2,:)
      else if( trim(caux).eq.'chg' ) then
        dlmp(idlmp,:) = auxtot(iauxof('chg'),:)
      else if( trim(caux).eq.'tei' ) then
        dlmp(idlmp,:) = auxtot(iauxof('tei'),:)
      else if( trim(caux).eq.'clr' ) then
        dlmp(idlmp,:) = auxtot(iauxof('clr'),:)
      endif
    enddo

    open(ionum,file=trim(cfname),status='replace')
    write(ionum,'(a)') 'ITEM: TIMESTEP'
    write(ionum,'(i10)') 0
    write(ionum,'(a)') 'ITEM: NUMBER OF ATOMS'
    write(ionum,'(i10)') ntot
    write(ionum,'(a)') 'ITEM: BOX BOUNDS xy xz yz'
    write(ionum,'(3f15.4)') xlo_bound, xhi_bound, xy
    write(ionum,'(3f15.4)') ylo_bound, yhi_bound, xz
    write(ionum,'(3f15.4)') zlo_bound, zhi_bound, yz
    write(ionum,'(a)',advance='no') 'ITEM: ATOMS id type x y z'
    do i=1,ndumpaux
      write(ionum,'(a)',advance='no') ' '//trim(cdumpauxarr(i))
    enddo
    write(ionum,*) ''
    do i=1,ntot
      write(ionum,'(i8)',advance='no') itotOf(tagtot(i))
      if( has_specorder ) then
        is = int(tagtot(i))
        csp = specorder(is)
        write(ionum,'(a4)',advance='no') trim(csp)
!!$      print *,'tag,i,csp = ',tagtot(i),itotOf(tagtot(i)),csp
      else
        write(ionum,'(i3)',advance='no') int(tagtot(i))
      endif
      write(ionum,'(3f12.5)',advance='no') dlmp(1:3,i)  ! pos
      write(ionum,'('//trim(cndlmp)//'es11.2e3)') dlmp(4:ndlmp,i)  ! except pos
    enddo

    close(ionum)
  end subroutine write_dump
!=======================================================================
  subroutine pmd2lammps(h,ntot,rtot,rlmp,vtot,vlmp, &
       xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz)
!
!     Convert pmd data format to LAMMPS format not only cell
!     but also atomic positions.
!     Parameters to be output:
!       xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz
!     LAMMPS cell should be defined as,
!       a = ( xhi-xlo,       0,       0 )
!       b = (      xy, yhi-hlo,       0 )
!       c = (      xz,      yz, zhi-zlo )
!     See, http://lammps.sandia.gov/doc/Section_howto.html, for detail.
!
    use vector,only: norm,dot,cross
    integer,intent(in):: ntot
    real(8),intent(in):: h(3,3),rtot(3,ntot),vtot(3,ntot)
    real(8),intent(out):: xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz &
         ,rlmp(3,ntot),vlmp(3,ntot)

    integer:: i,lxy,lxz,lyz
    real(8):: a0(3),b0(3),c0(3),a1(3),a2(3),a3(3) &
         ,b1(3),b2(3),b3(3),rt(3),vt(3),amat(3,3),bmat(3,3) &
         ,x,y,z,a23(3),a31(3),a12(3),vol,xyp
    real(8):: a,b,c,alpha,beta,gamma

    xlo = 0d0
    ylo = 0d0
    zlo = 0d0
    a0(1:3) = h(1:3,1)
    b0(1:3) = h(1:3,2)
    c0(1:3) = h(1:3,3)
    a = norm(a0)
    b = norm(b0)
    c = norm(c0)
    alpha = acos(dot(b0,c0)/b/c)
    beta  = acos(dot(a0,c0)/a/c)
    gamma = acos(dot(a0,b0)/a/b)
    xhi = a
    xy = b*cos(gamma)
    xz = c*cos(beta)
    yhi= sqrt(b*b -xy*xy)
    yz = (b*c*cos(alpha) -xy*xz)/yhi
    zhi= sqrt(c*c -xz*xz -yz*yz)
    x = xhi -xlo
    y = yhi -ylo
    z = zhi -zlo
    lxy = 0
    if( xy.gt.xhi/2 ) then
      xy = xy -xhi
      lxy = -1
    else if( xy.lt.-xhi/2 ) then
      xy = xy +xhi
      lxy = 1
    endif
    lxz = 0
    if( xz.gt.xhi/2 ) then
      xz = xz -xhi
      lxz = -1
    else if( xz.lt.-xhi/2 ) then
      xz = xz +xhi
      lxz = 1
    endif
    lyz = 0
    if( yz.gt.yhi/2 ) then
      yz = yz -yhi
      lyz = -1
    else if( yz.lt.-yhi/2 ) then
      yz = yz +yhi
      lyz = 1
    endif

    a1(1:3) = h(1:3,1)
    a2(1:3) = h(1:3,2)
    a3(1:3) = h(1:3,3)
    a23 = cross(a2,a3)
    a31 = cross(a3,a1)
    a12 = cross(a1,a2)
    vol = abs(dot(a1,a23))
    amat(1:3,1:3) = 0d0
    amat(1,1:3) = a23(1:3)
    amat(2,1:3) = a31(1:3)
    amat(3,1:3) = a12(1:3)
    b1(1:3) = (/ x, 0d0, 0d0 /)
    b2(1:3) = (/ xy,  y, 0d0 /)
    b3(1:3) = (/ xz, yz,   z /)
    bmat(1:3,1:3) = 0d0
    bmat(1:3,1) = b1(1:3)
    bmat(1:3,2) = b2(1:3)
    bmat(1:3,3) = b3(1:3)
    xyp = xy -lxy*x
    do i=1,ntot
      rlmp(1:3,i) = 0d0
      call shift_pos_for_lammps(rtot(1,i),rlmp(1,i),lxy,lxz,lyz &
           ,x,y,z,yz,xz,xy)
!.....Velocity as well
      vlmp(1:3,i) = vtot(1:3,i)
      vlmp(2,i) = vlmp(2,i) -lyz*vtot(3,i)
      vlmp(1,i) = vlmp(1,i) -lxz*vtot(3,i) +(vtot(2,i)*xyp -vlmp(2,i)*xy)/x
      rt = matmul(h,rlmp(1:3,i))
      vt = matmul(h,vlmp(1:3,i))
      rlmp(1:3,i) = matmul(bmat,matmul(amat,rt))/vol
      vlmp(1:3,i) = matmul(bmat,matmul(amat,vt))/vol
    enddo

    return
  end subroutine pmd2lammps
!=======================================================================
  subroutine shift_pos_for_lammps(r,rn,lxy,lxz,lyz,x,y,z,yz,xz,xy)
    use pmdvars,only: boundary
    real(8),intent(in):: r(3),x,y,z,yz,xz,xy
    integer,intent(in):: lxy,lxz,lyz
    real(8),intent(out):: rn(3)

    integer:: i
    real(8):: xyp

    xyp = xy -lxy*x
    rn(1:3) = r(1:3)
    rn(2) = rn(2) -lyz*r(3)
    rn(1) = rn(1) -lxz*r(3) +(r(2)*xyp -rn(2)*xy)/x
    do i=1,3
      if( boundary(i:i).eq.'p' ) then
        rn(i) = pbc(rn(i))
      endif
    enddo
    return
  end subroutine shift_pos_for_lammps
!=======================================================================
  function pbc(x)
    real(8),intent(in):: x
    real(8):: pbc

    if( x.lt.0d0 ) then
      pbc = x -int(x) +1d0
    else if( x.ge.1d0 ) then
      pbc = x -int(x)
    else
      pbc = x
    endif
    return
  end function pbc
!=======================================================================
  subroutine hmat2lammps(hmat,xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz)
!
!     Convert h-matrix to LAMMPS cell vectors.
!     Parameters to be output:
!       xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz
!     LAMMPS cell should be defined as,
!       a = ( xhi-xlo,       0,       0 )
!       b = (      xy, yhi-hlo,       0 )
!       c = (      xz,      yz, zhi-zlo )
!     See, http://lammps.sandia.gov/doc/Section_howto.html, for detail.
!
    use vector,only: norm,dot
    real(8),intent(in):: hmat(3,3)
    real(8),intent(out):: xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz

    real(8):: a0(3),b0(3),c0(3)
    real(8):: a,b,c,alpha,beta,gamma

    xlo = 0d0
    ylo = 0d0
    zlo = 0d0
    a0(1:3) = hmat(1:3,1)
    b0(1:3) = hmat(1:3,2)
    c0(1:3) = hmat(1:3,3)
    a = norm(a0)
    b = norm(b0)
    c = norm(c0)
    alpha = acos(dot(b0,c0)/b/c)
    beta  = acos(dot(a0,c0)/a/c)
    gamma = acos(dot(a0,b0)/a/b)
    xhi = a
    xy = b*cos(gamma)
    xz = c*cos(beta)
    yhi= sqrt(b*b -xy*xy)
    yz = (b*c*cos(alpha) -xy*xz)/yhi
    zhi= sqrt(c*c -xz*xz -yz*yz)
    return
  end subroutine hmat2lammps
!=======================================================================
  subroutine parse_option(cline)
!
!  Parse options from a comment line.
!  Lines starting from ! or # are treated as comment lines,
!  and options can be given at the comment lines.
!  The option words should be put after these comment characters with
!  one or more spaces between them for example,
!
!  specorder: Al Mg Si
!
!  Currently available options are:
!    - "specorder:", Species order. The number of species limited up to 9.
!
    use pmdvars,only: specorder,has_specorder,iprint
    use util, only: num_data
    include "./const.h"
    character(len=*),intent(in):: cline

    integer:: iopt1,isp,num
    real(8):: opt1, opt2
    character(len=10):: c1,copt
    logical:: lopt

    if( index(cline,'specorder:').ne.0 ) then
      num = num_data(trim(cline),' ')
      if( num.gt.11 ) stop 'ERROR: number of species exceeds the limit.'
      read(cline,*) c1, copt, specorder(1:num-2)
      if( iprint.ge.ipl_basic ) then
        print '(a)',' Species order read from pmdini option: '
        do isp=1,num-2
          print '(i5,": ",a4)',isp,trim(specorder(isp))
        enddo
      endif
      has_specorder = .true.
    endif
    
  end subroutine parse_option
!=======================================================================
  subroutine split_pair(strin,str1,str2)
!
!  Split the input string of a pair connected by hyphen, e.g.) Si-O,
!  as separated strings, str1=Si, str=O.
!
    character(len=*),intent(in):: strin
    character(len=*),intent(out):: str1,str2

    character(len=1),parameter:: delim='-'

    integer:: l,ind

    l = len_trim(strin)
    ind = index(strin,delim)
    str1 = strin(1:ind-1)
    str2 = strin(ind+1:l)
    return
  end subroutine split_pair
!=======================================================================
end module pmdio
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
