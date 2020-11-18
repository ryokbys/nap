module pmdio
!-----------------------------------------------------------------------
!                     Last modified: <2020-11-18 22:16:35 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
  implicit none
  save

  character(len=20),parameter:: cinpmd='in.pmd'

  integer:: ntot0,ntot,naux
!.....data of total system
  real(8),allocatable:: rtot(:,:),vtot(:,:),stot(:,:,:),epitot(:) &
       ,ekitot(:,:,:),tagtot(:),atot(:,:)
  real(8),allocatable:: chgtot(:),chitot(:),teitot(:),clrtot(:)
  real(8),allocatable:: auxtot(:,:)
  real(8):: hunit,h(3,3,0:1)

!.....max. num. of atoms in a node
  integer:: namax = 20000
!.....max. num. of boundary-particles
  integer:: nbmax = 10000
!.....max. num. of neighbors
  integer:: nnmax = 100

  integer:: nstp = 0
  integer:: minstp = 0
  integer:: nerg = 1000
  integer:: ifpmd= 2    ! 0:none, 1:pmd, 2:dump (default)
  integer:: npmd = 10
  integer:: ifsort= 1
  real(8):: dt = 1d0
  real(8):: vardt_len = 0.1d0  ! Length criterion for variable time-step
  real(8):: rc = 5.0d0
  real(8):: rc1nn = 2.5d0
  real(8):: rbuf= 0d0
  integer:: ifdmp= 0 ! 0:none, 1:damped-MD, 2:FIRE
  character(len=20):: cmin= ''
  real(8):: dmp  = 0.9d0
  real(8):: eps_conv = 1d-4
  integer:: n_conv = 1
!.....temperature
  character(len=20):: ctctl='none'
  integer:: iftctl= 0
  real(8):: tinit= -1d0
  real(8):: tfin = -1d0
  real(8):: ttgt(9)
  data ttgt / 300d0, 300d0, 300d0, 300d0, 300d0, 300d0, &
       300d0, 300d0, 300d0 /
  real(8):: trlx = 100d0
!.....Random seed
!  If positve, use (RSEED+MYID) as the seed for each process
!  If negative, use the same random seeds for all the parallel process
  real(8):: rseed = 12345d0
!.....Remove translational motion:
!     N< 0: not to remove translation
!     N==0: remove translation only at the beginning
!     N> 1: remove translation at the begining and every N step.
  integer:: nrmtrans = 0
!.....Coulomb system?
  integer:: ifcoulomb = 0
!.....temperature distribution on x-direction
  logical:: ltdst= .false.
  integer:: ntdst= 1
!.....shear stress
  real(8):: shrst = 0.0d0
!.....factors on each moving direction
  real(8):: fmv(3,0:9)
  data fmv &
       / 0d0, 0d0, 0d0, & ! 0
       1d0, 1d0, 1d0, & ! 1
       1d0, 1d0, 1d0, & ! 2
       1d0, 1d0, 1d0, & ! 3
       1d0, 1d0, 1d0, & ! 4
       1d0, 1d0, 1d0, & ! 5
       1d0, 1d0, 1d0, & ! 6
       1d0, 1d0, 1d0, & ! 7
       1d0, 1d0, 1d0, & ! 8
       1d0, 1d0, 1d0  & ! 9
       /
!.....whether compute stress or not
  logical:: lstrs0 = .true.
!.....barostat
  character(len=20):: cpctl='none'
  real(8):: ptgt   = 0d0
  real(8):: pini   = 0d0
  real(8):: pfin   = 0d0
  real(8):: srlx   = 100d0
  real(8):: stbeta = 1d0
  real(8):: strfin = 0.0d0
  real(8):: stnsr(3,3)
  real(8):: stgt(1:3,1:3)= 0d0
  logical:: lcellfix(1:3,1:3)= .false.
!.....charge optimize or variable charge
  logical:: lvc = .false.
!.....Charge setting
  character(len=20):: chgfix='input'
!.....print level
!  0:quiet, 1:normal,
!  >10:fitpot data
  integer:: iprint= 1

!.....Auxiliary data order for dump output
  logical:: ldumpaux_changed = .false.
  character(len=256):: cdumpaux = 'ekin epot sxx syy szz syz sxz sxy'
  integer:: ndumpaux
  character(len=6),allocatable:: cdumpauxarr(:)

  character(len=128):: cpmdini = 'pmdini'
  character(len=128):: cpmdfin = 'pmdfin'
  character(len=6):: ciofmt='ascii '
  character(len=20):: cforce='none'
!!$  integer:: numff = 0 ! number of force-fields
!!$  character(len=20),allocatable:: cffs(:)  ! force-fields
!.....max. num. of species
  integer,parameter:: nspmax= 9
!.....mass
  real(8):: am(1:nspmax)= 12.0d0
!.....charges
  real(8):: schg(1:nspmax)= 0d0
!.....species name
  character(len=3):: specorder(nspmax) = 'x'
  logical:: has_specorder = .false.

!.....Boundary condition: p = periodic, f = free, w = wall
  character(len=3):: boundary = 'ppp'

!.....PKA for radiation damage
  integer:: iatom_pka = -1
  real(8):: pka_energy = -1.d0 ! in eV
  real(8):: pka_theta = 0.d0  ! in degree
  real(8):: pka_phi = 0.d0    ! in degree

!.....Metadynamics
  logical:: lmetaD = .false. 
!.....Constraints
  logical:: lconst = .false. 
!.....Reduced force
  logical:: lrdcfrc = .false. 

!.....zload type: zload or shear
  character(len=128):: czload_type= 'none'
!.....top and bottom skin width in which atoms are fixed and/or controlled
  real(8):: zskin_width = 5.0d0
!.....Shear angle from x in degree, shear direction is on xy-plane
  real(8):: zshear_angle = 0d0

!.....Deformation
  character(len=128):: cdeform= 'none'
  real(8):: dhratio(3,3)

!.....Structure analysis: CNA, a-CNA
  character(len=128):: cstruct = 'none'
  integer:: istruct = 1

contains
!=======================================================================
  subroutine read_pmdtot_ascii(ionum,cfname)
    implicit none
    integer,intent(in):: ionum
    character(len=*),intent(in):: cfname

    integer:: ia,ib,l,i
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
    read(ionum,*) ntot0
    ntot = ntot0
    allocate(tagtot(ntot0),rtot(3,ntot0),vtot(3,ntot0),epitot(ntot0) &
         ,ekitot(3,3,ntot0),stot(3,3,ntot0),atot(3,ntot0))
    do i=1,ntot0
!!$      read(ionum,*) tagtot(i),rtot(1:3,i),vtot(1:3,i) &
!!$           ,ekitot(1,1,i),epitot(i) &
!!$           ,stot(1,1,i),stot(2,2,i),stot(3,3,i) &
!!$           ,stot(2,3,i),stot(3,1,i),stot(1,2,i)
      read(ionum,*) tagtot(i),rtot(1:3,i),vtot(1:3,i)
    enddo
    close(ionum)

  end subroutine read_pmdtot_ascii
!=======================================================================
  subroutine write_pmdtot_ascii(ionum,cfname)
    implicit none
    include './params_unit.h'
    integer,intent(in):: ionum
    character(len=*),intent(in) :: cfname

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
!!$      write(ionum,'(7es23.14e3,11es13.4e3)') tagtot(i) &
!!$           ,rtot(1:3,i) &
!!$           ,vtot(1:3,i) & !/dt
!!$           ,ekitot(1,1,i)+ekitot(2,2,i)+ekitot(3,3,i) &
!!$           ,epitot(i) &
!!$           ,stot(1,1,i),stot(2,2,i),stot(3,3,i) &
!!$           ,stot(2,3,i),stot(3,1,i),stot(1,2,i)
      write(ionum,'(7es23.14e3,11es13.4e3)') tagtot(i) &
           ,rtot(1:3,i),vtot(1:3,i)    ! dt
    enddo
    close(ionum)

  end subroutine write_pmdtot_ascii
!=======================================================================
  subroutine read_pmdtot_bin(ionum,cfname)
    implicit none
    integer,intent(in):: ionum
    character(len=*),intent(in):: cfname

    integer:: ia,ib,l,i

    open(ionum,file=trim(cfname),form='unformatted',status='old')
!-----natm: num. of particles in this node
    read(ionum) hunit
    read(ionum) (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
    h(1:3,1:3,0:1)= h(1:3,1:3,0:1)*hunit
    read(ionum) ntot0
    ntot = ntot0
    allocate(tagtot(ntot0),rtot(3,ntot0),atot(3,ntot0) &
         ,vtot(3,ntot0),epitot(ntot0) &
         ,ekitot(3,3,ntot0),stot(3,3,ntot0))
!!$    read(ionum) (tagtot(i),rtot(1:3,i),vtot(1:3,i) &
!!$         ,ekitot(1,1,i),epitot(i) &
!!$         ,stot(1,1,i),stot(2,2,i),stot(3,3,i) &
!!$         ,stot(2,3,i),stot(3,1,i),stot(1,2,i) &
!!$         ,i=1,ntot0)
    read(ionum) (tagtot(i),rtot(1:3,i),vtot(1:3,i) ,i=1,ntot0)
    close(ionum)

  end subroutine read_pmdtot_bin
!=======================================================================
  subroutine write_pmdtot_bin(ionum,cfname)
    implicit none
    include './params_unit.h'
    integer,intent(in):: ionum
    character(len=*),intent(in) :: cfname

    integer:: ia,ib,l,i

    open(ionum,file=cfname,form='unformatted' &
         ,status='replace')
    write(ionum) hunit
    write(ionum) (((h(ia,ib,l)/hunit,ia=1,3),ib=1,3),l=0,1)
    write(ionum) ntot
    do i=1,ntot
!!$      write(ionum) tagtot(i),rtot(1:3,i),vtot(1:3,i) & !/dt
!!$           ,ekitot(1,1,i)+ekitot(2,2,i)+ekitot(3,3,i) &
!!$           ,epitot(i) &
!!$           ,stot(1,1,i),stot(2,2,i),stot(3,3,i) &
!!$           ,stot(2,3,i),stot(3,1,i),stot(1,2,i)
      write(ionum) tagtot(i),rtot(1:3,i),vtot(1:3,i)
    enddo
    close(ionum)

  end subroutine write_pmdtot_bin
!=======================================================================
  subroutine write_dump(ionum,cfname)
!
!     Write atomic configuration in LAMMPS-dump format file.
!
    use util,only: itotOf
    implicit none 
    integer,intent(in):: ionum
    character(len=*),intent(in) :: cfname

    integer:: i,j,k,l,is
    real(8):: xi(3),ri(3),eki,epi,xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz, &
         xlo_bound,xhi_bound,ylo_bound,yhi_bound, &
         zlo_bound,zhi_bound,st(3,3)
!!$    integer,external:: itotOf
    real(8),allocatable,save:: rlmp(:,:),vlmp(:,:)
    character(len=3):: csp
    character(len=6):: caux

    real(8),parameter:: tiny = 1d-14

    if( .not. allocated(rlmp) ) allocate(rlmp(3,ntot),vlmp(3,ntot))
    if( size(rlmp).ne.3*ntot ) then
      deallocate(rlmp,vlmp)
      allocate(rlmp(3,ntot),vlmp(3,ntot))
    endif

    open(ionum,file=trim(cfname),status='replace')
    write(ionum,'(a)') 'ITEM: TIMESTEP'
    write(ionum,'(i10)') 0
    write(ionum,'(a)') 'ITEM: NUMBER OF ATOMS'
    write(ionum,'(i10)') ntot
    write(ionum,'(a)') 'ITEM: BOX BOUNDS xy xz yz'
    call pmd2lammps(h,ntot,rtot,rlmp,vtot,vlmp &
         ,xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz)
    xlo_bound = xlo +min(0d0, xy, xz, xy+xz)
    xhi_bound = xhi +max(0d0, xy, xz, xy+xz)
    ylo_bound = ylo +min(0d0, yz)
    yhi_bound = yhi +max(0d0, yz)
    zlo_bound = zlo
    zhi_bound = zhi
    write(ionum,'(3f15.4)') xlo_bound, xhi_bound, xy
    write(ionum,'(3f15.4)') ylo_bound, yhi_bound, xz
    write(ionum,'(3f15.4)') zlo_bound, zhi_bound, yz
    write(ionum,'(a)',advance='no') 'ITEM: ATOMS id type x y z'
    do i=1,ndumpaux
      write(ionum,'(a)',advance='no') ' '//trim(cdumpauxarr(i))
    enddo
    write(ionum,*) ''
    do i=1,ntot
!        xi(1:3)= rtot(1:3,i)
!        ri(1:3)= h(1:3,1,0)*xi(1) +h(1:3,2,0)*xi(2) +h(1:3,3,0)*xi(3)
      eki = ekitot(1,1,i) +ekitot(2,2,i) +ekitot(3,3,i)
      epi = epitot(i)
      st(1:3,1:3) = stot(1:3,1:3,i)
      if( eki.lt.tiny ) eki = 0d0
      if( abs(epi).lt.tiny ) epi = 0d0
      do l=1,3
        do k=l,3
          if( abs(st(l,k)).lt.tiny ) st(l,k) = 0d0
        enddo
      enddo
      write(ionum,'(i8)',advance='no') itotOf(tagtot(i))
      if( has_specorder ) then
        is = int(tagtot(i))
        csp = specorder(is)
        write(ionum,'(a4)',advance='no') trim(csp)
      else
        write(ionum,'(i3)',advance='no') int(tagtot(i))
      endif
      write(ionum,'(3f12.5)',advance='no') rlmp(1:3,i)
      if( ldumpaux_changed ) then
        do j=1,ndumpaux
          caux = cdumpauxarr(j)
          if( trim(caux).eq.'ekin'.or.trim(caux).eq.'eki' ) then
            write(ionum,'(es11.3)',advance='no') eki
          else if( trim(caux).eq.'epot'.or.trim(caux).eq.'epi' ) then
            write(ionum,'(es11.3)',advance='no') epi
          else if( trim(caux).eq.'vx' ) then
            write(ionum,'(es11.3)',advance='no') vlmp(1,i)
          else if( trim(caux).eq.'vy' ) then
            write(ionum,'(es11.3)',advance='no') vlmp(2,i)
          else if( trim(caux).eq.'vz' ) then
            write(ionum,'(es11.3)',advance='no') vlmp(3,i)
          else if( trim(caux).eq.'sxx' ) then
            write(ionum,'(es11.3)',advance='no') st(1,1)
          else if( trim(caux).eq.'syy' ) then
            write(ionum,'(es11.3)',advance='no') st(2,2)
          else if( trim(caux).eq.'szz' ) then
            write(ionum,'(es11.3)',advance='no') st(3,3)
          else if( trim(caux).eq.'syz' .or. trim(caux).eq.'szy' ) then
            write(ionum,'(es11.3)',advance='no') st(2,3)
          else if( trim(caux).eq.'sxz' .or. trim(caux).eq.'szx' ) then
            write(ionum,'(es11.3)',advance='no') st(1,3)
          else if( trim(caux).eq.'sxy' .or. trim(caux).eq.'syx' ) then
            write(ionum,'(es11.3)',advance='no') st(1,2)
          else if( trim(caux).eq.'chg' ) then
            write(ionum,'(f9.4)',advance='no') chgtot(i)
          else if( trim(caux).eq.'chi' ) then
            write(ionum,'(f9.2)',advance='no') chitot(i)
          else if( trim(caux).eq.'tei' ) then
            write(ionum,'(es11.3)',advance='no') teitot(i)
          else if( trim(caux).eq.'clr' ) then
            write(ionum,'(es11.3)',advance='no') clrtot(i)
          endif
        enddo
        write(ionum,*) ''
      else
        write(ionum,'(8es11.3)') eki,epi, &
             st(1,1),st(2,2),st(3,3), &
             st(2,3),st(1,3),st(1,2)
      end if
!!$      if( has_specorder ) then
!!$        is = int(tagtot(i))
!!$        csp = specorder(is)
!!$        write(ionum,'(i8,a4,3f12.5,8es11.3,f9.4,f9.2)') &
!!$             itotOf(tagtot(i)),trim(csp),rlmp(1:3,i),eki, &
!!$             epi, &
!!$             st(1,1),st(2,2),st(3,3), &
!!$             st(2,3),st(1,3),st(1,2), &
!!$             chgtot(i),chitot(i)
!!$      else
!!$        write(ionum,'(i8,i3,3f12.5,8es11.3,f9.4,f9.2)') &
!!$             itotOf(tagtot(i)),int(tagtot(i)),rlmp(1:3,i),eki, &
!!$             epi, &
!!$             st(1,1),st(2,2),st(3,3), &
!!$             st(2,3),st(1,3),st(1,2), &
!!$             chgtot(i),chitot(i)
!!$      endif
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
    implicit none
    integer,intent(in):: ntot
    real(8),intent(in):: h(3,3),rtot(3,ntot),vtot(3,ntot)
    real(8),intent(out):: xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz &
         ,rlmp(3,ntot),vlmp(3,ntot)

    integer:: i,lxy,lxz,lyz
    real(8):: a0(3),b0(3),c0(3),a1(3),a2(3),a3(3) &
         ,b1(3),b2(3),b3(3),rt(3),vt(3),amat(3,3),bmat(3,3) &
         ,x,y,z,a23(3),a31(3),a12(3),vol,xyp
    real(8):: a,b,c,alpha,beta,gamma
    real(8),external:: absv,sprod

    xlo = 0d0
    ylo = 0d0
    zlo = 0d0
    a0(1:3) = h(1:3,1)
    b0(1:3) = h(1:3,2)
    c0(1:3) = h(1:3,3)
    a = absv(3,a0)
    b = absv(3,b0)
    c = absv(3,c0)
    alpha = acos(sprod(3,b0,c0)/b/c)
    beta  = acos(sprod(3,a0,c0)/a/c)
    gamma = acos(sprod(3,a0,b0)/a/b)
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
    call vprod(a2,a3,a23)
    call vprod(a3,a1,a31)
    call vprod(a1,a2,a12)
    vol = abs( sprod(3,a1,a23) )
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
    implicit none
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
    implicit none
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
    implicit none
    real(8),intent(in):: hmat(3,3)
    real(8),intent(out):: xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz

    real(8):: a0(3),b0(3),c0(3)
    real(8):: a,b,c,alpha,beta,gamma
    real(8),external:: absv,sprod

    xlo = 0d0
    ylo = 0d0
    zlo = 0d0
    a0(1:3) = hmat(1:3,1)
    b0(1:3) = hmat(1:3,2)
    c0(1:3) = hmat(1:3,3)
    a = absv(3,a0)
    b = absv(3,b0)
    c = absv(3,c0)
    alpha = acos(sprod(3,b0,c0)/b/c)
    beta  = acos(sprod(3,a0,c0)/a/c)
    gamma = acos(sprod(3,a0,b0)/a/b)
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
    use util, only: num_data
    implicit none
    character(len=*),intent(in):: cline

    integer:: iopt1,isp,num
    real(8):: opt1, opt2
    character(len=10):: c1,copt
    logical:: lopt

    if( index(cline,'specorder:').ne.0 ) then
      num = num_data(trim(cline),' ')
      if( num.gt.11 ) stop 'ERROR: number of species exceeds the limit.'
      read(cline,*) c1, copt, specorder(1:num-2)
      if( iprint.gt.0 ) then
        print '(a)',' Species order read from pmdini option: '
        do isp=1,num-2
          print '(i5,": ",a4)',isp,trim(specorder(isp))
        enddo
      endif
      has_specorder = .true.
    endif
    
  end subroutine parse_option
!=======================================================================
  function csp2isp(csp,spcs)
!
!  Convert cspi to isp.
!  If not found, return -1.
!
    character(len=3),intent(in):: spcs(nspmax)
    character(len=*),intent(in):: csp
    integer:: csp2isp

    integer:: isp

    csp2isp = -1
    do isp=1,nspmax
      if( trim(csp).eq.trim(spcs(isp)) ) then
        csp2isp = isp
        return
      endif
    enddo
    return
  end function csp2isp
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
  subroutine make_cdumpauxarr()
!
!  Builld cdumpauxarr if ifpmd==2 (dump output)
!
    use util,only: num_data

    integer:: i

    ndumpaux = num_data(trim(cdumpaux),' ')
    if( allocated(cdumpauxarr) ) deallocate(cdumpauxarr)
    allocate(cdumpauxarr(ndumpaux))
    read(cdumpaux,*) (cdumpauxarr(i),i=1,ndumpaux)
    return
  end subroutine make_cdumpauxarr
!=======================================================================
  function get_vol(h)
    real(8),intent(in):: h(3,3)

    integer:: i,j,jm,jp,im,ip
    real(8):: sgm(3,3)
    real(8):: get_vol

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
    get_vol=h(1,1)*sgm(1,1)+h(2,1)*sgm(2,1)+h(3,1)*sgm(3,1)
    return
  end function get_vol
end module pmdio
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
