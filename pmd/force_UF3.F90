module UF3
!-----------------------------------------------------------------------
!                     Last modified: <2025-05-06 21:36:01 KOBAYASHI Ryo>
!-----------------------------------------------------------------------
!  Parallel implementation of Ultra-Fast Force-Field (UF3) for pmd
!    - 2024.09.02 by R.K., start to implement
!    - 2025.04.00 by R.K., implemented UF3L
!    - 2026.02.05 by R.K., start implementing UF3D
!-----------------------------------------------------------------------
  use pmdvars,only: nspmax
  use util,only: csp2isp, num_data
  use memory,only: accum_mem
  use vector,only: dot
  implicit none
  include "mpif.h"
  include './const.h'
  include "./params_unit.h"
  save

  character(len=128):: paramsdir = '.'
!.....parameter file name
  character(128),parameter:: cpfname= 'in.params.uf3'
  character(128),parameter:: cpfname_l= 'in.params.uf3l'
  character(128),parameter:: cpfname_d= 'in.params.uf3d'
  integer,parameter:: ioprms = 50

  logical:: lprms_read_uf3 = .false.
  logical:: lprms_read_uf3l = .false.
  logical:: lprms_read_uf3d = .false.
!!$  logical:: lprmset_uf3 = .false.
!!$  logical:: lprmset_uf3l = .false.
!!$  logical:: lprmset_uf3d = .false.

!.....uf2 parameters
  type prm2
    character(2):: cb, csi, csj, cknot
    integer:: isp,jsp
!  cb: NA, 2B or 3B
!  csi,csj,csk: species name
!  cknot: nk (non-uniform knot spacing) or uk (uniform knot spacing)
    integer:: nklead, nktrail
    integer:: nknot, ncoef
    real(8):: rc,rc2
    real(8),allocatable:: knots(:), coefs(:)
    real(8),allocatable:: gwe(:), gwf(:,:,:), gws(:,:)
  end type prm2

!.....uf3 parameters
  type prm3
    character(2):: cb, csi, csj, csk, cknot
    integer:: isp, jsp, ksp
!  cb: NA, 2B or 3B
!  csi,csj,csk: species name
!  cknot: nk (non-uniform knot spacing) or uk (uniform knot spacing)
    integer:: nklead, nktrail
    integer:: nknij, nknik, nknjk, ncfij, ncfik, ncfjk
    real(8):: rcij, rcik, rcjk, rcij2, rcik2, rcjk2
    real(8),allocatable:: knij(:), knik(:), knjk(:), coefs(:,:,:)
    real(8),allocatable:: gwe(:,:,:), gwf(:,:,:,:,:), gws(:,:,:,:)
  end type prm3

!.....uf3l 3B parameters
  type prm3l
    character(2):: cb, csi, csj, csk, cknot
    integer:: isp, jsp, ksp
!  cb: NA, 2B or 3B
!  csi,csj,csk: species name
!  cknot: nk (non-uniform knot spacing) or uk (uniform knot spacing)
    integer:: nklead, nktrail
    integer:: nknot, ncoef
    real(8):: rcij, rcik, rcij2, rcik2, gmj, gmk
    real(8),allocatable:: knots(:), coefs(:)
    real(8),allocatable:: gwe(:), gwf(:,:,:), gws(:,:)
  end type prm3l

!.....uf3d 3B parameters
  type prm3d
    character(2):: cb, csi, csj, csk, cknot
    integer:: isp, jsp, ksp
!  cb: NA, 2B or 3B
!  csi,csj,csk: species name
!  cknot: nk (non-uniform knot spacing) or uk (uniform knot spacing)
    integer:: nklead, nktrail
    integer:: nknij, nknik, nkncs, ncfij, ncfik, ncfcs
    real(8):: rcij, rcik, rcij2, rcik2
    real(8),allocatable:: knij(:), knik(:), kncs(:), cfij(:), cfik(:), cfcs(:)
    real(8),allocatable:: gwe(:), gwf(:,:,:), gws(:,:)
  end type prm3d

  integer:: n1b, n2b, n3b
  integer:: ncoef = 0
  real(8):: erg1s(nspmax), gerg1s(nspmax)
  integer,allocatable:: prm1s(:)
  type(prm2),allocatable:: prm2s(:)
  type(prm3),allocatable:: prm3s(:)
  type(prm3l),allocatable:: prm3ls(:)
  type(prm3d),allocatable:: prm3ds(:)
  logical:: has_trios = .false.
  logical:: has_solo = .false.
  real(8):: rcmax = 0.0d0
  real(8):: rc3max = 0.0d0
  real(8):: rc3max2 = 0.0d0

  real(8),allocatable:: aal2(:,:),aal3(:,:),strsl(:,:,:)
  integer,allocatable:: ls3b(:)

!.....Map of pairs (trios) to parameter set id
  integer:: interact2(nspmax,nspmax), interact3(nspmax,nspmax,nspmax)
!!$!.....Cutoffs
!!$  real(8):: rc2_3b(nspmax,nspmax)

!.....constants
  integer:: nelem,nexp,nsp
  integer,parameter:: ivoigt(3,3)= &
       reshape((/ 1, 6, 5,  6, 2, 4,  5, 4, 3 /),shape(ivoigt))

contains
  subroutine read_params_uf3(myid,mpi_world,iprint)
!
!  Read parameters of uf3 potential from in.params.uf3 file that is given
!  by uf3/lammps_plugin/scripts/generate_uf3_lammps_pots.py,
!  but it can accept 1-body energy as well.
!-----------------------------------------------------------------------
!  #UF3 POT UNTIS: metal DATE: 2024-09-12 17:12:36 AUTHOR: RK CITATION:
!  1B  W  erg
!  #
!  #UF3 POT UNITS: metal DATE: 2024-09-12 17:12:36 AUTHOR: RK CITATION:
!  2B W W 0 3 nk
!  5.5  22
!  0.001 0.001 0.001 0.001 0.3676 0.7342 ... 5.1334 5.5 5.5 5.5 5.5
!  18
!  36.82 36.82 26.69 ... -0.032 0 0 0
!  #
!  #UF3 POT UNITS: metal DATE: 2024-09-12 17:12:36 AUTHOR: RK CITATION:
!  3B W W W 0 3 nk
!  7.0 3.5 3.5 19 13 13
!  1.5 1.5 1.5 1.5 1.96 2.42 ... 6.54 7 7 7 7
!  1.5 1.5 1.5 1.5 1.83 2.17 ... 3.17 3.5 3.5 3.5 3.5
!  1.5 1.5 1.5 1.5 1.83 2.17 ... 3.17 3.5 3.5 3.5 3.5
!  9 9 15
!  0.0 0.0 0.0 0.0 ... 0.0 0.0 0.0
!  ...
!  #
!-----------------------------------------------------------------------
!  Each pair or trio is sandwitched with lines beginning with "#".
!  For trios, in the case of knot information the order is jk,ik,ij,
!  but for coefficient information the order is ij,ik,jk.
!  And there is a limitation in the current lammps implementation version
!  such that r_max_jk = 2*r_max_ij = 2_r_max_ik.
!
    implicit none
    integer,intent(in):: myid,mpi_world,iprint

    integer:: itmp,ierr,i,j,i1b,i2b,i3b,isp
    integer:: nklead, nktrail
    real(8):: etmp
!  nklead, nktrail: num of leading or trailing knots
    logical:: lexist
    character:: fname*128, cmode*4, cb*2, csi*2, csj*2, csk*2, &
         cknot*2, ctmp*128, cline*128
!  cmode: none or read

    if( myid == 0 ) then
      if( iprint >= ipl_basic ) print '(/,a)',' Read UF3 parameters...'
      fname = trim(paramsdir)//'/'//trim(cpfname)
!.....read parameters at the 1st call
      inquire(file=trim(fname),exist=lexist)
      if( .not. lexist ) then
        write(6,'(a)') '   [Error] '//trim(fname)//' does not exist !!!.'
!!$        call mpi_finalize(ierr)
        stop
      endif
      cmode = 'none'
      n1b = 0
      n2b = 0
      n3b = 0
      erg1s(:) = 0d0
      interact2(:,:) = -1
      interact3(:,:,:) = -1
      open(ioprms,file=trim(fname),status='old')
!.....1st, count number of 2B and 3B entries to allocate type objects
      do while(.true.)
        read(ioprms,'(a)',end=10) cline
        if( cline(1:1) == '#' ) then
          if( cline(1:3) == '#UF' ) then
            cmode = 'read'
            cb = 'NA'
          else
            cmode = 'none'
          endif
          cycle
        endif
        if( trim(cmode) == 'none' ) cycle
        if( cb == 'NA' ) then  ! if cb is not read yet
          read(cline,*,iostat=ierr) cb
          if( cb == '2B' ) then
            n2b = n2b + 1
          else if( cb == '3B' ) then
            n3b = n3b + 1
            has_trios = .true.
          else if( cb == '1B' ) then
            n1b = n1b + 1
            has_solo = .true.
          endif
! do nothing if cb is already set...
        endif
!        if( ierr /= 0 ) cycle
      end do  ! finished counting 2B & 3B entries
10    continue
      if( iprint >= ipl_debug ) then
        print *,'  n1b, n2b, n3b = ',n1b,n2b,n3b
        print *,'  has_solo= ',has_solo
        print *,'  has_trios= ',has_trios
      endif
!.....allocate lists of uf2prms and uf3prms
      if( allocated(prm1s) ) deallocate(prm1s)
      if( allocated(prm2s) ) deallocate(prm2s)
      if( allocated(prm3s) ) deallocate(prm3s)
      if( n1b /= 0 ) allocate(prm1s(n1b))
      if( n2b /= 0 ) allocate(prm2s(n2b))
      if( n3b /= 0 ) allocate(prm3s(n3b))
      rewind(ioprms)
      i1b = 0
      i2b = 0
      i3b = 0
      do while(.true.)
        read(ioprms,'(a)',end=20) cline
        if( cline(1:1) == '#') then
          if( cline(1:3) == '#UF') then
            cmode = 'read'
            cb = 'NA'
          else
            cmode = 'none'
          endif
          cycle
        endif
        if( trim(cmode) == 'none' ) cycle
        if( cb == 'NA' ) then  ! if cb is not read yet
          read(cline,*,iostat=ierr) cb
          if( cb == '2B' ) then
            backspace(ioprms)
            i2b = i2b +1
            call read_2b(prm2s(i2b),i2b)
            if( iprint >= ipl_basic ) call print_2b(prm2s(i2b))
          else if( cb == '3B' ) then
            backspace(ioprms)
            i3b = i3b +1
            call read_3b(prm3s(i3b),i3b)
            if( iprint >= ipl_basic ) call print_3b(prm3s(i3b))
          else if( cb == '1B' ) then
            i1b = i1b +1
            backspace(ioprms)
            read(ioprms,*) ctmp, csi, etmp
            ncoef = ncoef +1
            isp = csp2isp(csi)
            prm1s(i1b) = isp
            erg1s(isp) = etmp
            if( iprint >= ipl_basic ) then
              print '(a,f10.3)','   UF3 parameters of 1B for '//trim(csi) &
                   //': erg1 = ',etmp
            endif
          endif
        endif
!        if( ierr /= 0 ) cycle

      enddo  ! while(.true.)
20    continue  ! when the file reached the end
      if( iprint >= ipl_basic ) print '(/,a,i0)', &
           '   Total num of UF3 coefficients =  ',ncoef
    endif

    call bcast_uf3_params(mpi_world,myid)
    lprms_read_uf3 = .true.
    return
  end subroutine read_params_uf3
!=======================================================================
  subroutine read_params_uf3l(myid,mpi_world,iprint)
!
!  Read parameters of uf3l potential from in.params.uf3l file.
!  And the 3B format is much different from that of original uf3.
!-----------------------------------------------------------------------
!  #UF3 POT UNTIS: metal DATE: 2024-09-12 17:12:36 AUTHOR: RK CITATION:
!  1B  W  erg
!  #
!  #UF3 POT UNITS: metal DATE: 2024-09-12 17:12:36 AUTHOR: RK CITATION:
!  2B W W 0 3 nk
!  5.5  22
!  0.1 0.1 0.1 0.1 0.3676 0.7342 ... 5.1334 5.5 5.5 5.5 5.5
!  18
!  36.82 36.82 26.69 ... -0.032 0 0 0
!  #
!  #UF3 POT UNITS: metal DATE: 2024-09-12 17:12:36 AUTHOR: RK CITATION:
!  3B W W W 0 3 nk
!  3.5  3.0  19  1.0  1.0
!  -1.0 -1.0 -1.0 -1.0 -0.82 ... 0.82 1.0 1.0 1.0 1.0
!  15
!  2.3  2.3  2.1 ... 0.2 0 0 0
!  #
!-----------------------------------------------------------------------
!  Each pair or trio is sandwitched with lines beginning with "#".
!  For trios, the value range spans from -1.0 to 1.0 corresponding to cos.
!  The line in 3B "3.5  3.0  19  1.0  1.0" means that
!    - two 3B cutoffs 3.5 for ij and 3.0 for ik pairs
!    - 19 is the number of nodes
!    - following "1.0 1.0" are gamma for ij and for ik pairs.
!
    implicit none
    integer,intent(in):: myid,mpi_world,iprint

    integer:: itmp,ierr,i,j,i1b,i2b,i3b,isp
    integer:: nklead, nktrail
    real(8):: etmp
!  nklead, nktrail: num of leading or trailing knots
    logical:: lexist
    character:: fname*128, cmode*4, cb*2, csi*2, csj*2, csk*2, &
         cknot*2, ctmp*128, cline*128
!  cmode: none or read

    if( myid == 0 ) then
      if( iprint >= ipl_basic ) print '(/,a)',' Read UF3L parameters...'
      fname = trim(paramsdir)//'/'//trim(cpfname_l)
!.....read parameters at the 1st call
      inquire(file=trim(fname),exist=lexist)
      if( .not. lexist ) then
        write(6,'(a)') '   [Error] '//trim(fname)//' does not exist !!!.'
!!$        call mpi_finalize(ierr)
        stop
      endif
      cmode = 'none'
      n1b = 0
      n2b = 0
      n3b = 0
      erg1s(:) = 0d0
      interact2(:,:) = -1
      interact3(:,:,:) = -1
      open(ioprms,file=trim(fname),status='old')
!.....1st, count number of 2B and 3B entries to allocate type objects
      do while(.true.)
        read(ioprms,'(a)',end=10) cline
        if( cline(1:1) == '#' ) then
          if( cline(1:3) == '#UF' ) then
            cmode = 'read'
            cb = 'NA'
          else
            cmode = 'none'
          endif
          cycle
        endif
        if( trim(cmode) == 'none' ) cycle
        if( cb == 'NA' ) then  ! if cb is not read yet
          read(cline,*,iostat=ierr) cb
          if( cb == '2B' ) then
            n2b = n2b + 1
          else if( cb == '3B' ) then
            n3b = n3b + 1
            has_trios = .true.
          else if( cb == '1B' ) then
            n1b = n1b + 1
            has_solo = .true.
          endif
! do nothing if cb is already set...
        endif
!        if( ierr /= 0 ) cycle
      end do  ! finished counting 2B & 3B entries
10    continue
      if( iprint >= ipl_debug ) then
        print *,'  n1b, n2b, n3b = ',n1b,n2b,n3b
        print *,'  has_solo= ',has_solo
        print *,'  has_trios= ',has_trios
      endif
!.....allocate lists of uf2prms and uf3prms
      if( allocated(prm1s) ) deallocate(prm1s)
      if( allocated(prm2s) ) deallocate(prm2s)
      if( allocated(prm3ls) ) deallocate(prm3ls)
      if( n1b /= 0 ) allocate(prm1s(n1b))
      if( n2b /= 0 ) allocate(prm2s(n2b))
      if( n3b /= 0 ) allocate(prm3ls(n3b))
      rewind(ioprms)
      i1b = 0
      i2b = 0
      i3b = 0
      do while(.true.)
        read(ioprms,'(a)',end=20) cline
        if( cline(1:1) == '#') then
          if( cline(1:3) == '#UF') then
            cmode = 'read'
            cb = 'NA'
          else
            cmode = 'none'
          endif
          cycle
        endif
        if( trim(cmode) == 'none' ) cycle
        if( cb == 'NA' ) then  ! if cb is not read yet
          read(cline,*,iostat=ierr) cb
          if( cb == '2B' ) then
            backspace(ioprms)
            i2b = i2b +1
            call read_2b(prm2s(i2b),i2b)
            if( iprint >= ipl_basic ) call print_2b(prm2s(i2b))
          else if( cb == '3B' ) then
            backspace(ioprms)
            i3b = i3b +1
            call read_3bl(prm3ls(i3b),i3b)
            if( iprint >= ipl_basic ) call print_3bl(prm3ls(i3b))
          else if( cb == '1B' ) then
            i1b = i1b +1
            backspace(ioprms)
            read(ioprms,*) ctmp, csi, etmp
            ncoef = ncoef +1
            isp = csp2isp(csi)
            prm1s(i1b) = isp
            erg1s(isp) = etmp
            if( iprint >= ipl_basic ) then
              print '(a,f10.3)','   UF3 parameters of 1B for '//trim(csi) &
                   //': erg1 = ',etmp
            endif
          endif
        endif
!        if( ierr /= 0 ) cycle

      enddo  ! while(.true.)
20    continue  ! when the file reached the end
      if( iprint >= ipl_basic ) print '(/,a,i0)', &
           '   Total num of UF3L coefficients =  ',ncoef
    endif

    call bcast_uf3l_params(mpi_world,myid)
    lprms_read_uf3l = .true.
    return
  end subroutine read_params_uf3l
!=======================================================================
  subroutine read_params_uf3d(myid,mpi_world,iprint)
!
!  Read parameters of uf3d potential from in.params.uf3d file.
!  3B format is different from UF3 and UF3L,
!  rather in between these two.
!-----------------------------------------------------------------------
!  #UF3 POT DATE: 2024-09-12 17:12:36 AUTHOR: RK CITATION:
!  1B  W  erg
!  #
!  #UF3 POT DATE: 2024-09-12 17:12:36 AUTHOR: RK CITATION:
!  2B W W 0 3 nk
!  5.5  22
!  0.1 0.1 0.1 0.1 0.3676 0.7342 ... 5.1334 5.5 5.5 5.5 5.5
!  18
!  36.82 36.82 26.69 ... -0.032 0 0 0
!  #
!  #UF3 POT DATE: 2024-09-12 17:12:36 AUTHOR: RK CITATION:
!  3B W W W 0 3 nk
!  3.5  3.0  14  12  19
!  0.1 0.1 0.1 0.1 0.xxxx ... 0.yyyy 3.5 3.5 3.5 3.5
!  0.1 0.1 0.1 0.1 0.zzzz ... 0.wwww 3.0 3.0 3.0 3.0
!  -1.0 -1.0 -1.0 -1.0 -0.82 ... 0.82 1.0 1.0 1.0 1.0
!  10  8  15
!  10.3  8.3  2.0 ... 0.1 0 0 0
!  10.0  9.0  3.0 ... 0.1 0 0 0
!   2.3  2.3  2.1 ... 0.2 0 0 0
!  #
!-----------------------------------------------------------------------
!  Each pair or trio is sandwitched with lines beginning with "#".
!  For trios, the value range spans from -1.0 to 1.0 corresponding to cos.
!  The line in 3B "3.5  3.0  14  12  10" means that
!    - two 3B cutoffs 3.5 for ij and 3.0 for ik pairs
!    - 14 and 12 are the numbers of knots along ij and ik pairs
!    - following "10" is the number of knots for cos (-1.0 to 1.0)
!
    implicit none
    integer,intent(in):: myid,mpi_world,iprint

    integer:: itmp,ierr,i,j,i1b,i2b,i3b,isp
    integer:: nklead, nktrail
    real(8):: etmp
!  nklead, nktrail: num of leading or trailing knots
    logical:: lexist
    character:: fname*128, cmode*4, cb*2, csi*2, csj*2, csk*2, &
         cknot*2, ctmp*128, cline*128
!  cmode: none or read

    if( myid == 0 ) then
      if( iprint >= ipl_basic ) print '(/,a)',' Read UF3D parameters...'
      fname = trim(paramsdir)//'/'//trim(cpfname_d)
!.....read parameters at the 1st call
      inquire(file=trim(fname),exist=lexist)
      if( .not. lexist ) then
        write(6,'(a)') '   [Error] '//trim(fname)//' does not exist !!!.'
!!$        call mpi_finalize(ierr)
        stop
      endif
      cmode = 'none'
      n1b = 0
      n2b = 0
      n3b = 0
      erg1s(:) = 0d0
      interact2(:,:) = -1
      interact3(:,:,:) = -1
      open(ioprms,file=trim(fname),status='old')
!.....1st, count number of 2B and 3B entries to allocate type objects
      do while(.true.)
        read(ioprms,'(a)',end=10) cline
        if( cline(1:1) == '#' ) then
          if( cline(1:3) == '#UF' ) then
            cmode = 'read'
            cb = 'NA'
          else
            cmode = 'none'
          endif
          cycle
        endif
        if( trim(cmode) == 'none' ) cycle
        if( cb == 'NA' ) then  ! if cb is not read yet
          read(cline,*,iostat=ierr) cb
          if( cb == '2B' ) then
            n2b = n2b + 1
          else if( cb == '3B' ) then
            n3b = n3b + 1
            has_trios = .true.
          else if( cb == '1B' ) then
            n1b = n1b + 1
            has_solo = .true.
          endif
! do nothing if cb is already set...
        endif
!        if( ierr /= 0 ) cycle
      end do  ! finished counting 2B & 3B entries
10    continue
      if( iprint >= ipl_debug ) then
        print *,'  n1b, n2b, n3b = ',n1b,n2b,n3b
        print *,'  has_solo= ',has_solo
        print *,'  has_trios= ',has_trios
      endif
!.....allocate lists of uf2prms and uf3prms
      if( allocated(prm1s) ) deallocate(prm1s)
      if( allocated(prm2s) ) deallocate(prm2s)
      if( allocated(prm3ds) ) deallocate(prm3ds)
      if( n1b /= 0 ) allocate(prm1s(n1b))
      if( n2b /= 0 ) allocate(prm2s(n2b))
      if( n3b /= 0 ) allocate(prm3ds(n3b))
      rewind(ioprms)
      i1b = 0
      i2b = 0
      i3b = 0
      do while(.true.)
        read(ioprms,'(a)',end=20) cline
        if( cline(1:1) == '#') then
          if( cline(1:3) == '#UF') then
            cmode = 'read'
            cb = 'NA'
          else
            cmode = 'none'
          endif
          cycle
        endif
        if( trim(cmode) == 'none' ) cycle
        if( cb == 'NA' ) then  ! if cb is not read yet
          read(cline,*,iostat=ierr) cb
          if( cb == '2B' ) then
            backspace(ioprms)
            i2b = i2b +1
            call read_2b(prm2s(i2b),i2b)
            if( iprint >= ipl_basic ) call print_2b(prm2s(i2b))
          else if( cb == '3B' ) then
            backspace(ioprms)
            i3b = i3b +1
            call read_3bd(prm3ds(i3b),i3b)
            if( iprint >= ipl_basic ) call print_3bd(prm3ds(i3b))
          else if( cb == '1B' ) then
            i1b = i1b +1
            backspace(ioprms)
            read(ioprms,*) ctmp, csi, etmp
            ncoef = ncoef +1
            isp = csp2isp(csi)
            prm1s(i1b) = isp
            erg1s(isp) = etmp
            if( iprint >= ipl_basic ) then
              print '(a,f10.3)','   UF3 parameters of 1B for '//trim(csi) &
                   //': erg1 = ',etmp
            endif
          endif
        endif
!        if( ierr /= 0 ) cycle

      enddo  ! while(.true.)
20    continue  ! when the file reached the end
      if( iprint >= ipl_basic ) print '(/,a,i0)', &
           '   Total num of UF3d coefficients =  ',ncoef
    endif

    call bcast_uf3d_params(mpi_world,myid)
    lprms_read_uf3d = .true.
    return
  end subroutine read_params_uf3d
!=======================================================================
  subroutine read_2b(ps,i2b)
!
!  Read 2B part of file via ioprms.
!  Assuming that the starting line of IOPRMS is the 1st line after
!  the line starting with '#UF3'.
!
    type(prm2),intent(out):: ps
    integer,intent(in):: i2b
    integer:: i, isp, jsp

    read(ioprms,*) ps%cb, ps%csi, ps%csj, ps%nklead, ps%nktrail, ps%cknot
    if( ps%cb /= '2B' ) stop 'ERROR@read_2b: CB should be 2B.'
    if( ps%nktrail /= 3 ) stop 'ERROR@read_2b: nktrail must be 2 !'
    read(ioprms,*) ps%rc, ps%nknot
    ps%rc2 = ps%rc**2
    if( allocated(ps%knots) ) deallocate(ps%knots)
    allocate(ps%knots(ps%nknot))
    read(ioprms,*) (ps%knots(i), i=1,ps%nknot)
    read(ioprms,*) ps%ncoef
    if( allocated(ps%coefs) ) deallocate(ps%coefs)
    allocate(ps%coefs(ps%ncoef))
    read(ioprms,*) (ps%coefs(i), i=1,ps%ncoef)

    ncoef = ncoef + ps%ncoef
    isp = csp2isp(ps%csi)
    jsp = csp2isp(ps%csj)
    ps%isp = isp
    ps%jsp = jsp
    interact2(isp,jsp) = i2b
    interact2(jsp,isp) = i2b
    rcmax = max(ps%rc,rcmax)
  end subroutine read_2b
!=======================================================================
  subroutine read_3b(ps,i3b)
!
!  Read 3B part of file via ioprms
!  Assuming that the starting line of IOPRMS is the 1st line after
!  the line starting with '#UF3'.
!
    type(prm3),intent(out):: ps
    integer,intent(in):: i3b
    integer:: i,j,k,isp,jsp,ksp

    read(ioprms,*) ps%cb, ps%csi, ps%csj, ps%csk, ps%nklead, ps%nktrail, ps%cknot
    if( ps%cb /= '3B' ) stop 'ERROR: CB should be 3B.'
    read(ioprms,*) ps%rcjk, ps%rcij, ps%rcik, ps%nknjk, ps%nknij, ps%nknik
    ps%rcjk2 = ps%rcjk**2
    ps%rcik2 = ps%rcik**2
    ps%rcij2 = ps%rcij**2
    if( allocated(ps%knij) ) deallocate(ps%knij, ps%knik, ps%knjk)
    allocate(ps%knij(ps%nknij), ps%knik(ps%nknik), ps%knjk(ps%nknjk))
    read(ioprms,*) (ps%knjk(i), i=1,ps%nknjk)
    read(ioprms,*) (ps%knij(i), i=1,ps%nknij)
    read(ioprms,*) (ps%knik(i), i=1,ps%nknik)
    read(ioprms,*) ps%ncfij, ps%ncfik, ps%ncfjk
    if( allocated(ps%coefs) ) deallocate(ps%coefs)
    allocate(ps%coefs(ps%ncfjk, ps%ncfik, ps%ncfij))
    do i=1, ps%ncfij
      do j=1, ps%ncfik
        read(ioprms,*) (ps%coefs(k,j,i), k=1,ps%ncfjk)
      enddo
    enddo

    ncoef = ncoef +ps%ncfij *ps%ncfik *ps%ncfjk
    isp = csp2isp(ps%csi)
    jsp = csp2isp(ps%csj)
    ksp = csp2isp(ps%csk)
    ps%isp = isp
    ps%jsp = jsp
    ps%ksp = ksp
    interact3(isp,jsp,ksp) = i3b
    interact3(isp,ksp,jsp) = i3b
    rc3max = max(ps%rcij,rc3max)
    rc3max = max(ps%rcik,rc3max)
    rcmax  = max(rcmax, rc3max)
    rc3max2 = rc3max**2
!!$    rc2_3b(isp,jsp) = max(ps%rcij2, rc2_3b(isp,jsp))
!!$    rc2_3b(isp,ksp) = max(ps%rcik2, rc2_3b(isp,ksp))
!!$    rc2_3b(jsp,isp) = rc2_3b(isp,jsp)
!!$    rc2_3b(ksp,isp) = rc2_3b(isp,ksp)
  end subroutine read_3b
!=======================================================================
  subroutine read_3bl(ps,i3b)
!
!  Read 3B part of file via ioprms for uf3l potential.
!  Assuming that the starting line of IOPRMS is the 1st line after
!  the line starting with '#UF3'.
!
    type(prm3l),intent(out):: ps
    integer,intent(in):: i3b
    integer:: i,j,k,isp,jsp,ksp

    read(ioprms,*) ps%cb, ps%csi, ps%csj, ps%csk, ps%nklead, ps%nktrail, ps%cknot
    if( ps%cb /= '3B' ) stop 'ERROR@read_3bl: CB should be 3B.'
    read(ioprms,*) ps%rcij, ps%rcik, ps%nknot, ps%gmj, ps%gmk
    ps%rcij2 = ps%rcij**2
    ps%rcik2 = ps%rcik**2
    if( ps%csj.eq.ps%csk .and. abs(ps%gmj-ps%gmk).gt.1d-10 ) &
         stop 'ERROR@read_3bl: must be gmj==gmk for csj==csk in 3BL.'
    if( ps%csj.eq.ps%csk .and. abs(ps%rcij-ps%rcik).gt.1d-10 ) &
         stop 'ERROR@read_3bl: must be rcij==rcik for csj==csk in 3BL.'
    if( allocated(ps%knots) ) deallocate(ps%knots)
    allocate(ps%knots(ps%nknot))
    read(ioprms,*) (ps%knots(i), i=1,ps%nknot)
    read(ioprms,*) ps%ncoef
    if( allocated(ps%coefs) ) deallocate(ps%coefs)
    allocate(ps%coefs(ps%ncoef))
    read(ioprms,*) (ps%coefs(i), i=1,ps%ncoef)

    ncoef = ncoef +ps%ncoef
    isp = csp2isp(ps%csi)
    jsp = csp2isp(ps%csj)
    ksp = csp2isp(ps%csk)
    ps%isp = isp
    ps%jsp = jsp
    ps%ksp = ksp
    interact3(isp,jsp,ksp) = i3b
    interact3(isp,ksp,jsp) = i3b
    rc3max = max(ps%rcij, ps%rcik, rc3max)
    rcmax  = max(rcmax, rc3max)
    rc3max2 = rc3max**2
!!$    rc2_3b(isp,jsp) = max(ps%rc2, rc2_3b(isp,jsp))
!!$    rc2_3b(isp,ksp) = max(ps%rc2, rc2_3b(isp,ksp))
!!$    rc2_3b(jsp,isp) = rc2_3b(isp,jsp)
!!$    rc2_3b(ksp,isp) = rc2_3b(isp,ksp)
  end subroutine read_3bl
!=======================================================================
  subroutine read_3bd(ps,i3b)
!
!  Read 3B part of file via ioprms for uf3d potential.
!  Assuming that the starting line of IOPRMS is the 1st line after
!  the line starting with '#UF3'.
!
    type(prm3d),intent(out):: ps
    integer,intent(in):: i3b
    integer:: i,j,k,isp,jsp,ksp

    read(ioprms,*) ps%cb, ps%csi, ps%csj, ps%csk, ps%nklead, ps%nktrail, ps%cknot
    if( ps%cb /= '3B' ) stop 'ERROR@read_3bd: CB should be 3B.'
    read(ioprms,*) ps%rcij, ps%rcik, ps%nknij, ps%nknik, ps%nkncs
    ps%rcij2 = ps%rcij**2
    ps%rcik2 = ps%rcik**2
    if( ps%csj.eq.ps%csk .and. abs(ps%rcij-ps%rcik).gt.1d-10 ) &
         stop 'ERROR@read_3bd: must be rcij==rcik for csj==csk in 3BD.'
    if( ps%csj.eq.ps%csk .and. ps%nknij.ne.ps%nknik ) &
         stop 'ERROR@read_3bd: must be nknij==nknik for csj==csk in 3BD.'
    if( allocated(ps%knij) ) deallocate(ps%knij, ps%knik, ps%kncs)
    allocate(ps%knij(ps%nknij), ps%knik(ps%nknik), ps%kncs(ps%nkncs))
    read(ioprms,*) (ps%knij(i), i=1,ps%nknij)
    read(ioprms,*) (ps%knik(i), i=1,ps%nknik)
    read(ioprms,*) (ps%kncs(i), i=1,ps%nkncs)
    read(ioprms,*) ps%ncfij, ps%ncfik, ps%ncfcs
    if( allocated(ps%cfij) ) deallocate(ps%cfij, ps%cfik, ps%cfcs)
    allocate(ps%cfij(ps%ncfij), ps%cfik(ps%ncfik), ps%cfcs(ps%ncfcs))
    read(ioprms,*) (ps%cfij(i), i=1,ps%ncfij)
    read(ioprms,*) (ps%cfik(i), i=1,ps%ncfik)
    read(ioprms,*) (ps%cfcs(i), i=1,ps%ncfcs)

    ncoef = ncoef +ps%ncfij +ps%ncfik +ps%ncfcs
    isp = csp2isp(ps%csi)
    jsp = csp2isp(ps%csj)
    ksp = csp2isp(ps%csk)
    ps%isp = isp
    ps%jsp = jsp
    ps%ksp = ksp
    interact3(isp,jsp,ksp) = i3b
    interact3(isp,ksp,jsp) = i3b
    rc3max = max(ps%rcij, ps%rcik, rc3max)
    rcmax  = max(rcmax, rc3max)
    rc3max2 = rc3max**2
!!$    rc2_3b(isp,jsp) = max(ps%rc2, rc2_3b(isp,jsp))
!!$    rc2_3b(isp,ksp) = max(ps%rc2, rc2_3b(isp,ksp))
!!$    rc2_3b(jsp,isp) = rc2_3b(isp,jsp)
!!$    rc2_3b(ksp,isp) = rc2_3b(isp,ksp)
  end subroutine read_3bd
!=======================================================================
  subroutine bcast_uf3_params(mpi_world,myid)
!
!  Broadcast 2B & 3B parameters.
!
    integer,intent(in):: mpi_world, myid
    integer:: i, i2b, i3b, i1b, ierr
    character(2):: p1
    type(prm2):: p2
    type(prm3):: p3

    call mpi_bcast(n2b, 1, mpi_integer, 0, mpi_world, ierr)
    call mpi_bcast(n3b, 1, mpi_integer, 0, mpi_world, ierr)
    call mpi_bcast(n1b, 1, mpi_integer, 0, mpi_world, ierr)
    if( .not. allocated(prm2s) ) allocate(prm2s(n2b))
    call mpi_bcast(ncoef, 1, mpi_integer, 0, mpi_world, ierr)

    do i2b=1,n2b
      call mpi_bcast(prm2s(i2b)%cb,2,mpi_character,0,mpi_world,ierr)
      call mpi_bcast(prm2s(i2b)%csi,2,mpi_character,0,mpi_world,ierr)
      call mpi_bcast(prm2s(i2b)%csj,2,mpi_character,0,mpi_world,ierr)
      call mpi_bcast(prm2s(i2b)%isp,1,mpi_integer,0,mpi_world,ierr)
      call mpi_bcast(prm2s(i2b)%jsp,1,mpi_integer,0,mpi_world,ierr)
      call mpi_bcast(prm2s(i2b)%cknot,2,mpi_character,0,mpi_world,ierr)
      call mpi_bcast(prm2s(i2b)%nklead,1,mpi_integer,0,mpi_world,ierr)
      call mpi_bcast(prm2s(i2b)%nktrail,1,mpi_integer,0,mpi_world,ierr)
      call mpi_bcast(prm2s(i2b)%nknot,1,mpi_integer,0,mpi_world,ierr)
      call mpi_bcast(prm2s(i2b)%ncoef,1,mpi_integer,0,mpi_world,ierr)
      call mpi_bcast(prm2s(i2b)%rc,1,mpi_real8,0,mpi_world,ierr)
      if( .not.allocated(prm2s(i2b)%knots) ) &
           allocate(prm2s(i2b)%knots(prm2s(i2b)%nknot))
      if( .not.allocated(prm2s(i2b)%coefs) ) &
           allocate(prm2s(i2b)%coefs(prm2s(i2b)%ncoef))
      call mpi_bcast(prm2s(i2b)%knots,prm2s(i2b)%nknot,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(prm2s(i2b)%coefs,prm2s(i2b)%ncoef,mpi_real8,0,mpi_world,ierr)
    enddo

    call mpi_bcast(has_trios, 1, mpi_logical, 0,mpi_world,ierr)
    if( has_trios ) then
      if( .not. allocated(prm3s) ) allocate(prm3s(n3b))
      do i3b=1,n3b
        call mpi_bcast(prm3s(i3b)%cb,2,mpi_character,0,mpi_world,ierr)
        call mpi_bcast(prm3s(i3b)%csi,2,mpi_character,0,mpi_world,ierr)
        call mpi_bcast(prm3s(i3b)%csj,2,mpi_character,0,mpi_world,ierr)
        call mpi_bcast(prm3s(i3b)%csk,2,mpi_character,0,mpi_world,ierr)
        call mpi_bcast(prm3s(i3b)%isp,1,mpi_integer,0,mpi_world,ierr)
        call mpi_bcast(prm3s(i3b)%jsp,1,mpi_integer,0,mpi_world,ierr)
        call mpi_bcast(prm3s(i3b)%ksp,1,mpi_integer,0,mpi_world,ierr)
        call mpi_bcast(prm3s(i3b)%cknot,2,mpi_character,0,mpi_world,ierr)
        call mpi_bcast(prm3s(i3b)%nklead,1,mpi_integer,0,mpi_world,ierr)
        call mpi_bcast(prm3s(i3b)%nktrail,1,mpi_integer,0,mpi_world,ierr)
        call mpi_bcast(prm3s(i3b)%nknij,1,mpi_integer,0,mpi_world,ierr)
        call mpi_bcast(prm3s(i3b)%nknik,1,mpi_integer,0,mpi_world,ierr)
        call mpi_bcast(prm3s(i3b)%nknjk,1,mpi_integer,0,mpi_world,ierr)
        call mpi_bcast(prm3s(i3b)%ncfij,1,mpi_integer,0,mpi_world,ierr)
        call mpi_bcast(prm3s(i3b)%ncfik,1,mpi_integer,0,mpi_world,ierr)
        call mpi_bcast(prm3s(i3b)%ncfjk,1,mpi_integer,0,mpi_world,ierr)
        call mpi_bcast(prm3s(i3b)%rcij,1,mpi_real8,0,mpi_world,ierr)
        call mpi_bcast(prm3s(i3b)%rcik,1,mpi_real8,0,mpi_world,ierr)
        call mpi_bcast(prm3s(i3b)%rcjk,1,mpi_real8,0,mpi_world,ierr)
        if( .not.allocated(prm3s(i3b)%knij) ) &
             allocate(prm3s(i3b)%knij(prm3s(i3b)%nknij))
        if( .not.allocated(prm3s(i3b)%knik) ) &
             allocate(prm3s(i3b)%knik(prm3s(i3b)%nknik))
        if( .not.allocated(prm3s(i3b)%knjk) ) &
             allocate(prm3s(i3b)%knjk(prm3s(i3b)%nknjk))
        if( .not.allocated(prm3s(i3b)%coefs)) &
             allocate(prm3s(i3b)%coefs(prm3s(i3b)%ncfjk, &
             prm3s(i3b)%ncfik, prm3s(i3b)%ncfij))
        call mpi_bcast(prm3s(i3b)%knij,prm3s(i3b)%nknij, &
             mpi_real8,0,mpi_world,ierr)
        call mpi_bcast(prm3s(i3b)%knik,prm3s(i3b)%nknik, &
             mpi_real8,0,mpi_world,ierr)
        call mpi_bcast(prm3s(i3b)%knjk,prm3s(i3b)%nknjk, &
             mpi_real8,0,mpi_world,ierr)
        call mpi_bcast(prm3s(i3b)%coefs,&
             prm3s(i3b)%ncfjk*prm3s(i3b)%ncfik*prm3s(i3b)%ncfij,&
             mpi_real8,0,mpi_world,ierr)
      enddo
    endif  ! has_trios

    call mpi_bcast(has_solo, 1, mpi_logical, 0,mpi_world,ierr)
    if( has_solo ) then
      if( .not. allocated(prm1s) ) allocate(prm1s(n1b))
      call mpi_bcast(prm1s,n1b,mpi_integer,0,mpi_world,ierr)
      call mpi_bcast(erg1s,nspmax,mpi_real8,0,mpi_world,ierr)
    endif  ! has_solo

    call mpi_bcast(interact2, nspmax**2, mpi_integer, 0,mpi_world,ierr)
    call mpi_bcast(interact3, nspmax**3, mpi_integer, 0,mpi_world,ierr)

  end subroutine bcast_uf3_params
!=======================================================================
  subroutine bcast_uf3l_params(mpi_world,myid)
!
!  Broadcast 2B & 3B parameters.
!
    integer,intent(in):: mpi_world, myid
    integer:: i, i2b, i3b, i1b, ierr
    character(2):: p1
    type(prm2):: p2
    type(prm3l):: p3

    call mpi_bcast(n2b, 1, mpi_integer, 0, mpi_world, ierr)
    call mpi_bcast(n3b, 1, mpi_integer, 0, mpi_world, ierr)
    call mpi_bcast(n1b, 1, mpi_integer, 0, mpi_world, ierr)
    if( .not. allocated(prm2s) ) allocate(prm2s(n2b))
    call mpi_bcast(ncoef, 1, mpi_integer, 0, mpi_world, ierr)

    do i2b=1,n2b
      call mpi_bcast(prm2s(i2b)%cb,2,mpi_character,0,mpi_world,ierr)
      call mpi_bcast(prm2s(i2b)%csi,2,mpi_character,0,mpi_world,ierr)
      call mpi_bcast(prm2s(i2b)%csj,2,mpi_character,0,mpi_world,ierr)
      call mpi_bcast(prm2s(i2b)%isp,1,mpi_integer,0,mpi_world,ierr)
      call mpi_bcast(prm2s(i2b)%jsp,1,mpi_integer,0,mpi_world,ierr)
      call mpi_bcast(prm2s(i2b)%cknot,2,mpi_character,0,mpi_world,ierr)
      call mpi_bcast(prm2s(i2b)%nklead,1,mpi_integer,0,mpi_world,ierr)
      call mpi_bcast(prm2s(i2b)%nktrail,1,mpi_integer,0,mpi_world,ierr)
      call mpi_bcast(prm2s(i2b)%nknot,1,mpi_integer,0,mpi_world,ierr)
      call mpi_bcast(prm2s(i2b)%ncoef,1,mpi_integer,0,mpi_world,ierr)
      call mpi_bcast(prm2s(i2b)%rc,1,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(prm2s(i2b)%rc2,1,mpi_real8,0,mpi_world,ierr)
      if( .not.allocated(prm2s(i2b)%knots) ) &
           allocate(prm2s(i2b)%knots(prm2s(i2b)%nknot))
      if( .not.allocated(prm2s(i2b)%coefs) ) &
           allocate(prm2s(i2b)%coefs(prm2s(i2b)%ncoef))
      call mpi_bcast(prm2s(i2b)%knots,prm2s(i2b)%nknot,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(prm2s(i2b)%coefs,prm2s(i2b)%ncoef,mpi_real8,0,mpi_world,ierr)
    enddo

    call mpi_bcast(has_trios, 1, mpi_logical, 0,mpi_world,ierr)
    if( has_trios ) then
      call mpi_bcast(rc3max,1,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(rc3max2,1,mpi_real8,0,mpi_world,ierr)
      if( .not. allocated(prm3ls) ) allocate(prm3ls(n3b))
      do i3b=1,n3b
        call mpi_bcast(prm3ls(i3b)%cb,2,mpi_character,0,mpi_world,ierr)
        call mpi_bcast(prm3ls(i3b)%csi,2,mpi_character,0,mpi_world,ierr)
        call mpi_bcast(prm3ls(i3b)%csj,2,mpi_character,0,mpi_world,ierr)
        call mpi_bcast(prm3ls(i3b)%csk,2,mpi_character,0,mpi_world,ierr)
        call mpi_bcast(prm3ls(i3b)%cknot,2,mpi_character,0,mpi_world,ierr)
        call mpi_bcast(prm3ls(i3b)%isp,1,mpi_integer,0,mpi_world,ierr)
        call mpi_bcast(prm3ls(i3b)%jsp,1,mpi_integer,0,mpi_world,ierr)
        call mpi_bcast(prm3ls(i3b)%ksp,1,mpi_integer,0,mpi_world,ierr)
        call mpi_bcast(prm3ls(i3b)%nklead,1,mpi_integer,0,mpi_world,ierr)
        call mpi_bcast(prm3ls(i3b)%nktrail,1,mpi_integer,0,mpi_world,ierr)
        call mpi_bcast(prm3ls(i3b)%nknot,1,mpi_integer,0,mpi_world,ierr)
        call mpi_bcast(prm3ls(i3b)%ncoef,1,mpi_integer,0,mpi_world,ierr)
        call mpi_bcast(prm3ls(i3b)%rcij,1,mpi_real8,0,mpi_world,ierr)
        call mpi_bcast(prm3ls(i3b)%rcik,1,mpi_real8,0,mpi_world,ierr)
        call mpi_bcast(prm3ls(i3b)%rcij2,1,mpi_real8,0,mpi_world,ierr)
        call mpi_bcast(prm3ls(i3b)%rcik2,1,mpi_real8,0,mpi_world,ierr)
        call mpi_bcast(prm3ls(i3b)%gmj,1,mpi_real8,0,mpi_world,ierr)
        call mpi_bcast(prm3ls(i3b)%gmk,1,mpi_real8,0,mpi_world,ierr)
        if( .not.allocated(prm3ls(i3b)%knots) ) &
             allocate(prm3ls(i3b)%knots(prm3ls(i3b)%nknot))
        if( .not.allocated(prm3ls(i3b)%coefs)) &
             allocate(prm3ls(i3b)%coefs(prm3ls(i3b)%ncoef))
        call mpi_bcast(prm3ls(i3b)%knots,prm3ls(i3b)%nknot, &
             mpi_real8,0,mpi_world,ierr)
        call mpi_bcast(prm3ls(i3b)%coefs,prm3ls(i3b)%ncoef, &
             mpi_real8,0,mpi_world,ierr)
      enddo
    endif  ! has_trios

    call mpi_bcast(has_solo, 1, mpi_logical, 0,mpi_world,ierr)
    if( has_solo ) then
      if( .not. allocated(prm1s) ) allocate(prm1s(n1b))
      call mpi_bcast(prm1s,n1b,mpi_integer,0,mpi_world,ierr)
      call mpi_bcast(erg1s,nspmax,mpi_real8,0,mpi_world,ierr)
    endif  ! has_solo

    call mpi_bcast(interact2, nspmax**2, mpi_integer, 0,mpi_world,ierr)
    call mpi_bcast(interact3, nspmax**3, mpi_integer, 0,mpi_world,ierr)

  end subroutine bcast_uf3l_params
!=======================================================================
  subroutine bcast_uf3d_params(mpi_world,myid)
!
!  Broadcast 2B & 3B parameters.
!
    integer,intent(in):: mpi_world, myid
    integer:: i, i2b, i3b, i1b, ierr
    character(2):: p1
    type(prm2):: p2
    type(prm3d):: p3

    call mpi_bcast(n2b, 1, mpi_integer, 0, mpi_world, ierr)
    call mpi_bcast(n3b, 1, mpi_integer, 0, mpi_world, ierr)
    call mpi_bcast(n1b, 1, mpi_integer, 0, mpi_world, ierr)
    if( .not. allocated(prm2s) ) allocate(prm2s(n2b))
    call mpi_bcast(ncoef, 1, mpi_integer, 0, mpi_world, ierr)

    do i2b=1,n2b
      call mpi_bcast(prm2s(i2b)%cb,2,mpi_character,0,mpi_world,ierr)
      call mpi_bcast(prm2s(i2b)%csi,2,mpi_character,0,mpi_world,ierr)
      call mpi_bcast(prm2s(i2b)%csj,2,mpi_character,0,mpi_world,ierr)
      call mpi_bcast(prm2s(i2b)%isp,1,mpi_integer,0,mpi_world,ierr)
      call mpi_bcast(prm2s(i2b)%jsp,1,mpi_integer,0,mpi_world,ierr)
      call mpi_bcast(prm2s(i2b)%cknot,2,mpi_character,0,mpi_world,ierr)
      call mpi_bcast(prm2s(i2b)%nklead,1,mpi_integer,0,mpi_world,ierr)
      call mpi_bcast(prm2s(i2b)%nktrail,1,mpi_integer,0,mpi_world,ierr)
      call mpi_bcast(prm2s(i2b)%nknot,1,mpi_integer,0,mpi_world,ierr)
      call mpi_bcast(prm2s(i2b)%ncoef,1,mpi_integer,0,mpi_world,ierr)
      call mpi_bcast(prm2s(i2b)%rc,1,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(prm2s(i2b)%rc2,1,mpi_real8,0,mpi_world,ierr)
      if( .not.allocated(prm2s(i2b)%knots) ) &
           allocate(prm2s(i2b)%knots(prm2s(i2b)%nknot))
      if( .not.allocated(prm2s(i2b)%coefs) ) &
           allocate(prm2s(i2b)%coefs(prm2s(i2b)%ncoef))
      call mpi_bcast(prm2s(i2b)%knots,prm2s(i2b)%nknot,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(prm2s(i2b)%coefs,prm2s(i2b)%ncoef,mpi_real8,0,mpi_world,ierr)
    enddo

    call mpi_bcast(has_trios, 1, mpi_logical, 0,mpi_world,ierr)
    if( has_trios ) then
      call mpi_bcast(rc3max,1,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(rc3max2,1,mpi_real8,0,mpi_world,ierr)
      if( .not. allocated(prm3ds) ) allocate(prm3ds(n3b))
      do i3b=1,n3b
        call mpi_bcast(prm3ds(i3b)%cb,2,mpi_character,0,mpi_world,ierr)
        call mpi_bcast(prm3ds(i3b)%csi,2,mpi_character,0,mpi_world,ierr)
        call mpi_bcast(prm3ds(i3b)%csj,2,mpi_character,0,mpi_world,ierr)
        call mpi_bcast(prm3ds(i3b)%csk,2,mpi_character,0,mpi_world,ierr)
        call mpi_bcast(prm3ds(i3b)%cknot,2,mpi_character,0,mpi_world,ierr)
        call mpi_bcast(prm3ds(i3b)%isp,1,mpi_integer,0,mpi_world,ierr)
        call mpi_bcast(prm3ds(i3b)%jsp,1,mpi_integer,0,mpi_world,ierr)
        call mpi_bcast(prm3ds(i3b)%ksp,1,mpi_integer,0,mpi_world,ierr)
        call mpi_bcast(prm3ds(i3b)%nklead,1,mpi_integer,0,mpi_world,ierr)
        call mpi_bcast(prm3ds(i3b)%nktrail,1,mpi_integer,0,mpi_world,ierr)
        call mpi_bcast(prm3ds(i3b)%rcij,1,mpi_real8,0,mpi_world,ierr)
        call mpi_bcast(prm3ds(i3b)%rcik,1,mpi_real8,0,mpi_world,ierr)
        call mpi_bcast(prm3ds(i3b)%rcij2,1,mpi_real8,0,mpi_world,ierr)
        call mpi_bcast(prm3ds(i3b)%rcik2,1,mpi_real8,0,mpi_world,ierr)
        if( .not.allocated(prm3ds(i3b)%knij) ) &
             allocate(prm3ds(i3b)%knij(prm3ds(i3b)%nknij), &
             prm3ds(i3b)%knik(prm3ds(i3b)%nknik), &
             prm3ds(i3b)%kncs(prm3ds(i3b)%nkncs))
        if( .not.allocated(prm3ds(i3b)%cfij)) &
             allocate(prm3ds(i3b)%cfij(prm3ds(i3b)%ncfij), &
             prm3ds(i3b)%cfik(prm3ds(i3b)%ncfik), &
             prm3ds(i3b)%cfcs(prm3ds(i3b)%ncfcs) )
        call mpi_bcast(prm3ds(i3b)%knij,prm3ds(i3b)%nknij, &
             mpi_real8,0,mpi_world,ierr)
        call mpi_bcast(prm3ds(i3b)%knik,prm3ds(i3b)%nknik, &
             mpi_real8,0,mpi_world,ierr)
        call mpi_bcast(prm3ds(i3b)%kncs,prm3ds(i3b)%nkncs, &
             mpi_real8,0,mpi_world,ierr)
        call mpi_bcast(prm3ds(i3b)%cfij,prm3ds(i3b)%ncfij, &
             mpi_real8,0,mpi_world,ierr)
        call mpi_bcast(prm3ds(i3b)%cfik,prm3ds(i3b)%ncfik, &
             mpi_real8,0,mpi_world,ierr)
        call mpi_bcast(prm3ds(i3b)%cfcs,prm3ds(i3b)%ncfcs, &
             mpi_real8,0,mpi_world,ierr)
      enddo
    endif  ! has_trios

    call mpi_bcast(has_solo, 1, mpi_logical, 0,mpi_world,ierr)
    if( has_solo ) then
      if( .not. allocated(prm1s) ) allocate(prm1s(n1b))
      call mpi_bcast(prm1s,n1b,mpi_integer,0,mpi_world,ierr)
      call mpi_bcast(erg1s,nspmax,mpi_real8,0,mpi_world,ierr)
    endif  ! has_solo

    call mpi_bcast(interact2, nspmax**2, mpi_integer, 0,mpi_world,ierr)
    call mpi_bcast(interact3, nspmax**3, mpi_integer, 0,mpi_world,ierr)

  end subroutine bcast_uf3d_params
!=======================================================================
  subroutine force_uf3_tmp(namax,natm,tag,ra,nnmax,aa,strs,h,hi &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rcin,lspr &
       ,mpi_world,myid,epi,epot,lstrs,iprint,l1st)
!
!  UF3 implementation without using recursive function of b-spline.
!
!  TODO: More efficient B-spline implementation is available (see, lammps src/ML-UF3/pair_uf3.cpp).
!        But that is a bit complicated, use simple b_spl() routine for now.
!
    use util, only: itotOf
    implicit none
    integer,intent(in):: namax,natm,nnmax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_world,myid,lspr(0:nnmax,namax),nex(3)
    real(8),intent(in):: ra(3,namax),tag(namax) &
         ,h(3,3),hi(3,3),sv(3,6)
    real(8),intent(inout):: rcin
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: l1st
    logical:: lstrs

!.....local
    integer:: ia,ja,ka,jj,kk,l,is,nr2,n,nij3,inc,nik3,njk3,&
         nik,njk,nij,itot,jtot,ktot,i2b,i3b,js,ks,jsp,ksp,ierr, &
         ixyz,jxyz,lij,lik,ljk
    real(8):: epotl2,epotl3,epot2,epot3,tmp,tmp2,bij(-3:0),dbij(-3:0), &
         bij3(-3:0),dbij3(-3:0),bik3(-3:0),dbik3(-3:0),bjk3(-3:0), &
         dbjk3(-3:0),c2t,c3t,epotl1,epot1,fac3b,tmp3
    real(8):: xi(3),xj(3),xk(3),xij(3),xik(3),xjk(3),rij(3),rik(3),&
         rjk(3),dij2,dij,dik2,dik,djk2,djk,drijj(3),drikk(3),&
         drjkk(3),tmpij(3),tmpik(3),tmpjk(3)
    real(8),save:: rcin2 = -1d0

    type(prm2):: p2
    type(prm3):: p3

    if( rcin2 < 0d0 ) then  ! Probably it is the 1st call...
      if( rcin < rcmax ) then
        if( myid == 0 ) then
          write(6,'(1x,a)') "ERROR: Cutoff radius is not appropriate !!!"
          write(6,'(1x,a,f0.3)') "  rc should be longer than ", rcmax
        endif
        call mpi_finalize(ierr)
        stop
      endif
      rcin2 = rcin*rcin
    endif

    if( .not.allocated(aal2) ) then
      allocate(aal2(3,namax),aal3(3,namax),strsl(3,3,namax),ls3b(0:nnmax))
!!$      allocate(aal2(3,namax),aal3(3,namax),strsl(3,3,namax))
    endif
    if( size(aal2) < 3*namax ) then
      deallocate(aal2,aal3,strsl)
      allocate(aal2(3,namax),aal3(3,namax),strsl(3,3,namax))
    endif
    if( size(ls3b) < nnmax+1 ) then
      deallocate(ls3b)
      allocate(ls3b(0:nnmax))
    endif

    aal2(:,:) = 0d0
    aal3(:,:) = 0d0
    strsl(:,:,:) = 0d0
    epotl1 = 0d0
    epotl2 = 0d0
    epotl3 = 0d0
!$omp parallel
!$omp do private(ia,is,xi,jj,ja,js,i2b,p2,xj,xij,rij,dij2,dij, &
!$omp      drijj,nr2,n,bij,c2t,tmp,dbij,tmp2,ixyz,jxyz, &
!$omp      jsp,ksp,i3b,p3,nij3,inc,bij3,dbij3, &
!$omp      kk,ka,ks,xk,xik,xjk,rik,dik2,dik,rjk,djk2,djk, &
!$omp      drikk,drjkk,nik3,njk3,nik,bik3,dbik3,njk,bjk3,dbjk3, &
!$omp      l,nij,c3t,tmpij,tmpik,tmpjk) &
!$omp      reduction(+:epotl2,epotl3)
    do ia=1,natm
      is = int(tag(ia))
!!$      epi(ia) = epi(ia) +erg1s(is)
      epotl1 = epotl1 +erg1s(is)
      xi(1:3) = ra(1:3,ia)
      tmp3 = 0d0
      do jj=1,lspr(0,ia)
        ja = lspr(jj,ia)
!!$        if( ja <= ia ) cycle
        js = int(tag(ja))
!.....Pair terms
        i2b = interact2(is,js)
        if( i2b <= 0 ) cycle
        p2 = prm2s(i2b)
        xj(1:3) = ra(1:3,ja)
        xij(1:3) = xj(1:3) -xi(1:3)
        rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij2 = rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3)
        if( dij2 > p2%rc2 ) cycle
        dij = sqrt(dij2)
!.....NOTE: ls3b could be source of bugs and stop using ls3b.
!           To speed up by reducing pairs for 3-body, use other way
!           that is safer than this approach.
!!$!.....Make short-distance pair-list for 3-body term
!!$        if( has_trios ) then
!!$          if( dij2 < rc2_3b(is,js) ) then
!!$            ls3b(0) = ls3b(0) +1
!!$            ls3b(ls3b(0)) = ja
!!$          endif
!!$        endif
        drijj(1:3) = rij(1:3)/dij
        call b_spl(dij,p2%knots,p2%nknot,nr2,bij,dbij)
        do lij = -3,0
          n = nr2 +lij
          if( n < 1 .or. n > p2%nknot-4 ) cycle
          c2t = p2%coefs(n)
          tmp = c2t *bij(lij)
!.....Energy
          epi(ia) = epi(ia) +tmp
          epotl2 = epotl2 +tmp
!.....Forces
          tmp2 = c2t *dbij(lij)
          do ixyz=1,3
!$omp atomic
            aal2(ixyz,ia) = aal2(ixyz,ia) +drijj(ixyz)*tmp2
!$omp atomic
            aal2(ixyz,ja) = aal2(ixyz,ja) -drijj(ixyz)*tmp2
          enddo
!.....Stresses
          do ixyz=1,3
            do jxyz=1,3
!$omp atomic
              strsl(jxyz,ixyz,ia)= strsl(jxyz,ixyz,ia) &
                   -0.5d0 *tmp2*rij(ixyz)*drijj(jxyz)
!$omp atomic
              strsl(jxyz,ixyz,ja)= strsl(jxyz,ixyz,ja) &
                   -0.5d0 *tmp2*rij(ixyz)*drijj(jxyz)
            enddo
          enddo
        enddo  ! lij
      enddo  ! jj

!.....Trio part is separated from pair part,
!.....which may be slower because of double computation of dij,
!.....but this code is a bit simpler.
      if( .not.has_trios ) cycle
!!$!.....Compute ls3b before going to main 3b part
!!$      ls3b(0) = 0
!!$      do jj=1,lspr(0,ia)
!!$        ja = lspr(jj,ia)
!!$        xj(1:3) = ra(1:3,ja)
!!$        xij(1:3) = xj(1:3) -xi(1:3)
!!$        rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
!!$        dij2 = rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3)
!!$        if( dij2 > rc3max2 ) cycle
!!$        ls3b(0) = ls3b(0) +1
!!$        ls3b(ls3b(0)) = ja
!!$      enddo

      tmp3 = 0d0
!!$      do jsp=1,nspmax
!!$        do ksp=jsp,nspmax
!!$          i3b = interact3(is,jsp,ksp)
!!$          if( i3b <= 0 ) cycle
!!$          p3 = prm3s(i3b)
          do jj=1,lspr(0,ia)
            ja = lspr(jj,ia)
            js = int(tag(ja))
!!$            if( js /= jsp ) cycle
            jtot = itotOf(tag(ja))
            xj(1:3) = ra(1:3,ja)
            xij(1:3) = xj(1:3) -xi(1:3)
            rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
            dij2 = rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3)
!!$            if( dij2 > p3%rcij2 ) cycle
            if( dij2 > rc3max2 ) cycle
            dij = sqrt(dij2)
            drijj(1:3) = rij(1:3)/dij
!!$            call b_spl(dij, p3%knij, p3%nknij, nij3, bij3, dbij3)

            do kk=jj+1,lspr(0,ia)
              ka = lspr(kk,ia)
!!$              if( ka == ja ) cycle
              ktot = itotOf(tag(ka))
!.....Taking into account the double counting of symmetric terms
              fac3b = 1d0
!!$              if( jsp == ksp ) fac3b = 0.5d0
              ks = int(tag(ka))
!!$              if( ks /= ksp ) cycle
              i3b = interact3(is,js,ks)
              if( i3b <= 0 ) cycle
              p3 = prm3s(i3b)
              if( dij2 > p3%rcij2 ) cycle
              xk(1:3) = ra(1:3,ka)
              xik(1:3) = xk(1:3) -xi(1:3)
              xjk(1:3) = xk(1:3) -xj(1:3)
              rik(1:3) = h(1:3,1)*xik(1) +h(1:3,2)*xik(2) +h(1:3,3)*xik(3)
              dik2 = rik(1)*rik(1) +rik(2)*rik(2) +rik(3)*rik(3)
              if( dik2 > p3%rcik2 ) cycle
              dik = sqrt(dik2)
              rjk(1:3) = h(1:3,1)*xjk(1) +h(1:3,2)*xjk(2) +h(1:3,3)*xjk(3)
              djk2 = rjk(1)*rjk(1) +rjk(2)*rjk(2) +rjk(3)*rjk(3)
              if( djk2 > p3%rcjk2 ) cycle
              djk = sqrt(djk2)
              drikk(1:3) = rik(1:3)/dik
              drjkk(1:3) = rjk(1:3)/djk
!.....B-spline part
              call b_spl(dij, p3%knij, p3%nknij, nij3, bij3, dbij3)
              call b_spl(dik, p3%knik, p3%nknik, nik3, bik3, dbik3)
              call b_spl(djk, p3%knjk, p3%nknjk, njk3, bjk3, dbjk3)
              do lik = -3,0
                nik = nik3 +lik
                if( nik < 1 .or. nik > p3%nknik-4 ) cycle
                do ljk = -3,0
                  njk = njk3 +ljk
                  if( njk < 1 .or. njk > p3%nknjk-4 ) cycle
                  do lij = -3,0
                    nij = nij3 +lij
                    if( nij < 1 .or. nij > p3%nknij-4 ) cycle
!.....Energy
                    c3t = p3%coefs(njk,nik,nij) *fac3b
                    tmp = c3t*bij3(lij)*bik3(lik)*bjk3(ljk)
                    epi(ia) = epi(ia) +tmp
                    epotl3 = epotl3 +tmp
!.....Force
                    tmpij(1:3) = dbij3(lij)* bik3(lik)* bjk3(ljk) &
                         *drijj(1:3) *c3t
                    tmpik(1:3) =  bij3(lij)*dbik3(lik)* bjk3(ljk) &
                         *drikk(1:3) *c3t
                    tmpjk(1:3) =  bij3(lij)* bik3(lik)*dbjk3(ljk) &
                         *drjkk(1:3) *c3t
                    do ixyz=1,3
!$omp atomic
                      aal3(ixyz,ia)= aal3(ixyz,ia) &
                           +(tmpij(ixyz) +tmpik(ixyz))
!$omp atomic
                      aal3(ixyz,ja)= aal3(ixyz,ja) &
                           +(-tmpij(ixyz) +tmpjk(ixyz))
!$omp atomic
                      aal3(ixyz,ka)= aal3(ixyz,ka) &
                           +(-tmpik(ixyz) -tmpjk(ixyz))
                    enddo
!.....Stresses
                    do ixyz=1,3
                      do jxyz=1,3
!$omp atomic
                        strsl(jxyz,ixyz,ia)= strsl(jxyz,ixyz,ia) &
                             -0.5d0 *(rij(ixyz)*tmpij(jxyz) &
                             +rik(ixyz)*tmpik(jxyz))
!$omp atomic
                        strsl(jxyz,ixyz,ja)= strsl(jxyz,ixyz,ja) &
                             -0.5d0 *(rij(ixyz)*tmpij(jxyz) &
                             +rjk(ixyz)*tmpjk(jxyz))
!$omp atomic
                        strsl(jxyz,ixyz,ka)= strsl(jxyz,ixyz,ka) &
                             -0.5d0 *(rik(ixyz)*tmpik(jxyz) &
                             +rjk(ixyz)*tmpjk(jxyz))
                      enddo
                    enddo

                  enddo  ! lij
                enddo  ! ljk
              enddo  ! lik

            enddo  ! kk
          enddo  ! jj
!!$        enddo  ! ksp
!!$      enddo  ! jsp

    enddo ! ia
!$omp end do
!$omp end parallel

    call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,aal2,3)
    if( has_trios ) call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex, &
         lsrc,myparity,nn,mpi_world,aal3,3)
    aa(1:3,1:natm) = aa(1:3,1:natm) +aal2(1:3,1:natm) +aal3(1:3,1:natm)

    call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,strsl,9)
    strs(1:3,1:3,1:natm) = strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)

!-----gather epot
    epot1 = 0d0
    epot2 = 0d0
    epot3 = 0d0
    call mpi_allreduce(epotl2,epot2,1,mpi_real8,mpi_sum,mpi_world,ierr)
    if( has_solo ) call mpi_allreduce(epotl1,epot1,1,mpi_real8, &
         mpi_sum,mpi_world,ierr)
    if( has_trios ) call mpi_allreduce(epotl3,epot3,1,mpi_real8, &
         mpi_sum,mpi_world,ierr)
    epot= epot +epot1 +epot2 +epot3
    if( myid == 0 .and. iprint > 2 ) &
         print '(a,3es12.4)',' force_uf3 epot1,epot2,epot3 = ', &
         epot1,epot2,epot3

    return
  end subroutine force_uf3_tmp
!=======================================================================
  subroutine force_uf3(namax,natm,tag,ra,nnmax,aa,strs,h,hi &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rcin,lspr &
       ,mpi_world,myid,epi,epot,lstrs,iprint,l1st)
!
!  UF3 implementation without using non-recursive function of b-spline.
!
!  TODO: More efficient B-spline implementation is available (see, lammps src/ML-UF3/pair_uf3.cpp).
!        But that is a bit complicated, use simple b_spl() routine for now.
!
    use util, only: itotOf
    implicit none
    integer,intent(in):: namax,natm,nnmax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_world,myid,lspr(0:nnmax,namax),nex(3)
    real(8),intent(in):: ra(3,namax),tag(namax) &
         ,h(3,3),hi(3,3),sv(3,6)
    real(8),intent(inout):: rcin
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: l1st
    logical:: lstrs

!.....local
    integer:: ia,ja,ka,jj,kk,l,is,nr2,n,nij3,inc,nik3,njk3,&
         nik,njk,nij,itot,jtot,ktot,i2b,i3b,js,ks,jsp,ksp,ierr, &
         ixyz,jxyz,lij,lik,ljk
    real(8):: epotl2,epotl3,epot2,epot3,tmp,tmp2,bij(-3:0),dbij(-3:0), &
         bij3(-3:0),dbij3(-3:0),bik3(-3:0),dbik3(-3:0),bjk3(-3:0), &
         dbjk3(-3:0),c2t,c3t,epotl1,epot1,fac3b,tmp3
    real(8):: xi(3),xj(3),xk(3),xij(3),xik(3),xjk(3),rij(3),rik(3),&
         rjk(3),dij2,dij,dik2,dik,djk2,djk,drijj(3),drikk(3),&
         drjkk(3),tmpij(3),tmpik(3),tmpjk(3)
    real(8),save:: rcin2 = -1d0

    type(prm2):: p2
    type(prm3):: p3

    if( rcin2 < 0d0 ) then  ! Probably it is the 1st call...
      if( rcin < rcmax ) then
        if( myid == 0 ) then
          write(6,'(1x,a)') "ERROR: Cutoff radius is not appropriate !!!"
          write(6,'(1x,a,f0.3)') "  rc should be longer than ", rcmax
        endif
        call mpi_finalize(ierr)
        stop
      endif
      rcin2 = rcin*rcin
    endif

    if( .not.allocated(aal2) ) then
      allocate(aal2(3,namax),aal3(3,namax),strsl(3,3,namax),ls3b(0:nnmax))
!!$      allocate(aal2(3,namax),aal3(3,namax),strsl(3,3,namax))
    endif
    if( size(aal2) < 3*namax ) then
      deallocate(aal2,aal3,strsl)
      allocate(aal2(3,namax),aal3(3,namax),strsl(3,3,namax))
    endif
    if( size(ls3b) < nnmax+1 ) then
      deallocate(ls3b)
      allocate(ls3b(0:nnmax))
    endif

    aal2(:,:) = 0d0
    aal3(:,:) = 0d0
    strsl(:,:,:) = 0d0
    epotl1 = 0d0
    epotl2 = 0d0
    epotl3 = 0d0
!$omp parallel
!$omp do private(ia,is,xi,jj,ja,js,i2b,p2,xj,xij,rij,dij2,dij, &
!$omp      drijj,nr2,n,bij,c2t,tmp,dbij,tmp2,ixyz,jxyz, &
!$omp      jsp,ksp,i3b,p3,nij3,inc,bij3,dbij3, &
!$omp      kk,ka,ks,xk,xik,xjk,rik,dik2,dik,rjk,djk2,djk, &
!$omp      drikk,drjkk,nik3,njk3,nik,bik3,dbik3,njk,bjk3,dbjk3, &
!$omp      l,nij,c3t,tmpij,tmpik,tmpjk) &
!$omp      reduction(+:epotl2,epotl3)
    do ia=1,natm
      is = int(tag(ia))
!!$      epi(ia) = epi(ia) +erg1s(is)
      epotl1 = epotl1 +erg1s(is)
      xi(1:3) = ra(1:3,ia)
      tmp3 = 0d0
      do jj=1,lspr(0,ia)
        ja = lspr(jj,ia)
!!$        if( ja <= ia ) cycle
        js = int(tag(ja))
!.....Pair terms
        i2b = interact2(is,js)
        if( i2b <= 0 ) cycle
        p2 = prm2s(i2b)
        xj(1:3) = ra(1:3,ja)
        xij(1:3) = xj(1:3) -xi(1:3)
        rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij2 = rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3)
        if( dij2 > p2%rc2 ) cycle
        dij = sqrt(dij2)
!.....NOTE: ls3b could be source of bugs and stop using ls3b.
!           To speed up by reducing pairs for 3-body, use other way
!           that is safer than this approach.
!!$!.....Make short-distance pair-list for 3-body term
!!$        if( has_trios ) then
!!$          if( dij2 < rc2_3b(is,js) ) then
!!$            ls3b(0) = ls3b(0) +1
!!$            ls3b(ls3b(0)) = ja
!!$          endif
!!$        endif
        drijj(1:3) = rij(1:3)/dij
        call b_spl(dij,p2%knots,p2%nknot,nr2,bij,dbij)
        do lij = -3,0
          n = nr2 +lij
          if( n < 1 .or. n > p2%nknot-4 ) cycle
          c2t = p2%coefs(n)
          tmp = c2t *bij(lij)
!.....Energy
          epi(ia) = epi(ia) +tmp
          epotl2 = epotl2 +tmp
!.....Forces
          tmp2 = c2t *dbij(lij)
          do ixyz=1,3
!$omp atomic
            aal2(ixyz,ia) = aal2(ixyz,ia) +drijj(ixyz)*tmp2
!$omp atomic
            aal2(ixyz,ja) = aal2(ixyz,ja) -drijj(ixyz)*tmp2
          enddo
!.....Stresses
          do ixyz=1,3
            do jxyz=1,3
!$omp atomic
              strsl(jxyz,ixyz,ia)= strsl(jxyz,ixyz,ia) &
                   -0.5d0 *tmp2*rij(ixyz)*drijj(jxyz)
!$omp atomic
              strsl(jxyz,ixyz,ja)= strsl(jxyz,ixyz,ja) &
                   -0.5d0 *tmp2*rij(ixyz)*drijj(jxyz)
            enddo
          enddo
        enddo  ! lij
      enddo  ! jj

!.....Trio part is separated from pair part,
!.....which may be slower because of double computation of dij,
!.....but this code is a bit simpler.
      if( .not.has_trios ) cycle
!!$!.....Compute ls3b before going to main 3b part
!!$      ls3b(0) = 0
!!$      do jj=1,lspr(0,ia)
!!$        ja = lspr(jj,ia)
!!$        xj(1:3) = ra(1:3,ja)
!!$        xij(1:3) = xj(1:3) -xi(1:3)
!!$        rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
!!$        dij2 = rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3)
!!$        if( dij2 > rc3max2 ) cycle
!!$        ls3b(0) = ls3b(0) +1
!!$        ls3b(ls3b(0)) = ja
!!$      enddo

      tmp3 = 0d0
      do jsp=1,nspmax
        do ksp=jsp,nspmax
          i3b = interact3(is,jsp,ksp)
          if( i3b <= 0 ) cycle
          p3 = prm3s(i3b)
          do jj=1,lspr(0,ia)
            ja = lspr(jj,ia)
            js = int(tag(ja))
            if( js /= jsp ) cycle
            jtot = itotOf(tag(ja))
            xj(1:3) = ra(1:3,ja)
            xij(1:3) = xj(1:3) -xi(1:3)
            rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
            dij2 = rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3)
            if( dij2 > p3%rcij2 ) cycle
            dij = sqrt(dij2)
            drijj(1:3) = rij(1:3)/dij
            call b_spl(dij, p3%knij, p3%nknij, nij3, bij3, dbij3)

            do kk=1,lspr(0,ia)
              ka = lspr(kk,ia)
              if( ka == ja ) cycle
              ktot = itotOf(tag(ka))
!.....Taking into account the double counting of symmetric terms
              fac3b = 1d0
              if( jsp == ksp ) fac3b = 0.5d0
              ks = int(tag(ka))
              if( ks /= ksp ) cycle
              xk(1:3) = ra(1:3,ka)
              xik(1:3) = xk(1:3) -xi(1:3)
              xjk(1:3) = xk(1:3) -xj(1:3)
              rik(1:3) = h(1:3,1)*xik(1) +h(1:3,2)*xik(2) +h(1:3,3)*xik(3)
              dik2 = rik(1)*rik(1) +rik(2)*rik(2) +rik(3)*rik(3)
              if( dik2 > p3%rcik2 ) cycle
              dik = sqrt(dik2)
              rjk(1:3) = h(1:3,1)*xjk(1) +h(1:3,2)*xjk(2) +h(1:3,3)*xjk(3)
              djk2 = rjk(1)*rjk(1) +rjk(2)*rjk(2) +rjk(3)*rjk(3)
              if( djk2 > p3%rcjk2 ) cycle
              djk = sqrt(djk2)
              drikk(1:3) = rik(1:3)/dik
              drjkk(1:3) = rjk(1:3)/djk
!.....B-spline part
              call b_spl(dik, p3%knik, p3%nknik, nik3, bik3, dbik3)
              call b_spl(djk, p3%knjk, p3%nknjk, njk3, bjk3, dbjk3)
              do lik = -3,0
                nik = nik3 +lik
                if( nik < 1 .or. nik > p3%nknik-4 ) cycle
                do ljk = -3,0
                  njk = njk3 +ljk
                  if( njk < 1 .or. njk > p3%nknjk-4 ) cycle
                  do lij = -3,0
                    nij = nij3 +lij
                    if( nij < 1 .or. nij > p3%nknij-4 ) cycle
!.....Energy
                    c3t = p3%coefs(njk,nik,nij) *fac3b
                    tmp = c3t*bij3(lij)*bik3(lik)*bjk3(ljk)
                    tmp3 = tmp3 + tmp
                    epi(ia) = epi(ia) +tmp
                    epotl3 = epotl3 +tmp
!.....Force
                    tmpij(1:3) = dbij3(lij)* bik3(lik)* bjk3(ljk) &
                         *drijj(1:3) *c3t
                    tmpik(1:3) =  bij3(lij)*dbik3(lik)* bjk3(ljk) &
                         *drikk(1:3) *c3t
                    tmpjk(1:3) =  bij3(lij)* bik3(lik)*dbjk3(ljk) &
                         *drjkk(1:3) *c3t
                    do ixyz=1,3
!$omp atomic
                      aal3(ixyz,ia)= aal3(ixyz,ia) &
                           +(tmpij(ixyz) +tmpik(ixyz))
!$omp atomic
                      aal3(ixyz,ja)= aal3(ixyz,ja) &
                           +(-tmpij(ixyz) +tmpjk(ixyz))
!$omp atomic
                      aal3(ixyz,ka)= aal3(ixyz,ka) &
                           +(-tmpik(ixyz) -tmpjk(ixyz))
                    enddo
!.....Stresses
                    do ixyz=1,3
                      do jxyz=1,3
!$omp atomic
                        strsl(jxyz,ixyz,ia)= strsl(jxyz,ixyz,ia) &
                             -0.5d0 *(rij(ixyz)*tmpij(jxyz) &
                             +rik(ixyz)*tmpik(jxyz))
!$omp atomic
                        strsl(jxyz,ixyz,ja)= strsl(jxyz,ixyz,ja) &
                             -0.5d0 *(rij(ixyz)*tmpij(jxyz) &
                             +rjk(ixyz)*tmpjk(jxyz))
!$omp atomic
                        strsl(jxyz,ixyz,ka)= strsl(jxyz,ixyz,ka) &
                             -0.5d0 *(rik(ixyz)*tmpik(jxyz) &
                             +rjk(ixyz)*tmpjk(jxyz))
                      enddo
                    enddo

                  enddo  ! lij
                enddo  ! ljk
              enddo  ! lik

            enddo  ! kk
          enddo  ! jj
        enddo  ! ksp
      enddo  ! jsp

    enddo ! ia
!$omp end do
!$omp end parallel

    call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,aal2,3)
    if( has_trios ) call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex, &
         lsrc,myparity,nn,mpi_world,aal3,3)
    aa(1:3,1:natm) = aa(1:3,1:natm) +aal2(1:3,1:natm) +aal3(1:3,1:natm)

    call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,strsl,9)
    strs(1:3,1:3,1:natm) = strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)

!-----gather epot
    epot1 = 0d0
    epot2 = 0d0
    epot3 = 0d0
    call mpi_allreduce(epotl2,epot2,1,mpi_real8,mpi_sum,mpi_world,ierr)
    if( has_solo ) call mpi_allreduce(epotl1,epot1,1,mpi_real8, &
         mpi_sum,mpi_world,ierr)
    if( has_trios ) call mpi_allreduce(epotl3,epot3,1,mpi_real8, &
         mpi_sum,mpi_world,ierr)
    epot= epot +epot1 +epot2 +epot3
    if( myid == 0 .and. iprint > 2 ) &
         print '(a,3es12.4)',' force_uf3 epot1,epot2,epot3 = ', &
         epot1,epot2,epot3

    return
  end subroutine force_uf3
!=======================================================================
  subroutine force_uf3l(namax,natm,tag,ra,nnmax,aa,strs,h,hi &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rcin,lspr &
       ,mpi_world,myid,epi,epot,lstrs,iprint,l1st)
!
!  Light-weight UF3 implementation without using non-recursive function of b-spline.
!
    use util, only: itotOf
#ifdef IMPULSE
    use impulse,only: ftaul, itot_impls, tau_impls
#endif
    implicit none
    real(8),parameter:: tiny = 1d-8
    integer,intent(in):: namax,natm,nnmax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_world,myid,lspr(0:nnmax,namax),nex(3)
    real(8),intent(in):: ra(3,namax),tag(namax) &
         ,h(3,3),hi(3,3),sv(3,6)
    real(8),intent(inout):: rcin
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: l1st, lstrs

!.....local
    integer:: ia,ja,ka,jj,kk,l,is,nr2,n,inc,&
         nij,itot,jtot,ktot,i2b,i3b,js,ks,jsp,ksp,ierr, &
         ixyz,jxyz,lcs,lij,ncs
    real(8):: epotl2,epotl3,epot2,epot3,tmp,tmp2,bij(-3:0),dbij(-3:0), &
         bcs(-3:0),dbcs(-3:0),c2t,c3t,epotl1,epot1,fac3b
    real(8):: xi(3),xj(3),xk(3),xij(3),xik(3),xjk(3),rij(3),rik(3),&
         rjk(3),dij2,dij,dik2,dik,drijj(3),drikk(3),&
         drjkk(3),diji,diki,drijc,drikc,dv3csn,dv3rij,dv3rik,sumcb,sumcdb
    real(8):: dcsnj(3),dcsnk(3),dcsni(3),tmpj(3),tmpk(3),gmj,gmk,csn,vexp
    real(8):: rcij, rcik, rcij2, rcik2
    real(8),save:: rcin2 = -1d0

#ifdef CONTRIB
    real(8):: epot_LiLa, epot_NbO, epot_X, epot_Xp, epot_Xm, epot_elem, &
         epot_Li, epot_La
#endif

    type(prm2):: p2
    type(prm3l):: p3

    if( rcin2 < 0d0 ) then  ! Probably it is the 1st call...
      if( rcin < rcmax ) then
        if( myid == 0 ) then
          write(6,'(1x,a)') "ERROR: Cutoff radius is not appropriate !!!"
          write(6,'(1x,a,f0.3)') "  rc should be longer than ", rcmax
        endif
        call mpi_finalize(ierr)
        stop
      endif
      rcin2 = rcin*rcin
    endif

    if( .not.allocated(aal2) ) then
      allocate(aal2(3,namax),aal3(3,namax),strsl(3,3,namax),ls3b(0:nnmax))
!!$      allocate(aal2(3,namax),aal3(3,namax),strsl(3,3,namax))
    endif
    if( size(aal2) < 3*namax ) then
      deallocate(aal2,aal3,strsl)
      allocate(aal2(3,namax),aal3(3,namax),strsl(3,3,namax))
    endif
    if( size(ls3b) < nnmax+1 ) then
      deallocate(ls3b)
      allocate(ls3b(0:nnmax))
    endif

    aal2(:,:) = 0d0
    aal3(:,:) = 0d0
    strsl(:,:,:) = 0d0
    epotl1 = 0d0
    epotl2 = 0d0
    epotl3 = 0d0
#ifdef IMPULSE
    ftaul(:) = 0d0
#endif
#ifdef CONTRIB
    epot_LiLa = 0d0
    epot_Li = 0d0
    epot_La = 0d0
    epot_NbO = 0d0
    epot_X = 0d0
    epot_Xp = 0d0
    epot_Xm = 0d0
    epot_elem = 0d0
#endif
    do ia=1,natm
      is = int(tag(ia))
      itot = itotOf(tag(ia))
      epi(ia) = epi(ia) +erg1s(is)
      epotl1 = epotl1 +erg1s(is)
#ifdef CONTRIB
      epot_elem = epot_elem +erg1s(is)
#endif
      xi(1:3) = ra(1:3,ia)
      do jj=1,lspr(0,ia)
        ja = lspr(jj,ia)
!!$        if( ja <= ia ) cycle
        js = int(tag(ja))
        jtot = itotOf(tag(ja))
!.....Pair terms
        i2b = interact2(is,js)
        if( i2b <= 0 ) cycle
        p2 = prm2s(i2b)
        xj(1:3) = ra(1:3,ja)
        xij(1:3) = xj(1:3) -xi(1:3)
        rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij2 = rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3)
        if( dij2 > p2%rc2 ) cycle
        dij = sqrt(dij2)
!.....NOTE: ls3b could be source of bugs and stop using ls3b.
!           To speed up by reducing pairs for 3-body, use other way
!           that is safer than this approach.
!!$!.....Make short-distance pair-list for 3-body term
!!$        if( has_trios ) then
!!$          if( dij2 < rc2_3b(is,js) ) then
!!$            ls3b(0) = ls3b(0) +1
!!$            ls3b(ls3b(0)) = ja
!!$          endif
!!$        endif
        drijj(1:3) = rij(1:3)/dij
        call b_spl(dij,p2%knots,p2%nknot,nr2,bij,dbij)
        do lij = -3,0
          n = nr2 +lij
          if( n < 1 .or. n > p2%nknot-4 ) cycle
          c2t = p2%coefs(n)
          tmp = c2t *bij(lij)
!.....Energy
#ifdef CONTRIB
          if( is==1 .and. js==1 ) then
            epot_Li = epot_Li +tmp
          else if( is==2 .and. js==2 ) then
            epot_La = epot_La +tmp
          else if( (is==1 .and. js==2) .or. (is==2.and.js==1) ) then
            epot_LiLa = epot_LiLa +tmp
          else if( (is==3 .or. is==4) .and. (js==3 .or. js==4) ) then
            epi(ia) = epi(ia) +tmp
            epot_NbO = epot_NbO +tmp
          else
            if( is==2 ) epi(ia) = epi(ia) +tmp
            epot_X = epot_X +tmp
            if( tmp > 0d0 ) then
              epot_Xp = epot_Xp +tmp
            else
              epot_Xm = epot_Xm +tmp
            endif
          endif
#else
          epi(ia) = epi(ia) +tmp
#endif
          epotl2 = epotl2 +tmp
!.....Forces
          tmp2 = c2t *dbij(lij)
          do ixyz=1,3
            aal2(ixyz,ia) = aal2(ixyz,ia) +drijj(ixyz)*tmp2
            aal2(ixyz,ja) = aal2(ixyz,ja) -drijj(ixyz)*tmp2
          enddo
#ifdef IMPULSE
          if( itot == itot_impls ) then
            ftaul(js) = ftaul(js) +tmp2 &
                 *(drijj(1)*tau_impls(1) &
                 +drijj(2)*tau_impls(2) &
                 +drijj(3)*tau_impls(3) )
          endif
          if( jtot == itot_impls ) then
            ftaul(is) = ftaul(is) -tmp2 &
                 *(drijj(1)*tau_impls(1) &
                 +drijj(2)*tau_impls(2) &
                 +drijj(3)*tau_impls(3) )
          endif
#endif
!.....Stresses
          do ixyz=1,3
            do jxyz=1,3
              strsl(jxyz,ixyz,ia)= strsl(jxyz,ixyz,ia) &
                   -0.5d0 *tmp2*rij(ixyz)*drijj(jxyz)
              strsl(jxyz,ixyz,ja)= strsl(jxyz,ixyz,ja) &
                   -0.5d0 *tmp2*rij(ixyz)*drijj(jxyz)
            enddo
          enddo
        enddo  ! lij
      enddo  ! jj

!.....Trio part is separated from pair part,
!.....which may be slower because of double computation of dij,
!.....but this code is a bit simpler.
      if( .not.has_trios ) cycle

      do jj=1,lspr(0,ia)-1
        ja = lspr(jj,ia)
        js = int(tag(ja))
        jtot = itotOf(tag(ja))
        xj(1:3)= ra(1:3,ja)
        xij(1:3)= xj(1:3) - xi(1:3)
        rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij2 = rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3)
        if( dij2 > rc3max2 ) cycle
        dij = sqrt(dij2)
        diji = 1d0/dij
        drijj(1:3) = rij(1:3)*diji
        do kk=jj+1,lspr(0,ia)
          ka = lspr(kk,ia)
          ks = int(tag(ka))
          i3b = interact3(is,js,ks)
          if( i3b <= 0 ) cycle
          p3 = prm3ls(i3b)
          rcij = p3%rcij
          rcij2 = p3%rcij2
          rcik = p3%rcik
          rcik2 = p3%rcik2
          if( js==ks ) then
            rcij = (rcij+rcik)/2
            rcik = rcij
            rcij2= (rcij2+rcik2)/2
            rcik2= rcij2
          endif
          if(  dij2 > rcij2 ) exit
          ktot = itotOf(tag(ka))
!!$          fac3b = 1d0
!!$          if( js == ks ) fac3b = 0.5d0
          xk(1:3) = ra(1:3,ka)
          xik(1:3) = xk(1:3) -xi(1:3)
          rik(1:3) = h(1:3,1)*xik(1) +h(1:3,2)*xik(2) +h(1:3,3)*xik(3)
          dik2 = rik(1)*rik(1) +rik(2)*rik(2) +rik(3)*rik(3)
          if( dik2 > rcik2 ) cycle
          dik = sqrt(dik2)
          diki = 1d0/dik
          drikk(1:3) = rik(1:3)*diki
          drijc = 1d0/(dij-rcij)
          drikc = 1d0/(dik-rcik)
!.....parameters
          gmj = p3%gmj
          gmk = p3%gmk
          if( js==ks ) then
            gmj = (gmj+gmk) /2
            gmk = gmj
          endif
!.....common terms
          vexp = dexp(gmj*drijc +gmk*drikc)
          csn = (rij(1)*rik(1) +rij(2)*rik(2) +rij(3)*rik(3)) *(diji*diki)
          csn = max(min(csn, 1d0-tiny), -1d0+tiny)
          call b_spl(-csn, p3%knots, p3%nknot, ncs, bcs, dbcs)
          sumcb = 0d0
          sumcdb= 0d0
          do lcs = -3,0
            n = ncs +lcs
            if( n < 1 .or. n > p3%ncoef ) cycle
            c3t = p3%coefs(n) !*fac3b
            sumcb = sumcb +c3t*bcs(lcs)
            sumcdb= sumcdb +c3t*dbcs(lcs)
          enddo  ! lcs
!.....energy
          tmp = vexp *sumcb
#ifdef CONTRIB
          if( (is==1.or.is==2) .and. (js==1.or.js==2) .and. &
               (ks==1.or.ks==2) ) then
            epot_LiLa = epot_LiLa +tmp
          else if( (is==3.or.is==4) .and. (js==3.or.js==4) .and. &
               (ks==3.or.ks==4) ) then
            epot_NbO = epot_NbO +tmp
            epi(ia) = epi(ia) +tmp
          else
            epot_X = epot_X +tmp
            if( is==2 ) epi(ia) = epi(ia) +tmp
            if( tmp > 0d0 ) then
              epot_Xp = epot_Xp +tmp
            else
              epot_Xm = epot_Xm +tmp
            endif
          endif
#else
          epi(ia) = epi(ia) +tmp
#endif
          epotl3 = epotl3 +tmp
!.....force
          dv3rij = -drijc*drijc *gmj *tmp
          dv3rik = -drikc*drikc *gmk *tmp
          dv3csn = -vexp *sumcdb
          dcsnj(1:3)= (-rij(1:3)*csn*(diji*diji) +rik(1:3)*(diji*diki))
          dcsnk(1:3)= (-rik(1:3)*csn*(diki*diki) +rij(1:3)*(diji*diki))
!!$          dcsni(1:3)= -dcsnj(1:3) -dcsnk(1:3)
          tmpj(1:3)= dv3rij*drijj(1:3) +dv3csn*dcsnj(1:3)
          tmpk(1:3)= dv3rik*drikk(1:3) +dv3csn*dcsnk(1:3)
          do ixyz=1,3
            aal3(ixyz,ia) = aal3(ixyz,ia) +tmpj(ixyz) +tmpk(ixyz)
            aal3(ixyz,ja) = aal3(ixyz,ja) -tmpj(ixyz)
            aal3(ixyz,ka) = aal3(ixyz,ka) -tmpk(ixyz)
          enddo
#ifdef IMPULSE
          if( itot == itot_impls ) then
            ftaul(js) = ftaul(js) +tmpj(1)*tau_impls(1) &
                 +tmpj(2)*tau_impls(2) &
                 +tmpj(3)*tau_impls(3)
            ftaul(ks) = ftaul(ks) +tmpk(1)*tau_impls(1) &
                 +tmpk(2)*tau_impls(2) &
                 +tmpk(3)*tau_impls(3)
          endif
          if( jtot == itot_impls ) then
            ftaul(is) = ftaul(is) -tmpj(1)*tau_impls(1) &
                 -tmpj(2)*tau_impls(2) &
                 -tmpj(3)*tau_impls(3)
          endif
          if( ktot == itot_impls ) then
            ftaul(is) = ftaul(is) -tmpk(1)*tau_impls(1) &
                 -tmpk(2)*tau_impls(2) &
                 -tmpk(3)*tau_impls(3)
          endif
#endif
!.....stress
          do jxyz=1,3
            do ixyz=1,3
              strsl(ixyz,jxyz,ia)= strsl(ixyz,jxyz,ia) + &
                   (-0.5d0*rij(jxyz)*tmpj(ixyz) &
                   -0.5d0*rik(jxyz)*tmpk(ixyz))
              strsl(ixyz,jxyz,ja)=strsl(ixyz,jxyz,ja) &
                   -0.5d0*rij(jxyz)*tmpj(ixyz)
              strsl(ixyz,jxyz,ka)=strsl(ixyz,jxyz,ka) &
                   -0.5d0*rik(jxyz)*tmpk(ixyz)
            enddo
          enddo ! jxyz

        enddo  ! kk

      enddo  ! jj

    enddo ! ia

    call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,aal2,3)
    if( has_trios ) call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex, &
         lsrc,myparity,nn,mpi_world,aal3,3)
    aa(1:3,1:natm) = aa(1:3,1:natm) +aal2(1:3,1:natm) +aal3(1:3,1:natm)

    call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,strsl,9)
    strs(1:3,1:3,1:natm) = strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)

!-----gather epot
    epot1 = 0d0
    epot2 = 0d0
    epot3 = 0d0
    call mpi_allreduce(epotl2,epot2,1,mpi_real8,mpi_sum,mpi_world,ierr)
    if( has_solo ) call mpi_allreduce(epotl1,epot1,1,mpi_real8, &
         mpi_sum,mpi_world,ierr)
    if( has_trios ) call mpi_allreduce(epotl3,epot3,1,mpi_real8, &
         mpi_sum,mpi_world,ierr)
    epot= epot +epot1 +epot2 +epot3
    if( myid == 0 .and. iprint > 2 ) &
         print '(a,3es12.4)',' force_uf3l epot1,epot2,epot3 = ', &
         epot1,epot2,epot3

#ifdef CONTRIB
    if( myid == 0 ) print '(a,6es12.4)', ' epot_{elem,Li,La,LiLa,NbO,X} =', &
         epot_elem, epot_Li, epot_La, epot_LiLa, epot_NbO, epot_X
#endif

    return
  end subroutine force_uf3l
!=======================================================================
  subroutine force_uf3d(namax,natm,tag,ra,nnmax,aa,strs,h,hi &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rcin,lspr &
       ,mpi_world,myid,epi,epot,lstrs,iprint,l1st)
!
!  Decomposed UF3 implementation without using non-recursive function of b-spline.
!
    use util, only: itotOf
#ifdef IMPULSE
    use impulse,only: ftaul, itot_impls, tau_impls
#endif
    implicit none
    real(8),parameter:: tiny = 1d-8
    integer,intent(in):: namax,natm,nnmax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_world,myid,lspr(0:nnmax,namax),nex(3)
    real(8),intent(in):: ra(3,namax),tag(namax) &
         ,h(3,3),hi(3,3),sv(3,6)
    real(8),intent(inout):: rcin
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: l1st, lstrs

!.....local
    integer:: ia,ja,ka,jj,kk,l,is,nr2,n,inc,&
         nij,itot,jtot,ktot,i2b,i3b,js,ks,jsp,ksp,ierr, &
         ixyz,jxyz,lcs,lij,ncs,nij3,nik3,lik,nik
    real(8):: epotl2,epotl3,epot2,epot3,tmp,tmp2,bij(-3:0),dbij(-3:0), &
         bcs(-3:0),dbcs(-3:0),c2t,c3t,epotl1,epot1,fac3b, &
         bij3(-3:0),dbij3(-3:0),bik3(-3:0),dbik3(-3:0)
    real(8):: xi(3),xj(3),xk(3),xij(3),xik(3),xjk(3),rij(3),rik(3),&
         rjk(3),dij2,dij,dik2,dik,drijj(3),drikk(3),&
         drjkk(3),diji,diki,drijc,drikc,dv3csn,dv3rij,dv3rik, &
         sumcb,sumcdb,sumcbij,sumcdbij,sumcbik,sumcdbik, &
         c3ij,c3ik
    real(8):: dcsnj(3),dcsnk(3),dcsni(3),tmpj(3),tmpk(3),csn
    real(8):: rcij, rcik, rcij2, rcik2
    real(8),save:: rcin2 = -1d0

    type(prm2):: p2
    type(prm3d):: p3

    if( rcin2 < 0d0 ) then  ! Probably it is the 1st call...
      if( rcin < rcmax ) then
        if( myid == 0 ) then
          write(6,'(1x,a)') "ERROR: Cutoff radius is not appropriate !!!"
          write(6,'(1x,a,f0.3)') "  rc should be longer than ", rcmax
        endif
        call mpi_finalize(ierr)
        stop
      endif
      rcin2 = rcin*rcin
    endif

    if( .not.allocated(aal2) ) then
      allocate(aal2(3,namax),aal3(3,namax),strsl(3,3,namax),ls3b(0:nnmax))
!!$      allocate(aal2(3,namax),aal3(3,namax),strsl(3,3,namax))
    endif
    if( size(aal2) < 3*namax ) then
      deallocate(aal2,aal3,strsl)
      allocate(aal2(3,namax),aal3(3,namax),strsl(3,3,namax))
    endif
    if( size(ls3b) < nnmax+1 ) then
      deallocate(ls3b)
      allocate(ls3b(0:nnmax))
    endif

    aal2(:,:) = 0d0
    aal3(:,:) = 0d0
    strsl(:,:,:) = 0d0
    epotl1 = 0d0
    epotl2 = 0d0
    epotl3 = 0d0
#ifdef IMPULSE
    ftaul(:) = 0d0
#endif
    do ia=1,natm
      is = int(tag(ia))
      itot = itotOf(tag(ia))
      epi(ia) = epi(ia) +erg1s(is)
      epotl1 = epotl1 +erg1s(is)
      xi(1:3) = ra(1:3,ia)
      do jj=1,lspr(0,ia)
        ja = lspr(jj,ia)
!!$        if( ja <= ia ) cycle
        js = int(tag(ja))
        jtot = itotOf(tag(ja))
!.....Pair terms
        i2b = interact2(is,js)
        if( i2b <= 0 ) cycle
        p2 = prm2s(i2b)
        xj(1:3) = ra(1:3,ja)
        xij(1:3) = xj(1:3) -xi(1:3)
        rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij2 = rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3)
        if( dij2 > p2%rc2 ) cycle
        dij = sqrt(dij2)
!.....NOTE: ls3b could be source of bugs and stop using ls3b.
!           To speed up by reducing pairs for 3-body, use other way
!           that is safer than this approach.
!!$!.....Make short-distance pair-list for 3-body term
!!$        if( has_trios ) then
!!$          if( dij2 < rc2_3b(is,js) ) then
!!$            ls3b(0) = ls3b(0) +1
!!$            ls3b(ls3b(0)) = ja
!!$          endif
!!$        endif
        drijj(1:3) = rij(1:3)/dij
        call b_spl(dij,p2%knots,p2%nknot,nr2,bij,dbij)
        do lij = -3,0
          n = nr2 +lij
          if( n < 1 .or. n > p2%nknot-4 ) cycle
          c2t = p2%coefs(n)
          tmp = c2t *bij(lij)
!.....Energy
          epi(ia) = epi(ia) +tmp
          epotl2 = epotl2 +tmp
!.....Forces
          tmp2 = c2t *dbij(lij)
          do ixyz=1,3
            aal2(ixyz,ia) = aal2(ixyz,ia) +drijj(ixyz)*tmp2
            aal2(ixyz,ja) = aal2(ixyz,ja) -drijj(ixyz)*tmp2
          enddo
#ifdef IMPULSE
          if( itot == itot_impls ) then
            ftaul(js) = ftaul(js) +tmp2 &
                 *(drijj(1)*tau_impls(1) &
                 +drijj(2)*tau_impls(2) &
                 +drijj(3)*tau_impls(3) )
          endif
          if( jtot == itot_impls ) then
            ftaul(is) = ftaul(is) -tmp2 &
                 *(drijj(1)*tau_impls(1) &
                 +drijj(2)*tau_impls(2) &
                 +drijj(3)*tau_impls(3) )
          endif
#endif
!.....Stresses
          do ixyz=1,3
            do jxyz=1,3
              strsl(jxyz,ixyz,ia)= strsl(jxyz,ixyz,ia) &
                   -0.5d0 *tmp2*rij(ixyz)*drijj(jxyz)
              strsl(jxyz,ixyz,ja)= strsl(jxyz,ixyz,ja) &
                   -0.5d0 *tmp2*rij(ixyz)*drijj(jxyz)
            enddo
          enddo
        enddo  ! lij
      enddo  ! jj

!.....Trio part is separated from pair part,
!.....which may be slower because of double computation of dij,
!.....but this code is a bit simpler.
      if( .not.has_trios ) cycle

      do jj=1,lspr(0,ia)-1
        ja = lspr(jj,ia)
        js = int(tag(ja))
        jtot = itotOf(tag(ja))
        xj(1:3)= ra(1:3,ja)
        xij(1:3)= xj(1:3) - xi(1:3)
        rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij2 = rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3)
        if( dij2 > rc3max2 ) cycle
        dij = sqrt(dij2)
        diji = 1d0/dij
        drijj(1:3) = rij(1:3)*diji
        
        do kk=jj+1,lspr(0,ia)
          ka = lspr(kk,ia)
          ks = int(tag(ka))
          i3b = interact3(is,js,ks)
          if( i3b <= 0 ) cycle
          p3 = prm3ds(i3b)
          rcij = p3%rcij
          rcij2 = p3%rcij2
          rcik = p3%rcik
          rcik2 = p3%rcik2
          if( js==ks ) then
            rcij = (rcij+rcik)/2
            rcik = rcij
            rcij2= (rcij2+rcik2)/2
            rcik2= rcij2
          endif
          if(  dij2 > rcij2 ) exit
          ktot = itotOf(tag(ka))
!!$          fac3b = 1d0
!!$          if( js == ks ) fac3b = 0.5d0
          xk(1:3) = ra(1:3,ka)
          xik(1:3) = xk(1:3) -xi(1:3)
          rik(1:3) = h(1:3,1)*xik(1) +h(1:3,2)*xik(2) +h(1:3,3)*xik(3)
          dik2 = rik(1)*rik(1) +rik(2)*rik(2) +rik(3)*rik(3)
          if( dik2 > rcik2 ) cycle
          dik = sqrt(dik2)
          diki = 1d0/dik
          drikk(1:3) = rik(1:3)*diki
          drijc = 1d0/(dij-rcij)
          drikc = 1d0/(dik-rcik)
!.....r_ij term
          call b_spl(dij, p3%knij, p3%nknij, nij3, bij3, dbij3)
          sumcbij = 0d0
          sumcdbij= 0d0
          do lij = -3,0
            nij = nij3 +lij
            if( nij < 1 .or. nij > p3%nknij-4 ) cycle
            c3ij = p3%cfij(nij)
            sumcbij = sumcbij +c3ij*bij3(lij)
            sumcdbij= sumcdbij +c3ij*dbij3(lij)
          enddo
!.....r_ik term
          call b_spl(dik, p3%knik, p3%nknik, nik3, bik3, dbik3)
          sumcbik = 0d0
          sumcdbik= 0d0
          do lik = -3,0
            nik = nik3 +lik
            if( nik < 1 .or. nik > p3%nknik-4 ) cycle
            c3ik = p3%cfik(nik)
            sumcbik = sumcbik +c3ik*bik3(lik)
            sumcdbik= sumcdbik +c3ik*dbik3(lik)
          enddo
!.....Cos term
          csn = (rij(1)*rik(1) +rij(2)*rik(2) +rij(3)*rik(3)) *(diji*diki)
          csn = max(min(csn, 1d0-tiny), -1d0+tiny)
          call b_spl(-csn, p3%kncs, p3%nkncs, ncs, bcs, dbcs)
          sumcb = 0d0
          sumcdb= 0d0
          do lcs = -3,0
            n = ncs +lcs
            if( n < 1 .or. n > p3%ncfcs ) cycle
            c3t = p3%cfcs(n) !*fac3b
            sumcb = sumcb +c3t*bcs(lcs)
            sumcdb= sumcdb +c3t*dbcs(lcs)
          enddo  ! lcs
!.....energy
          tmp = sumcbij * sumcbik * sumcb
          epi(ia) = epi(ia) +tmp
          epotl3 = epotl3 +tmp
!.....force
!!$          dv3rij = -drijc*drijc *gmj *tmp
!!$          dv3rik = -drikc*drikc *gmk *tmp
!!$          dv3csn = -vexp *sumcdb
          dv3rij = sumcdbij *sumcbik *sumcb
          dv3rik = sumcdbik *sumcbij *sumcb
          dv3csn = -sumcdb *sumcbij *sumcbik
          dcsnj(1:3)= (-rij(1:3)*csn*(diji*diji) +rik(1:3)*(diji*diki))
          dcsnk(1:3)= (-rik(1:3)*csn*(diki*diki) +rij(1:3)*(diji*diki))
!!$          dcsni(1:3)= -dcsnj(1:3) -dcsnk(1:3)
          tmpj(1:3)= dv3rij*drijj(1:3) +dv3csn*dcsnj(1:3)
          tmpk(1:3)= dv3rik*drikk(1:3) +dv3csn*dcsnk(1:3)
          do ixyz=1,3
            aal3(ixyz,ia) = aal3(ixyz,ia) +tmpj(ixyz) +tmpk(ixyz)
            aal3(ixyz,ja) = aal3(ixyz,ja) -tmpj(ixyz)
            aal3(ixyz,ka) = aal3(ixyz,ka) -tmpk(ixyz)
          enddo
#ifdef IMPULSE
          if( itot == itot_impls ) then
            ftaul(js) = ftaul(js) +tmpj(1)*tau_impls(1) &
                 +tmpj(2)*tau_impls(2) &
                 +tmpj(3)*tau_impls(3)
            ftaul(ks) = ftaul(ks) +tmpk(1)*tau_impls(1) &
                 +tmpk(2)*tau_impls(2) &
                 +tmpk(3)*tau_impls(3)
          endif
          if( jtot == itot_impls ) then
            ftaul(is) = ftaul(is) -tmpj(1)*tau_impls(1) &
                 -tmpj(2)*tau_impls(2) &
                 -tmpj(3)*tau_impls(3)
          endif
          if( ktot == itot_impls ) then
            ftaul(is) = ftaul(is) -tmpk(1)*tau_impls(1) &
                 -tmpk(2)*tau_impls(2) &
                 -tmpk(3)*tau_impls(3)
          endif
#endif
!.....stress
          do jxyz=1,3
            do ixyz=1,3
              strsl(ixyz,jxyz,ia)= strsl(ixyz,jxyz,ia) + &
                   (-0.5d0*rij(jxyz)*tmpj(ixyz) &
                   -0.5d0*rik(jxyz)*tmpk(ixyz))
              strsl(ixyz,jxyz,ja)=strsl(ixyz,jxyz,ja) &
                   -0.5d0*rij(jxyz)*tmpj(ixyz)
              strsl(ixyz,jxyz,ka)=strsl(ixyz,jxyz,ka) &
                   -0.5d0*rik(jxyz)*tmpk(ixyz)
            enddo
          enddo ! jxyz

        enddo  ! kk

      enddo  ! jj

    enddo ! ia

    call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,aal2,3)
    if( has_trios ) call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex, &
         lsrc,myparity,nn,mpi_world,aal3,3)
    aa(1:3,1:natm) = aa(1:3,1:natm) +aal2(1:3,1:natm) +aal3(1:3,1:natm)

    call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,strsl,9)
    strs(1:3,1:3,1:natm) = strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)

!-----gather epot
    epot1 = 0d0
    epot2 = 0d0
    epot3 = 0d0
    call mpi_allreduce(epotl2,epot2,1,mpi_real8,mpi_sum,mpi_world,ierr)
    if( has_solo ) call mpi_allreduce(epotl1,epot1,1,mpi_real8, &
         mpi_sum,mpi_world,ierr)
    if( has_trios ) call mpi_allreduce(epotl3,epot3,1,mpi_real8, &
         mpi_sum,mpi_world,ierr)
    epot= epot +epot1 +epot2 +epot3
    if( myid == 0 .and. iprint > 2 ) &
         print '(a,3es12.4)',' force_uf3d epot1,epot2,epot3 = ', &
         epot1,epot2,epot3

    return
  end subroutine force_uf3d
!=======================================================================
  subroutine force_uf3_rec(namax,natm,tag,ra,nnmax,aa,strs,h,hi &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rcin,lspr &
       ,mpi_world,myid,epi,epot,lstrs,iprint,l1st)
!
!  Recursive implementation of force_uf3.
!
    use util, only: itotOf
    implicit none
    integer,intent(in):: namax,natm,nnmax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_world,myid,lspr(0:nnmax,namax),nex(3)
    real(8),intent(in):: ra(3,namax),tag(namax) &
         ,h(3,3),hi(3,3),sv(3,6)
    real(8),intent(inout):: rcin
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: l1st
    logical:: lstrs

!.....local
    integer:: ia,ja,ka,jj,kk,l,is,nr2,n,nij3,inc,nik3,njk3,&
         nik,njk,nij,itot,jtot,ktot,i2b,i3b,js,ks,jsp,ksp,ierr, &
         ixyz,jxyz,lij,lik,ljk
    real(8):: epotl2,epotl3,epot2,epot3,tmp,tmp2,bij,dbij, &
         bij3(4),dbij3(4),bik3,dbik3,bjk3,dbjk3, &
         c2t,c3t
    real(8):: xi(3),xj(3),xk(3),xij(3),xik(3),xjk(3),rij(3),rik(3),&
         rjk(3),dij2,dij,dik2,dik,djk2,djk,drijj(3),drikk(3),&
         drjkk(3),tmpij(3),tmpik(3),tmpjk(3)
    real(8),save:: rcin2

    type(prm2):: p2
    type(prm3):: p3

    if( l1st ) then
      if( allocated(aal2) ) deallocate(aal2,strsl)
      allocate(aal2(3,namax),aal3(3,namax),strsl(3,3,namax))
      if( rcin < rcmax ) then
        if( myid == 0 ) then
          write(6,'(1x,a)') "ERROR: Cutoff radius is not appropriate !!!"
          write(6,'(1x,a,f0.3)') "  rc should be longer than ", rcmax
        endif
        call mpi_finalize(ierr)
        stop
      endif
      rcin2 = rcin*rcin
    endif

    if( size(aal2) < 3*namax ) then
      deallocate(aal2,aal3,strsl)
      allocate(aal2(3,namax),aal3(3,namax),strsl(3,3,namax))
    endif

    aal2(:,:) = 0d0
    aal3(:,:) = 0d0
    strsl(1:3,1:3,1:namax) = 0d0
    epotl2 = 0d0
    epotl3 = 0d0
!$omp parallel
!$omp do private(ia,is,xi,jj,ja,js,i2b,p2,xj,xij,rij,dij2,dij,drijj,nr2, &
!$omp      n,bij,c2t,tmp,dbij,tmp2,ixyz,jxyz, &
!$omp      jsp,ksp,i3b,p3,nij3,inc,bij3,dbij3, &
!$omp      kk,ka,ks,xk,xik,xjk,rik,dik2,dik,rjk,djk2,djk, &
!$omp      drikk,drjkk,nik3,njk3,nik,bik3,dbik3,njk,bjk3,dbjk3, &
!$omp      l,nij,c3t,tmpij,tmpik,tmpjk) &
!$omp      reduction(+:epotl2,epotl3)
    do ia=1,natm
      is = int(tag(ia))
      xi(1:3) = ra(1:3,ia)
      do jj=1,lspr(0,ia)
        ja = lspr(jj,ia)
!!$        if( ja <= ia ) cycle
        js = int(tag(ja))
!.....Pair terms
        i2b = interact2(is,js)
        if( i2b <= 0 ) cycle
        p2 = prm2s(i2b)
        xj(1:3) = ra(1:3,ja)
        xij(1:3) = xj(1:3) -xi(1:3)
        rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij2 = rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3)
        if( dij2 > p2%rc2 ) cycle
        dij = sqrt(dij2)
        drijj(1:3) = rij(1:3)/dij
!.....Look for nr2 from knots data
        nr2 = knot_index(dij, p2%nknot, p2%knots)
        do n=max(1,nr2-3),min(nr2,p2%nknot-4)
          bij = b_spl_rec(n,3,dij,p2%knots,p2%nknot)
          c2t = p2%coefs(n)
          tmp = c2t *bij
          epi(ia) = epi(ia) +tmp
          epotl2 = epotl2 + tmp
!.....Forces
          dbij = db_spl(n,3,dij,p2%knots,p2%nknot)
          tmp2 = dbij*c2t
          do ixyz=1,3
!$omp atomic
            aal2(ixyz,ia) = aal2(ixyz,ia) +drijj(ixyz)*tmp2
!$omp atomic
            aal2(ixyz,ja) = aal2(ixyz,ja) -drijj(ixyz)*tmp2
          enddo
!.....Stresses
          do ixyz=1,3
            do jxyz=1,3
!$omp atomic
              strsl(jxyz,ixyz,ia)= strsl(jxyz,ixyz,ia) &
                   -0.5d0 *tmp2*rij(ixyz)*drijj(jxyz)
!$omp atomic
              strsl(jxyz,ixyz,ja)= strsl(jxyz,ixyz,ja) &
                   -0.5d0 *tmp2*rij(ixyz)*drijj(jxyz)
            enddo
          enddo
        enddo
      enddo

!.....Trio part is separated from pair part,
!.....which may be slower because of double computation of dij,
!.....but this code is a bit simpler.
      if( .not.has_trios ) cycle
      do jsp=1,nspmax
        do ksp=jsp,nspmax
          i3b = interact3(is,jsp,ksp)
          if( i3b <= 0 ) cycle
          p3 = prm3s(i3b)
          do jj=1,lspr(0,ia)
            ja = lspr(jj,ia)
            js = int(tag(ja))
            if( js /= jsp ) cycle
            xj(1:3) = ra(1:3,ja)
            xij(1:3) = xj(1:3) -xi(1:3)
            rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
            dij2 = rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3)
            if( dij2 > p3%rcij2 ) cycle
            dij = sqrt(dij2)
            drijj(1:3) = rij(1:3)/dij
!!$            call b_spl(dij, p3%knij, p3%nknij, nij3, bij3, dbij3)
            nij3 = knot_index(dij, p3%nknij, p3%knij)
            inc = 0
            bij3(:) = 0d0
            dbij3(:) = 0d0
            do n=nij3-3,nij3
              inc = inc +1
              if( n < 1 .or. n > p3%nknij-4 ) cycle
              bij3(inc) = b_spl_rec(n,3,dij, p3%knij, p3%nknij)
              dbij3(inc) = db_spl(n,3,dij, p3%knij, p3%nknij)
            enddo

            do kk=1,lspr(0,ia)
              ka = lspr(kk,ia)
              if( jsp == ksp .and. ka <= ja ) cycle
              ks = int(tag(ka))
              if( ks /= ksp ) cycle
              xk(1:3) = ra(1:3,ka)
              xik(1:3) = xk(1:3) -xi(1:3)
              xjk(1:3) = xk(1:3) -xj(1:3)
              rik(1:3) = h(1:3,1)*xik(1) +h(1:3,2)*xik(2) +h(1:3,3)*xik(3)
              dik2 = rik(1)*rik(1) +rik(2)*rik(2) +rik(3)*rik(3)
              if( dik2 > p3%rcik2 ) cycle
              dik = sqrt(dik2)
              rjk(1:3) = h(1:3,1)*xjk(1) +h(1:3,2)*xjk(2) +h(1:3,3)*xjk(3)
              djk2 = rjk(1)*rjk(1) +rjk(2)*rjk(2) +rjk(3)*rjk(3)
              if( djk2 > p3%rcjk2 ) cycle
              djk = sqrt(djk2)
              nik3 = knot_index(dik, p3%nknik, p3%knik)
              njk3 = knot_index(djk, p3%nknjk, p3%knjk)
!.....B-spline part
              do nik=nik3-3,nik3
                if( nik < 1 .or. nik > p3%nknik-4 ) cycle
                bik3 = b_spl_rec(nik,3,dik, p3%knik, p3%nknik)
                dbik3 = db_spl(nik,3,dik, p3%knik, p3%nknik)
                do njk=njk3-3,njk3
                  if( njk < 1 .or. njk > p3%nknjk-4 ) cycle
                  bjk3 = b_spl_rec(njk,3,djk, p3%knjk, p3%nknjk)
                  dbjk3 = db_spl(njk,3,djk, p3%knjk, p3%nknjk)
                  l = 0
                  do nij=nij3-3,nij3
                    l = l +1
                    if( nij < 1 .or. nij > p3%nknij-4 ) cycle
!.....Energy
                    c3t = p3%coefs(njk,nik,nij)
                    tmp = c3t*bij3(l)*bik3*bjk3
                    epi(ia) = epi(ia) +tmp
                    epotl3 = epotl3 +tmp
!.....Force
                    tmpij(1:3) = dbij3(l)*bik3*bjk3*drijj(1:3)
                    tmpik(1:3) = bij3(l)*dbik3*bjk3*drikk(1:3)
                    tmpjk(1:3) = bij3(l)*bik3*dbjk3*drjkk(1:3)
                    do ixyz=1,3
!$omp atomic
                      aal3(ixyz,ia)= aal3(ixyz,ia) +c3t*(tmpij(ixyz) +tmpik(ixyz))
!$omp atomic
                      aal3(ixyz,ja)= aal3(ixyz,ja) +c3t*(-tmpij(ixyz) +tmpjk(ixyz))
!$omp atomic
                      aal3(ixyz,ka)= aal3(ixyz,ka) +c3t*(-tmpik(ixyz) -tmpjk(ixyz))
                    enddo
!.....Stresses
                    do jxyz=1,3
                      do ixyz=1,3
!$omp atomic
                        strsl(ixyz,jxyz,ia)= strsl(ixyz,jxyz,ia) &
                             -0.5d0 *c3t *(xij(jxyz)*tmpij(ixyz) +xik(jxyz)*tmpik(ixyz))
!$omp atomic
                        strsl(ixyz,jxyz,ja)= strsl(ixyz,jxyz,ja) &
                             -0.5d0 *c3t *(xij(jxyz)*tmpij(ixyz) +xjk(jxyz)*tmpjk(ixyz))
!$omp atomic
                        strsl(ixyz,jxyz,ka)= strsl(ixyz,jxyz,ka) &
                             -0.5d0 *c3t *(xik(jxyz)*tmpik(ixyz) +xjk(jxyz)*tmpjk(ixyz))
                      enddo
                    enddo
                  enddo ! nij
                enddo ! njk
              enddo ! nik

            enddo  ! kk
          enddo  ! jj
        enddo  ! ksp
      enddo  ! jsp

    enddo ! ia
!$omp end do
!$omp end parallel

    call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,aal2,3)
    if( has_trios ) call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,aal3,3)
    aa(1:3,1:natm) = aa(1:3,1:natm) +aal2(1:3,1:natm) +aal3(1:3,1:natm)

    call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,strsl,9)
    strs(1:3,1:3,1:natm) = strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)

!-----gather epot
    epot2 = 0d0
    epot3 = 0d0
    call mpi_allreduce(epotl2,epot2,1,mpi_real8,mpi_sum,mpi_world,ierr)
    if( has_trios ) call mpi_allreduce(epotl3,epot3,1,mpi_real8,mpi_sum,mpi_world,ierr)
    epot= epot +epot2 +epot3
    if( myid == 0 .and. iprint > 2 ) &
         print '(a,2es12.4)',' force_uf3 epot2,epot3 = ',epot2,epot3

    return
  end subroutine force_uf3_rec
!=======================================================================
  subroutine gradw_uf3(namax,natm,tag,ra,nnmax,h,rcin,lspr, &
       iprint,ndimp,gwe,gwf,gws,lematch,lfmatch,lsmatch,iprm0, &
       lgrad_done,nfcal,lfrc_eval)
!
!  Gradient of UF3 wrt weights.
!  Note: This routine is always called in single run,
!  thus no need of parallel implementation.
!  - iprm0: The starting point -1 in parameter array for this FF.
!
    use util, only: itotOf
    implicit none
    integer,intent(in):: namax,natm,nnmax,iprint,iprm0
    integer,intent(in):: lspr(0:nnmax,namax)
    real(8),intent(in):: ra(3,namax),tag(namax),h(3,3)
    real(8),intent(inout):: rcin
    integer,intent(in):: ndimp
    integer,intent(in):: nfcal
    logical,intent(in):: lfrc_eval(natm)
    real(8),intent(inout):: gwe(ndimp),gwf(3,ndimp,nfcal),gws(6,ndimp)
    logical,intent(in):: lematch,lfmatch,lsmatch,lgrad_done

!.....local
    integer:: i,j,k,ia,ja,ka,jj,kk,l,is,nr2,n,nij3,inc,nik3,njk3,&
         nik,njk,nij,itot,jtot,ktot,i1b,i2b,i3b,js,ks,jsp,ksp,ierr, &
         ixyz,jxyz,lij,lik,ljk,ip,jra,kra,ij,ik,jk
    integer:: ifcal,jfcal,kfcal
    real(8):: epotl2,epotl3,epot2,epot3,tmp,tmp2,bij(-3:0),dbij(-3:0), &
         bij3(-3:0),dbij3(-3:0),bik3(-3:0),dbik3(-3:0),bjk3(-3:0), &
         dbjk3(-3:0),c2t,c3t,epotl1,epot1,ttmp,fac3b
    real(8):: xi(3),xj(3),xk(3),xij(3),xik(3),xjk(3),rij(3),rik(3),&
         rjk(3),dij2,dij,dik2,dik,djk2,djk,drijj(3),drikk(3),&
         drjkk(3),tmpij(3),tmpik(3),tmpjk(3)
    real(8),save:: rcin2 = -1d0
    integer,save,allocatable:: ia2ifcal(:)

    type(prm2):: p2
    type(prm3):: p3

    if( lgrad_done ) goto 10

    if( .not.allocated(ia2ifcal) ) allocate(ia2ifcal(namax))
    if( rcin2 < 0d0 ) then
      if( rcin < rcmax ) then
        call mpi_finalize(ierr)
        stop
      endif
      rcin2 = rcin*rcin
    endif

    if( .not.allocated(prm2s(1)%gwe) &
         .or. size(prm2s(1)%gwf).lt.prm2s(1)%ncoef*3*nfcal) then
      tmp = 0d0
      do i2b=1,n2b
        p2 = prm2s(i2b)
        if( lematch ) then
          if( allocated(prm2s(i2b)%gwe) ) deallocate(prm2s(i2b)%gwe)
          allocate(prm2s(i2b)%gwe(p2%ncoef))
          tmp = tmp +size(prm2s(i2b)%gwe)*8d0
        endif
        if( lfmatch ) then
          if( allocated(prm2s(i2b)%gwf) ) deallocate(prm2s(i2b)%gwf)
          allocate(prm2s(i2b)%gwf(3,p2%ncoef,nfcal))
          tmp = tmp +size(prm2s(i2b)%gwf)*8d0
        endif
        if( lsmatch ) then
          if( allocated(prm2s(i2b)%gws) ) deallocate(prm2s(i2b)%gws)
          allocate(prm2s(i2b)%gws(6,p2%ncoef))
          tmp = tmp +size(prm2s(i2b)%gws)*8d0
        endif
      enddo
      do i3b=1,n3b
        p3 = prm3s(i3b)
        if( lematch ) then
          if( allocated(prm3s(i3b)%gwe) ) deallocate(prm3s(i3b)%gwe)
          allocate(prm3s(i3b)%gwe(p3%ncfjk,p3%ncfik,p3%ncfij))
          tmp = tmp +size(prm3s(i3b)%gwe)*8d0
        endif
        if( lfmatch ) then
          if( allocated(prm3s(i3b)%gwf) ) deallocate(prm3s(i3b)%gwf)
          allocate(prm3s(i3b)%gwf(3,p3%ncfjk,p3%ncfik,p3%ncfij,nfcal))
          tmp = tmp +size(prm3s(i3b)%gwf)*8d0
        endif
        if( lsmatch ) then
          if( allocated(prm3s(i3b)%gws) ) deallocate(prm3s(i3b)%gws)
          allocate(prm3s(i3b)%gws(6,p3%ncfjk,p3%ncfik,p3%ncfij))
          tmp = tmp +size(prm3s(i3b)%gws)*8d0
        endif
      enddo
    endif

!!$    if( size(ls3b) < nnmax+1 ) then
!!$      deallocate(ls3b)
!!$      allocate(ls3b(0:nnmax))
!!$    endif
    if( size(ia2ifcal) < namax ) then
      deallocate(ia2ifcal)
      allocate(ia2ifcal(namax))
    endif
!.....Construct ia2ifcal
    ifcal = 0
    ia2ifcal(:) = 0
    do ia=1,natm
      if( lfrc_eval(ia) ) then
        ifcal = ifcal +1
        ia2ifcal(ia) = ifcal
      endif
    enddo

    if( lematch ) gerg1s(:) = 0d0
    do i2b=1,n2b
      if( lematch ) prm2s(i2b)%gwe(:) = 0d0
      if( lfmatch ) prm2s(i2b)%gwf(:,:,:) = 0d0
      if( lsmatch ) prm2s(i2b)%gws(:,:) = 0d0
    enddo
    do i3b=1,n3b
      if( lematch ) prm3s(i3b)%gwe(:,:,:) = 0d0
      if( lfmatch ) prm3s(i3b)%gwf(:,:,:,:,:) = 0d0
      if( lsmatch ) prm3s(i3b)%gws(:,:,:,:) = 0d0
    enddo

    ttmp = mpi_wtime()

    do ia=1,natm
      is = int(tag(ia))
!!$      epi(ia) = epi(ia) +erg1s(is)
!!$      epotl1 = epotl1 +erg1s(is)
      gerg1s(is) = gerg1s(is) +1d0
      xi(1:3) = ra(1:3,ia)
      do jj=1,lspr(0,ia)
        ja = lspr(jj,ia)
!!$        if( ja <= ia ) cycle
        js = int(tag(ja))
        jra = itotOf(tag(ja))
!.....Pair terms
        i2b = interact2(is,js)
        if( i2b <= 0 ) cycle
        p2 = prm2s(i2b)
        xj(1:3) = ra(1:3,ja)
        xij(1:3) = xj(1:3) -xi(1:3)
        rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij2 = rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3)
        if( dij2 > p2%rc2 ) cycle
        dij = sqrt(dij2)
!!$!.....Pair-list for 3-body cannot be used for gwf
!!$!     because it uses ia2ifcal
!!$        if( has_trios ) then
!!$          if( dij2 < rc2_3b(is,js) ) then
!!$            ls3b(0) = ls3b(0) +1
!!$            ls3b(ls3b(0)) = ja
!!$          endif
!!$        endif
        drijj(1:3) = rij(1:3)/dij
        call b_spl(dij,p2%knots,p2%nknot,nr2,bij,dbij)
        do lij = -3,0
          n = nr2 +lij
          if( n < 1 .or. n > p2%nknot-4 ) cycle
!!$          c2t = p2%coefs(n)
!.....Energy
          if( lematch ) prm2s(i2b)%gwe(n) = prm2s(i2b)%gwe(n) +bij(lij)
!.....Forces
          if( lfmatch ) then
            tmp2 = dbij(lij)
            ifcal = ia2ifcal(ia)
            jfcal = ia2ifcal(jra)
            if( ifcal.ne.0 ) prm2s(i2b)%gwf(:,n,ifcal) &
                 = prm2s(i2b)%gwf(:,n,ifcal) +drijj(:)*tmp2
            if( jfcal.ne.0 ) prm2s(i2b)%gwf(:,n,jfcal) &
                 = prm2s(i2b)%gwf(:,n,jfcal) -drijj(:)*tmp2
!!$            prm2s(i2b)%gwf(:,n,ia)  = prm2s(i2b)%gwf(:,n,ia) +drijj(:)*tmp2
!!$            prm2s(i2b)%gwf(:,n,jra) = prm2s(i2b)%gwf(:,n,jra) -drijj(:)*tmp2
          endif
!.....Stresses
          if( lsmatch ) then
            tmp2 = dbij(lij)
            prm2s(i2b)%gws(1,n)=prm2s(i2b)%gws(1,n) -rij(1)*drijj(1)*tmp2
            prm2s(i2b)%gws(2,n)=prm2s(i2b)%gws(2,n) -rij(2)*drijj(2)*tmp2
            prm2s(i2b)%gws(3,n)=prm2s(i2b)%gws(3,n) -rij(3)*drijj(3)*tmp2
            prm2s(i2b)%gws(4,n)=prm2s(i2b)%gws(4,n) -rij(2)*drijj(3)*tmp2
            prm2s(i2b)%gws(5,n)=prm2s(i2b)%gws(5,n) -rij(1)*drijj(3)*tmp2
            prm2s(i2b)%gws(6,n)=prm2s(i2b)%gws(6,n) -rij(1)*drijj(2)*tmp2
          endif
        enddo  ! lij
      enddo

!.....Trio part is separated from pair part,
!.....which may be slower because of double computation of dij,
!.....but this code is a bit simpler.
      if( .not.has_trios ) cycle
      do jsp=1,nspmax
        do ksp=jsp,nspmax
          i3b = interact3(is,jsp,ksp)
          if( i3b <= 0 ) cycle
          p3 = prm3s(i3b)
          do jj=1,lspr(0,ia)
            ja = lspr(jj,ia)
            js = int(tag(ja))
            if( js /= jsp ) cycle
            jra = itotOf(tag(ja))
            xj(1:3) = ra(1:3,ja)
            xij(1:3) = xj(1:3) -xi(1:3)
            rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
            dij2 = rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3)
            if( dij2 > p3%rcij2 ) cycle
            dij = sqrt(dij2)
            drijj(1:3) = rij(1:3)/dij
            call b_spl(dij, p3%knij, p3%nknij, nij3, bij3, dbij3)

            do kk=1,lspr(0,ia)
              ka = lspr(kk,ia)
              if( ka == ja ) cycle
              kra = itotOf(tag(ka))
!.....Taking into account the double counting of symmetric terms
              fac3b = 1d0
              if( jsp == ksp ) fac3b = 0.5d0
              ks = int(tag(ka))
              if( ks /= ksp ) cycle
              xk(1:3) = ra(1:3,ka)
              xik(1:3) = xk(1:3) -xi(1:3)
              xjk(1:3) = xk(1:3) -xj(1:3)
              rik(1:3) = h(1:3,1)*xik(1) +h(1:3,2)*xik(2) +h(1:3,3)*xik(3)
              dik2 = rik(1)*rik(1) +rik(2)*rik(2) +rik(3)*rik(3)
              if( dik2 > p3%rcik2 ) cycle
              dik = sqrt(dik2)
              rjk(1:3) = h(1:3,1)*xjk(1) +h(1:3,2)*xjk(2) +h(1:3,3)*xjk(3)
              djk2 = rjk(1)*rjk(1) +rjk(2)*rjk(2) +rjk(3)*rjk(3)
              if( djk2 > p3%rcjk2 ) cycle
              djk = sqrt(djk2)
              drikk(1:3) = rik(1:3)/dik
              drjkk(1:3) = rjk(1:3)/djk
!.....B-spline part
              call b_spl(dik, p3%knik, p3%nknik, nik3, bik3, dbik3)
              call b_spl(djk, p3%knjk, p3%nknjk, njk3, bjk3, dbjk3)
              do lik = -3,0
                nik = nik3 +lik
                if( nik < 1 .or. nik > p3%nknik-4 ) cycle
                do ljk = -3,0
                  njk = njk3 +ljk
                  if( njk < 1 .or. njk > p3%nknjk-4 ) cycle
                  do lij = -3,0
                    nij = nij3 +lij
                    if( nij < 1 .or. nij > p3%nknij-4 ) cycle
!.....Energy
!!$                    c3t = p3%coefs(nij,nik,njk)
                    if( lematch ) then
                      tmp = bij3(lij)*bik3(lik)*bjk3(ljk)*fac3b
                      prm3s(i3b)%gwe(njk,nik,nij) = &
                           prm3s(i3b)%gwe(njk,nik,nij) +tmp
                    endif
!.....Force
                    tmpij(1:3)= dbij3(lij)* bik3(lik)* bjk3(ljk) &
                         *drijj(1:3)*fac3b
                    tmpik(1:3)=  bij3(lij)*dbik3(lik)* bjk3(ljk) &
                         *drikk(1:3)*fac3b
                    tmpjk(1:3)=  bij3(lij)* bik3(lik)*dbjk3(ljk) &
                         *drjkk(1:3)*fac3b
                    if( lfmatch ) then
                      ifcal = ia2ifcal(ia)
                      jfcal = ia2ifcal(jra)
                      kfcal = ia2ifcal(kra)
                      if( ifcal.ne.0 ) prm3s(i3b)%gwf(:,njk,nik,nij,ifcal) &
                           = prm3s(i3b)%gwf(:,njk,nik,nij,ifcal) &
                           +(tmpij(:) +tmpik(:))
                      if( jfcal.ne.0 ) prm3s(i3b)%gwf(:,njk,nik,nij,jfcal) &
                           = prm3s(i3b)%gwf(:,njk,nik,nij,jfcal) &
                           +(-tmpij(:)+tmpjk(:))
                      if( kfcal.ne.0 ) prm3s(i3b)%gwf(:,njk,nik,nij,kfcal) &
                           = prm3s(i3b)%gwf(:,njk,nik,nij,kfcal) &
                           +(-tmpik(:)-tmpjk(:))
                    endif
!.....Stresses
                    if( lsmatch ) then
                      prm3s(i3b)%gws(1,njk,nik,nij) = &
                           prm3s(i3b)%gws(1,njk,nik,nij) &
                           -(rij(1)*tmpij(1) +rik(1)*tmpik(1) &
                           +rjk(1)*tmpjk(1))
                      prm3s(i3b)%gws(2,njk,nik,nij) = &
                           prm3s(i3b)%gws(2,njk,nik,nij) &
                           -(rij(2)*tmpij(2) +rik(2)*tmpik(2) &
                           +rjk(2)*tmpjk(2))
                      prm3s(i3b)%gws(3,njk,nik,nij) = &
                           prm3s(i3b)%gws(3,njk,nik,nij) &
                           -(rij(3)*tmpij(3) +rik(3)*tmpik(3) &
                           +rjk(3)*tmpjk(3))
                      prm3s(i3b)%gws(4,njk,nik,nij) = &
                           prm3s(i3b)%gws(4,njk,nik,nij) &
                           -(rij(2)*tmpij(3) +rik(2)*tmpik(3) &
                           +rjk(2)*tmpjk(3))
                      prm3s(i3b)%gws(5,njk,nik,nij) = &
                           prm3s(i3b)%gws(5,njk,nik,nij) &
                           -(rij(1)*tmpij(3) +rik(1)*tmpik(3) &
                           +rjk(1)*tmpjk(3))
                      prm3s(i3b)%gws(6,njk,nik,nij) = &
                           prm3s(i3b)%gws(6,njk,nik,nij) &
                           -(rij(1)*tmpij(2) +rik(1)*tmpik(2) &
                           +rjk(1)*tmpjk(2))
                    endif

                  enddo  ! lij
                enddo  ! ljk
              enddo  ! lik

            enddo  ! kk
          enddo  ! jj
        enddo  ! ksp
      enddo  ! jsp
    enddo ! ia

!.....Skip to here if grad calc is already done before.
10  continue

!.....Tidy up gradient arrays
    if( lematch ) then  ! energy matching
      ip = iprm0
      do i1b=1,n1b
        ip = ip +1
        gwe(ip) = gwe(ip) +gerg1s(i1b)
      enddo
      do i2b=1,n2b
        p2 = prm2s(i2b)
        do i=1,p2%ncoef
          ip = ip +1
          gwe(ip) = gwe(ip) +p2%gwe(i)
        enddo
      enddo
      do i3b=1,n3b
        p3 = prm3s(i3b)
        do ij=1,p3%ncfij
          do ik=1,p3%ncfik
            do jk=1,p3%ncfjk
              ip = ip +1
!!$              if( (ip>=14520 .and. ip<=14530)  )&
!!$                  print '(a,4i4,i7,4es11.3)','i3b,ij,ik,jk,ip,gwe(ik,ij),gwe(ij,ik)=',&
!!$                   i3b,ij,ik,jk,ip,p3%gwe(jk,ik,ij),p3%gwe(jk,ij,ik)
!.....if jsp==ksp, coef should be symmetric such that c_lmn=c_mln.
!!$              if( p3%jsp==p3%ksp .and. ij.ne.ik ) then
!!$                gwe(ip) = gwe(ip) +(p3%gwe(jk,ik,ij)+p3%gwe(jk,ij,ik)) !/2
!!$              else
                gwe(ip) = gwe(ip) +p3%gwe(jk,ik,ij)
!!$              endif
            enddo
          enddo
        enddo
      enddo  ! i3b
    endif

    if( lfmatch ) then  ! force matching
      do ia=1,natm
        if( .not. lfrc_eval(ia) ) cycle
        ifcal = ia2ifcal(ia)
        is = int(tag(ia))
        ip = iprm0 +n1b  ! no contrib. from solo term to forces
        do i2b=1,n2b
          p2 = prm2s(i2b)
          do i=1,p2%ncoef
            ip = ip +1
            gwf(1:3,ip,ifcal) = gwf(1:3,ip,ifcal)  +p2%gwf(1:3,i,ifcal)
          enddo
        enddo  ! i2b
        do i3b=1,n3b
          p3 = prm3s(i3b)
          do ij=1,p3%ncfij
            do ik=1,p3%ncfik
              do jk=1,p3%ncfjk
                ip = ip +1
!.....if jsp==ksp, coef should be symmetric such that c_lmn=c_mln.
!!$                if( p3%jsp==p3%ksp .and. ij.ne.ik ) then
!!$                  gwf(1:3,ip,ifcal) = gwf(1:3,ip,ifcal) &
!!$                       +(p3%gwf(1:3,jk,ik,ij,ifcal) &
!!$                       +p3%gwf(1:3,jk,ij,ik,ifcal)) !/2
!!$                else
                  gwf(1:3,ip,ifcal) = gwf(1:3,ip,ifcal) &
                       +p3%gwf(1:3,jk,ik,ij,ifcal)
!!$                endif
              enddo
            enddo
          enddo
        enddo
      enddo ! ia
    endif

    if( lsmatch ) then  ! stress matching
      ip = iprm0 +n1b  ! no contrib. from solo term to forces
      do i2b=1,n2b
        p2 = prm2s(i2b)
        do i=1,p2%ncoef
          ip = ip +1
          gws(1:6,ip) = gws(1:6,ip) +p2%gws(1:6,i)
        enddo
      enddo  ! i2b
      do i3b=1,n3b
        p3 = prm3s(i3b)
        do ij=1,p3%ncfij
          do ik=1,p3%ncfik
            do jk=1,p3%ncfjk
              ip = ip +1
!.....if jsp==ksp, coef should be symmetric such that c_lmn=c_mln.
!!$              if( p3%jsp==p3%ksp .and. ij.ne.ik ) then
!!$                gws(1:6,ip) = gws(1:6,ip) &
!!$                       +(p3%gws(1:6,jk,ik,ij) &
!!$                       +p3%gws(1:6,jk,ij,ik)) !/2
!!$              else
                gws(1:6,ip) = gws(1:6,ip) +p3%gws(1:6,jk,ik,ij)
!!$              endif
            enddo
          enddo
        enddo
      enddo  ! i3b
    endif

    return
  end subroutine gradw_uf3
!=======================================================================
  subroutine gradw_uf3l(namax,natm,tag,ra,nnmax,h,rcin,lspr, &
       iprint,ndimp,gwe,gwf,gws,lematch,lfmatch,lsmatch,iprm0, &
       lgrad_done,nfcal,lfrc_eval)
!
!  Gradient of UF3L wrt weights.
!  Note: This routine is always called in single run,
!  thus no need of parallel implementation.
!  - iprm0: The starting point -1 in parameter array for this FF.
!
    use util, only: itotOf
    implicit none
    integer,intent(in):: namax,natm,nnmax,iprint,iprm0
    integer,intent(in):: lspr(0:nnmax,namax)
    real(8),intent(in):: ra(3,namax),tag(namax),h(3,3)
    real(8),intent(inout):: rcin
    integer,intent(in):: ndimp
    integer,intent(in):: nfcal
    logical,intent(in):: lfrc_eval(natm)
    real(8),intent(inout):: gwe(ndimp),gwf(3,ndimp,nfcal),gws(6,ndimp)
    logical,intent(in):: lematch,lfmatch,lsmatch,lgrad_done

!.....local
    real(8),parameter:: tiny = 1d-8
    integer:: i,ia,ja,ka,jj,kk,l,is,nr2,n,inc,&
         nij,itot,i1b,i2b,i3b,js,ks,jsp,ksp,ierr, &
         ixyz,jxyz,lcs,lij,ncs
    integer:: ifcal,jfcal,kfcal,ip,iv,jra,kra
    real(8):: epotl2,epotl3,epot2,epot3,tmp,tmp2,bij(-3:0),dbij(-3:0), &
         bcs(-3:0),dbcs(-3:0),c2t,c3t,epotl1,epot1,fac3b,tmp3
    real(8):: xi(3),xj(3),xk(3),xij(3),xik(3),xjk(3),rij(3),rik(3),&
         rjk(3),dij2,dij,dik2,dik,drijj(3),drikk(3),&
         drjkk(3),diji,diki,drijc,drikc,dv3csn,dv3rij,dv3rik,sumcb,sumcdb
    real(8):: dcsnj(3),dcsnk(3),dcsni(3),tmpj(3),tmpk(3),gmj,gmk,csn,vexp
    real(8):: dv3rijdgj,dv3rijdgk,dv3rijc, dv3rikdgj,dv3rikdgk,dv3rikc, &
         dv3csndgj,dv3csndgk,dv3csnc, &
         dv3rijdcj,dv3rijdck,dv3rikdcj,dv3rikdck,dv3csndcj,dv3csndck
    real(8):: rcij,rcik,rcij2,rcik2
    real(8),save:: rcin2 = -1d0
    integer,save,allocatable:: ia2ifcal(:)

    type(prm2):: p2
    type(prm3l):: p3

!!$    if( lgrad_done ) goto 10

    if( .not.allocated(ia2ifcal) ) allocate(ia2ifcal(namax))
    if( rcin2 < 0d0 ) then
      if( rcin < rcmax ) then
        call mpi_finalize(ierr)
        stop
      endif
      rcin2 = rcin*rcin
    endif

    if( .not.allocated(prm2s(1)%gwe) &
         .or. size(prm2s(1)%gwf).lt.prm2s(1)%ncoef*3*nfcal) then
      tmp = 0d0
      do i2b=1,n2b
        p2 = prm2s(i2b)
        if( lematch ) then
          if( allocated(prm2s(i2b)%gwe) ) deallocate(prm2s(i2b)%gwe)
          allocate(prm2s(i2b)%gwe(p2%ncoef))
          tmp = tmp +size(prm2s(i2b)%gwe)*8d0
        endif
        if( lfmatch ) then
          if( allocated(prm2s(i2b)%gwf) ) deallocate(prm2s(i2b)%gwf)
          allocate(prm2s(i2b)%gwf(3,p2%ncoef,nfcal))
          tmp = tmp +size(prm2s(i2b)%gwf)*8d0
        endif
        if( lsmatch ) then
          if( allocated(prm2s(i2b)%gws) ) deallocate(prm2s(i2b)%gws)
          allocate(prm2s(i2b)%gws(6,p2%ncoef))
          tmp = tmp +size(prm2s(i2b)%gws)*8d0
        endif
      enddo
      do i3b=1,n3b
        p3 = prm3ls(i3b)
        if( lematch ) then
          if( allocated(prm3ls(i3b)%gwe) ) deallocate(prm3ls(i3b)%gwe)
          allocate(prm3ls(i3b)%gwe(p3%ncoef+4))  ! ncoef +4 (rcij,rik,gmj & gmk)
          tmp = tmp +size(prm3ls(i3b)%gwe)*8d0
        endif
        if( lfmatch ) then
          if( allocated(prm3ls(i3b)%gwf) ) deallocate(prm3ls(i3b)%gwf)
          allocate(prm3ls(i3b)%gwf(3,p3%ncoef+4,nfcal))  ! ncoef +4
          tmp = tmp +size(prm3ls(i3b)%gwf)*8d0
        endif
        if( lsmatch ) then
          if( allocated(prm3ls(i3b)%gws) ) deallocate(prm3ls(i3b)%gws)
          allocate(prm3ls(i3b)%gws(6,p3%ncoef+4))  ! ncoef +4
          tmp = tmp +size(prm3ls(i3b)%gws)*8d0
        endif
      enddo
    endif

!!$    if( size(ls3b) < nnmax+1 ) then
!!$      deallocate(ls3b)
!!$      allocate(ls3b(0:nnmax))
!!$    endif
    if( size(ia2ifcal) < namax ) then
      deallocate(ia2ifcal)
      allocate(ia2ifcal(namax))
    endif
!.....Construct ia2ifcal
    ifcal = 0
    ia2ifcal(:) = 0
    do ia=1,natm
      if( lfrc_eval(ia) ) then
        ifcal = ifcal +1
        ia2ifcal(ia) = ifcal
      endif
    enddo

    if( lematch ) gerg1s(:) = 0d0
    do i2b=1,n2b
      if( lematch ) prm2s(i2b)%gwe(:) = 0d0
      if( lfmatch ) prm2s(i2b)%gwf(:,:,:) = 0d0
      if( lsmatch ) prm2s(i2b)%gws(:,:) = 0d0
    enddo
    do i3b=1,n3b
      if( lematch ) prm3ls(i3b)%gwe(:) = 0d0
      if( lfmatch ) prm3ls(i3b)%gwf(:,:,:) = 0d0
      if( lsmatch ) prm3ls(i3b)%gws(:,:) = 0d0
    enddo

    do ia=1,natm
      is = int(tag(ia))
      gerg1s(is) = gerg1s(is) +1d0
      xi(1:3) = ra(1:3,ia)
      do jj=1,lspr(0,ia)
        ja = lspr(jj,ia)
        js = int(tag(ja))
        jra = itotOf(tag(ja))
!.....Pair terms
        i2b = interact2(is,js)
        if( i2b <= 0 ) cycle
        p2 = prm2s(i2b)
        xj(1:3) = ra(1:3,ja)
        xij(1:3) = xj(1:3) -xi(1:3)
        rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij2 = rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3)
        if( dij2 > p2%rc2 ) cycle
        dij = sqrt(dij2)
        drijj(1:3) = rij(1:3)/dij
        call b_spl(dij,p2%knots,p2%nknot,nr2,bij,dbij)
        do lij = -3,0
          n = nr2 +lij
          if( n < 1 .or. n > p2%nknot-4 ) cycle
!.....Energy
          if( lematch ) prm2s(i2b)%gwe(n) = prm2s(i2b)%gwe(n) +bij(lij)
!.....Forces
          if( lfmatch ) then
            tmp2 = dbij(lij)
            ifcal = ia2ifcal(ia)
            jfcal = ia2ifcal(jra)
            if( ifcal.ne.0 ) prm2s(i2b)%gwf(:,n,ifcal) &
                 = prm2s(i2b)%gwf(:,n,ifcal) +drijj(:)*tmp2
            if( jfcal.ne.0 ) prm2s(i2b)%gwf(:,n,jfcal) &
                 = prm2s(i2b)%gwf(:,n,jfcal) -drijj(:)*tmp2
          endif
!.....Stresses
          if( lsmatch ) then
            tmp2 = dbij(lij)
            prm2s(i2b)%gws(1,n)=prm2s(i2b)%gws(1,n) -rij(1)*drijj(1)*tmp2
            prm2s(i2b)%gws(2,n)=prm2s(i2b)%gws(2,n) -rij(2)*drijj(2)*tmp2
            prm2s(i2b)%gws(3,n)=prm2s(i2b)%gws(3,n) -rij(3)*drijj(3)*tmp2
            prm2s(i2b)%gws(4,n)=prm2s(i2b)%gws(4,n) -rij(2)*drijj(3)*tmp2
            prm2s(i2b)%gws(5,n)=prm2s(i2b)%gws(5,n) -rij(1)*drijj(3)*tmp2
            prm2s(i2b)%gws(6,n)=prm2s(i2b)%gws(6,n) -rij(1)*drijj(2)*tmp2
          endif
        enddo  ! lij
      enddo

!.....Trio part is separated from pair part,
!.....which may be slower because of double computation of dij,
!.....but this code is a bit simpler.
      if( .not.has_trios ) cycle
      tmp3 = 0d0
      do jj=1,lspr(0,ia)
        ja = lspr(jj,ia)
        js = int(tag(ja))
        jra = itotOf(tag(ja))
        xj(1:3)= ra(1:3,ja)
        xij(1:3)= xj(1:3) - xi(1:3)
        rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij2 = rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3)
        dij = sqrt(dij2)
        diji = 1d0/dij
        drijj(1:3) = rij(1:3)*diji
        do kk=jj+1,lspr(0,ia)
          ka = lspr(kk,ia)
          ks = int(tag(ka))
          i3b = interact3(is,js,ks)
          if( i3b <= 0 ) cycle
          p3 = prm3ls(i3b)
          rcij = p3%rcij
          rcij2 = p3%rcij2
          rcik = p3%rcik
          rcik2 = p3%rcik2
          if( js==ks ) then
            rcij = (rcij+rcik)/2
            rcik = rcij
            rcij2= (rcij2+rcik2)/2
            rcik2= rcij2
          endif
          if(  dij2 > rcij2 ) exit
          kra = itotOf(tag(ka))
          xk(1:3) = ra(1:3,ka)
          xik(1:3) = xk(1:3) -xi(1:3)
          rik(1:3) = h(1:3,1)*xik(1) +h(1:3,2)*xik(2) +h(1:3,3)*xik(3)
          dik2 = rik(1)*rik(1) +rik(2)*rik(2) +rik(3)*rik(3)
          if( dik2 > rcik2 ) cycle
          dik = sqrt(dik2)
          diki = 1d0/dik
          drikk(1:3) = rik(1:3)*diki
          drijc = 1d0/(dij-rcij)
          drikc = 1d0/(dik-rcik)
!.....parameters
          gmj = p3%gmj
          gmk = p3%gmk
          if( js==ks ) then
            gmj = (gmj+gmk) /2
            gmk = gmj
          endif
!.....common terms
          vexp = dexp(gmj*drijc +gmk*drikc)
          csn = (rij(1)*rik(1) +rij(2)*rik(2) +rij(3)*rik(3)) *(diji*diki)
          csn = max(min(csn, 1d0-tiny), -1d0+tiny)
          call b_spl(-csn, p3%knots, p3%nknot, ncs, bcs, dbcs)
          sumcb = 0d0
          sumcdb= 0d0
          do lcs = -3,0
            n = ncs +lcs
            if( n < 1 .or. n > p3%ncoef ) cycle
            c3t = p3%coefs(n) !*fac3b
            sumcb = sumcb +c3t*bcs(lcs)
            sumcdb= sumcdb +c3t*dbcs(lcs)
          enddo  ! lcs
!.....energy
          tmp = vexp *sumcb
          if( lematch ) then
!.....deriv. wrt rcij and rcik
            prm3ls(i3b)%gwe(1)= prm3ls(i3b)%gwe(1) +gmj*drijc**2*tmp
            prm3ls(i3b)%gwe(2)= prm3ls(i3b)%gwe(2) +gmk*drikc**2*tmp
!.....deriv. wrt gmj and gmk
            prm3ls(i3b)%gwe(3)= prm3ls(i3b)%gwe(3) +drijc*tmp
            prm3ls(i3b)%gwe(4)= prm3ls(i3b)%gwe(4) +drikc*tmp
            do lcs=-3,0
              n = ncs +lcs
              if( n < 1 .or. n > p3%ncoef ) cycle
              prm3ls(i3b)%gwe(4+n)= prm3ls(i3b)%gwe(4+n) +vexp*bcs(lcs)
            enddo
          endif
!.....force
          if( .not.lfmatch .and. .not.lsmatch ) cycle
          dv3rij = -drijc*drijc *gmj *tmp
          dv3rik = -drikc*drikc *gmk *tmp
          dv3csn = -vexp *sumcdb
          dcsnj(1:3)= (-rij(1:3)*csn*(diji*diji) +rik(1:3)*(diji*diki))
          dcsnk(1:3)= (-rik(1:3)*csn*(diki*diki) +rij(1:3)*(diji*diki))
!!$          dcsni(1:3)= -dcsnj(1:3) -dcsnk(1:3)
          dv3rijdcj = -drijc**3 *gmj *tmp *(2d0 + gmj*drijc)
          dv3rijdck = -gmj*drijc**2 *gmk*drikc**2 *tmp
          dv3rijdgj= -drijc*drijc*(1d0+gmj*drijc)*tmp
          dv3rijdgk= -drijc*drijc*(    gmk*drikc)*tmp
          dv3rijc = -drijc*drijc*gmj*vexp !*bcs(lcs)
          dv3rikdck = -drikc**3 *gmk *tmp *(2d0 + gmk*drikc)
          dv3rikdcj = -gmk*drikc**2 *gmj*drijc**2 *tmp
          dv3rikdgj= -drikc*drikc*(    gmj*drijc)*tmp
          dv3rikdgk= -drikc*drikc*(1d0+gmk*drikc)*tmp
          dv3rikc = -drikc*drikc*gmk*vexp !*bcs(lcs)
          dv3csndcj = dv3csn *drijc**2 *gmj
          dv3csndck = dv3csn *drikc**2 *gmk
          dv3csndgj = dv3csn *drijc
          dv3csndgk = dv3csn *drikc
          dv3csnc  = -vexp !*dbcs(lcs)
          if( lfmatch ) then
            ifcal = ia2ifcal(ia)
            jfcal = ia2ifcal(jra)
            kfcal = ia2ifcal(kra)
            if( ifcal.ne.0 ) then
!.....deriv. wrt rcij
              prm3ls(i3b)%gwf(:,1,ifcal)= prm3ls(i3b)%gwf(:,1,ifcal) &
                   +drijj(:)*dv3rijdcj +dcsnj(:)*dv3csndcj &
                   +drikk(:)*dv3rikdcj +dcsnk(:)*dv3csndcj
!.....deriv. wrt rcik
              prm3ls(i3b)%gwf(:,2,ifcal)= prm3ls(i3b)%gwf(:,2,ifcal) &
                   +drijj(:)*dv3rijdck +dcsnj(:)*dv3csndck &
                   +drikk(:)*dv3rikdck +dcsnk(:)*dv3csndck
!.....deriv. wrt gmj
              prm3ls(i3b)%gwf(:,3,ifcal)= prm3ls(i3b)%gwf(:,3,ifcal) &
                   +drijj(:)*dv3rijdgj +dcsnj(:)*dv3csndgj &
                   +drikk(:)*dv3rikdgj +dcsnk(:)*dv3csndgj
!.....deriv. wrt gmk
              prm3ls(i3b)%gwf(:,4,ifcal)= prm3ls(i3b)%gwf(:,4,ifcal) &
                   +drijj(:)*dv3rijdgk +dcsnj(:)*dv3csndgk &
                   +drikk(:)*dv3rikdgk +dcsnk(:)*dv3csndgk
!.....deriv. wrt c3t
              do lcs=-3,0
                n = ncs +lcs
                if( n < 1 .or. n > p3%ncoef ) cycle
                prm3ls(i3b)%gwf(:,4+n,ifcal)= prm3ls(i3b)%gwf(:,4+n,ifcal) &
                     +drijj(:)*dv3rijc*bcs(lcs) +dcsnj(:)*dv3csnc*dbcs(lcs) &
                     +drikk(:)*dv3rikc*bcs(lcs) +dcsnk(:)*dv3csnc*dbcs(lcs)
              enddo
            endif
            if( jfcal.ne.0 ) then
              prm3ls(i3b)%gwf(:,1,jfcal)= prm3ls(i3b)%gwf(:,1,jfcal) &
                   -drijj(:)*dv3rijdcj -dcsnj(:)*dv3csndcj
              prm3ls(i3b)%gwf(:,2,jfcal)= prm3ls(i3b)%gwf(:,2,jfcal) &
                   -drijj(:)*dv3rijdck -dcsnj(:)*dv3csndck
              prm3ls(i3b)%gwf(:,3,jfcal)= prm3ls(i3b)%gwf(:,3,jfcal) &
                   -drijj(:)*dv3rijdgj -dcsnj(:)*dv3csndgj
              prm3ls(i3b)%gwf(:,4,jfcal)= prm3ls(i3b)%gwf(:,4,jfcal) &
                   -drijj(:)*dv3rijdgk -dcsnj(:)*dv3csndgk
              do lcs=-3,0
                n = ncs +lcs
                if( n < 1 .or. n > p3%ncoef ) cycle
                prm3ls(i3b)%gwf(:,4+n,jfcal)= prm3ls(i3b)%gwf(:,4+n,jfcal) &
                     -drijj(:)*dv3rijc*bcs(lcs) -dcsnj(:)*dv3csnc*dbcs(lcs)
              enddo
            endif
            if( kfcal.ne.0 ) then
              prm3ls(i3b)%gwf(:,1,kfcal)= prm3ls(i3b)%gwf(:,1,kfcal) &
                   -drikk(:)*dv3rikdcj -dcsnk(:)*dv3csndcj
              prm3ls(i3b)%gwf(:,2,kfcal)= prm3ls(i3b)%gwf(:,2,kfcal) &
                   -drikk(:)*dv3rikdck -dcsnk(:)*dv3csndck
              prm3ls(i3b)%gwf(:,3,kfcal)= prm3ls(i3b)%gwf(:,3,kfcal) &
                   -drikk(:)*dv3rikdgj -dcsnk(:)*dv3csndgj
              prm3ls(i3b)%gwf(:,4,kfcal)= prm3ls(i3b)%gwf(:,4,kfcal) &
                   -drikk(:)*dv3rikdgk -dcsnk(:)*dv3csndgk
              do lcs=-3,0
                n = ncs +lcs
                if( n < 1 .or. n > p3%ncoef ) cycle
                prm3ls(i3b)%gwf(:,4+n,kfcal)= prm3ls(i3b)%gwf(:,4+n,kfcal) &
                     -drikk(:)*dv3rikc*bcs(lcs) -dcsnk(:)*dv3csnc*dbcs(lcs)
              enddo
            endif
          endif  ! lfmatch
!.....stress
          if( lsmatch ) then
            do ixyz=1,3
              do jxyz=ixyz,3
                iv = ivoigt(ixyz,jxyz)
                prm3ls(i3b)%gws(iv,1)= prm3ls(i3b)%gws(iv,1) &
                     -rij(ixyz)*(drijj(jxyz)*dv3rijdcj+dcsnj(jxyz)*dv3csndcj) &
                     -rik(ixyz)*(drikk(jxyz)*dv3rikdcj+dcsnk(jxyz)*dv3csndcj)
                prm3ls(i3b)%gws(iv,2)= prm3ls(i3b)%gws(iv,2) &
                     -rij(ixyz)*(drijj(jxyz)*dv3rijdck+dcsnj(jxyz)*dv3csndck) &
                     -rik(ixyz)*(drikk(jxyz)*dv3rikdck+dcsnk(jxyz)*dv3csndck)
                prm3ls(i3b)%gws(iv,3)= prm3ls(i3b)%gws(iv,3) &
                     -rij(ixyz)*(drijj(jxyz)*dv3rijdgj+dcsnj(jxyz)*dv3csndgj) &
                     -rik(ixyz)*(drikk(jxyz)*dv3rikdgj+dcsnk(jxyz)*dv3csndgj)
                prm3ls(i3b)%gws(iv,4)= prm3ls(i3b)%gws(iv,4) &
                     -rij(ixyz)*(drijj(jxyz)*dv3rijdgk+dcsnj(jxyz)*dv3csndgk) &
                     -rik(ixyz)*(drikk(jxyz)*dv3rikdgk+dcsnk(jxyz)*dv3csndgk)
                do lcs=-3,0
                  n = ncs +lcs
                  if( n < 1 .or. n > p3%ncoef ) cycle
                  prm3ls(i3b)%gws(iv,4+n)= prm3ls(i3b)%gws(iv,4+n) &
                       -rij(ixyz)*(drijj(jxyz)*dv3rijc*bcs(lcs) +dcsnj(jxyz)*dv3csnc*dbcs(lcs)) &
                       -rik(ixyz)*(drikk(jxyz)*dv3rikc*bcs(lcs) +dcsnk(jxyz)*dv3csnc*dbcs(lcs))
                enddo ! lcs
              enddo
            enddo ! ixyz
          endif  ! lsmatch

        enddo  ! kk
      enddo  ! jj
    enddo ! ia

!!$!.....Skip to here if grad calc is already done before.
!!$10  continue

!.....Tidy up gradient arrays
    if( lematch ) then  ! energy matching
      ip = iprm0
      do i1b=1,n1b
        ip = ip +1
        gwe(ip) = gwe(ip) +gerg1s(i1b)
      enddo
      do i2b=1,n2b
        p2 = prm2s(i2b)
        do i=1,p2%ncoef
          ip = ip +1
          gwe(ip) = gwe(ip) +p2%gwe(i)
        enddo
      enddo
      do i3b=1,n3b
        p3 = prm3ls(i3b)
        if( p3%jsp==p3%ksp ) then
          p3%gwe(1) = (p3%gwe(1)+p3%gwe(2))/2
          p3%gwe(2) = p3%gwe(1)
          p3%gwe(3) = (p3%gwe(3)+p3%gwe(4))/2
          p3%gwe(4) = p3%gwe(3)
        endif
        ip = ip +1
        gwe(ip) = gwe(ip) +p3%gwe(1)  ! rcij
        ip = ip +1
        gwe(ip) = gwe(ip) +p3%gwe(2)  ! rcik
        ip = ip +1
        gwe(ip) = gwe(ip) +p3%gwe(3)  ! gmj
        ip = ip +1
        gwe(ip) = gwe(ip) +p3%gwe(4)  ! gmk
        do i=1,p3%ncoef
          ip = ip +1
          gwe(ip) = gwe(ip) +p3%gwe(4+i)
        enddo
      enddo  ! i3b
    endif

    if( lfmatch ) then  ! force matching
      do ia=1,natm
        if( .not. lfrc_eval(ia) ) cycle
        ifcal = ia2ifcal(ia)
        is = int(tag(ia))
        ip = iprm0 +n1b  ! no contrib. from solo term to forces
        do i2b=1,n2b
          p2 = prm2s(i2b)
          do i=1,p2%ncoef
            ip = ip +1
            gwf(1:3,ip,ifcal) = gwf(1:3,ip,ifcal)  +p2%gwf(1:3,i,ifcal)
          enddo
        enddo  ! i2b
        do i3b=1,n3b
          p3 = prm3ls(i3b)
          if( p3%jsp==p3%ksp ) then
            p3%gwf(1:3,1,ifcal) = (p3%gwf(1:3,1,ifcal)+p3%gwf(1:3,2,ifcal))/2
            p3%gwf(1:3,2,ifcal) = p3%gwf(1:3,1,ifcal)
            p3%gwf(1:3,3,ifcal) = (p3%gwf(1:3,3,ifcal)+p3%gwf(1:3,4,ifcal))/2
            p3%gwf(1:3,4,ifcal) = p3%gwf(1:3,3,ifcal)
          endif
          ip = ip +1
          gwf(1:3,ip,ifcal) = gwf(1:3,ip,ifcal) +p3%gwf(1:3,1,ifcal)
          ip = ip +1
          gwf(1:3,ip,ifcal) = gwf(1:3,ip,ifcal) +p3%gwf(1:3,2,ifcal)
          ip = ip +1
          gwf(1:3,ip,ifcal) = gwf(1:3,ip,ifcal) +p3%gwf(1:3,3,ifcal)
          ip = ip +1
          gwf(1:3,ip,ifcal) = gwf(1:3,ip,ifcal) +p3%gwf(1:3,4,ifcal)
          do i=1,p3%ncoef
            ip = ip +1
            gwf(1:3,ip,ifcal) = gwf(1:3,ip,ifcal) +p3%gwf(1:3,4+i,ifcal)
          enddo
        enddo
      enddo ! ia
    endif

    if( lsmatch ) then  ! stress matching
      ip = iprm0 +n1b  ! no contrib. from solo term to forces
      do i2b=1,n2b
        p2 = prm2s(i2b)
        do i=1,p2%ncoef
          ip = ip +1
          gws(1:6,ip) = gws(1:6,ip) +p2%gws(1:6,i)
        enddo
      enddo  ! i2b
      do i3b=1,n3b
        p3 = prm3ls(i3b)
        if( p3%jsp==p3%ksp ) then
          p3%gws(1:6,1) = (p3%gws(1:6,1)+p3%gws(1:6,2))/2
          p3%gws(1:6,2) = p3%gws(1:6,1)
          p3%gws(1:6,3) = (p3%gws(1:6,3)+p3%gws(1:6,4))/2
          p3%gws(1:6,4) = p3%gws(1:6,3)
        endif
        ip = ip +1
        gws(1:6,ip) = gws(1:6,ip) +p3%gws(1:6,1)
        ip = ip +1
        gws(1:6,ip) = gws(1:6,ip) +p3%gws(1:6,2)
        ip = ip +1
        gws(1:6,ip) = gws(1:6,ip) +p3%gws(1:6,3)
        ip = ip +1
        gws(1:6,ip) = gws(1:6,ip) +p3%gws(1:6,4)
        do i=1,p3%ncoef
          ip = ip +1
          gws(1:6,ip) = gws(1:6,ip) +p3%gws(1:6,4+i)
        enddo
      enddo  ! i3b
    endif

    return
  end subroutine gradw_uf3l
!=======================================================================
  subroutine b_spl(r,ts,nmax,nr,b,db,ddb)
!
!  Non-recursive implementation of B-spline function at R.
!  Assuming maximum order (d) = 3.
!
!  Args:
!    r --- position to be evaluated
!    ts --- a list of knots {t_n}
!    nmax --- length of ts
!
!  Return:
!    nr: index n in ts of position r
!    b: B array
!    db: dB array (derivative of B)
!    ddb: ddB array (2nd derivative of B) (optional)
!
    integer,intent(in):: nmax
    real(8),intent(in):: r,ts(nmax)
    integer,intent(out):: nr
    real(8),intent(out):: b(-3:0), db(-3:0)
    real(8),intent(out),optional:: ddb(-3:0)
!.....local variables
    real(8):: btmp(-3:+1,0:3)  ! Temporal B(n,d) array with n in (-3,+1)
    real(8):: dbtmp(-3:0), ddbtmp(-3:0)
    integer:: id, n, l
    real(8):: tn0,tn1,tn2,tn3,tn4,dt1,dt2,tmp1,tmp2
    real(8),parameter:: teps = 1d-8  ! epsilon for neighboring knot distance

!...index in the knot where ts(nr) <= r < ts(nr+1)
    nr = knot_index(r,nmax,ts)

!.....Compute B(n,d)
    btmp(:,:) = 0d0
    btmp(0,0) = 1d0
    do id = 1,3
      do l = -id,0
        n = nr +l
        tn0 = ts(n)
        tn1 = ts(n+id)
        dt1 = tn1 -tn0
        tmp1 = 0d0
        if( abs(dt1) > teps ) tmp1 = (r-tn0)/dt1 *btmp(l,id-1)
        tn2 = ts(n+1)
        tn3 = ts(n+id+1)
        dt2 = tn3 -tn2
        tmp2 = 0d0
        if( abs(dt2) > teps ) tmp2 = (tn3-r)/dt2 *btmp(l+1,id-1)
        btmp(l,id) = tmp1 + tmp2
      enddo
    enddo

!.....Compute dB(n) where d=3 is fixed
    dbtmp(:) = 0d0
    do l = -3,0
      n = nr +l
      tn0 = ts(n)
      tn1 = ts(n+1)
      tn3 = ts(n+3)
      tn4 = ts(n+1+3)
      tmp1 = 0d0
      if( abs(tn3-tn0) > teps ) tmp1 = btmp(l,2)/(tn3 -tn0)
      tmp2 = 0d0
      if( abs(tn4-tn1) > teps ) tmp2 = btmp(l+1,2)/(tn4 -tn1)
      dbtmp(l) = 3d0 *(tmp1 -tmp2)
    enddo

    b(-3:0) = btmp(-3:0,3)
    db(-3:0) = dbtmp(-3:0)

    if( .not.present(ddb) ) return
!.....Optional: 2nd derivative
    dbtmp(:) = 0d0
!.....Compute dB_{n,d-1} and dB_{n+1,d-1} first
    do l= -3,0
      n = nr +l
      tn0 = ts(n)
      tn1 = ts(n+1)
      tn2 = ts(n+2)
      tn3 = ts(n+3)
      if( abs(tn2-tn0) > teps ) tmp1 = btmp(l,1)/(tn2 -tn0)
      if( abs(tn3-tn1) > teps ) tmp2 = btmp(l+1,1)/(tn3 -tn1)
      dbtmp(l) = 2d0 *(tmp1 -tmp2)
    enddo
!.....Then, compute ddB_{n,d} using dBs
    ddbtmp(:) = 0d0
    do l= -3,0
      n = nr +l
      tn0 = ts(n)
      tn1 = ts(n+1)
      tn3 = ts(n+3)
      tn4 = ts(n+1+3)
      if( abs(tn3-tn0) > teps ) tmp1 = dbtmp(l)/(tn3 -tn0)
      if( abs(tn4-tn1) > teps ) tmp2 = dbtmp(l+1)/(tn4 -tn1)
      ddbtmp(l) = 3d0 *(tmp1 -tmp2)
    enddo
    ddb(-3:0) = ddbtmp(-3:0)

    return
  end subroutine b_spl
!=======================================================================
  recursive function b_spl_rec(n,d,r,ts,nmax) result(val)
!
!  TODO: check the efficiency of recursive func
!
!  Recursive implementation of B-spline function with N and D as indices
!  and R as an argument.
!  TS --- a list of {t_n}
!
    integer,intent(in):: n,d,nmax
    real(8),intent(in):: r,ts(nmax)
    real(8):: val
    real(8):: denom1,denom2

    if( d == 0 ) then
      if( r >= ts(n) .and. r < ts(n+1) ) then
        val = 1d0
        return
      else
        val = 0d0
        return
      endif
    else
      val = 0d0
      denom1 = ts(n+d)-ts(n)
      if( abs(denom1) > 1d-8 ) then
        val = val +(r-ts(n))/denom1 *b_spl_rec(n,d-1,r,ts,nmax)
      endif
      denom2 = ts(n+d+1)-ts(n+1)
      if( abs(denom2) > 1d-8 ) then
        val = val +(ts(n+d+1)-r)/denom2 *b_spl_rec(n+1,d-1,r,ts,nmax)
      endif
    endif
    return
  end function b_spl_rec
!=======================================================================
  function db_spl(n,d,r,ts,nmax)
    integer,intent(in):: n,d,nmax
    real(8),intent(in):: r,ts(nmax)
    real(8):: db_spl
    real(8):: denom1, denom2

    db_spl = 0d0
    denom1 = ts(n+d)-ts(n)
    if( denom1 > 1d-8 ) then
      db_spl = db_spl +d*b_spl_rec(n,d-1,r,ts,nmax) /denom1
    endif
    denom2 = ts(n+d+1)-ts(n+1)
    if( denom2 > 1d-8 ) then
      db_spl = db_spl -d*b_spl_rec(n+1,d-1,r,ts,nmax)/denom2
    endif
    return
  end function db_spl
!=======================================================================
  function knot_index(r, nknot, knots) result(n)
    real(8),intent(in):: r
    integer,intent(in):: nknot
    real(8),intent(in):: knots(nknot)
    integer:: n, i

!.....TODO: use faster algorithm
    n = 0
    do i=1,nknot
      if( r > knots(i) ) then
        n = i
      else
        return
      endif
    enddo
    return
  end function knot_index
!=======================================================================
  subroutine set_paramsdir_uf3(dname)
!
!  Accessor routine to set paramsdir.
!
    implicit none
    character(len=*),intent(in):: dname

    paramsdir = trim(dname)
    return
  end subroutine set_paramsdir_uf3
!=======================================================================
  subroutine set_params_uf3(ndimp,params)
!
!  Accesor routine to set uf3 parameters from outside.
!  It is supposed to be called from fitpot in a seriral process.
!  This will compare num of parameters given from outside
!  and num of coefficients already read in read_params_uf3().
!  Also it takes into account the symmetry of 3-body term
!  when species of j and k are identical.
!
    integer,intent(in):: ndimp
    real(8),intent(in):: params(ndimp)

    integer:: i1b,i2b,i3b,ncoef,ic,icfij,icfik,icfjk,inc,itmp
    type(prm2):: p2
    type(prm3):: p3

    if( .not. lprms_read_uf3 ) then
      print *,'ERROR(set_params_uf3): read_params_uf3 has not been called yet.'
      stop
    endif

!!$    if( .not. has_solo ) then
!!$      print *,'ERROR(set_params_uf3): .not.has_solo which should not happen.'
!!$    endif
!!$    if( .not. has_trios ) then
!!$      print *,'ERROR(set_params_uf3): .not.has_trio which should not happen.'
!!$    endif

!.....Count num of coeffs in force_uf3
    ncoef = 0
    do i1b=1,n1b
      ncoef = ncoef +1
    enddo
    do i2b=1,n2b
      p2 = prm2s(i2b)
      ncoef = ncoef +p2%ncoef
    enddo
    do i3b=1,n3b
      p3 = prm3s(i3b)
      ncoef = ncoef +p3%ncfij *p3%ncfik *p3%ncfjk
    enddo
    if( ncoef.ne.ndimp ) then
      print *,'ERROR(set_params_uf3): ncoef != ndimp'
      print *,'    This error may be caused when in.vars.fitpot and ' &
           //'in.params.uf3 are not consistent.'
      stop
    endif

!.....Replace coefficients with params given from outside.
    inc = 0
    do i1b=1,n1b
      inc = inc +1
      erg1s(i1b) = params(inc)
    enddo
    do i2b=1,n2b
      do ic=1,prm2s(i2b)%ncoef
        inc = inc +1
        prm2s(i2b)%coefs(ic) = params(inc)
      enddo
    enddo
    do i3b=1,n3b
      itmp = inc
      do icfij=1,prm3s(i3b)%ncfij
        do icfik=1,prm3s(i3b)%ncfik
          do icfjk=1,prm3s(i3b)%ncfjk
            inc = inc +1
            prm3s(i3b)%coefs(icfjk,icfik,icfij) = params(inc)
          enddo
        enddo
      enddo
    enddo

    return
  end subroutine set_params_uf3
!=======================================================================
  subroutine set_params_uf3l(ndimp,params)
!
!  Accesor routine to set uf3l parameters from outside.
!  It is supposed to be called from fitpot in a seriral process.
!  This will compare num of parameters given from outside
!  and num of coefficients already read in read_params_uf3().
!
    integer,intent(in):: ndimp
    real(8),intent(in):: params(ndimp)

    integer:: i1b,i2b,i3b,ncoef,ic,icfij,icfik,icfjk,inc,itmp, &
         nklead, nktrail
    real(8):: rcij, rcik
    type(prm2):: p2
    type(prm3l):: p3

    if( .not. lprms_read_uf3l ) then
      print *,'ERROR(set_params_uf3l): read_params_uf3l has not been called yet.'
      stop
    endif

!.....Count num of coeffs in force_uf3
    ncoef = 0
    do i1b=1,n1b
      ncoef = ncoef +1
    enddo
    do i2b=1,n2b
      p2 = prm2s(i2b)
      ncoef = ncoef +p2%ncoef
    enddo
    do i3b=1,n3b
      p3 = prm3ls(i3b)
      ncoef = ncoef +p3%ncoef +4  ! ncoef+4 (rcij,rcik,gmj,gmk)
    enddo
    if( ncoef.ne.ndimp ) then
      print *,'ERROR(set_params_uf3l): ncoef != ndimp'
      print *,'    This error may be caused when in.vars.fitpot and ' &
           //'in.params.uf3l are not consistent.'
      stop
    endif

!.....Replace coefficients with params given from outside.
    inc = 0
    do i1b=1,n1b
      inc = inc +1
      erg1s(i1b) = params(inc)
    enddo
    do i2b=1,n2b
      ncoef = prm2s(i2b)%ncoef
      nklead = prm2s(i2b)%nklead
      nktrail = prm2s(i2b)%nktrail
      do ic=1,ncoef
        inc = inc +1
        prm2s(i2b)%coefs(ic) = params(inc)
      enddo
    enddo
    rc3max = 0d0
    do i3b=1,n3b
      inc = inc +1
      rcij = params(inc)
      prm3ls(i3b)%rcij = rcij
      prm3ls(i3b)%rcij2= rcij**2
      inc = inc +1
      rcik = params(inc)
      prm3ls(i3b)%rcik = rcik
      prm3ls(i3b)%rcik2= rcik**2
      rc3max = max(rcij, rcik, rc3max)
      rc3max2 = rc3max**2
      inc = inc +1
      prm3ls(i3b)%gmj = params(inc)
      inc = inc +1
      prm3ls(i3b)%gmk = params(inc)
      ncoef = prm3ls(i3b)%ncoef
      do ic=1,ncoef
        inc = inc +1
        prm3ls(i3b)%coefs(ic) = params(inc)
      enddo
    enddo

    return
  end subroutine set_params_uf3l
!=======================================================================
  subroutine symmetrize_params_uf3(ndimp,params)
!
!  Accesor routine to set uf3 parameters from outside.
!  Make the 3-body parameters symmetric when species of j and k are identical.
!
    integer,intent(in):: ndimp
    real(8),intent(inout):: params(ndimp)

    integer:: i1b,i2b,i3b,ncoef,ic,icfij,icfik,icfjk,inc,itmp
    type(prm2):: p2
    type(prm3):: p3

    if( .not. lprms_read_uf3 ) then
      print *,'ERROR(set_params_uf3): read_params_uf3 has not been called yet.'
      stop
    endif

!.....Count num of coeffs in force_uf3
    ncoef = 0
    do i1b=1,n1b
      ncoef = ncoef +1
    enddo
    do i2b=1,n2b
      p2 = prm2s(i2b)
      ncoef = ncoef +p2%ncoef
    enddo
    do i3b=1,n3b
      p3 = prm3s(i3b)
      ncoef = ncoef +p3%ncfij *p3%ncfik *p3%ncfjk
    enddo
    if( ncoef.ne.ndimp ) then
      print *,'ERROR(set_params_uf3): ncoef != ndimp'
      print *,'    This error may be caused when in.vars.fitpot and ' &
           //'in.params.uf3 are not consistent.'
      stop
    endif

!.....Replace coefficients with params given from outside.
    inc = 0
    do i1b=1,n1b
      inc = inc +1
!.....pass
!!$      erg1s(i1b) = params(inc)
    enddo
    do i2b=1,n2b
      do ic=1,prm2s(i2b)%ncoef
        inc = inc +1
!.....pass
!!$        prm2s(i2b)%coefs(ic) = params(inc)
      enddo
    enddo
    do i3b=1,n3b
      p3 = prm3s(i3b)
      do icfij=1,p3%ncfij
        do icfik=1,p3%ncfik
          do icfjk=1,p3%ncfjk
            inc = inc +1
            if( p3%jsp == p3%ksp ) then
              params(inc) = (p3%coefs(icfjk,icfik,icfij) &
                   +p3%coefs(icfjk,icfij,icfik))/2
            else
              params(inc) = p3%coefs(icfjk,icfik,icfij)
            endif
          enddo
        enddo
      enddo
    enddo

    return
  end subroutine symmetrize_params_uf3
!=======================================================================
  subroutine calc_penalty_uf3(ndimp,params_in,pwgt2b,pwgt2bd, &
       pwgt2bs,pwgt3b,pwgt3bd,repul_radii,penalty)
!
!  Accesor routine to set uf3 parameters from outside.
!  It is supposed to be called from fitpot in a seriral process.
!
    integer,intent(in):: ndimp
    real(8),intent(in):: params_in(ndimp)
    real(8),intent(in):: pwgt2b,pwgt2bd,pwgt2bs,pwgt3b,pwgt3bd, &
         repul_radii(nspmax,nspmax)
    real(8),intent(out):: penalty

    integer:: i1b,i2b,i3b,ncoef,ic,icfij,icfik,icfjk,inc,k,nr2,isp,jsp
    type(prm2):: p2
    type(prm3):: p3
    real(8):: tmp,p2b,p2bd,p3b,p3bd,p2bs,dc1,dc2,rc
    logical:: ledge
    real(8),parameter:: tiny = 1d-8

    inc = 0
    do i1b=1,n1b
      inc = inc +1
      erg1s(i1b) = params_in(inc)
    enddo
!.....2-body
    p2b = 0d0
    p2bd= 0d0
    p2bs= 0d0
    do i2b=1,n2b
      p2 = prm2s(i2b)
!.....variables only for short-distance penalty
      isp = p2%isp
      jsp = p2%jsp
      rc = repul_radii(isp,jsp)
      nr2 = knot_index(rc,p2%nknot,p2%knots)

      do ic=1,p2%ncoef
        inc = inc +1
        prm2s(i2b)%coefs(ic) = params_in(inc)  ! replace coefs (p2 cannot be used here)
        if( .not.(abs(pwgt2bs)>1d-14 .and. ic<=nr2-2) ) p2b = p2b +prm2s(i2b)%coefs(ic)**2
      enddo
      if( abs(pwgt2bs)>1d-14 ) then
        do ic=1,nr2-3
          dc1 = prm2s(i2b)%coefs(ic)-prm2s(i2b)%coefs(ic+1)
          if( dc1 < 0d0 ) p2bs = p2bs +dc1**2
          if( ic-2 < 1 ) cycle
          dc2 = (prm2s(i2b)%coefs(ic-2) -3d0*prm2s(i2b)%coefs(ic-1) &
               +3d0*prm2s(i2b)%coefs(ic) -prm2s(i2b)%coefs(ic+1))
          if( dc2 < 0d0 ) p2bs = p2bs +dc2**2
        enddo
      endif
      if( abs(pwgt2bd).lt.1d-14 ) cycle
      do ic=2,p2%ncoef-1
        if( abs(pwgt2bs)>1d-14 .and. ic<=nr2-2 ) cycle
        tmp = 0d0
        do k=-1,1,2
          tmp = tmp +(prm2s(i2b)%coefs(ic+k) -prm2s(i2b)%coefs(ic))
        enddo
        p2bd = p2bd +tmp**2
      enddo

    enddo
    p2b = p2b*pwgt2b
    p2bd = p2bd*pwgt2bd
    p2bs = p2bs*pwgt2bs

!.....3-body
    p3b = 0d0
    p3bd= 0d0
    if( abs(pwgt3b).gt.1d-14 .or. abs(pwgt3bd).gt.1d-14 ) then
      do i3b=1,n3b
        do icfij=1,prm3s(i3b)%ncfij
          do icfik=1,prm3s(i3b)%ncfik
            do icfjk=1,prm3s(i3b)%ncfjk
              inc = inc +1
              prm3s(i3b)%coefs(icfjk,icfik,icfij) = params_in(inc)
              p3b = p3b +prm3s(i3b)%coefs(icfjk,icfik,icfij)**2
            enddo
          enddo
        enddo
        if( abs(pwgt3bd).lt.1d-14 ) cycle
        p3 = prm3s(i3b)
        do icfij=1,p3%ncfij
          do icfik=1,p3%ncfik
            do icfjk=1,p3%ncfjk
!.....ij
              if( icfij-1 > 0 .and. icfij+1 <= p3%ncfij ) p3bd = p3bd + &
                   (p3%coefs(icfjk,icfik,icfij-1) &
                   -2d0*p3%coefs(icfjk,icfik,icfij) &
                   +p3%coefs(icfjk,icfik,icfij+1))**2
!.....ik
              if( icfik-1 > 0 .and. icfik+1 <= p3%ncfik ) p3bd = p3bd + &
                   (p3%coefs(icfjk,icfik-1,icfij) &
                   -2d0*p3%coefs(icfjk,icfik,icfij) &
                   +p3%coefs(icfjk,icfik+1,icfij))**2
!.....jk
              if( icfjk-1 > 0 .and. icfjk+1 <= p3%ncfjk ) p3bd = p3bd + &
                   (p3%coefs(icfjk-1,icfik,icfij) &
                   -2d0*p3%coefs(icfjk,icfik,icfij) &
                   +p3%coefs(icfjk+1,icfik,icfij))**2
            enddo
          enddo
        enddo
      enddo
    endif
    p3b = p3b *pwgt3b
    p3bd= p3bd*pwgt3bd
    penalty = p2b +p2bd +p2bs +p3b +p3bd

    return
  end subroutine calc_penalty_uf3
!=======================================================================
  subroutine calc_penalty_grad_uf3(ndimp,params_in,pwgt2b,pwgt2bd, &
       pwgt2bs,pwgt3b,pwgt3bd,repul_radii,grad)
!
!  Accesor routine to get gradient of uf3 ridge penalty.
!  It is supposed to be called from fitpot in a seriral process.
!
    integer,intent(in):: ndimp
    real(8),intent(in):: params_in(ndimp)
    real(8),intent(in):: pwgt2b,pwgt2bd,pwgt2bs,pwgt3b,pwgt3bd, &
         repul_radii(nspmax,nspmax)
    real(8),intent(out):: grad(ndimp)

    integer:: i1b,i2b,i3b,ncoef,ic,icfij,icfik,icfjk,inc,k,ip, &
         nr2,isp,jsp
    type(prm2):: p2
    type(prm3):: p3
    real(8):: tmp,dc1,dc2,rc
    logical:: ledge
    real(8),save,allocatable:: gp2b(:),gp2bd(:),gp2bs(:),gp3b(:),gp3bd(:)
    integer,save:: nc2max, ncfijmax, ncfikmax, ncfjkmax
    integer,save,allocatable:: ic2ip(:),ic3ip(:,:,:)

    if( .not.allocated(gp2b) ) then
      allocate(gp2b(ndimp),gp2bd(ndimp),gp2bs(ndimp), &
           gp3b(ndimp),gp3bd(ndimp))
      nc2max = 0
      do i2b=1,n2b
        nc2max = max(nc2max, prm2s(i2b)%ncoef)
      enddo
      allocate(ic2ip(nc2max))
      ncfijmax = 0
      ncfikmax = 0
      ncfjkmax = 0
      do i3b=1,n3b
        ncfijmax = max(ncfijmax, prm3s(i3b)%ncfij)
        ncfikmax = max(ncfikmax, prm3s(i3b)%ncfik)
        ncfjkmax = max(ncfjkmax, prm3s(i3b)%ncfjk)
      enddo
      allocate(ic3ip(ncfjkmax,ncfikmax,ncfijmax))
    endif

    inc = 0
    do i1b=1,n1b
      inc = inc +1
      erg1s(i1b) = params_in(inc)
    enddo
!.....2-body
    gp2b(:) = 0d0
    gp2bd(:)= 0d0
    gp2bs(:)= 0d0
    do i2b=1,n2b
      p2 = prm2s(i2b)
!.....variables only for short-distance penalty
      isp = p2%isp
      jsp = p2%jsp
      rc = repul_radii(isp,jsp)
      nr2 = knot_index(rc,p2%nknot,p2%knots)

      ic2ip(:) = 0
      do ic=1,prm2s(i2b)%ncoef
        inc = inc +1
        prm2s(i2b)%coefs(ic) = params_in(inc)
        ic2ip(ic) = inc
        if( .not. (abs(pwgt2bs)>1d-14 .and. ic<=nr2-2) ) gp2b(inc) = &
             gp2b(inc) +2d0*pwgt2b*prm2s(i2b)%coefs(ic)
      enddo
      if( abs(pwgt2bs)>1d-14 ) then
        do ic=1,nr2-3
          dc1 = prm2s(i2b)%coefs(ic)-prm2s(i2b)%coefs(ic+1)
          ip = ic2ip(ic)
          if( dc1 < 0d0 ) then
            gp2bs(ip  )= gp2bs(ip  ) +2d0*dc1 *pwgt2bs
            gp2bs(ip+1)= gp2bs(ip+1) -2d0*dc1 *pwgt2bs
          endif
          if( ic-2 < 1 ) cycle
          dc2 = (prm2s(i2b)%coefs(ic-2) -3d0*prm2s(i2b)%coefs(ic-1) &
               +3d0*prm2s(i2b)%coefs(ic) -prm2s(i2b)%coefs(ic+1))
          if( dc2 < 0d0 ) then
            gp2bs(ip-2)= gp2bs(ip-2) +2d0*dc2 *pwgt2bs
            gp2bs(ip-1)= gp2bs(ip-1) -6d0*dc2 *pwgt2bs
            gp2bs(ip  )= gp2bs(ip  ) +6d0*dc2 *pwgt2bs
            gp2bs(ip+1)= gp2bs(ip+1) -2d0*dc2 *pwgt2bs
          endif
        enddo
      endif
      if( abs(pwgt2bd).lt.1d-14 ) cycle
      do ic=1,p2%ncoef
        if( abs(pwgt2bs)>1d-14 .and. ic<=nr2-2 ) cycle
        ip = ic2ip(ic)
        if( ic-2 > 0 ) gp2bd(ip) = gp2bd(ip) +2d0*pwgt2bd &
             *(p2%coefs(ic-2) &
             -2d0*p2%coefs(ic-1) &
             +p2%coefs(ic))
        if( ic-1 > 0 .and. ic+1 <= p2%ncoef ) &
             gp2bd(ip) = gp2bd(ip) +2d0*pwgt2bd*(-2d0) &
             *(p2%coefs(ic-1) &
             -2d0*p2%coefs(ic) &
             +p2%coefs(ic+1))
        if( ic+2 <= p2%ncoef ) gp2bd(ip) = gp2bd(ip) +2d0*pwgt2bd &
             *(p2%coefs(ic) &
             -2d0*p2%coefs(ic+1) &
             +p2%coefs(ic+2))
        enddo
      enddo

!.....3-body
    gp3b(:) = 0d0
    gp3bd(:)= 0d0
    if( abs(pwgt3b).gt.1d-14 .or. abs(pwgt3bd).gt.1d-14 ) then
      do i3b=1,n3b
        ic3ip(:,:,:) = 0
        do icfij=1,prm3s(i3b)%ncfij
          do icfik=1,prm3s(i3b)%ncfik
            do icfjk=1,prm3s(i3b)%ncfjk
              inc = inc +1
              prm3s(i3b)%coefs(icfjk,icfik,icfij) = params_in(inc)
              gp3b(inc) = gp3b(inc) +2d0*pwgt3b*prm3s(i3b)%coefs(icfjk,icfik,icfij)
              ic3ip(icfjk,icfik,icfij) = inc
            enddo
          enddo
        enddo
        if( abs(pwgt3bd).lt.1d-14 ) cycle
        p3 = prm3s(i3b)
        do icfij=1,p3%ncfij
          do icfik=1,p3%ncfik
            do icfjk=1,p3%ncfjk
              ip = ic3ip(icfjk,icfik,icfij)
!.....ij
              if( icfij-2 > 0 ) gp3bd(ip) = gp3bd(ip) +2d0*pwgt3bd &
                   *(   p3%coefs(icfjk,icfik,icfij-2) &
                   -2d0*p3%coefs(icfjk,icfik,icfij-1) &
                   +    p3%coefs(icfjk,icfik,icfij))
              if( icfij-1 > 0 .and. icfij+1 <= p3%ncfij ) &
                   gp3bd(ip) = gp3bd(ip) +2d0*pwgt3bd*(-2d0) &
                   *(   p3%coefs(icfjk,icfik,icfij-1) &
                   -2d0*p3%coefs(icfjk,icfik,icfij) &
                   +    p3%coefs(icfjk,icfik,icfij+1))
              if( icfij+2 <= p3%ncfij ) gp3bd(ip) = gp3bd(ip) +2d0*pwgt3bd &
                   *(   p3%coefs(icfjk,icfik,icfij) &
                   -2d0*p3%coefs(icfjk,icfik,icfij+1) &
                   +    p3%coefs(icfjk,icfik,icfij+2))
!.....ik
              if( icfik-2 > 0 ) gp3bd(ip) = gp3bd(ip) +2d0*pwgt3bd &
                   *(   p3%coefs(icfjk,icfik-2,icfij) &
                   -2d0*p3%coefs(icfjk,icfik-1,icfij) &
                   +    p3%coefs(icfjk,icfik,  icfij))
              if( icfik-1 > 0 .and. icfik+1 <= p3%ncfik ) &
                   gp3bd(ip) = gp3bd(ip) +2d0*pwgt3bd*(-2d0) &
                   *(   p3%coefs(icfjk,icfik-1,icfij) &
                   -2d0*p3%coefs(icfjk,icfik,  icfij) &
                   +    p3%coefs(icfjk,icfik+1,icfij))
              if( icfik+2 <= p3%ncfik ) gp3bd(ip) = gp3bd(ip) +2d0*pwgt3bd &
                   *(   p3%coefs(icfjk,icfik,  icfij) &
                   -2d0*p3%coefs(icfjk,icfik+1,icfij) &
                   +    p3%coefs(icfjk,icfik+2,icfij))
!.....jk
              if( icfjk-2 > 0 ) gp3bd(ip) = gp3bd(ip) +2d0*pwgt3bd &
                   *(   p3%coefs(icfjk-2,icfik,icfij) &
                   -2d0*p3%coefs(icfjk-1,icfik,icfij) &
                   +    p3%coefs(icfjk,  icfik,icfij))
              if( icfjk-1 > 0 .and. icfjk+1 <= p3%ncfjk ) &
                   gp3bd(ip) = gp3bd(ip) +2d0*pwgt3bd*(-2d0) &
                   *(   p3%coefs(icfjk-1,icfik,icfij) &
                   -2d0*p3%coefs(icfjk,  icfik,icfij) &
                   +    p3%coefs(icfjk+1,icfik,icfij))
              if( icfjk+2 <= p3%ncfjk ) gp3bd(ip) = gp3bd(ip) +2d0*pwgt3bd &
                   *(   p3%coefs(icfjk,  icfik,icfij) &
                   -2d0*p3%coefs(icfjk+1,icfik,icfij) &
                   +    p3%coefs(icfjk+2,icfik,icfij))
            enddo  ! icfjk
          enddo  ! icfik
        enddo  ! icfij
      enddo  ! i3b
    endif
    grad(:) = gp2b(:) +gp2bd(:) +gp2bs(:) +gp3b(:) +gp3bd(:)

    return
  end subroutine calc_penalty_grad_uf3
!=======================================================================
  subroutine calc_penalty_uf3l(ndimp,params_in,pwgt2b,pwgt2bd, &
       pwgt2bs,pwgt3b,pwgt3bd,repul_radii,penalty)
!
!  Accesor routine to set uf3l parameters from outside.
!  It is supposed to be called from fitpot in a seriral process.
!
    integer,intent(in):: ndimp
    real(8),intent(in):: params_in(ndimp)
    real(8),intent(in):: pwgt2b,pwgt2bd,pwgt2bs,pwgt3b,pwgt3bd, &
         repul_radii(nspmax,nspmax)
    real(8),intent(out):: penalty

    integer:: i1b,i2b,i3b,ncoef,ic,icfij,icfik,icfjk,inc,k,nr2,isp,jsp, &
         nklead, nktrail
    type(prm2):: p2
    type(prm3l):: p3
    real(8):: tmp,p2b,p2bd,p3b,p3bd,p2bs,dc1,dc2,rc
    logical:: ledge
    real(8),parameter:: tiny = 1d-8

    inc = 0
    do i1b=1,n1b
      inc = inc +1
      erg1s(i1b) = params_in(inc)
    enddo
!.....2-body
    p2b = 0d0
    p2bd= 0d0
    p2bs= 0d0
    do i2b=1,n2b
      p2 = prm2s(i2b)
!.....variables only for short-distance penalty
      isp = p2%isp
      jsp = p2%jsp
      nklead = p2%nklead
      nktrail= p2%nktrail
      rc = repul_radii(isp,jsp)
      nr2 = knot_index(rc,p2%nknot,p2%knots)

      do ic=1,p2%ncoef
        inc = inc +1
        prm2s(i2b)%coefs(ic) = params_in(inc)  ! replace coefs (p2 cannot be used here)
        if( .not.(abs(pwgt2bs)>1d-14 .and. ic<=nr2-2) ) p2b = p2b +prm2s(i2b)%coefs(ic)**2
      enddo
      if( abs(pwgt2bs)>1d-14 ) then
        do ic=1,nr2-3
          dc1 = prm2s(i2b)%coefs(ic)-prm2s(i2b)%coefs(ic+1)
          if( dc1 < 0d0 ) p2bs = p2bs +dc1**2
          if( ic-2 < 1 ) cycle
          dc2 = (prm2s(i2b)%coefs(ic-2) -3d0*prm2s(i2b)%coefs(ic-1) &
               +3d0*prm2s(i2b)%coefs(ic) -prm2s(i2b)%coefs(ic+1))
          if( dc2 < 0d0 ) p2bs = p2bs +dc2**2
        enddo
      endif
      if( abs(pwgt2bd).lt.1d-14 ) cycle
      do ic=2,p2%ncoef-1
        if( abs(pwgt2bs)>1d-14 .and. ic<=nr2-2 ) cycle
        tmp = 0d0
        do k=-1,1,2
          tmp = tmp +(prm2s(i2b)%coefs(ic+k) -prm2s(i2b)%coefs(ic))
        enddo
        p2bd = p2bd +tmp**2
      enddo

    enddo
    p2b = p2b*pwgt2b
    p2bd = p2bd*pwgt2bd
    p2bs = p2bs*pwgt2bs

!.....3-body
    p3b = 0d0
    p3bd= 0d0
    if( abs(pwgt3b).gt.1d-14 .or. abs(pwgt3bd).gt.1d-14 ) then
      do i3b=1,n3b
        p3 = prm3ls(i3b)
        inc = inc +1
        prm3ls(i3b)%gmj = params_in(inc)
        p3b = p3b + prm3ls(i3b)%gmj**2
        inc = inc +1
        prm3ls(i3b)%gmk = params_in(inc)
        p3b = p3b + prm3ls(i3b)%gmk**2
        do ic=1,p3%ncoef
          inc = inc +1
          prm3ls(i3b)%coefs(ic) = params_in(inc)
          p3b = p3b +prm3ls(i3b)%coefs(ic)**2
        enddo

!.....penalty to derivative
        if( abs(pwgt3bd).lt.1d-14 ) cycle
        do ic=1,p3%ncoef
          if( ic-1 > 0 .and. ic+1 <= p3%ncoef ) p3bd = p3bd + &
               (prm3ls(i3b)%coefs(ic-1) -2d0*prm3ls(i3b)%coefs(ic) &
               +prm3ls(i3b)%coefs(ic+1))**2
        enddo
      enddo
    endif
    p3b = p3b *pwgt3b
    p3bd= p3bd*pwgt3bd
    penalty = p2b +p2bd +p2bs +p3b +p3bd

    return
  end subroutine calc_penalty_uf3l
!=======================================================================
  subroutine calc_penalty_grad_uf3l(ndimp,params_in,pwgt2b,pwgt2bd, &
       pwgt2bs,pwgt3b,pwgt3bd,repul_radii,grad)
!
!  Accesor routine to get gradient of uf3 ridge penalty.
!  It is supposed to be called from fitpot in a seriral process.
!
    integer,intent(in):: ndimp
    real(8),intent(in):: params_in(ndimp)
    real(8),intent(in):: pwgt2b,pwgt2bd,pwgt2bs,pwgt3b,pwgt3bd, &
         repul_radii(nspmax,nspmax)
    real(8),intent(out):: grad(ndimp)

    integer:: i1b,i2b,i3b,ncoef,ic,icfij,icfik,icfjk,inc,k,ip, &
         nr2,isp,jsp, nklead, nktrail
    type(prm2):: p2
    type(prm3l):: p3
    real(8):: tmp,dc1,dc2,rc
    logical:: ledge
    real(8),save,allocatable:: gp2b(:),gp2bd(:),gp2bs(:),gp3b(:),gp3bd(:)
    integer,save:: nc2max, nc3max
    integer,save,allocatable:: ic2ip(:),ic3ip(:)

    if( .not.allocated(gp2b) ) then
      allocate(gp2b(ndimp),gp2bd(ndimp),gp2bs(ndimp), &
           gp3b(ndimp),gp3bd(ndimp))
      nc2max = 0
      do i2b=1,n2b
        nc2max = max(nc2max, prm2s(i2b)%ncoef)
      enddo
      allocate(ic2ip(nc2max))
      nc3max = 0
      do i3b=1,n3b
        nc3max = max(nc3max, prm3ls(i3b)%ncoef)
      enddo
      allocate(ic3ip(nc3max))
    endif

    inc = 0
    do i1b=1,n1b
      inc = inc +1
      erg1s(i1b) = params_in(inc)
    enddo
!.....2-body
    gp2b(:) = 0d0
    gp2bd(:)= 0d0
    gp2bs(:)= 0d0
    do i2b=1,n2b
      p2 = prm2s(i2b)
!.....variables only for short-distance penalty
      isp = p2%isp
      jsp = p2%jsp
      rc = repul_radii(isp,jsp)
      nr2 = knot_index(rc,p2%nknot,p2%knots)

      ic2ip(:) = 0
      do ic=1,p2%ncoef
        inc = inc +1
        prm2s(i2b)%coefs(ic) = params_in(inc)
        ic2ip(ic) = inc
        if( .not. (abs(pwgt2bs)>1d-14 .and. ic<=nr2-2) ) gp2b(inc) = &
             gp2b(inc) +2d0*pwgt2b*prm2s(i2b)%coefs(ic)
      enddo
      if( abs(pwgt2bs)>1d-14 ) then
        do ic=1,nr2-3
          dc1 = prm2s(i2b)%coefs(ic)-prm2s(i2b)%coefs(ic+1)
          ip = ic2ip(ic)
          if( dc1 < 0d0 ) then
            gp2bs(ip  )= gp2bs(ip  ) +2d0*dc1 *pwgt2bs
            gp2bs(ip+1)= gp2bs(ip+1) -2d0*dc1 *pwgt2bs
          endif
          if( ic-2 < 1 ) cycle
          dc2 = (prm2s(i2b)%coefs(ic-2) -3d0*prm2s(i2b)%coefs(ic-1) &
               +3d0*prm2s(i2b)%coefs(ic) -prm2s(i2b)%coefs(ic+1))
          if( dc2 < 0d0 ) then
            gp2bs(ip-2)= gp2bs(ip-2) +2d0*dc2 *pwgt2bs
            gp2bs(ip-1)= gp2bs(ip-1) -6d0*dc2 *pwgt2bs
            gp2bs(ip  )= gp2bs(ip  ) +6d0*dc2 *pwgt2bs
            gp2bs(ip+1)= gp2bs(ip+1) -2d0*dc2 *pwgt2bs
          endif
        enddo
      endif
      if( abs(pwgt2bd).lt.1d-14 ) cycle
      do ic=1,p2%ncoef
        if( abs(pwgt2bs)>1d-14 .and. ic<=nr2-2 ) cycle
        ip = ic2ip(ic)
        if( ic-2 > 0 ) gp2bd(ip) = gp2bd(ip) +2d0*pwgt2bd &
             *(p2%coefs(ic-2) &
             -2d0*p2%coefs(ic-1) &
             +p2%coefs(ic))
        if( ic-1 > 0 .and. ic+1 <= p2%ncoef ) &
             gp2bd(ip) = gp2bd(ip) +2d0*pwgt2bd*(-2d0) &
             *(p2%coefs(ic-1) &
             -2d0*p2%coefs(ic) &
             +p2%coefs(ic+1))
        if( ic+2 <= p2%ncoef ) gp2bd(ip) = gp2bd(ip) +2d0*pwgt2bd &
             *(p2%coefs(ic) &
             -2d0*p2%coefs(ic+1) &
             +p2%coefs(ic+2))
      enddo
    enddo

!.....3-body
    gp3b(:) = 0d0
    gp3bd(:)= 0d0
    if( abs(pwgt3b).gt.1d-14 .or. abs(pwgt3bd).gt.1d-14 ) then
      do i3b=1,n3b
        ncoef = prm3ls(i3b)%ncoef
        ic3ip(:) = 0
        inc = inc +1
        prm3ls(i3b)%gmj = params_in(inc)
        gp3b(inc) = gp3b(inc) +2d0*pwgt3b*prm3ls(i3b)%rcij
        inc = inc +1
        prm3ls(i3b)%gmk = params_in(inc)
        gp3b(inc) = gp3b(inc) +2d0*pwgt3b*prm3ls(i3b)%rcik
        inc = inc +1
        prm3ls(i3b)%gmj = params_in(inc)
        gp3b(inc) = gp3b(inc) +2d0*pwgt3b*prm3ls(i3b)%gmj
        inc = inc +1
        prm3ls(i3b)%gmk = params_in(inc)
        gp3b(inc) = gp3b(inc) +2d0*pwgt3b*prm3ls(i3b)%gmk
        do ic=1, ncoef
          inc = inc +1
          prm3ls(i3b)%coefs(ic) = params_in(inc)
          gp3b(inc) = gp3b(inc) +2d0*pwgt3b*prm3ls(i3b)%coefs(ic)
          ic3ip(ic) = inc
        enddo

        if( abs(pwgt3bd).lt.1d-14 ) cycle
        p3 = prm3ls(i3b)
        do ic=1,ncoef
          ip = ic3ip(ic)
!.....ij
          if( ic-2 > 0 ) gp3bd(ip) = gp3bd(ip) +2d0*pwgt3bd &
               *(   p3%coefs(ic-2) &
               -2d0*p3%coefs(ic-1) &
               +    p3%coefs(ic))
          if( ic-1 > 0 .and. ic+1 <= p3%ncoef ) &
               gp3bd(ip) = gp3bd(ip) +2d0*pwgt3bd*(-2d0) &
               *(   p3%coefs(ic-1) &
               -2d0*p3%coefs(ic) &
               +    p3%coefs(ic+1))
          if( ic+2 <= p3%ncoef ) gp3bd(ip) = gp3bd(ip) +2d0*pwgt3bd &
               *(   p3%coefs(ic) &
               -2d0*p3%coefs(ic+1) &
               +    p3%coefs(ic+2))
        enddo  ! ic
      enddo  ! i3b
    endif
    grad(:) = gp2b(:) +gp2bd(:) +gp2bs(:) +gp3b(:) +gp3bd(:)

    return
  end subroutine calc_penalty_grad_uf3l
!=======================================================================
  subroutine calc_short_lossfunc(npnts,radii,drepul,floss)
!
!  Compute loss function for short-distance repulsion correction.
!
    integer,intent(in):: npnts
    real(8),intent(in):: radii(nspmax,nspmax),drepul(npnts,nspmax,nspmax)
    real(8),intent(out):: floss

    integer:: i2b,isp,jsp,ir,lij,n,nr2
    real(8):: fli,ri,tmp,c2t,bij(-3:0),dbij(-3:0)
    type(prm2):: p2

    floss = 0d0
    do i2b=1,n2b
      p2 = prm2s(i2b)
      isp = csp2isp(p2%csi)
      jsp = csp2isp(p2%csj)
      fli = 0d0
      do ir=1,npnts
!.....r-point to be evaluated as mid-point of the section
        ri = radii(isp,jsp)/npnts *(dble(ir)-0.5d0)
        call b_spl(ri,p2%knots,p2%nknot,nr2,bij,dbij)
        tmp = 0d0
        do lij = -3,0
          n = nr2 +lij
          if( n < 1 .or. n > p2%nknot-4 ) cycle
          c2t = p2%coefs(n)
          tmp = tmp +c2t *dbij(lij)
        enddo
        fli = fli +(tmp -drepul(ir,isp,jsp))**2
      enddo
      floss = floss +fli/npnts
    enddo
    floss = floss/n2b *0.5d0
    return
  end subroutine calc_short_lossfunc
!=======================================================================
  subroutine calc_short_lossgrad(npnts,radii,drepul,ndimp,gloss)
!
!  Compute loss function gradient for short-distance repulsion correction.
!
    integer,intent(in):: npnts,ndimp
    real(8),intent(in):: radii(nspmax,nspmax),drepul(npnts,nspmax,nspmax)
    real(8),intent(out):: gloss(ndimp)

    integer:: i2b,isp,jsp,ir,lij,n,i1b,inc,ic,ip,nr2
    real(8):: fli,ri,tmp,c2t,bij(-3:0),dbij(-3:0)
    type(prm2):: p2
    integer,save:: nc2max
    integer,save,allocatable:: ic2ip(:)

    if( .not.allocated(ic2ip) ) then
      nc2max = 0
      do i2b=1,n2b
        nc2max = max(nc2max, prm2s(i2b)%ncoef)
      enddo
      allocate(ic2ip(nc2max))
    endif

    inc = 0
    do i1b=1,n1b
      inc = inc +1
    enddo

    gloss(:) = 0d0
    do i2b=1,n2b
      p2 = prm2s(i2b)
      isp = csp2isp(p2%csi)
      jsp = csp2isp(p2%csj)
      ic2ip(:) = 0
      do ic=1,p2%ncoef
        inc = inc +1
        ic2ip(ic) = inc
      enddo

      do ir=1,npnts
!.....r-point to be evaluated as mid-point of the section
        ri = radii(isp,jsp)/npnts *(dble(ir)-0.5d0)
        call b_spl(ri,p2%knots,p2%nknot,nr2,bij,dbij)
        tmp = 0d0
        do lij = -3,0
          n = nr2 +lij
          if( n < 1 .or. n > p2%nknot-4 ) cycle
          c2t = p2%coefs(n)
          tmp = tmp +c2t *dbij(lij)
        enddo  ! lij
        fli = tmp -drepul(ir,isp,jsp)
        do lij = -3,0
          n = nr2 +lij
          if( n < 1 .or. n > p2%nknot-4 ) cycle
          ip = ic2ip(n)
          gloss(ip) = gloss(ip) +fli*dbij(lij)/npnts/n2b
        enddo  ! lij
      enddo  ! ir
    enddo  ! i2b
    return

  end subroutine calc_short_lossgrad
!=======================================================================
  subroutine uf3_short_correction(ndimp,params,nsp,radii,ldcover)
!
!  Modify some coefficients to correct short-range repulsive potential.
!
    integer,intent(in):: ndimp,nsp
    real(8),intent(inout):: params(ndimp),radii(nspmax,nspmax)
    logical,intent(in):: ldcover(ndimp)

    integer:: inc,i1b,i2b,isp,jsp,ic,nr2,ip,min_ic_data,l,n
    real(8):: ri,tgt,bij(-3:0),dbij(-3:0),ddbij(-3:0),dr,rc, &
         val,dval,ddval,c
    type(prm2):: p2
    integer,save:: nc2max
    integer,allocatable,save:: ic2ip(:)
    integer,parameter:: nbuf = 3

    if( .not.allocated(ic2ip) ) then
      nc2max = 0
      do i2b=1,n2b
        nc2max = max(nc2max, prm2s(i2b)%ncoef)
      enddo
      allocate(ic2ip(nc2max))
    endif

!.....Replace coefficients with params given from outside
!     so that easily capture the relationshiop between neighboring parameters.
    inc = 0
    do i1b=1,n1b
      inc = inc +1
      erg1s(i1b) = params(inc)
    enddo
    do i2b=1,n2b
      ic2ip(:) = 0
      do ic=1,prm2s(i2b)%ncoef
        inc = inc +1
        prm2s(i2b)%coefs(ic) = params(inc)
        ic2ip(ic) = inc
      enddo
      min_ic_data = 0
      do ic=prm2s(i2b)%ncoef,1,-1
        ip = ic2ip(ic)
        if( ldcover(ip) ) min_ic_data = ic
      enddo
!!$      print '(a,2i5,2f10.3)','i2b,min_ic_data,r,rr=',i2b,min_ic_data,&
!!$           prm2s(i2b)%knots(min_ic_data),prm2s(i2b)%knots(min_ic_data+4)

!.....Short-range correction
      p2 = prm2s(i2b)
      isp = csp2isp(p2%csi)
      jsp = csp2isp(p2%csj)
      if( isp < 1 .or. isp > nsp .or. jsp < 1 .or. jsp > nsp ) cycle
      ri = radii(isp,jsp)
      if( ri < p2%knots(min_ic_data) ) ri = (p2%knots(min_ic_data)+p2%knots(min_ic_data+1))/2
      nr2 = knot_index(ri, p2%nknot, p2%knots)
      call b_spl(ri,p2%knots,p2%nknot,nr2,bij,dbij,ddbij)
      val  = 0d0
      dval = 0d0
      ddval= 0d0
      do l= -3,0
        n = nr2 +l
        c = p2%coefs(n)
        val  = val   +c*bij(l)
        dval = dval  +c*dbij(l)
        ddval= ddval +c*ddbij(l)
      enddo
!.....There must be some condition to dval and/or ddval
!     to achieve decent repulsion to the pair...
      dval = min(dval, 0d0)
      ddval= max(ddval, 1d0)
!!$      print '(a,3i3,f6.3,i3,3es12.3)','i2b,isp,jsp,ri,nr2,v,dv,ddv=', &
!!$           i2b,isp,jsp,ri,nr2,val,dval,ddval
!.....Correct coefs using above information
      do ic= min_ic_data-1,1,-1
        rc = p2%knots(ic+2)  ! peak position of B_{ic,3} as t_{ic+2}
        dr = rc -ri
        tgt = val +dval*dr +0.5d0*ddval*dr**2
        ip = ic2ip(ic)
        params(ip) = max(tgt,0d0)
!!$        print '(a,i3,4es12.3)','  ic,ri,rc,dr,tgt=',ic,ri,rc,dr,tgt
      enddo

!!$      ri = radii(isp,jsp)
!!$      nr2 = knot_index(ri, p2%nknot, p2%knots)
!!$      do ic=nr2-nbuf,1,-1
!!$        if( ic < 1 ) exit
!!$      do ic=min_ic_data-1,1,-1
!!$        ip = ic2ip(ic)
!!$        tgt = 3d0*(p2%coefs(ic+1)-p2%coefs(ic+2)) +p2%coefs(ic+3)
!!$        params(ip) = max(params(ip+1),tgt,0d0)
!!$        print '(a,3i5,2es12.3)', '   i2b,ic,ip,tgt,prm=',&
!!$             i2b,ic,ip,tgt,params(ip)
!!$      enddo
    enddo
    return
  end subroutine uf3_short_correction
!=======================================================================
  subroutine print_1b()
    integer:: i
    print '(/,a)','   UF3 parameters of 1B:'
    print '(a,9(2x,f0.3))','   ',(erg1s(i),i=1,nspmax)
  end subroutine print_1b
!=======================================================================
  subroutine print_2b(ps)
    type(prm2),intent(in):: ps
    integer:: i
    character:: c*2
    print '(/,a)','   UF3 parameters of 2B for '&
         //trim(ps%csi)//'-'//trim(ps%csj)
    print '(5(a,1x))','     cb,csi,csj,cknot = ',ps%cb,ps%csi,ps%csj,ps%cknot
    print '(a,2i3)','     nklead,nktrail = ',ps%nklead,ps%nktrail
    print '(a,2i3)','     nknot,ncoef = ',ps%nknot, ps%ncoef
    print '(a,f6.3)','     rc = ',ps%rc
    write(c,'(i2)') size(ps%knots)
    print '(a,'//c//'(1x,f7.3))','     knots =',(ps%knots(i),i=1,size(ps%knots))
    write(c,'(i2)') size(ps%coefs)
    print '(a,'//trim(c)//'(1x,es11.2))','     coefs =',(ps%coefs(i),i=1,size(ps%coefs))
  end subroutine print_2b
!=======================================================================
  subroutine print_3b(ps)
    type(prm3),intent(in):: ps
    integer:: i
    character:: c*2
    print '(/,a)','   UF3 parameters of 3B for '&
         //trim(ps%csi)//'-'//trim(ps%csj)//'-'//trim(ps%csk)
    print '(6(a,1x))','     cb,csi,csj,csk,cknot = ',ps%cb,ps%csi,ps%csj,ps%csk,ps%cknot
    print '(a,2i3)','     nklead,nktrail = ',ps%nklead,ps%nktrail
    print '(a,6i3)','     nknij,nknik,nknjk,ncfij,ncfik,ncfjk = ', &
         ps%nknij, ps%nknik, ps%nknjk, ps%ncfij, ps%ncfik, ps%ncfjk
    print '(a,3f6.3)','     rcij,rcik,rcjk = ',ps%rcij,ps%rcik,ps%rcjk
    write(c,'(i2)') ps%nknij
    print '(a,'//c//'(1x,f7.3))','     knij =',(ps%knij(i),i=1,ps%nknij)
    write(c,'(i2)') ps%nknik
    print '(a,'//c//'(1x,f7.3))','     knik =',(ps%knik(i),i=1,ps%nknik)
    write(c,'(i2)') ps%nknjk
    print '(a,'//c//'(1x,f7.3))','     knjk =',(ps%knjk(i),i=1,ps%nknjk)
!!$    write(c,'(i2)') size(ps%coefs)
!!$    print '(a,'//trim(c)//'es11.2)','     coefs =',(ps%coefs(i),i=1,size(ps%coefs))
  end subroutine print_3b
!=======================================================================
  subroutine print_3bl(ps)
    type(prm3l),intent(in):: ps
    integer:: i
    character:: c*2
    print '(/,a)','   UF3L parameters of 3B for '&
         //trim(ps%csi)//'-'//trim(ps%csj)//'-'//trim(ps%csk)
    print '(6(a,1x))','     cb,csi,csj,csk,cknot = ',ps%cb,ps%csi,ps%csj,ps%csk,ps%cknot
    print '(a,2i3)','     nklead,nktrail = ',ps%nklead,ps%nktrail
    print '(a,6i3)','     nknot,ncoef = ', ps%nknot, ps%ncoef
    print '(a,3f6.3)','     rcij,rcik = ',ps%rcij, ps%rcik
    write(c,'(i2)') size(ps%knots)
    print '(a,'//c//'(1x,f7.3))','     knots =',(ps%knots(i),i=1,size(ps%knots))
    write(c,'(i2)') size(ps%coefs)
    print '(a,'//trim(c)//'es11.2)','     coefs =',(ps%coefs(i),i=1,size(ps%coefs))
  end subroutine print_3bl
!=======================================================================
  subroutine print_3bd(ps)
    type(prm3d),intent(in):: ps
    integer:: i
    character:: c*2
    print '(/,a)','   UF3d parameters of 3B for '&
         //trim(ps%csi)//'-'//trim(ps%csj)//'-'//trim(ps%csk)
    print '(6(a,1x))','     cb,csi,csj,csk,cknot = ',ps%cb,ps%csi,ps%csj,ps%csk,ps%cknot
    print '(a,2i3)','     nklead,nktrail = ',ps%nklead,ps%nktrail
    print '(a,6i3)','     nknij,nknik,nkncs = ', ps%nknij, ps%nknik, ps%nkncs
    print '(a,6i3)','     ncfij,ncfik,ncfcs = ', ps%ncfij, ps%ncfik, ps%ncfcs
    print '(a,3f6.3)','     rcij,rcik = ',ps%rcij, ps%rcik
    write(c,'(i2)') size(ps%knij)
    print '(a,'//c//'(1x,f7.3))','     knij =',(ps%knij(i),i=1,size(ps%knij))
    write(c,'(i2)') size(ps%knik)
    print '(a,'//c//'(1x,f7.3))','     knik =',(ps%knik(i),i=1,size(ps%knik))
    write(c,'(i2)') size(ps%kncs)
    print '(a,'//c//'(1x,f7.3))','     kncs =',(ps%kncs(i),i=1,size(ps%kncs))
    write(c,'(i2)') size(ps%cfij)
    print '(a,'//trim(c)//'es11.2)','     cfij =',(ps%cfij(i),i=1,size(ps%cfij))
    write(c,'(i2)') size(ps%cfik)
    print '(a,'//trim(c)//'es11.2)','     cfik =',(ps%cfik(i),i=1,size(ps%cfik))
    write(c,'(i2)') size(ps%cfij)
    print '(a,'//trim(c)//'es11.2)','     cfcs =',(ps%cfcs(i),i=1,size(ps%cfcs))

  end subroutine print_3bd
!=======================================================================
  function get_mem_uf3() result(dmem)
!
!  Compute and return memory usage in this module.
!
    real(8):: dmem
    integer:: i2b, i3b
    type(prm2):: p2
    type(prm3):: p3

    dmem = 0d0

    dmem = dmem +8d0*(size(aal2) +size(aal3) +size(strsl))
!!$    dmem = dmem +4d0*size(ls3b)

    do i2b=1,n2b
      p2 = prm2s(i2b)
      dmem = dmem +8d0*(size(p2%knots) +size(p2%coefs))
      if( allocated(p2%gwe) ) dmem = dmem +8d0*size(p2%gwe)
      if( allocated(p2%gwf) ) dmem = dmem +8d0*size(p2%gwf)
      if( allocated(p2%gws) ) dmem = dmem +8d0*size(p2%gws)
    enddo
    do i3b=1,n3b
      p3 = prm3s(i3b)
      dmem = dmem +8d0*( size(p3%knij) +size(p3%knik) +size(p3%knjk) &
           +size(p3%coefs) )
      if( allocated(p3%gwe) ) dmem = dmem +8d0*size(p3%gwe)
      if( allocated(p3%gwf) ) dmem = dmem +8d0*size(p3%gwf)
      if( allocated(p3%gws) ) dmem = dmem +8d0*size(p3%gws)
    enddo
    return
  end function get_mem_uf3
!=======================================================================
  function get_mem_uf3l() result(dmem)
!
!  Compute and return memory usage in this module.
!
    real(8):: dmem
    integer:: i2b, i3b
    type(prm2):: p2
    type(prm3l):: p3

    dmem = 0d0

    dmem = dmem +8d0*(size(aal2) +size(aal3) +size(strsl))
!!$    dmem = dmem +4d0*size(ls3b)

    do i2b=1,n2b
      p2 = prm2s(i2b)
      dmem = dmem +8d0*(size(p2%knots) +size(p2%coefs))
      if( allocated(p2%gwe) ) dmem = dmem +8d0*size(p2%gwe)
      if( allocated(p2%gwf) ) dmem = dmem +8d0*size(p2%gwf)
      if( allocated(p2%gws) ) dmem = dmem +8d0*size(p2%gws)
    enddo
    do i3b=1,n3b
      p3 = prm3ls(i3b)
      dmem = dmem +8d0*( size(p3%knots) +size(p3%coefs) )
      if( allocated(p3%gwe) ) dmem = dmem +8d0*size(p3%gwe)
      if( allocated(p3%gwf) ) dmem = dmem +8d0*size(p3%gwf)
      if( allocated(p3%gws) ) dmem = dmem +8d0*size(p3%gws)
    enddo
    return
  end function get_mem_uf3l
!=======================================================================
  subroutine dealloc_gwx_uf3()
!
!  Release memories for gw{e,f,s} for the purpose of efficiency.
!
    integer:: i2b,i3b
    do i2b=1,n2b
      if( allocated(prm2s(i2b)%gwe) ) deallocate(prm2s(i2b)%gwe)
      if( allocated(prm2s(i2b)%gwf) ) deallocate(prm2s(i2b)%gwf)
      if( allocated(prm2s(i2b)%gws) ) deallocate(prm2s(i2b)%gws)
    enddo
    do i3b=1,n3b
      if( allocated(prm3s(i3b)%gwe) ) deallocate(prm3s(i3b)%gwe)
      if( allocated(prm3s(i3b)%gwf) ) deallocate(prm3s(i3b)%gwf)
      if( allocated(prm3s(i3b)%gws) ) deallocate(prm3s(i3b)%gws)
    enddo
  end subroutine dealloc_gwx_uf3
end module UF3
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd lib"
!     End:
