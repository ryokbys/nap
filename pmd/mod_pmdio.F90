module pmdio
!-----------------------------------------------------------------------
!                     Last modified: <2025-05-03 18:27:59 KOBAYASHI Ryo>
!-----------------------------------------------------------------------
  use mod_precision
  use util, only: num_data
  implicit none
  save

contains
!=======================================================================
  function tag_encode(isp, ifmv, igrp, itot) result(tag_r)
!  Encode 4 tag components into a real(8) tag value (legacy format).
!  Always uses real(8) regardless of rp to preserve precision.
    integer,intent(in):: isp, ifmv, igrp(4), itot
    real(8):: tag_r
    tag_r = real(isp,8) &
          + real(ifmv,8) * 0.1d0 &
          + real(igrp(1),8) * 0.01d0 &
          + real(igrp(2),8) * 0.001d0 &
          + real(igrp(3),8) * 0.0001d0 &
          + real(igrp(4),8) * 0.00001d0 &
          + real(itot,8)    * 1.0d-14
  end function tag_encode
!=======================================================================
  subroutine tag_decode(tag_r, isp, ifmv, igrp, itot)
!  Decode a real(8) tag value into 4 tag components (legacy format).
    real(8),intent(in):: tag_r
    integer,intent(out):: isp, ifmv, igrp(4), itot
    real(8):: tmp
    isp     = int(tag_r)
    ifmv    = int(mod(tag_r * 10.0d0, 10.0d0))
    igrp(1) = int(mod(tag_r * 1.0d2,  10.0d0))
    igrp(2) = int(mod(tag_r * 1.0d3,  10.0d0))
    igrp(3) = int(mod(tag_r * 1.0d4,  10.0d0))
    igrp(4) = int(mod(tag_r * 1.0d5,  10.0d0))
    tmp     = tag_r * 1.0d5
    tmp     = tmp - int(tmp)
    itot    = nint(tmp * 1.0d9)
  end subroutine tag_decode
!=======================================================================
  function get_ntot_ascii(ionum,cfname) result(ntot)
    integer,intent(in):: ionum
    character(len=*),intent(in):: cfname
    integer:: ntot

    integer:: ia,ib,l,i
    real(rp):: hunit,h(3,3,0:1)
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
    real(rp):: h(3,3,0:1),hunit
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
  subroutine read_pmdtot_ascii(ionum,cfname,ntot,hunit,h,tagtot_isp, &
       tagtot_ifmv,tagtot_igrp,tagtot_itot,rtot,vtot)
    use pmdvars,only: ngrpmax
    integer,intent(in):: ionum,ntot
    character(len=*),intent(in):: cfname
    real(rp),intent(out):: hunit,h(3,3,0:1)
    integer,intent(out):: tagtot_isp(ntot),tagtot_ifmv(ntot)
    integer,intent(out):: tagtot_igrp(ngrpmax,ntot),tagtot_itot(ntot)
    real(rp),intent(out):: rtot(3,ntot),vtot(3,ntot)

    integer:: ia,ib,l,i,itmp,num
    real(8):: tmp_tag
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
!!$    read(ionum,*) (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
!.....H-matrix IO format changed since 2024-03-07
!.....thus check the format at the 1st line of H-matrix.
    read(ionum,'(a)') ctmp
    num = num_data(trim(ctmp),' ')
    if( num .ne. 6 ) then
      print *,' WARNIGN: The pmdini format seems to be old,'
      print *,'          see the document, http://ryokbys.web.nitech.ac.jp/contents/nap_docs/pmd-file.html'
      print *,'          and check the pmdini file carefully.'
!.....Backward compatibility...
      backspace(ionum)
      read(ionum,*) (h(ia,1,0),ia=1,3)
      read(ionum,*) (h(ia,2,0),ia=1,3)
      read(ionum,*) (h(ia,3,0),ia=1,3)
      read(ionum,*) (h(ia,1,1),ia=1,3)
      read(ionum,*) (h(ia,2,1),ia=1,3)
      read(ionum,*) (h(ia,3,1),ia=1,3)
    else
      backspace(ionum)
      read(ionum,*) ((h(ia,1,l),ia=1,3),l=0,1)
      read(ionum,*) ((h(ia,2,l),ia=1,3),l=0,1)
      read(ionum,*) ((h(ia,3,l),ia=1,3),l=0,1)
    endif
    h(1:3,1:3,0:1)= h(1:3,1:3,0:1)*hunit
    read(ionum,*) itmp
    if( itmp.ne.ntot ) then
      print *,' ERROR: itmp.ne.ntot'
      stop
    endif
    do i=1,ntot
      read(ionum,*) tmp_tag,rtot(1:3,i),vtot(1:3,i)
      call tag_decode(tmp_tag, tagtot_isp(i), tagtot_ifmv(i), tagtot_igrp(:,i), tagtot_itot(i))
      rtot(1,i) = pbc(rtot(1,i))
      rtot(2,i) = pbc(rtot(2,i))
      rtot(3,i) = pbc(rtot(3,i))
    enddo
    close(ionum)

  end subroutine read_pmdtot_ascii
!=======================================================================
  subroutine write_pmdtot_ascii(ionum,cfname,ntot,hunit,h,tagtot_isp, &
       tagtot_ifmv,tagtot_igrp,tagtot_itot,rtot,vtot,atot,epot,ekin,stnsr,lforce,istp)
    use pmdvars,only: has_specorder,specorder,lcomb_pos,ngrpmax
    include './params_unit.h'
    integer,intent(in):: ionum,ntot,istp
    character(len=*),intent(in) :: cfname
    real(rp),intent(in):: hunit,h(3,3,0:1)
    integer,intent(in):: tagtot_isp(ntot),tagtot_ifmv(ntot)
    integer,intent(in):: tagtot_igrp(ngrpmax,ntot),tagtot_itot(ntot)
    real(rp),intent(in):: rtot(3,ntot),vtot(3,ntot),atot(3,ntot)
    real(rp),intent(in):: epot,ekin,stnsr(3,3)
    logical,intent(in):: lforce

    integer:: ia,ib,l,i,msp,num
    real(rp):: atmp(3)
    real(8):: tmp_tag
    character(len=128):: cftmp
    logical:: lopen = .false.
    logical:: lclose = .false.
    logical,save:: l1st = .true.

    if( l1st ) then
      inquire(file=trim(cfname), number=num, opened=lopen)
      if( lcomb_pos .and. .not.lopen ) open(ionum,file=trim(cfname),status='replace')
      lclose = .false.
      l1st = .false.
    endif

    inquire(ionum, name=cftmp, number=num, opened=lopen)
!!$    print *,'name,number,cfname,opened = ',trim(cftmp),num,trim(cfname),lopen
    if( .not.lcomb_pos ) then
      open(ionum,file=trim(cfname),status='replace')
      lclose = .true.
    else if( .not.lopen .or. &
         (lopen .and. trim(cfname).ne.trim(cftmp)) ) then ! the unit number is used by other file
      if( lopen ) close(ionum)
      open(ionum,file=trim(cfname),status='replace')
      lclose = .true.
    endif
!.....Since 240307, there must be at least one comment line at the top of a configuration.
    write(ionum,'(a)') '#'
    write(ionum,'(a,2x,i0)') '#  timestep: ',istp
    if( has_specorder ) then
      msp = 0
      do i=1,ntot
        msp = max(msp,tagtot_isp(i))
      enddo
      write(ionum,'(a,9(2x,a))') '#  specorder: ',(trim(specorder(i)),i=1,msp)
    endif
    write(ionum,'(a,es14.6)') '#  potential_energy: ', epot
    write(ionum,'(a,es14.6)') '#  kinetic_energy:   ', ekin
!.....Positive stress as compressive, negative as tensile
    write(ionum,'(a,6es11.3)') '#  stress:   ',  &
           stnsr(1,1), stnsr(2,2), stnsr(3,3), &
           stnsr(3,2), stnsr(1,3), stnsr(1,2)
    if(lforce) write(ionum,'(a,l1)') '#  auxiliary_data:  fx fy fz'
    write(ionum,'(a)') '#'
    write(ionum,'('//rpfmt3//')') hunit
!!$    write(ionum,'(3es23.14e3)') (((h(ia,ib,l)/hunit,ia=1,3) &
!!$         ,ib=1,3),l=0,1)
!.....H-matrix IO format changed since 2024-03-07
    write(ionum,'(3'//rpfmt3//', 3es12.3e3)') ((h(ia,1,l)/hunit,ia=1,3),l=0,1)
    write(ionum,'(3'//rpfmt3//', 3es12.3e3)') ((h(ia,2,l)/hunit,ia=1,3),l=0,1)
    write(ionum,'(3'//rpfmt3//', 3es12.3e3)') ((h(ia,3,l)/hunit,ia=1,3),l=0,1)
    write(ionum,'(i10)') ntot
!.....All the length values (r,v,a) are scaled by h-mat in pmd format
    if( lforce ) then ! write forces in [eV/A/A] (scaled by h-mat)
      do i=1,ntot
        tmp_tag = tag_encode(tagtot_isp(i),tagtot_ifmv(i),tagtot_igrp(:,i),tagtot_itot(i))
        write(ionum,'(f17.14, 6'//rpfmt//', 11es12.4)') tmp_tag &
             ,rtot(1:3,i) ,vtot(1:3,i) ,atot(1:3,i)    ! dt
      enddo
    else
      do i=1,ntot
        tmp_tag = tag_encode(tagtot_isp(i),tagtot_ifmv(i),tagtot_igrp(:,i),tagtot_itot(i))
        write(ionum,'(f17.14, 6'//rpfmt//', 11es12.4)') tmp_tag &
             ,rtot(1:3,i),vtot(1:3,i)    ! dt
      enddo
    endif
    if( lclose ) close(ionum)
    return
  end subroutine write_pmdtot_ascii
!=======================================================================
  subroutine read_pmdtot_bin(ionum,cfname,ntot,hunit,h,tagtot_isp, &
       tagtot_ifmv,tagtot_igrp,tagtot_itot,rtot,vtot)
    use pmdvars,only: specorder,ngrpmax
    integer,intent(in):: ionum,ntot
    character(len=*),intent(in):: cfname
    real(rp),intent(out):: hunit,h(3,3,0:1)
    integer,intent(out):: tagtot_isp(ntot),tagtot_ifmv(ntot)
    integer,intent(out):: tagtot_igrp(ngrpmax,ntot),tagtot_itot(ntot)
    real(rp),intent(out):: rtot(3,ntot),vtot(3,ntot)

    integer:: ia,ib,l,i,msp,itmp
    real(8),allocatable:: tmp_tags(:)

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
    allocate(tmp_tags(ntot))
    read(ionum) tmp_tags(1:ntot)
    read(ionum) rtot(1:3,1:ntot)
    read(ionum) vtot(1:3,1:ntot)
    close(ionum)
    do i=1,ntot
      call tag_decode(tmp_tags(i), tagtot_isp(i), tagtot_ifmv(i), tagtot_igrp(:,i), tagtot_itot(i))
      rtot(1,i) = pbc(rtot(1,i))
      rtot(2,i) = pbc(rtot(2,i))
      rtot(3,i) = pbc(rtot(3,i))
    enddo
    deallocate(tmp_tags)

  end subroutine read_pmdtot_bin
!=======================================================================
  subroutine write_pmdtot_bin(ionum,cfname,ntot,hunit,h,tagtot_isp, &
       tagtot_ifmv,tagtot_igrp,tagtot_itot,rtot,vtot)
    use pmdvars,only: specorder,ngrpmax
    include './params_unit.h'
    integer,intent(in):: ionum,ntot
    character(len=*),intent(in) :: cfname
    real(rp),intent(in):: hunit,h(3,3,0:1)
    integer,intent(in):: tagtot_isp(ntot),tagtot_ifmv(ntot)
    integer,intent(in):: tagtot_igrp(ngrpmax,ntot),tagtot_itot(ntot)
    real(rp),intent(in):: rtot(3,ntot),vtot(3,ntot)

    integer:: ia,ib,l,i,msp
    real(8),allocatable:: tmp_tags(:)

    open(ionum,file=cfname,form='unformatted',status='replace')
    msp = 0
    do ia=1,ntot
      msp = max(msp,tagtot_isp(ia))
    enddo
    write(ionum) msp
    write(ionum) (specorder(i),i=1,msp)
    write(ionum) hunit
    write(ionum) (((h(ia,ib,l)/hunit,ia=1,3),ib=1,3),l=0,1)
    write(ionum) ntot
    allocate(tmp_tags(ntot))
    do i=1,ntot
      tmp_tags(i) = tag_encode(tagtot_isp(i),tagtot_ifmv(i),tagtot_igrp(:,i),tagtot_itot(i))
    enddo
    write(ionum) tmp_tags(1:ntot)
    deallocate(tmp_tags)
    write(ionum) rtot(1:3,1:ntot)
    write(ionum) vtot(1:3,1:ntot)
    close(ionum)

  end subroutine write_pmdtot_bin
!=======================================================================
  subroutine write_dump(ionum,cfname,ntot,hunit,h,tagtot_isp,tagtot_itot, &
       rtot,vtot,atot,stot,ekitot,epitot,naux,auxtot,istp)
!
!     Write atomic configuration in LAMMPS-dump format file.
!
    use pmdvars,only: ndumpaux,cdumpauxarr,specorder,has_specorder,&
         iaux_chg,iaux_tei,iaux_clr,iaux_edesc,lcomb_pos
    use util,only: iauxof
    use time,only: accum_time
    implicit none
    include "mpif.h"
    integer,intent(in):: ionum,ntot,naux,istp
    integer,intent(in):: tagtot_isp(ntot),tagtot_itot(ntot)
    character(len=*),intent(in) :: cfname
    real(rp),intent(in):: hunit,h(3,3,0:1)
    real(rp),intent(in):: rtot(3,ntot),vtot(3,ntot), &
         atot(3,ntot),stot(3,3,ntot),ekitot(3,3,ntot), &
         epitot(ntot),auxtot(naux,ntot)

    integer:: i,j,k,l,is,idlmp
    real(rp):: xi(3),ri(3),xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz, &
         xlo_bound,xhi_bound,ylo_bound,yhi_bound, &
         zlo_bound,zhi_bound,st(3,3)
!!$    integer,external:: itotOf
!!$    real(rp),allocatable,save:: rlmp(:,:),vlmp(:,:)
    integer,save:: ndlmp,ndim
    real(rp),allocatable,save:: dlmp(:,:)
    character(len=3),save:: cndlmp
    character(len=3):: csp
    character(len=6):: caux
    real(rp),parameter:: tiny = 1e-14_rp
    logical,save:: l1st = .true.

    if( l1st ) then
      ndlmp = 3 +ndumpaux
      ndim = max(ndlmp,6)
!.....If vx,vy,vz are included in dumpauxarr, take that into account
!!$      if( idumpauxof('vx').gt.0 ) ndlmp = ndlmp -1
!!$      if( idumpauxof('vy').gt.0 ) ndlmp = ndlmp -1
!!$      if( idumpauxof('vz').gt.0 ) ndlmp = ndlmp -1
      write(cndlmp,'(i0)') ndlmp
      if( allocated(dlmp) .and. size(dlmp).ne.ndim*ntot) deallocate(dlmp)
      allocate(dlmp(ndim,ntot))
      if( lcomb_pos ) open(ionum,file=trim(cfname),status='replace')
      l1st = .false.
    endif

    if( .not.allocated(dlmp) ) then
      allocate(dlmp(ndim,ntot))
    else if( size(dlmp).ne.ndim*ntot ) then
      deallocate(dlmp)
      allocate(dlmp(ndim,ntot))
    endif

    call pmd2lammps(h,ntot,rtot,dlmp(1:3,:),vtot,dlmp(4:6,:) &
         ,xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz)
    xlo_bound = xlo +min(0.0_rp, xy, xz, xy+xz)
    xhi_bound = xhi +max(0.0_rp, xy, xz, xy+xz)
    ylo_bound = ylo +min(0.0_rp, yz)
    yhi_bound = yhi +max(0.0_rp, yz)
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
      else if( trim(caux).eq.'fx' ) then
        dlmp(idlmp,:) =h(1,1,0)*atot(1,:)+h(1,2,0)*atot(2,:)+h(1,3,0)*atot(3,:)
      else if( trim(caux).eq.'fy' ) then
        dlmp(idlmp,:) =h(2,1,0)*atot(1,:)+h(2,2,0)*atot(2,:)+h(2,3,0)*atot(3,:)
      else if( trim(caux).eq.'fz' ) then
        dlmp(idlmp,:) =h(3,1,0)*atot(1,:)+h(3,2,0)*atot(2,:)+h(3,3,0)*atot(3,:)
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
        dlmp(idlmp,:) = auxtot(iaux_chg,:)
      else if( trim(caux).eq.'tei' ) then
        dlmp(idlmp,:) = auxtot(iaux_tei,:)
      else if( trim(caux).eq.'clr' ) then
        dlmp(idlmp,:) = auxtot(iaux_clr,:)
      else if( trim(caux).eq.'edesc' ) then
        dlmp(idlmp,:) = auxtot(iaux_edesc,:)
      endif
    enddo

    if( .not. lcomb_pos ) open(ionum,file=trim(cfname),status='replace')
    write(ionum,'(a)') 'ITEM: TIMESTEP'
    write(ionum,'(3x,i0)') istp
    write(ionum,'(a)') 'ITEM: NUMBER OF ATOMS'
    write(ionum,'(3x,i0)') ntot
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
      write(ionum,'(i8)',advance='no') tagtot_itot(i)
      if( has_specorder ) then
        is = tagtot_isp(i)
        csp = specorder(is)
        write(ionum,'(a4)',advance='no') trim(csp)
      else
        write(ionum,'(i3)',advance='no') tagtot_isp(i)
      endif
      write(ionum,'(3f12.5)',advance='no') dlmp(1:3,i)  ! pos
      write(ionum,'('//trim(cndlmp)//'es11.2e3)') dlmp(4:ndlmp,i)  ! except pos
    enddo

    if( .not. lcomb_pos ) close(ionum)
  end subroutine write_dump
!=======================================================================
  subroutine write_extxyz(ionum,cfname,ntot,hunit,h,tagtot_isp, &
       rtot,vtot,atot,stot,ekitot,epitot,epot,ekin,stnsr,istp)
!
!  Format of the extxyz is like the following:
!  ---
!  8
!  Lattice="5.44 0.0 0.0 0.0 5.44 0.0 0.0 0.0 5.44" Properties=species:S:1:pos:R:3:vel:R:3:frc:R:3
!  Si  0.00000000  0.00000000  0.00000000  0.00e-00  0.00e-00  0.00e-00  1.6215e-03 -6.4788e-03  2.6939e-05
!  Si  1.36000000  1.36000000  1.36000000  0.00e-00  0.00e-00  0.00e-00 -9.4438e-05 -5.7187e-04 -2.6944e-04
!  ...
!  ---
    use pmdvars,only: has_specorder,specorder,lcomb_pos
    use util, only: basename
    include './params_unit.h'
    integer,intent(in):: ionum,ntot,istp
    integer,intent(in):: tagtot_isp(ntot)
    character(len=*),intent(in) :: cfname
    real(rp),intent(in):: hunit,h(3,3,0:1)
    real(rp),intent(in):: rtot(3,ntot),vtot(3,ntot), &
         atot(3,ntot),stot(3,3,ntot),ekitot(3,3,ntot),epitot(ntot)
    real(rp),intent(in):: epot,ekin,stnsr(3,3)

    integer:: ia,ja,ib,l,i,msp,num,is
    real(rp):: atmp(3),ri(3),vi(3),ai(3),epi,eki
    character(len=3):: csp
    character(len=128):: cftmp,str
    character(len=256):: cpath
    logical:: lopen = .false.
    logical:: lclose = .false.
    logical,save:: l1st = .true.

    if( l1st ) then
      inquire(file=trim(cfname), number=num, opened=lopen)
      if( .not.lopen ) open(ionum,file=trim(cfname),status='replace')
      lclose = .false.
      l1st = .false.
    endif

    inquire(ionum, name=cpath, number=num, opened=lopen)
    cftmp = trim(basename(cpath))
    if( .not.lopen .or. &
         (lopen .and. trim(cfname).ne.trim(cftmp)) ) then ! the unit number is used by other file
      if( lopen ) close(ionum)
      open(ionum,file=trim(cfname),status='replace')
      lclose = .true.
    endif
!!$    print *,'name,number,cfname,opened,closed = ',trim(cftmp),num,trim(cfname),lopen,lclose

!===== 1st line: number of atoms
    write(ionum,'(2x,i0)') ntot

!===== Property line starts
    write(ionum,'(a)',advance='no') 'Lattice="'
    do ia=1,3
      do ja=1,3
        write(str,'(f20.5)') h(ja,ia,0)
        write(ionum,'(1x,a)',advance='no') trim(adjustl(str))
      enddo
    enddo
    write(ionum,'(a)',advance='no') '"'

    write(ionum,'(a)',advance='no') ' Properties=species:S:1:pos:R:3:vel:R:3:forces:R:3:epot:R:1:ekin:R:1 energy='
    write(str,'(f20.6)') epot
    write(ionum,'(a)',advance='no') trim(adjustl(str))
    write(ionum,'(a)',advance='no') ' stress="'
    do ia=1,3
      do ja=1,3
        write(str,'(f20.5)') stnsr(ia,ja)/(-gpa2up)
        write(ionum,'(1x,a)',advance='no') trim(adjustl(str))
      enddo
    enddo
    write(ionum,'(a)') '"'
!===== Property line ends

!===== Atom information starts
    do i=1,ntot
      is = tagtot_isp(i)
      csp = specorder(is)
      ri(1:3)= h(1:3,1,0)*rtot(1,i) +h(1:3,2,0)*rtot(2,i) +h(1:3,3,0)*rtot(3,i)
      vi(1:3)= h(1:3,1,0)*vtot(1,i) +h(1:3,2,0)*vtot(2,i) +h(1:3,3,0)*vtot(3,i)
      ai(1:3)= h(1:3,1,0)*atot(1,i) +h(1:3,2,0)*atot(2,i) +h(1:3,3,0)*atot(3,i)
      eki = ekitot(1,1,i) +ekitot(2,2,i) +ekitot(3,3,i)
      epi = epitot(i)
      write(ionum,'(1x,a3,3(1x,f12.5),3(1x,es11.3),3(1x,es11.3),2(1x,f8.4))') &
           trim(csp), &
           ri(1:3), vi(1:3), ai(1:3),epi,eki
    enddo
!===== Atom information ends

    if( lclose ) close(ionum)
    return
  end subroutine write_extxyz
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
    use pmdvars,only: boundary
    integer,intent(in):: ntot
    real(rp),intent(in):: h(3,3),rtot(3,ntot),vtot(3,ntot)
    real(rp),intent(out):: xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz &
         ,rlmp(3,ntot),vlmp(3,ntot)

    integer:: i,lxy,lxz,lyz,ixyz
    real(rp):: a0(3),b0(3),c0(3),a1(3),a2(3),a3(3) &
         ,b1(3),b2(3),b3(3),rt(3),vt(3),amat(3,3),bmat(3,3) &
         ,x,y,z,a23(3),a31(3),a12(3),vol,xyp
    real(rp):: a,b,c,alpha,beta,gamma

    xlo = 0.0_rp
    ylo = 0.0_rp
    zlo = 0.0_rp
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
    amat(1:3,1:3) = 0.0_rp
    amat(1,1:3) = a23(1:3)
    amat(2,1:3) = a31(1:3)
    amat(3,1:3) = a12(1:3)
    b1(1:3) = (/ x, 0.0_rp, 0.0_rp /)
    b2(1:3) = (/ xy,  y, 0.0_rp /)
    b3(1:3) = (/ xz, yz,   z /)
    bmat(1:3,1:3) = 0.0_rp
    bmat(1:3,1) = b1(1:3)
    bmat(1:3,2) = b2(1:3)
    bmat(1:3,3) = b3(1:3)
    xyp = xy -lxy*x
    do i=1,ntot
      rlmp(1:3,i) = 0.0_rp
!!$      call shift_pos_for_lammps(rtot(1,i),rlmp(1,i),lxy,lxz,lyz &
!!$           ,x,y,z,yz,xz,xy)
!.....Shift positions
      rlmp(1:3,i) = rtot(1:3,i)
      rlmp(2,i) = rlmp(2,i) -lyz*rtot(3,i)
      rlmp(1,i) = rlmp(1,i) -lxz*rtot(3,i) +(rtot(2,i)*xyp -rlmp(2,i)*xy)/x
      do ixyz=1,3
        if( boundary(ixyz:ixyz).eq.'p' ) rlmp(ixyz,i) = pbc(rlmp(ixyz,i))
      enddo
!.....Velocity is in real unit (A/fs)
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
    real(rp),intent(in):: r(3),x,y,z,yz,xz,xy
    integer,intent(in):: lxy,lxz,lyz
    real(rp),intent(out):: rn(3)

    integer:: i
    real(rp):: xyp

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
    real(rp),intent(in):: x
    real(rp):: pbc

    if( x.lt.0.0_rp ) then
      pbc = x -int(x) +1.0_rp
    else if( x.ge.1.0_rp ) then
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
    real(rp),intent(in):: hmat(3,3)
    real(rp),intent(out):: xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz

    real(rp):: a0(3),b0(3),c0(3)
    real(rp):: a,b,c,alpha,beta,gamma

    xlo = 0.0_rp
    ylo = 0.0_rp
    zlo = 0.0_rp
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
!  "# specorder: Al Mg Si"
!
!  Currently available options are:
!    - "specorder:", Species order. The number of species limited up to 9.
!
    use pmdvars,only: specorder,has_specorder,iprint,has_forces
    include "./const.h"
    character(len=*),intent(in):: cline

    integer:: iopt1,isp,num
    real(rp):: opt1, opt2
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
!     compile-command: "make pmd lib"
!     End:
