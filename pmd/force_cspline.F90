module cspline
!-----------------------------------------------------------------------
!  Parallel implementation of cubic spline force field.
!-----------------------------------------------------------------------
  implicit none
  save

  real(8),parameter:: pi = 3.14159265358979d0

  character(len=128):: paramsdir = '.'
!.....parameter file name
  character(128),parameter:: cpfname= 'in.params.cspline'
  integer,parameter:: ionum = 50
  integer:: nspmax
  integer,allocatable:: idpair(:,:)
  integer,allocatable:: idtriplet(:,:,:)
  
  integer:: nspl = 0 ! # of splines

  type spline
    character(128):: ctype
    integer:: isp,jsp
    integer:: ksp = -1
    integer:: npnts = 0
    real(8):: rcut
    real(8),allocatable:: pnts(:),vals(:)
    real(8),allocatable:: coefs(:,:)
  end type spline
  type(spline),allocatable:: spls(:)

  real(8):: epot_cspln
contains
!=======================================================================
  subroutine force_cspline(namax,natm,tag,ra,nnmax,aa,strs,h,hi,tcom &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,lstrs,iprint,l1st)
    include "mpif.h"
    include "./params_unit.h"
    integer,intent(in):: namax,natm,nnmax,nismax,lspr(0:nnmax,namax)&
         ,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_md_world,myid_md,nex(3)
    real(8),intent(in):: ra(3,namax),h(3,3),hi(3,3),sv(3,6) &
         ,rc,tag(namax)
    real(8),intent(inout):: tcom
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: lstrs,l1st

    integer:: ispln,jj,ia,ja,kk,ka,ierr,is,js,ks,ispl,jxyz
    real(8):: epotl,xi(3),xj(3),xij(3),rij(3),drij(3) &
         ,xk(3),xik(3),rik(3),rjk(3),drik(3),dcsdj(3),dcsdk(3),tmpj(3) &
         ,tmpk(3),csn,dfcij,dfcik,dgcs,dgdij,dgdik,dij,dij2,diji &
         ,dik,diki,dik2,epotl2,epotl3,fcij,fcik,gijk,phi,dphi,rc2,rct &
         ,rct2,tmp,eta3,texpij,texpik,texpjk,tmpjk(3),djk2,djk,fcjk,dfcjk &
         ,drjk(3),dgdjk,xjk(3)
    real(8),save:: rcmax,rcmax2
    character(len=128):: ctype
    real(8),save,allocatable:: aa2(:,:),aa3(:,:),strsl(:,:,:)
    logical,save:: l3b,l2b
    logical:: l3b_exist
    type(spline):: spl

!.....Only at the 1st call
    if( l1st ) then
      nspmax = nismax
      if( allocated(idpair) ) deallocate(idpair,idtriplet)
      allocate(idpair(nspmax,nspmax),idtriplet(nspmax,nspmax,nspmax))
      call read_params_cspline(myid_md,mpi_md_world,iprint)
      if( iprint.gt.1 .and. myid_md.eq.0 ) then
        call write_cspline_curves()
      endif
!.....Check rc
      rcmax = 0d0
      do ispl=1,nspl
        spl = spls(ispl)
        rcmax = max(rcmax,spl%rcut)
      enddo
      rcmax2 = rcmax*rcmax
      if( rc.lt.rcmax ) then
        if( myid_md.eq.0 ) then
          print *,'ERROR: Cutoff radius is not appropriate in force_cspline.'
          print *,'  rc, rcmax = ',rc,rcmax
        endif
        call mpi_finalize(ierr)
        stop 1
      endif
!.....Check whether 2-body terms exist
      l2b = .false.
      do is=1,nismax
        do js=1,nismax
          if( idpair(is,js).gt.0 ) then
            l2b = .true.
            goto 10
          endif
        enddo
      enddo
10    continue
!.....Check whether 3-body terms exist
      l3b = .false.
      do is=1,nismax
        do js=1,nismax
          do ks=1,nismax
            if( idtriplet(is,js,ks).gt.0 ) then
              l3b = .true.
              goto 11
            endif
          enddo
        enddo
      enddo
11    continue

      if( allocated(strsl) ) deallocate(strsl)
      allocate(strsl(3,3,namax))
      if( allocated(aa2) ) deallocate(aa2,aa3)
      allocate(aa2(3,namax),aa3(3,namax))
    endif

    if( size(strsl).lt.3*3*namax ) then
      deallocate(strsl)
      allocate(strsl(3,3,namax))
    endif
    if( size(aa2).lt.3*namax ) then
      deallocate(aa2,aa3)
      allocate(aa2(3,namax),aa3(3,namax))
    endif
    strsl(1:3,1:3,1:namax) = 0d0

    epotl= 0d0
    epi(1:natm+nb)= 0d0

!.....2-body term
    epotl2 = 0d0
    aa2(1:3,1:natm+nb) = 0d0
    if( .not. l2b ) goto 20
    do ia=1,natm
      xi(1:3)= ra(1:3,ia)
      is = int(tag(ia))
      do jj=1,lspr(0,ia)
        ja=lspr(jj,ia)
        if( ja.eq.0 ) exit
        if( ja.le.ia ) cycle
        js = int(tag(ja))
        ispl = idpair(is,js)
        if( ispl.le.0 ) cycle
        spl = spls(ispl)
        rc2 = spl%rcut
        xj(1:3)= ra(1:3,ja)
        xij(1:3)= xj(1:3)-xi(1:3)
        rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij2= rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3)
        if( dij2.gt.rc2*rc2 ) cycle
        dij = sqrt(dij2)
        diji= 1d0/dij
        call eval_spline(dij,spl%npnts,spl%pnts,spl%vals,spl%coefs,phi,dphi)
!.....Potential
        epi(ia)= epi(ia) +phi*0.5d0
        epotl2 = epotl2 +phi*0.5d0
        if( ja.le.natm ) then
          epi(ja)= epi(ja) +phi*0.5d0
          epotl2 = epotl2 +phi*0.5d0
        endif
!.....Force
        drij(1:3)= -rij(1:3)*diji
        aa2(1:3,ia)= aa2(1:3,ia) -dphi*drij(1:3)
        aa2(1:3,ja)= aa2(1:3,ja) +dphi*drij(1:3)
!.....Stress
        do jxyz=1,3
          strsl(1:3,jxyz,ia)= strsl(1:3,jxyz,ia) &
               -0.5d0*rij(jxyz)*(-dphi*drij(1:3))
          strsl(1:3,jxyz,ja)= strsl(1:3,jxyz,ja) &
               -0.5d0*rij(jxyz)*(-dphi*drij(1:3))
        enddo
      enddo
    enddo
20  continue

!.....3-body term
    epotl3 = 0d0
    aa3(1:3,1:natm+nb) = 0d0
    if( .not. l3b ) goto 21
    do ia=1,natm
      xi(1:3)= ra(1:3,ia)
      is = int(tag(ia))
      do jj=1,lspr(0,ia)-1
        ja= lspr(jj,ia)
        js= int(tag(ja))
        l3b_exist = .false.
        do ks=1,nspmax
          if( idtriplet(is,js,ks).gt.0 ) then
            l3b_exist = .true.
            exit
          endif
        enddo
        if( .not. l3b_exist ) cycle
        xj(1:3)= ra(1:3,ja)
        xij(1:3)= xj(1:3)-xi(1:3)
        rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij2= rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3)
        if( dij2.gt.rcmax2 ) cycle
        dij= sqrt(dij2)
        do kk=jj+1,lspr(0,ia)
          ka= lspr(kk,ia)
          ks= int(tag(ka))
          ispl = idtriplet(is,js,ks)
          if( ispl.le.0 ) cycle
          spl = spls(ispl)
          rct = spl%rcut
          rct2 = rct*rct
          if( dij2.gt.rct2 ) cycle
          xk(1:3)= ra(1:3,ka)
          xik(1:3)= xk(1:3)-xi(1:3)
          rik(1:3)= h(1:3,1)*xik(1) +h(1:3,2)*xik(2) +h(1:3,3)*xik(3)
          dik2= rik(1)*rik(1) +rik(2)*rik(2) +rik(3)*rik(3)
          if( dik2.gt.rct2 ) cycle
          dik = sqrt(dik2)
          csn = (rij(1)*rik(1) +rij(2)*rik(2) +rij(3)*rik(3))/(dij*dik)
          call eval_spline(csn,spl%npnts,spl%pnts,spl%vals,spl%coefs,phi,dphi)
          fcij = fc1(dij,0d0,rct)
          dfcij = dfc1(dij,0d0,rct)
          fcik = fc1(dik,0d0,rct)
          dfcik = dfc1(dik,0d0,rct)
          ctype = spl%ctype
          if( trim(ctype).eq.'angular'&
               .or. trim(ctype).eq.'angular1' ) then
            eta3 = 0.5d0 /(rct/2)**2
            texpij = exp(-eta3 *dij2)
            texpik = exp(-eta3 *dik2)
            tmp = phi *texpij *texpik
            gijk = tmp *fcij *fcik
!.....Potential
            epi(ia)= epi(ia) +gijk
            epotl3 = epotl3  +gijk
!.....Force
            drij(1:3) = -rij(1:3)/dij
            drik(1:3) = -rik(1:3)/dik
            dgdij = dfcij*fcik*tmp +tmp*(-2d0*eta3*dij)*fcij*fcik
            dgdik = dfcik*fcij*tmp +tmp*(-2d0*eta3*dik)*fcij*fcik
            dgcs = dphi *texpij*texpik *fcij*fcik
            dcsdj(1:3)= rik(1:3)/dij/dik -rij(1:3)*csn/dij**2
            dcsdk(1:3)= rij(1:3)/dij/dik -rik(1:3)*csn/dik**2
!!$          dcsdi(1:3)= -dcsdj(1:3)  -dcsdk(1:3)
            tmpj(1:3)= dgcs*dcsdj(1:3) -drij(1:3)*dgdij
            tmpk(1:3)= dgcs*dcsdk(1:3) -drik(1:3)*dgdik
!!$            print *,'dgdij,dgdik,dgcs=',dgdij,dgdik,dgcs
!!$            print *,'tmpj(:)=',tmpj(:)
!!$            print *,'tmpk(:)=',tmpk(:)
            aa3(1:3,ja)= aa3(1:3,ja) -tmpj(1:3)
            aa3(1:3,ka)= aa3(1:3,ka) -tmpk(1:3)
            aa3(1:3,ia)= aa3(1:3,ia) +tmpj(1:3)+tmpk(1:3)
!.....Stress
            if( lstrs ) then
              do jxyz=1,3
                strsl(1:3,jxyz,ia)= strsl(1:3,jxyz,ia) &
                     -0.5d0*rij(jxyz)*tmpj(1:3) -0.5d0*rik(jxyz)*tmpk(1:3)
                strsl(1:3,jxyz,ja)= strsl(1:3,jxyz,ja) &
                     -0.5d0*rij(jxyz)*tmpj(1:3)
                strsl(1:3,jxyz,ka)= strsl(1:3,jxyz,ka) &
                     -0.5d0*rik(jxyz)*tmpk(1:3)
              enddo
            endif
          else if( trim(ctype).eq.'angular2' ) then
            xjk(1:3)= xk(1:3)-xj(1:3)
            rjk(1:3)= h(1:3,1)*xjk(1) +h(1:3,2)*xjk(2) +h(1:3,3)*xjk(3)
            djk2= rjk(1)*rjk(1) +rjk(2)*rjk(2) +rjk(3)*rjk(3)
            if( djk2.gt.rct2 ) cycle
            djk= sqrt(djk2)
            fcjk = fc1(djk,0d0,rct)
            dfcjk = dfc1(djk,0d0,rct)
            eta3 = 0.5d0 /(rct/2)**2
            texpij = exp(-eta3 *dij2)
            texpik = exp(-eta3 *dik2)
            texpjk = exp(-eta3 *djk2)
            tmp = phi *texpij *texpik *texpjk
            gijk = tmp *fcij *fcik *fcjk
!.....Potential
            epi(ia)= epi(ia) +gijk
            epotl3 = epotl3  +gijk
!.....Force
            drij(1:3) = -rij(1:3)/dij
            drik(1:3) = -rik(1:3)/dik
            drjk(1:3) = -rjk(1:3)/djk
            dgdij = dfcij*fcik*fcjk*tmp +tmp*(-2d0*eta3*dij)*fcij*fcik*fcjk
            dgdik = dfcik*fcij*fcjk*tmp +tmp*(-2d0*eta3*dik)*fcij*fcik*fcjk
            dgdjk = dfcjk*fcij*fcik*tmp +tmp*(-2d0*eta3*djk)*fcij*fcik*fcjk
            dgcs = dphi *texpij*texpik*texpjk *fcij*fcik*fcjk
            dcsdj(1:3)= rik(1:3)/dij/dik -rij(1:3)*csn/dij**2
            dcsdk(1:3)= rij(1:3)/dij/dik -rik(1:3)*csn/dik**2
!!$          dcsdi(1:3)= -dcsdj(1:3)  -dcsdk(1:3)
            tmpj(1:3)= dgcs*dcsdj(1:3) -drij(1:3)*dgdij
            tmpk(1:3)= dgcs*dcsdk(1:3) -drik(1:3)*dgdik
            tmpjk(1:3)= drjk(1:3)*dgdjk
            aa3(1:3,ja)= aa3(1:3,ja) -tmpj(1:3) -tmpjk(1:3)
            aa3(1:3,ka)= aa3(1:3,ka) -tmpk(1:3) +tmpjk(1:3)
            aa3(1:3,ia)= aa3(1:3,ia) +tmpj(1:3)+tmpk(1:3)
!.....Stress
            if( lstrs ) then
              do jxyz=1,3
                strsl(1:3,jxyz,ia)= strsl(1:3,jxyz,ia) &
                     -0.5d0*rij(jxyz)*tmpj(1:3) -0.5d0*rik(jxyz)*tmpk(1:3)
                strsl(1:3,jxyz,ja)= strsl(1:3,jxyz,ja) &
                     -0.5d0*rij(jxyz)*tmpj(1:3) -0.5d0*rjk(jxyz)*tmpjk(1:3)
                strsl(1:3,jxyz,ka)= strsl(1:3,jxyz,ka) &
                     -0.5d0*rik(jxyz)*tmpk(1:3) -0.5d0*rjk(jxyz)*tmpjk(1:3)
              enddo
            endif
          else if( trim(ctype).eq.'angular3' &
               .or. trim(ctype).eq.'angular4' ) then
            gijk = phi *fcij *fcik
!.....Potential
            epi(ia)= epi(ia) +gijk
            epotl3 = epotl3  +gijk
!.....Force
            drij(1:3) = -rij(1:3)/dij
            drik(1:3) = -rik(1:3)/dik
            dgdij = dfcij*fcik*tmp
            dgdik = dfcik*fcij*tmp
            dgcs = dphi*fcij*fcik
            dcsdj(1:3)= rik(1:3)/dij/dik -rij(1:3)*csn/dij**2
            dcsdk(1:3)= rij(1:3)/dij/dik -rik(1:3)*csn/dik**2
!!$          dcsdi(1:3)= -dcsdj(1:3)  -dcsdk(1:3)
            tmpj(1:3)= dgcs*dcsdj(1:3) -drij(1:3)*dgdij
            tmpk(1:3)= dgcs*dcsdk(1:3) -drik(1:3)*dgdik
            aa3(1:3,ja)= aa3(1:3,ja) -tmpj(1:3)
            aa3(1:3,ka)= aa3(1:3,ka) -tmpk(1:3)
            aa3(1:3,ia)= aa3(1:3,ia) +tmpj(1:3)+tmpk(1:3)
!.....Stress
            if( lstrs ) then
              do jxyz=1,3
                strsl(1:3,jxyz,ia)= strsl(1:3,jxyz,ia) &
                     -0.5d0*rij(jxyz)*tmpj(1:3) -0.5d0*rik(jxyz)*tmpk(1:3)
                strsl(1:3,jxyz,ja)= strsl(1:3,jxyz,ja) &
                     -0.5d0*rij(jxyz)*tmpj(1:3)
                strsl(1:3,jxyz,ka)= strsl(1:3,jxyz,ka) &
                     -0.5d0*rik(jxyz)*tmpk(1:3)
              enddo
            endif
          endif
        enddo
      enddo
    enddo
21  continue

    call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_md_world,strsl,9)
    call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_md_world,aa3,3)
    call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_md_world,epi,1)

!.....Sum up forces and stress
    aa(1:3,1:natm)= aa(1:3,1:natm) +aa2(1:3,1:natm) +aa3(1:3,1:natm)
    strs(1:3,1:3,1:natm)= strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)

    epotl = epotl2 +epotl3
    call mpi_allreduce(epotl,epot_cspln,1,mpi_real8,mpi_sum &
         ,mpi_md_world,ierr)
    epot = epot +epot_cspln
    return
  end subroutine force_cspline
!=======================================================================
  subroutine eval_spline(r,npnts,pnts,vals,coefs,spl,dspl)
    real(8),intent(in):: r
    integer,intent(in):: npnts
    real(8),intent(in):: pnts(npnts),vals(npnts),coefs(4,npnts-1)
    real(8),intent(inout):: spl,dspl

    integer:: i
    real(8):: a(4),r2,rt

    if( r.ge.pnts(npnts) ) then ! outside the right edge of the range
      a(:) = coefs(:,npnts-1)
      rt = pnts(npnts) ! right-most data point
      r2 = rt*rt
      dspl = a(2) +2d0*a(3)*rt +3d0*a(4)*r2 ! gradient of the right-most data point
      spl = vals(npnts) +dspl*(r -rt) ! extrapolation
    else if( r.lt.pnts(1) ) then ! outside left edge of the range
      a(:) = coefs(:,1)
      rt = pnts(1)  ! left-most data point
      r2 = rt*rt
      dspl = a(2) +2d0*a(3)*rt +3d0*a(4)*r2
      spl = vals(1) +dspl*(r -rt)
    else
      do i=2,npnts
        if( r.lt.pnts(i) ) exit
      enddo
      a(:) = coefs(:,i-1)
      r2 = r*r
      spl = a(1) +a(2)*r +a(3)*r2 +a(4)*r*r2
      dspl = a(2) +2d0*a(3)*r +3d0*a(4)*r2
    endif
    return
  end subroutine eval_spline
!=======================================================================
  subroutine spl2d(rij,npnts,pnts,coefs,spl,dspl)
    real(8),intent(in):: rij
    integer,intent(in):: npnts
    real(8),intent(in):: pnts(npnts),coefs(4,npnts-1)
    real(8),intent(inout):: spl,dspl

    integer:: i,ip
    real(8):: a(4),rij2

    do i=2,npnts
      if( rij.lt.pnts(i) ) exit
    enddo
    ip = i -1
    a(:) = coefs(:,ip)
    rij2 = rij*rij
    spl = a(1) +a(2)*rij +a(3)*rij2 +a(4)*rij*rij2
    dspl = a(2) +2d0*a(3)*rij +3d0*a(4)*rij2
    return
  end subroutine spl2d
!=======================================================================
  subroutine spl3d(csn,npnts,pnts,coefs,spl,dspl)
    real(8),intent(in):: csn
    integer,intent(in):: npnts
    real(8),intent(in):: pnts(npnts),coefs(4,npnts-1)
    real(8),intent(inout):: spl,dspl

    integer:: i,ip
    real(8):: a(4),csn2

    do i=2,npnts
      if( csn.lt.pnts(i) ) exit
    enddo
    ip = i -1
    a(:) = coefs(:,ip)
    csn2 = csn*csn
    spl = a(1) +a(2)*csn +a(3)*csn2 +a(4)*csn*csn2
    dspl = a(2) +2d0*a(3)*csn +3d0*a(4)*csn2
    return
  end subroutine spl3d
!=======================================================================
  function fc1(r,rin,rout)
    implicit none
    real(8),intent(in):: r,rin,rout
    real(8):: fc1

    if( r.le.rin ) then
      fc1= 1d0
    else if( r.gt.rin .and. r.le.rout ) then
      fc1= 0.5d0 *(cos((r-rin)/(rout-rin)*pi)+1d0)
    else
      fc1= 0d0
    endif
    return
  end function fc1
!=======================================================================
  function dfc1(r,rin,rout)
    implicit none
    real(8),intent(in):: r,rin,rout
    real(8):: dfc1

    if( r.le.rin ) then
      Dfc1= 0d0
    else if( r.gt.rin .and. r.le.rout ) then
      dfc1= -0.5d0*pi/(rout-rin) *sin((r-rin)/(rout-rin)*pi)
    else
      dfc1= 0d0
    endif
    return
  end function dfc1
!=======================================================================
  subroutine set_paramsdir_cspline(dname)
!
!  Accessor routine to set paramsdir.
!
    implicit none
    character(len=*),intent(in):: dname

    paramsdir = trim(dname)
    return
  end subroutine set_paramsdir_cspline
!=======================================================================
  subroutine read_params_cspline(myid,mpi_world,iprint)
!
!  Read knot positions and potential values of each knot.
!  Coefficients are computed from these values.
!
    include 'mpif.h'
    integer,intent(in):: myid,mpi_world,iprint

    integer:: i,ierr,isp,jsp,ksp,ndat,npnts,ispl
    real(8):: rcut
    logical:: lexist
    character(len=128):: fname,ctmp,ctmp2,ctype
    integer,external:: num_data
    
    
    if( myid.eq.0 ) then
      fname = trim(paramsdir)//'/'//trim(cpfname)
      inquire(file=trim(fname),exist=lexist)
      if( .not. lexist ) then
        write(6,'(a)') ' [Error] '//trim(cpfname)//' does not exist !!!.'
        write(6,'(a)') '   The cspline potential needs '//trim(cpfname)//'.'
!!$        call mpi_finalize(ierr)
        stop
      endif

      open(ionum,file=trim(fname),status='old')
1     read(ionum,'(a)',end=10) ctmp
      if( ctmp(1:1).eq.'!' .or. ctmp(1:1).eq.'#' ) goto 1 ! skip comment line
      ndat = num_data(trim(ctmp),' ')
      if( ndat.lt.1 ) goto 1 ! skip blank line
      read(ctmp,*) nspl
      if( allocated(spls)) deallocate(spls)
      allocate(spls(nspl))
      ispl = 0
      do while(.true.)
        read(ionum,'(a)',end=10) ctmp
        if( ctmp(1:1).eq.'!' .or. ctmp(1:1).eq.'#' ) cycle ! skip comment line
        ndat = num_data(trim(ctmp),' ')
        if( ndat.lt.1 ) cycle  ! skip blank line
        if( ndat.gt.2 ) then  ! 1st line of a spline data
          backspace(ionum)
          read(ionum,*) ctmp2
          if( trim(ctmp2).eq.'radial' ) then  ! radial spline data
            ispl = ispl +1
!!$            if( ispl.gt.1 ) stop 1
            if( ispl.gt.nspl ) then
              print *,'ERROR: ispl > nspl'
              stop 1
            endif
            backspace(ionum)
            read(ionum,*) ctype, isp, jsp, rcut, npnts
            spls(ispl)%ctype = ctype
            spls(ispl)%isp = isp
            spls(ispl)%jsp = jsp
            spls(ispl)%rcut = rcut
            spls(ispl)%npnts = npnts
            call alloc_spline(spls(ispl))
            do i=1,npnts
              read(ionum,*,err=99) spls(ispl)%pnts(i), spls(ispl)%vals(i)
            enddo
!.....Natural(left) and clamped(right) BCs for radial
!!$            call comp_coefs(npnts,spls(ispl)%pnts,spls(ispl)%vals &
!!$                 ,spls(i)%coefs,'natural','clamped')
            call comp_coefs(spls(ispl),'natural','clamped',iprint)
!!$            if( ispl.eq.1 ) then
!!$              print *,' ispl=',ispl
!!$              do i=1,npnts-1
!!$                print *,'i,pnt,a(1:4)=',i,spls(ispl)%pnts(i),spls(ispl)%coefs(1:4,i)
!!$              enddo
!!$            endif
          else if( ctmp2(1:7).eq.'angular' ) then  ! angular spline data
            ispl = ispl +1
            if( ispl.gt.nspl ) then
              print *,'ERROR: ispl > nspl'
              stop 1
            endif
            backspace(ionum)
            read(ionum,*) ctype, isp, jsp, ksp, rcut, npnts
            spls(ispl)%ctype = ctype
            spls(ispl)%isp = isp
            spls(ispl)%jsp = jsp
            spls(ispl)%ksp = ksp
            spls(ispl)%rcut = rcut
            spls(ispl)%npnts = npnts
            call alloc_spline(spls(ispl))
            do i=1,npnts
              read(ionum,*,err=99) spls(ispl)%pnts(i), spls(ispl)%vals(i)
            enddo
!.....Natural and natural BCs for angular with cosine points
!!$            call comp_coefs(npnts,spls(ispl)%pnts,spls(ispl)%vals &
!!$                 ,spls(ispl)%coefs,'natural','natural')
            call comp_coefs(spls(ispl),'natural','natural',iprint)
          endif
        endif
      enddo
10    close(ionum)
    endif

    call bcast_splines(myid,mpi_world,ierr)
!.....Assign pair/triplet to spline-ID
    idpair(:,:) = -1
    idtriplet(:,:,:) = -1
    do i=1,nspl
      isp = spls(i)%isp
      jsp = spls(i)%jsp
      ksp = spls(i)%ksp
      if( ksp .lt. 0 ) then  ! radial
        idpair(isp,jsp) = i
        idpair(jsp,isp) = i
      else  ! angular
        idtriplet(isp,jsp,ksp) = i
        idtriplet(isp,ksp,jsp) = i
      endif
    enddo
    return
    
99  continue
    print *,'ERROR: while reading spline data...'
    stop 1
  end subroutine read_params_cspline
!=======================================================================
  subroutine write_cspline_curves()
!
!  Write out.cspline file that contains spline curve data for each pair.
!
    integer,parameter:: ndpnts = 100
    integer:: is,js,ks,ispl,i
    real(8):: rmin,rmax,dr,ri,phi,dphi,da
    type(spline):: spl
    character(len=128),parameter:: cfrad = 'out.cspline.radial'
    character(len=128),parameter:: cfang = 'out.cspline.angular'
    integer,parameter:: iorad = 51
    integer,parameter:: ioang = 52

    rmin = 1d+10
    rmax = -1d0
    do is=1,nspmax
      do js=is,nspmax
        ispl = idpair(is,js)
        if( ispl.le.0 ) cycle
        spl = spls(ispl)
!!$        print *,'is,js,ispl=',is,js,ispl
!!$        do i=1,spl%npnts
!!$          if( i.eq.spl%npnts ) then
!!$            print *,'i,pnts,vals=',i,spl%pnts(i),spl%vals(i)
!!$          else
!!$            print *,'i,pnts,vals=',i,spl%pnts(i),spl%vals(i),spl%coefs(1:4,i)
!!$          endif
!!$        enddo
!!$        print *,'rmin,rmax=',spl%pnts(1),spl%pnts(spl%npnts)
        rmin = min(rmin,spl%pnts(1))
        rmax = max(rmax,spl%pnts(spl%npnts))
      enddo
    enddo
    print *,'total rmin,rmax=',rmin,rmax

    open(iorad,file=trim(cfrad),status='replace')
    dr = (rmax-rmin)/(ndpnts-1)
    do i=1,ndpnts
      ri = rmin +dr*(i-1)
      write(iorad,'(2x,f7.3)',advance='no') ri
      do is=1,nspmax
        do js=is,nspmax
          ispl = idpair(is,js)
          if( ispl.le.0 ) cycle
          spl = spls(ispl)
          call eval_spline(ri,spl%npnts,spl%pnts,spl%vals,spl%coefs,phi,dphi)
          write(iorad,'(2x,es12.3e3)',advance='no') phi
        enddo
      enddo
      write(iorad,*) ''
    enddo
    close(iorad)

    open(ioang,file=trim(cfang),status='replace')
    da = (1.0-(-1.0))/(ndpnts-1)
    do i=1,ndpnts
      ri = -1d0 +da*(i-1)
      write(ioang,'(2x,f7.3)',advance='no') ri
      do is=1,nspmax
        do js=1,nspmax
          do ks=1,nspmax
            ispl = idtriplet(is,js,ks)
            if( ispl.le.0 ) cycle
!!$            print *,'is,js,ks,ispl=',is,js,ks,ispl
            spl = spls(ispl)
            call eval_spline(ri,spl%npnts,spl%pnts,spl%vals,spl%coefs,phi,dphi)
            write(ioang,'(2x,es12.3e3)',advance='no') phi
          enddo
        enddo
      enddo
      write(ioang,*) ''
    enddo
    close(ioang)
    
  end subroutine write_cspline_curves
!=======================================================================
  subroutine alloc_spline(spl)
    type(spline),intent(inout):: spl

    integer:: np
    
    if( spl%npnts .le. 0 ) return

    if( allocated(spl%pnts) ) return

    np = spl%npnts
    allocate(spl%pnts(np), spl%vals(np), spl%coefs(4,np-1))
    return
  end subroutine alloc_spline
!=======================================================================
  subroutine bcast_splines(myid,mpi_world,iprint)
    include 'mpif.h'
    integer,intent(in):: myid,mpi_world,iprint

    integer:: i,ierr,npnts

    call mpi_bcast(nspl,1,mpi_integer,0,mpi_world,ierr)
    if( myid.ne.0 ) then
      if( allocated(spls) ) deallocate(spls)
      allocate(spls(nspl))
    endif
    do i=1,nspl
      call mpi_bcast(spls(i)%ctype,128,mpi_character,0,mpi_world,ierr)
      call mpi_bcast(spls(i)%isp,1,mpi_integer,0,mpi_world,ierr)
      call mpi_bcast(spls(i)%jsp,1,mpi_integer,0,mpi_world,ierr)
      call mpi_bcast(spls(i)%ksp,1,mpi_integer,0,mpi_world,ierr)
      call mpi_bcast(spls(i)%npnts,1,mpi_integer,0,mpi_world,ierr)
      call mpi_bcast(spls(i)%rcut,1,mpi_real8,0,mpi_world,ierr)
      if( myid.ne.0 ) call alloc_spline(spls(i))
      call mpi_barrier(mpi_world,ierr)
      npnts = spls(i)%npnts
      call mpi_bcast(spls(i)%pnts,npnts,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(spls(i)%vals,npnts,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(spls(i)%coefs,4*(npnts-1),mpi_real8,0,mpi_world,ierr)
    enddo
  end subroutine bcast_splines
!=======================================================================
  subroutine comp_coefs(spl,bcl,bcr,iprint)
!
!  Compute 4 coefficients for each section.
!  Left and right boundary conditions are specified in bcl and bcr.
!  bcl, bcr should either be:
!    - 'natural': 2nd derivative = 0.0
!    - 'clamped': 1st derivative = 0.0
!
    type(spline),intent(inout):: spl
!!$    integer,intent(in):: npnts
!!$    real(8),intent(in):: pnts(npnts),vals(npnts)
!!$    real(8),intent(out):: coefs(4,npnts-1)
    character(len=7),intent(in):: bcl,bcr
    integer,intent(in):: iprint

    integer:: i,j,npnts,n,ndim,ibase
    real(8),allocatable:: vb(:),vx(:),amat(:,:),amati(:,:),dat(:,:)
    real(8),parameter:: eps = 1e-10

    npnts = spl%npnts
    n = npnts -1
    ndim = 4*n
    allocate(dat(2,npnts),vb(ndim),amat(ndim,ndim),vx(ndim)&
         ,amati(ndim,ndim))
!!$    if( iprint.gt.1 ) print *,'Left & right BC for c-spline: ' &
!!$         ,trim(bcl),' ',trim(bcr)
    do i=1,npnts
      dat(1,i) = spl%pnts(i)
      dat(2,i) = spl%vals(i)
!!$      if( iprint.gt.1 ) print *,'i,pnt,val=',i,dat(1:2,i)
    enddo

!.....Create vector b
    vb(1:ndim)= 0d0
    do i=1,n
      vb(2*i-1)= dat(2,i)
      vb(2*i  )= dat(2,i+1)
    enddo

!.....Create matrix A
    amat(:,:) = 0d0
!.....0th order
    do i=1,n
      amat(2*i-1,4*i-3)= 1d0
      amat(2*i-1,4*i-2)= dat(1,i)
      amat(2*i-1,4*i-1)= dat(1,i)**2
      amat(2*i-1,4*i  )= dat(1,i)**3
      amat(2*i  ,4*i-3)= 1d0
      amat(2*i  ,4*i-2)= dat(1,i+1)
      amat(2*i  ,4*i-1)= dat(1,i+1)**2
      amat(2*i  ,4*i  )= dat(1,i+1)**3
    enddo
!.....1st order derivative
    ibase= 2*n
    do i=1,n-1
      amat(ibase+i,4*i-2)= 1d0
      amat(ibase+i,4*i-1)= 2d0*dat(1,i+1)
      amat(ibase+i,4*i  )= 3d0*dat(1,i+1)**2
      amat(ibase+i,4*(i+1)-2)= -1d0
      amat(ibase+i,4*(i+1)-1)= -2d0*dat(1,i+1)
      amat(ibase+i,4*(i+1)  )= -3d0*dat(1,i+1)**2
    enddo
!.....2nd order derivative
    ibase= 2*n +n-1
    do i=1,n-1
      amat(ibase+i,4*i-1)= 2d0
      amat(ibase+i,4*i  )= 6d0*dat(1,i+1)
      amat(ibase+i,4*(i+1)-1)= -2d0
      amat(ibase+i,4*(i+1)  )= -6d0*dat(1,i+1)
    enddo

    ibase= 2*n +2*(n-1)  ! ibase for both left and right edges
!.....Left BC
    if( bcl.eq.'clamped' ) then
      amat(ibase+1,2)= 1d0
      amat(ibase+1,3)= 2d0*dat(1,1)
      amat(ibase+1,4)= 3d0*dat(1,1)**2
    else ! default: 'natural'
      amat(ibase+1,3)= 2d0
      amat(ibase+1,4)= 6d0*dat(1,1)
    endif
!.....Right BC
    if( bcr.eq.'clamped' ) then
      amat(ibase+2,ndim-2)= 1d0
      amat(ibase+2,ndim-1)= 2d0*dat(1,npnts)
      amat(ibase+2,ndim  )= 3d0*dat(1,npnts)**2
    else ! default: 'natural'
      amat(ibase+2,ndim-1)= 2d0
      amat(ibase+2,ndim  )= 6d0*dat(1,npnts)
    endif

!.....Get inverse of matrix A
!!$    amati(1:ndim,1:ndim)= 0d0
!!$    do i=1,ndim
!!$      amat(i,i) = amat(i,i) +eps
!!$    enddo
    call ludc_inv(ndim,amat,amati)

!.....Compute Ai*b to get x
    vx(1:ndim)= 0d0
    do i=1,ndim
      do j=1,ndim
        vx(i)=vx(i) +amati(i,j)*vb(j)
      enddo
    enddo
    do i=1,npnts-1
      spl%coefs(1:4,i) = vx(4*i-3:4*i)
    enddo
    deallocate(vb,amat,amati,vx,dat)
  end subroutine comp_coefs
end module cspline
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
