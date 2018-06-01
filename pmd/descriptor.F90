module descriptor
!=======================================================================
! Descriptor module
!=======================================================================
  implicit none
  save
  character(len=128):: paramsdir = '.'
  
  character(128),parameter:: cpfname = 'in.params.desc'

  real(8),parameter:: pi = 3.14159265358979d0

  integer:: nsf
  integer,allocatable:: itype(:)
  real(8),allocatable:: cnst(:,:),rcs(:),rcs2(:)
!.....symmetry function values and their derivatives
  real(8),allocatable:: gsf(:,:),dgsf(:,:,:,:)
  logical:: lupdate_gsf = .true.
!.....start and end points of symmetry function IDs for each species pair
  integer,allocatable:: iaddr2(:,:,:),iaddr3(:,:,:,:)
!.....symmetry function IDs for each pair
  integer(2),allocatable:: igsf(:,:,:)
!.....function types and num of constatns for types
  integer,parameter:: max_ncnst= 2
  integer:: ncnst_type(200)
  integer:: ncomb_type(200)

!.....cutoff radii for 2- or 3-body
  real(8):: rc2 = 4.0d0
  real(8):: rc3 = 3.0d0
  real(8):: rcw2 = 0.0d0
  real(8):: rcw3 = 0.0d0
  real(8):: rcmax,rcmax2

  integer:: nalmax,nnlmax,nal,nnl

contains
  subroutine init_desc()
!
!  Initialize descriptor module
!
!.....initialize some
    ncnst_type(1) = 2   ! Gaussian
    ncnst_type(2) = 1   ! cosine
    ncnst_type(3) = 1   ! polynomial
    ncnst_type(4) = 2   ! Morse
    ncnst_type(101) = 1 ! angular
    ncomb_type(1:100) = 2    ! pair
    ncomb_type(101:200) = 3  ! triplet

  end subroutine init_desc
!=======================================================================
  subroutine make_gsf_arrays(l1st,namax,natm,tag,nnmax,lspr &
       ,myid,mpi_world,iprint)
!
!  Make or update gsf arrays
!
    include "mpif.h"
    logical,intent(in):: l1st
    integer,intent(in):: namax,natm,nnmax,lspr(0:nnmax,namax),iprint
    integer,intent(in):: myid,mpi_world
    real(8),intent(in):: tag(namax)

    integer:: i,nnltmp,ierr
    logical,save:: lrealloc = .false.

    if( .not. lupdate_gsf ) return

    if( l1st ) then
    
!  To reduce the memory usage, compute num of atoms and num of neighbors,
!  and add some margin for those numbers because they can change during
!  the simulation.
      nal = int(natm*1.1)
      nnltmp = 0
      do i=1,natm
        nnltmp = max(nnltmp,lspr(0,i))
      enddo
      nnl = int(nnltmp*1.1)
      if( nal .gt. namax ) then
        write(6,'(a)') ' [Error] nal .gt.namax'
        write(6,'(a,3i10)') '   myid,nal,namax = ',myid,nal,namax
        stop
      endif
      if( nnl.gt.nnmax ) then
        write(6,'(a)') ' [Error] nnl.gt.nnmax'
        write(6,'(a,3i10)') '   myid,nnl,nnmax = ',myid,nnl,nnmax
        stop
      endif
      if( myid.ge.0 ) then
        call mpi_reduce(nal,nalmax,1,mpi_integer,mpi_max,0,mpi_world,ierr)
        call mpi_reduce(nnl,nnlmax,1,mpi_integer,mpi_max,0,mpi_world,ierr)
      else
        nalmax = nal
        nnlmax = nnl
      endif
      if( myid.le.0 .and. iprint.ne.0 ) then
        write(6,'(a,2i10)') ' max num of (local atoms *1.1) = ',nalmax
        write(6,'(a,2i10)') ' max num of (neighbors *1.1)   = ',nnlmax
        write(6,'(a,f10.3,a)') ' gsf size  = ', &
             dble(nsf*nalmax*8)/1000/1000,' MB'
        write(6,'(a,f10.3,a)') ' dgsf size = ', &
             dble(3*nsf*(nnlmax+1)*nalmax*8)/1000/1000,' MB'
        write(6,'(a,f10.3,a)') ' igsf size = ', &
             dble(nsf*(nnlmax+1)*nalmax*2)/1000/1000,' MB'
      endif
      if( allocated(gsf) ) deallocate(gsf,dgsf,igsf)
      allocate( gsf(nsf,nal),dgsf(3,nsf,0:nnl,nal) &
           ,igsf(nsf,0:nnl,nal) )

      lrealloc = .false.

    endif

!  Since natm and nn can change every step of MD,
!  if natm/nnltmp becomes nal/nnl, they should be updated and
!  gsf/dgsf as well.
    if( natm.gt.nal ) then
      nal = int(natm*1.1)
      if( nal .gt. namax ) then
        write(6,'(a)') ' [Error] nal .gt.namax'
        write(6,'(a,3i10)') '   myid,nal,namax = ',myid,nal,namax
        stop
      else
        write(6,*) ' Since natm.gt.nal, nal was updated at myid =',myid
      endif
      lrealloc=.true.
    endif
    nnltmp = 0
    do i=1,natm
      nnltmp = max(nnltmp,lspr(0,i))
    enddo
    if( nnltmp.gt.nnl ) then
      nnl = int(nnltmp*1.1)
      if( nnlmax.gt.nnmax ) then
        write(6,'(a)') ' [Error] nnl.gt.nnmax'
        write(6,'(a,3i10)') '   myid,nnl,nnmax = ',myid,nnl,nnmax
        stop
      else
        write(6,*) ' Since nnltmp.gt.nnl, nnl was updated at myid =',myid
      endif
      lrealloc=.true.
    endif
    if( allocated(dgsf).and.lrealloc ) then
      deallocate( gsf,dgsf,igsf )
      allocate( gsf(nsf,nal),dgsf(3,nsf,0:nnl,nal) &
           ,igsf(nsf,0:nnl,nal))
      lrealloc=.false.
    endif
    
    return
  end subroutine make_gsf_arrays
!=======================================================================
  subroutine calc_desc(namax,natm,nb,nnmax,h,tag,ra,lspr,rc &
       ,myid,mpi_world,l1st,iprint)
!
!  Evaluate descriptors (symmetry functions)
!  and their derivatives wrt positions for multi-species system.
!
!  Cutoff radii from outside are written as rc2o and rc3o,
!  which could be different from rcs given from params file.
!  Smaller values would be adopted.
!
    implicit none
    include "mpif.h"
    integer,intent(in):: namax,natm,nb,nnmax,lspr(0:nnmax,namax)
    integer,intent(in):: myid,mpi_world,iprint
    real(8),intent(in):: h(3,3),tag(namax),ra(3,namax),rc
    logical,intent(in):: l1st 

    integer:: isf,ia,jj,ja,kk,ka,is,js,ks,ierr,i
    real(8):: xi(3),xj(3),xij(3),rij(3),dij,dij2,fcij,eta,rs,texp,driji(3), &
         dfcij,drijj(3),dgdr,xk(3),xik(3),rik(3),dik,fcik,dfcik, &
         driki(3),drikk(3),almbd,spijk,cs,t1,t2,dgdij,dgdik,dgcs, &
         dcsdj(3),dcsdk(3),dcsdi(3),tcos,tpoly,a1,a2,tmorse,dik2,tmp
!!$    real(8),save:: rc22,rc32,rcs2,rcs3,eta3
    
    real(8):: texpij,texpik,eta3

    if( .not.lupdate_gsf ) return

    if( l1st ) then
!.....Check all the rcs and compare them with rc
      if( rc.lt.rcmax ) then
!!$      if( rc.lt.rc2 .or. rc.lt.rc3 ) then
        if( myid.eq.0 ) then
          print *,'ERROR: cutoff radius rc is smaller than rcmax.'
          print *,'  rc,rcmax = ',rc,rcmax
        endif
        call mpi_finalize(ierr)
        stop
      endif
!!$!.....Define squares
!!$      rc22 = rc2*rc2
!!$      rc32 = rc3*rc3
!!$      rcs2 = rc2*rcw2
!!$      rcs3 = rc3*rcw3
!!$      eta3 = 0.5d0 /(rc3/2)**2
    endif

    gsf(1:nsf,1:nal)= 0d0
    dgsf(1:3,1:nsf,0:nnl,1:nal)= 0d0
    igsf(1:nsf,0:nnl,1:nal) = 0
    do ia=1,natm
      xi(1:3)= ra(1:3,ia)
      is= int(tag(ia))
      do jj=1,lspr(0,ia)
        ja= lspr(jj,ia)
        if( ja.eq.ia ) cycle
        xj(1:3)= ra(1:3,ja)
        xij(1:3)= xj(1:3)-xi(1:3)
        rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij2= rij(1)**2 +rij(2)**2 +rij(3)**2
        if( dij2.ge.rcmax2 ) cycle
        dij = sqrt(dij2)
        js= int(tag(ja))
        driji(1:3)= -rij(1:3)/dij
        drijj(1:3)= -driji(1:3)
!!$        fcij= fc0(dij,rc2)
!!$        dfcij= dfc0(dij,rc2)
        do isf=iaddr2(1,is,js),iaddr2(2,is,js)
          if( dij.ge.rcs(isf) ) cycle
          if( itype(isf).eq.1 ) then ! Gaussian
            fcij= fc1(dij,rcs(isf))
            dfcij= dfc1(dij,rcs(isf))
            eta= cnst(1,isf)
            rs=  cnst(2,isf)
            !.....function value
            texp= exp(-eta*(dij-rs)**2)
            gsf(isf,ia)= gsf(isf,ia) +texp*fcij
            !.....derivative
            ! dgsf(ixyz,isf,jj,ia): derivative of isf-th basis of atom-ia
            ! by ixyz coordinate of atom-jj.
            ! jj=0 means derivative by atom-ia.
            dgdr= -2d0*eta*(dij-rs)*texp*fcij +texp*dfcij
            dgsf(1:3,isf,0,ia)= dgsf(1:3,isf,0,ia) +driji(1:3)*dgdr
            dgsf(1:3,isf,jj,ia)= dgsf(1:3,isf,jj,ia) +drijj(1:3)*dgdr
            igsf(isf,0,ia) = 1
            igsf(isf,jj,ia) = 1
          else if( itype(isf).eq.2 ) then ! cosine
            fcij= fc1(dij,rcs(isf))
            dfcij= dfc1(dij,rcs(isf))
            a1= cnst(1,isf)
            !.....func value
            tcos= (1d0+cos(dij*a1))
            gsf(isf,ia)= gsf(isf,ia) +tcos*fcij
            !.....derivative
            dgdr= -a1*sin(dij*a1)*fcij +tcos*dfcij
            dgsf(1:3,isf,0,ia)= dgsf(1:3,isf,0,ia) +driji(1:3)*dgdr
            dgsf(1:3,isf,jj,ia)= dgsf(1:3,isf,jj,ia) +drijj(1:3)*dgdr
            igsf(isf,0,ia) = 1
            igsf(isf,jj,ia) = 1
          else if( itype(isf).eq.3 ) then ! polynomial
            fcij= fc1(dij,rcs(isf))
            dfcij= dfc1(dij,rcs(isf))
            a1= cnst(1,isf)
            !.....func value
            tpoly= 1d0*dij**(-a1)
            gsf(isf,ia)= gsf(isf,ia) +tpoly*fcij
            !.....derivative
            dgdr= -a1*dij**(-a1-1d0)*fcij +tpoly*dfcij
            dgsf(1:3,isf,0,ia)= dgsf(1:3,isf,0,ia) +driji(1:3)*dgdr
            dgsf(1:3,isf,jj,ia)= dgsf(1:3,isf,jj,ia) +drijj(1:3)*dgdr
            igsf(isf,0,ia) = 1
            igsf(isf,jj,ia) = 1
          else if( itype(isf).eq.4 ) then ! Morse-type
            fcij= fc1(dij,rcs(isf))
            dfcij= dfc1(dij,rcs(isf))
            a1= cnst(1,isf)
            a2= cnst(2,isf)
            !.....func value
            texp= exp(-a1*(dij-a2))
            tmorse= ((1d0-texp)**2 -1d0)
            gsf(isf,ia)= gsf(isf,ia) +tmorse*fcij
            !.....derivative
            dgdr= 2d0*a1*(1d0-texp)*texp*fcij +tmorse*dfcij
            dgsf(1:3,isf,0,ia)= dgsf(1:3,isf,0,ia) +driji(1:3)*dgdr
            dgsf(1:3,isf,jj,ia)= dgsf(1:3,isf,jj,ia) +drijj(1:3)*dgdr
            igsf(isf,0,ia) = 1
            igsf(isf,jj,ia) = 1
          endif
        enddo

!!$        fcij= fc0(dij,rc3)
!!$        dfcij= dfc0(dij,rc3)
!!$        texpij = exp(-eta3*dij2)
        do kk=1,lspr(0,ia)
          ka= lspr(kk,ia)
          ks= int(tag(ka))
          if( iaddr3(1,is,js,ks).lt.0 ) cycle
          if( ka.eq.ia .or. ka.le.ja ) cycle
          xk(1:3)= ra(1:3,ka)
          xik(1:3)= xk(1:3)-xi(1:3)
          rik(1:3)= h(1:3,1)*xik(1) +h(1:3,2)*xik(2) +h(1:3,3)*xik(3)
          dik2= rik(1)**2 +rik(2)**2 +rik(3)**2
!!$          if( dik2.ge.rc32 ) cycle
          dik= sqrt(dik2)
!!$          fcik= fc0(dik,rc3)
!!$          dfcik= dfc0(dik,rc3)
!!$          texpik= exp(-eta3*dik2)
          do isf=iaddr3(1,is,js,ks),iaddr3(2,is,js,ks)
            if( dij.ge.rcs(isf) .or. dik.ge.rcs(isf) ) cycle
!.....fcij's can be computed after rcs is determined
            fcij= fc1(dij,rcs(isf))
            dfcij= dfc1(dij,rcs(isf))
            fcik= fc1(dik,rcs(isf))
            dfcik= dfc1(dik,rcs(isf))
            almbd= cnst(1,isf)
            t2= (abs(almbd)+1d0)**2
            driki(1:3)= -rik(1:3)/dik
            drikk(1:3)= -driki(1:3)
            !.....function value
            spijk= rij(1)*rik(1) +rij(2)*rik(2) +rij(3)*rik(3)
            cs= spijk/dij/dik
            t1= (almbd +cs)**2
            eta3 = 0.5d0 /(rcs(isf)/2)**2
            texpij = exp(-eta3*dij2)
            texpik = exp(-eta3*dik2)
            tmp = t1/t2 *texpij *texpik
            gsf(isf,ia)= gsf(isf,ia) +tmp*fcij*fcik 
!!$            gsf(isf,ia)= gsf(isf,ia) +t1/t2 *fcij*fcik
            !.....derivative
            dgdij= dfcij *fcik *tmp &
                 +tmp *(-2d0*eta3*dij) *fcij*fcik 
            dgdik= fcij *dfcik *tmp &
                 +tmp *(-2d0*eta3*dik) *fcij*fcik 
!!$            dgdij= dfcij *fcik *t1/t2
!!$            dgdik= fcij *dfcik *t1/t2
            dgsf(1:3,isf,0,ia)= dgsf(1:3,isf,0,ia) &
                 +dgdij*driji(1:3) +dgdik*driki(1:3)
            dgsf(1:3,isf,jj,ia)= dgsf(1:3,isf,jj,ia) +dgdij*drijj(1:3)
            dgsf(1:3,isf,kk,ia)= dgsf(1:3,isf,kk,ia) +dgdik*drikk(1:3)
            dgcs= 2d0*(almbd+cs)/t2 *fcij*fcik *texpij*texpik
            dcsdj(1:3)= rik(1:3)/dij/dik -rij(1:3)*cs/dij**2
            dcsdk(1:3)= rij(1:3)/dij/dik -rik(1:3)*cs/dik**2
            dcsdi(1:3)= -dcsdj(1:3) -dcsdk(1:3)
            dgsf(1:3,isf,0,ia)= dgsf(1:3,isf,0,ia) +dgcs*dcsdi(1:3)
            dgsf(1:3,isf,jj,ia)= dgsf(1:3,isf,jj,ia) +dgcs*dcsdj(1:3)
            dgsf(1:3,isf,kk,ia)= dgsf(1:3,isf,kk,ia) +dgcs*dcsdk(1:3)
            igsf(isf,0,ia) = 1
            igsf(isf,jj,ia) = 1
            igsf(isf,kk,ia) = 1
          enddo
        enddo
20      continue
      enddo
    enddo

  end subroutine calc_desc
!=======================================================================
  function fc1(r,rc)
    implicit none
    real(8),intent(in):: r,rc
    real(8):: fc1
    real(8),parameter:: rcs = 0d0

    if( r.le.rcs ) then
      fc1= 1d0
    else if( r.gt.rcs .and. r.le.rc ) then
      fc1= 0.5d0 *(cos((r-rcs)/(rc-rcs)*pi)+1d0)
    else
      fc1= 0d0
    endif
    return
  end function fc1
!=======================================================================
  function dfc1(r,rc)
    implicit none
    real(8),intent(in):: r,rc
    real(8):: dfc1
    real(8),parameter:: rcs = 0d0

    if( r.le.rcs ) then
      dfc1= 0d0
    else if( r.gt.rcs .and. r.le.rc ) then
      dfc1= -pi/2/(rc-rcs) *sin((r-rcs)/(rc-rcs)*pi)
    else
      dfc1= 0d0
    endif
    return
  end function dfc1
!=======================================================================
  subroutine read_params_desc(myid,mpi_world,iprint)
    implicit none
    include 'mpif.h'

    integer,intent(in):: myid,mpi_world,iprint
!!$    real(8),intent(in):: rcin

    integer:: ierr,i,j,k,nc,ncoeff,nsp &
         ,ihl0,ihl1,ihl2,icmb(3),nsf1,nsf2,iap,jap,kap,ndat
    logical:: lexist
    character:: ctmp*128,fname*128
    integer,external:: num_data

    if( myid.eq.0 ) then

!.....read constants at the 1st call
      fname = trim(paramsdir)//'/'//trim(cpfname)
      inquire(file=trim(fname),exist=lexist)
      if( .not. lexist ) then
        if( myid.eq.0 ) then
          write(6,'(a)') ' [Error] '//trim(fname)//' does not exist !!!.'
        endif
        call mpi_finalize(ierr)
        stop
      endif
      open(51,file=trim(fname),status='old')
!.....num of symmetry functions, num of node in 1st hidden layer
10    read(51,'(a)') ctmp
      if( ctmp(1:1).eq.'!' .or. ctmp(1:1).eq.'#' ) then
!        call parse_option(ctmp,iprint,ierr)
        goto 10
      else
        backspace(51)
      endif
!!$!.....Check rc's given by in.params.desc and by in.pmd
!!$      if( rc2.gt.rcin .or. rc3.gt.rcin ) then
!!$        if( myid.eq.0 ) then
!!$          print *,'ERROR: Cutoff radius in in.pmd shorter than that in in.params.desc'
!!$          print *,'  rcin,rc2,rc3 = ',rcin,rc2,rc3
!!$        endif
!!$        stop
!!$      endif
!.....Read numbers of species and symmetry functions
      read(51,*) nsp, nsf
    endif

!.....Bcast nsp and nsf before allocating arrays
    call mpi_bcast(nsp,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(nsf,1,mpi_integer,0,mpi_world,ierr)

!.....Allocate arrays of lenths, nsp and/or nsf
    if( .not.allocated(itype) ) then
      allocate(itype(nsf),cnst(max_ncnst,nsf),rcs(nsf),rcs2(nsf))
      allocate(iaddr2(2,nsp,nsp),iaddr3(2,nsp,nsp,nsp))
    endif
    
    if( myid.eq.0 ) then
      iaddr2(1:2,1:nsp,1:nsp)= -1
      iaddr3(1:2,1:nsp,1:nsp,1:nsp)= -1
      nsf1= 0
      nsf2= 0
      iap= 0
      jap= 0
      kap= 0
      do i=1,nsf
        read(51,*) itype(i),(icmb(k),k=1,ncomb_type(itype(i))) &
             ,rcs(i),(cnst(j,i),j=1,ncnst_type(itype(i)))
        if( itype(i).le.100 ) then
          if( icmb(1).ne.iap .or. icmb(2).ne.jap ) then
            iaddr2(1,icmb(1),icmb(2))= i
            iaddr2(1,icmb(2),icmb(1))= i
          endif
          iaddr2(2,icmb(1),icmb(2))= i
          iaddr2(2,icmb(2),icmb(1))= i
          nsf1= nsf1 +1
          iap= icmb(1)
          jap= icmb(2)
        else if( itype(i).le.200 ) then
          if( icmb(1).ne.iap .or. icmb(2).ne.jap .or. &
               icmb(3).ne.kap ) then
            iaddr3(1,icmb(1),icmb(2),icmb(3))= i
            iaddr3(1,icmb(1),icmb(3),icmb(2))= i
          endif
          iaddr3(2,icmb(1),icmb(2),icmb(3))= i
          iaddr3(2,icmb(1),icmb(3),icmb(2))= i
          nsf2= nsf2 +1
          iap= icmb(1)
          jap= icmb(2)
          kap= icmb(3)
        endif
      enddo
      if( nsf.ne.nsf1+nsf2 ) then
        print *,'[Error] nsf.ne.nsf1+nsf2 !!!'
!        call mpi_finalize(ierr)
        stop
      endif
      close(51)
    endif

!!$    call mpi_bcast(interact,msp*msp,mpi_logical,0,mpi_world,ierr)
    call mpi_bcast(itype,nsf,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(cnst,max_ncnst*nsf,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(iaddr2,2*nsp*nsp,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(iaddr3,2*nsp*nsp*nsp,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(rcs,nsf,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(rc2,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(rc3,1,mpi_real8,0,mpi_world,ierr)

!.....Compute maximum rcut in all descriptors
    rcmax = 0d0
    do i=1,nsf
      rcmax = max(rcmax,rcs(i))
      rcs2(i) = rcs(i)**2
    enddo
    rcmax2 = rcmax**2
    
    return
  end subroutine read_params_desc
!=======================================================================
  subroutine parse_option(cline,iprint,ierr)
!
!  Parse options from a comment line.
!  Lines starting from ! or # are treated as comment lines,
!  but options can be given at the comment lines.
!  The option words should be put after these comment characters with
!  one or more spaces between them for example,
!
!  rc2:  4.0
!
!  Currently available options are:
!    - "rc2:" or "rc3:" cutoff radii for 2- or 3-body.
!
    implicit none
    character(len=*),intent(in):: cline
    integer,intent(in):: iprint
    integer,intent(out):: ierr

    real(8):: ropt
    character(len=10):: c1,copt
    logical:: lopt
    integer,external:: num_data

    ierr = 0
    if( index(cline,'rc').ne.0 ) then
      read(cline,*) c1,copt,ropt
      if( trim(copt).eq.'rc2:' ) then
        rc2 = ropt
      else if( trim(copt).eq.'rc3:' ) then
        rc3 = ropt
      else
        print *, 'Error: copt is not "rc2:" or "rc3:" !!!'
        ierr = 2
      endif
!!$      rcw2 = 0d0
!!$      rcw3 = 0d0
    endif
    
  end subroutine parse_option
!=======================================================================
  subroutine set_paramsdir_desc(dname)
!
!  Accessor routine to set paramsdir.
!
    implicit none
    character(len=*),intent(in):: dname

    paramsdir = trim(dname)
    return
  end subroutine set_paramsdir_desc
!=======================================================================
  subroutine write_descs(ionum,natm,namax,nnmax,lspr,tag)
!
!   Write out descriptor data (gsf,dgsf,igsf).
!   Buffer atom indices are replaced to resident atom ones.
!
    implicit none
    integer,intent(in):: ionum
    integer,intent(in):: natm,namax,nnmax,lspr(0:nnmax,namax)
    real(8),intent(in):: tag(namax)

    integer:: ia,jj,ja,jra,isf,ihl0
    real(8),allocatable:: dgsfo(:,:,:,:)
    integer(2),allocatable:: igsfo(:,:,:)
    integer,external:: itotOf

    open(ionum,file='out.desc.gsf',status='replace',form='unformatted')
    write(ionum) nsf
    do ia=1,natm
      write(ionum) (gsf(ihl0,ia),ihl0=1,nsf)
    enddo
    close(ionum)

    allocate(dgsfo(3,natm,nsf,natm),igsfo(nsf,natm,natm))
!.....reduce d(i)gsf data of buffer atoms to those of resident atoms
    dgsfo(1:3,1:natm,1:nsf,1:natm)= 0d0
    igsfo(1:nsf,1:natm,1:natm)= 0
    do ia=1,natm
      do jj=0,lspr(0,ia)
        if( jj.eq.0 ) then
          ja= ia
        else
          ja= lspr(jj,ia)
        endif
        jra= itotOf(tag(ja))
        do isf=1,nsf
          dgsfo(1:3,jra,isf,ia)= dgsfo(1:3,jra,isf,ia) &
               +dgsf(1:3,isf,jj,ia)
          igsfo(isf,jra,ia) = igsf(isf,jj,ia)
        enddo
      enddo
    enddo
!.....write
    open(ionum+1,file='out.desc.dgsf',status='replace',form='unformatted')
    do ia=1,natm
      do isf=1,nsf
        write(ionum+1) (dgsfo(1:3,jra,isf,ia),jra=1,natm)
      enddo
    enddo
    close(ionum+1)
!.....write igsf to ionum+1
    open(ionum+2,file='out.desc.igsf',status='replace',form='unformatted')
    do ia=1,natm
      do ja=1,natm
        write(ionum+2) (igsfo(isf,ja,ia),isf=1,nsf)
      enddo
    enddo
    close(ionum+2)

    deallocate(dgsfo,igsfo)
  end subroutine write_descs
!=======================================================================
  subroutine get_ints(nsfo,nalo,nnlo)
!
!  Access to some integers.
!
    integer,intent(out):: nsfo,nalo,nnlo

    nsfo = nsf
    nalo = nal
    nnlo = nnl
    return
  end subroutine get_ints
!=======================================================================
  subroutine get_descs(nsfo,nalo,nnlo,gsfo,dgsfo,igsfo)
!
!  Access to descriptors from outside
!
    integer,intent(in):: nsfo,nalo,nnlo
    real(8),intent(out):: gsfo(nsfo,nalo),dgsfo(3,nsfo,0:nnlo,nalo)&
         ,igsfo(nsfo,0:nnlo,nalo)

    gsfo(:,:) = gsf(:,:)
    dgsfo(:,:,:,:) = dgsf(:,:,:,:)
    igsfo(:,:,:) = igsf(:,:,:)
    return
  end subroutine get_descs
!=======================================================================
  subroutine set_descs(nsfo,nalo,nnlo,gsfo,dgsfo,igsfo)
!
!  Set descriptors from outside, which is called from fitpot.
!
    integer,intent(in):: nsfo,nalo,nnlo
    real(8),intent(in):: gsfo(nsfo,nalo),dgsfo(3,nsfo,0:nnlo,nalo)&
         ,igsfo(nsfo,0:nnlo,nalo)

    integer:: isf

    if( nsf.ne.nsfo .or. nal.ne.nalo .or. nnl.ne.nnlo )  then
!!$      print *,'ERROR: nsf or nal or nnl is different.'
!!$      print *,'nsf,nsfo=',nsf,nsfo
!!$      print *,'nal,nalo=',nal,nalo
!!$      print *,'nnl,nnlo=',nnl,nnlo
      if( allocated(gsf) ) deallocate(gsf,dgsf,igsf)
      nsf = nsfo
      nal = nalo
      nnl = nnlo
      allocate( gsf(nsf,nal),dgsf(3,nsf,0:nnl,nal) &
           ,igsf(nsf,0:nnl,nal) )
!!$      print *,'WARNING: gsfs are reallocated...'
    endif
    gsf(:,:) = gsfo(:,:)
    dgsf(:,:,:,:) = dgsfo(:,:,:,:)
    igsf(:,:,:) = igsfo(:,:,:)
!!$    print *,'gsf(:,1) after set:'
!!$    do isf=1,nsf
!!$      print *,'isf,gsf(isf,1)=',isf,gsf(isf,1)
!!$    enddo
    return
  end subroutine set_descs
end module descriptor
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
