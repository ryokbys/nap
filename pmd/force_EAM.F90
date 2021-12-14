module EAM
!-----------------------------------------------------------------------
!                     Last modified: <2021-11-27 13:53:59 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
!  Parallel implementation of the EAM pontential.
!-----------------------------------------------------------------------
  use pmdvars, only: nspmax
  use util,only: csp2isp
  implicit none
  include "./const.h"
  save
  character(len=128):: paramsdir = '.'
  character(len=128),parameter:: paramsfname = 'in.params.EAM'

  real(8),external:: fcut1,dfcut1

  integer,parameter:: ioprms = 20
  
!.....Number of species
  integer:: nsp = 0

!!$.....Default parameter set for Al
!!$  real(8):: ea_a  = 0.763905d0
!!$  real(8):: ea_b  = 0.075016d0
!!$  real(8):: ea_c  = 0.159472d0
!!$  real(8):: ea_re = 3.389513d-10 /ang
!!$  real(8):: ea_al = 1.755162d+10 *ang
!!$  real(8):: ea_bt = 2.003449d+10 *ang
!!$  real(8):: ea_xi = 0.147699d0
!!$  real(8):: am_al = 26.9815d0

  real(8),allocatable:: ea_a(:),ea_xi(:)
  real(8),allocatable,dimension(:,:):: ea_b,ea_c,ea_re, &
       ea_alp,ea_beta,ea_rcin,ea_rcout
  logical,allocatable:: pair_interact(:,:),ea_interact(:)

  logical:: lprmset_EAM = .false.

!.....Parameters from fitpot
  integer:: nprms
  real(8),allocatable:: prms(:)

!.....Types of forms of potential terms
  character(len=128),allocatable:: type_rho(:,:),type_frho(:),type_phi(:,:)
  character(len=128):: default_type_rho = 'exp1'
  character(len=128):: default_type_frho = 'sqrt1'
  character(len=128):: default_type_phi = 'SM'

contains
  subroutine init_EAM()
!
!  Allocate and initialize parameters to be used.
!
    if( .not.allocated(ea_a) ) then
      allocate(ea_a(nspmax), ea_b(nspmax,nspmax), ea_c(nspmax,nspmax), &
           ea_re(nspmax,nspmax), ea_alp(nspmax,nspmax), ea_beta(nspmax,nspmax),&
           ea_xi(nspmax), pair_interact(nspmax,nspmax), ea_interact(nspmax), &
           ea_rcin(nspmax,nspmax), ea_rcout(nspmax,nspmax), &
           type_rho(nspmax,nspmax),type_frho(nspmax),type_phi(nspmax,nspmax))
    endif
    
  end subroutine init_EAM
!=======================================================================
  subroutine read_params_EAM(myid_md,mpi_md_world,iprint,specorder)
!
!  Read parameters from file
!
    include 'mpif.h'
    include './const.h'
    integer,intent(in):: myid_md,mpi_md_world,iprint
    character(len=3),intent(in):: specorder(nspmax)

    integer:: isp,jsp,ierr
    real(8):: a,b,c,re,alp,beta,xi,rc,rcin,rcout
    character(len=128):: cline,fname,c1,c2,c3
    character(len=3):: cspi,cspj

    if( myid_md.eq.0 ) then
      ea_a(:) = 0d0
      ea_xi(:) = 0d0
      ea_b(:,:) = 0d0
      ea_c(:,:) = 0d0
      ea_re(:,:) = 0d0
      ea_alp(:,:) = 0d0
      ea_beta(:,:) = 0d0
      ea_rcin(:,:) = 0d0
      ea_rcout(:,:) = 0d0
      ea_interact(:) = .false.
      pair_interact(:,:) = .false.
      fname = trim(paramsdir)//'/'//trim(paramsfname)
      open(ioprms,file=trim(fname),status='old')
      do while(.true.)
        read(ioprms,'(a)',end=10) cline
!.....Check comment lines
        if( cline(1:1).eq.'#' .or. cline(1:1).eq.'!' ) cycle
        if( len_trim(cline).eq.0 ) cycle
        backspace(ioprms)
!.....The 1st entry specifies ATOMIC or PAIR parameters
        read(ioprms,*) c1
        if( trim(c1).eq.'atomic' ) then
          backspace(ioprms)
          read(ioprms,*) c1,cspi,c2,a,xi,rcin,rcout
          isp = csp2isp(trim(cspi))
          if( isp.gt.nspmax ) then
            print *,'Warning @read_params_EAM: isp is greater than nspmax,' &
                 //' skip reading the line'
            cycle
          endif
          if( isp.le.0 ) cycle
          ea_interact(isp) = .true.
          ea_a(isp) = a
          ea_xi(isp) = xi
          ea_rcin(isp,isp) = rcin
          ea_rcout(isp,isp) = rcout
          type_frho(isp) = trim(c2)
!.....Otherwise the entry is for pair parameter
        else if( trim(c1).eq.'pair' ) then
          backspace(ioprms)
          read(ioprms,*) c1,cspi,cspj,c2,c3,b,c,re,alp,beta
          isp = csp2isp(trim(cspi))
          jsp = csp2isp(trim(cspj))
          if( isp.gt.nspmax .or. jsp.gt.nspmax ) then
            print *,'Warning @read_params_EAM: isp/jsp is greater than nspmax,'&
                 //' skip reading the line.'
            cycle
          endif
          if( isp.le.0 .or. jsp.le.0 ) cycle
          type_rho(isp,jsp) = trim(c2)
          type_phi(isp,jsp) = trim(c3)
          ea_b(isp,jsp) = b
          ea_c(isp,jsp) = c
          ea_re(isp,jsp) = re
          ea_alp(isp,jsp) = alp
          ea_beta(isp,jsp) = beta
          pair_interact(isp,jsp) = .true.
!.....Symmetrize parameters
          type_rho(jsp,isp) = trim(c2)
          type_phi(jsp,isp) = trim(c3)
          ea_b(jsp,isp) = b
          ea_c(jsp,isp) = c
          ea_re(jsp,isp) = re
          ea_alp(jsp,isp) = alp
          ea_beta(jsp,isp) = beta
          pair_interact(jsp,isp) = .true.
        endif
      enddo
10    continue
      close(ioprms)

!.....Define rcin and rcout for pairs as averages bewteen species
      do isp=1,nspmax
        do jsp=1,nspmax
          ea_rcin(isp,jsp) = (ea_rcin(isp,isp)+ea_rcin(jsp,jsp))/2
          ea_rcin(jsp,isp) = ea_rcin(isp,jsp)
          ea_rcout(isp,jsp) = (ea_rcout(isp,isp)+ea_rcout(jsp,jsp))/2
          ea_rcout(jsp,isp) = ea_rcout(isp,jsp)
        enddo
      enddo
      
      if( iprint.ge.ipl_basic ) then
        print *,''
        write(6,'(a)') ' EAM parameters read from file '//trim(fname) &
             //':'
        do isp=1,nspmax
          if( ea_interact(isp) ) then
            cspi = trim(specorder(isp))
            write(6,'(a8,2(2x,a),4f8.3)') 'atomic',trim(cspi),trim(type_frho(isp)), &
                 ea_a(isp),ea_xi(isp), ea_rcin(isp,isp), ea_rcout(isp,isp)
          endif
        enddo
        do isp=1,nspmax
          cspi = trim(specorder(isp))
          do jsp=isp,nspmax
            if( pair_interact(isp,jsp) ) then
              cspj = trim(specorder(jsp))
              write(6,'(a8,4(2x,a),5f8.3)') 'pair',trim(cspi),trim(cspj), &
                   trim(type_rho(isp,jsp)), trim(type_phi(isp,jsp)), &
                   ea_b(isp,jsp),ea_c(isp,jsp), &
                   ea_re(isp,jsp),ea_alp(isp,jsp),ea_beta(isp,jsp)
            endif
          enddo
        enddo
        print *,''
      endif
    endif  ! myid_md == 0

    call mpi_bcast(ea_a,nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(ea_xi,nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(ea_interact,nspmax,mpi_logical,0,mpi_md_world,ierr)
    call mpi_bcast(type_frho,128*nspmax,mpi_character,0,mpi_md_world,ierr)

    call mpi_bcast(ea_b,nspmax*nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(ea_c,nspmax*nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(ea_re,nspmax*nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(ea_alp,nspmax*nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(ea_beta,nspmax*nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(ea_rcin,nspmax*nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(ea_rcout,nspmax*nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(pair_interact,nspmax*nspmax,mpi_logical,0,mpi_md_world,ierr)
    call mpi_bcast(type_rho,128*nspmax*nspmax,mpi_character,0,mpi_md_world,ierr)
    call mpi_bcast(type_phi,128*nspmax*nspmax,mpi_character,0,mpi_md_world,ierr)
    return
  end subroutine read_params_EAM
!=======================================================================
  subroutine force_EAM(namax,natm,tag,ra,nnmax,aa,strs,h,hi &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr,d2lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,lstrs,iprint,l1st)
!-----------------------------------------------------------------------
!  Parallel implementation of EAM poetntial of SM paper.
!    - smoothing is applied to both 2- and many-body terms
!    - rho of boundary atoms are sent to the neighbor nodes
!    - only force on i is calculated, not necessary to send-back
!-----------------------------------------------------------------------
!  See only EAM part of Streitz and Mintmire, PRB 50(16), 11996 (1994)
!-----------------------------------------------------------------------
    implicit none
    include "mpif.h"
    include "./params_unit.h"
    integer,intent(in):: namax,natm,nnmax,nismax,lspr(0:nnmax,namax)&
         ,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_md_world,myid_md,nex(3)
    real(8),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),sv(3,6) &
         ,rc,tag(namax),d2lspr(nnmax,namax)
    real(8),intent(inout):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: l1st
    logical:: lstrs

    integer:: i,j,k,l,m,n,ierr,is,js,ixyz,jxyz
    real(8):: xij(3),rij,rcin,rcout,rc2,dfi,dfj,drho,drdxi(3),drdxj(3),r,at(3)
    real(8):: x,y,z,xi(3),epotl,epott,tmp,dtmp
    real(8),allocatable,save:: rho(:)
    real(8),allocatable,save:: strsl(:,:,:)

    if( l1st ) then
      if( allocated(rho) ) deallocate(rho)
      allocate(rho(namax+nbmax))
      if( allocated(strsl) ) deallocate(strsl)
      allocate(strsl(3,3,namax))
!.....Check cutoff radius
      do is=1,nspmax
        do js=is,nspmax
          if( .not.pair_interact(is,js) ) cycle
          if( rc.lt.ea_rcout(is,js) ) then
            if( myid_md.eq.0 ) then
              write(6,*) ' Error: rc is smaller than one of EAM rcouts'
            endif
            call mpi_finalize(ierr)
            stop
          endif
        enddo
      enddo
    endif

    if( size(strsl).lt.3*3*namax ) then
      deallocate(strsl)
      allocate(strsl(3,3,namax))
    endif
    if( size(rho).lt.namax+nbmax ) then
      deallocate(rho)
      allocate(rho(namax+nbmax))
    endif
    
    epotl= 0d0
    rho(1:namax)= 0d0
    strsl(1:3,1:3,1:namax) = 0d0

!-----rho(i)
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is = int(tag(i))
      if( .not. ea_interact(is) ) cycle
      do k=1,lspr(0,i)
        j=lspr(k,i)
        if(j.eq.0) exit
        js = int(tag(j))
        rcout= ea_rcout(is,js)
        rc2 = rcout*rcout
        if( d2lspr(k,i).gt.rc2 ) cycle
        rcin = ea_rcin(is,js)
        x= ra(1,j) -xi(1)
        y= ra(2,j) -xi(2)
        z= ra(3,j) -xi(3)
        xij(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
        rij=sqrt(xij(1)*xij(1)+ xij(2)*xij(2) +xij(3)*xij(3))
        rho(i) = rho(i) +rhoij(is,js,rij,rcin,rcout,type_rho(is,js))
      enddo
    enddo

!-----copy rho of boundary atoms
    call copy_dba_fwd(namax,natm,nb,nbmax,lsb,nex,&
         lsrc,myparity,nn,sv,mpi_md_world,rho,1)

!-----dE/dr_i
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is = int(tag(i))
      if( ea_interact(is) ) dfi = dfrho(is,rho(i),type_frho(is))
      do k=1,lspr(0,i)
        j=lspr(k,i)
        if(j.eq.0) exit
        if(j.le.i) cycle
        js = int(tag(j))
        rcout = ea_rcout(is,js)
        rc2 = rcout*rcout
        if( d2lspr(k,i).ge.rc2 ) cycle
        rcin = ea_rcin(is,js)
        x= ra(1,j) -xi(1)
        y= ra(2,j) -xi(2)
        z= ra(3,j) -xi(3)
        xij(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
        rij=sqrt(xij(1)**2+ xij(2)**2 +xij(3)**2)
        drdxi(1:3)= -xij(1:3)/rij
        tmp = 0.5d0 *phi(is,js,rij,rcin,rcout,type_phi(is,js))
        epi(i)= epi(i) +tmp
        epi(j)= epi(j) +tmp
        if(j.le.natm) then
          epotl=epotl +tmp +tmp
        else
          epotl=epotl +tmp
        endif
        dtmp = dphi(is,js,rij,rcin,rcout,type_phi(is,js))
        aa(1:3,i)=aa(1:3,i) -dtmp*drdxi(1:3)
        aa(1:3,j)=aa(1:3,j) +dtmp*drdxi(1:3)
!.....Atomic stress for 2-body terms
        if( lstrs ) then
          do ixyz=1,3
            do jxyz=1,3
              strsl(jxyz,ixyz,i)=strsl(jxyz,ixyz,i) &
                   -0.5d0*dtmp*xij(ixyz)*(-drdxi(jxyz))
              strsl(jxyz,ixyz,j)=strsl(jxyz,ixyz,j) &
                   -0.5d0*dtmp*xij(ixyz)*(-drdxi(jxyz))
            enddo
          enddo
        endif
!.....Embedded term
        if( ea_interact(js) ) dfj = dfrho(js,rho(j),type_frho(js))
        if( .not. (ea_interact(is).or.ea_interact(js)) ) cycle
        drho = drhoij(is,js,rij,rcin,rcout,type_rho(is,js))
        tmp = (dfi+dfj)*drho
        aa(1:3,i)=aa(1:3,i) -tmp*drdxi(1:3)
        aa(1:3,j)=aa(1:3,j) +tmp*drdxi(1:3)
!.....Atomic stress of many-body contributions
        if( lstrs ) then
          do ixyz=1,3
            do jxyz=1,3
              strsl(jxyz,ixyz,i)=strsl(jxyz,ixyz,i) &
                   -0.5d0*tmp*xij(ixyz)*(-drdxi(jxyz))
              strsl(jxyz,ixyz,j)=strsl(jxyz,ixyz,j) &
                   -0.5d0*tmp*xij(ixyz)*(-drdxi(jxyz))
            enddo
          enddo
        endif
      enddo
      if( .not. ea_interact(is) ) cycle
      tmp = frho(is,rho(i),type_frho(is))
      epi(i)=epi(i) +tmp
      epotl=epotl +tmp
    enddo

    if( lstrs ) then
!!$      call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
!!$           ,nn,mpi_md_world,strsl,9)
      strs(1:3,1:3,1:natm)= strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)
    endif

!-----gather epot
    call mpi_allreduce(epotl,epott,1,mpi_real8,mpi_sum,mpi_md_world,ierr)
    epot= epot +epott
    if( myid_md.eq.0 .and. iprint.ge.ipl_info ) &
         write(6,'(a,es15.7)') ' epot EAM = ',epott
    return
  end subroutine force_EAM
!=======================================================================
  function rhoij(is,js,rij,rcin,rcout,ctype)
!
! Calculate rho_ij from rij and rcij.
! The type of rho form is given by TYPE_RHO global variable.
!
    real(8),intent(in):: rij,rcin,rcout
    integer,intent(in):: is,js
    character(len=*),intent(in):: ctype
    real(8):: rhoij
    real(8),external:: fcut1

    rhoij = 0d0
    if( trim(ctype).eq.'exp1' ) then
      rhoij = ea_xi(is)*exp(-ea_beta(is,js)*(rij-ea_re(is,js))) &
           *fcut1(rij,rcin,rcout)
    endif
    return
  end function rhoij
!=======================================================================
  function drhoij(is,js,rij,rcin,rcout,ctype)
    integer,intent(in):: is,js
    real(8),intent(in):: rij,rcin,rcout
    character(len=*),intent(in):: ctype 
    real(8):: drhoij
    real(8):: r

    drhoij = 0d0
    if( trim(ctype).eq.'exp1' ) then
      r = rij -ea_re(is,js)
      drhoij= -ea_xi(is)*ea_beta(is,js)*exp(-ea_beta(is,js)*r)&
           *fcut1(rij,rcin,rcout) &
           +ea_xi(is)*exp(-ea_beta(is,js)*r)*dfcut1(rij,rcin,rcout)
    endif
    return
  end function drhoij
!=======================================================================
  function frho(is,rho,ctype)
!
! Calc F[rho] with given rho
! The type of F[.] form is given by TYPE_FRHO global variable.
!
    integer,intent(in):: is
    real(8),intent(in):: rho
    character(len=*),intent(in):: ctype 
    
    real(8):: frho

    frho = 0d0
    if( trim(ctype).eq.'sqrt1' ) then
      frho = -ea_a(is)*sqrt(rho/ea_xi(is))
    endif
    return
  end function frho
!=======================================================================
  function dfrho(is,rho,ctype)
    integer,intent(in):: is
    real(8),intent(in):: rho
    character(len=*),intent(in):: ctype 
    real(8):: dfrho

    dfrho = 0d0
    if( trim(ctype).eq.'sqrt1') then
      dfrho = -ea_a(is) /ea_xi(is) *0.5d0 /sqrt(rho/ea_xi(is))
    endif
    return
  end function dfrho
!=======================================================================
  function phi(is,js,rij,rcin,rcout,ctype)
    integer,intent(in):: is,js
    real(8),intent(in):: rij,rcin,rcout
    character(len=*),intent(in):: ctype 
    real(8):: phi
    real(8):: r

    phi = 0d0
    if( trim(ctype).eq.'SM' ) then
      r = rij -ea_re(is,js)
      phi = 2d0*ea_b(is,js)*exp(-0.5d0*ea_beta(is,js)*r) &
           -ea_c(is,js)*(1d0+ea_alp(is,js)*r)*exp(-ea_alp(is,js)*r)
      phi = phi*fcut1(rij,rcin,rcout)
!!$    else if( trim(ctype).eq.'Bonny' ) then
!!$      
    endif
    return
  end function phi
!=======================================================================
  function dphi(is,js,rij,rcin,rcout,ctype)
    integer,intent(in):: is,js
    real(8),intent(in):: rij,rcin,rcout
    character(len=*),intent(in):: ctype 
    real(8):: dphi,tmp
    real(8):: r

    dphi = 0d0
    if( trim(ctype).eq.'SM' ) then
      r = rij -ea_re(is,js)
!!$      dphi = 2d0*ea_b(is,js)*exp(-0.5d0*ea_beta(is,js)*r) &
!!$           -ea_c(is,js)*(1d0+ea_alp(is,js)*r)*exp(-ea_alp(is,js)*r)
      tmp = 2d0*ea_b(is,js)*exp(-0.5d0*ea_beta(is,js)*r) &
           -ea_c(is,js)*(1d0+ea_alp(is,js)*r)*exp(-ea_alp(is,js)*r)
      dphi= -ea_beta(is,js)*ea_b(is,js) &
           *exp(-0.5d0*ea_beta(is,js)*r)*fcut1(rij,rcin,rcout)  &
           + ea_c(is,js)*ea_alp(is,js)*ea_alp(is,js)*r &
           *exp(-ea_alp(is,js)*r)*fcut1(rij,rcin,rcout) &
           + tmp*dfcut1(rij,rcin,rcout)
      
    endif
    return
  end function dphi
!=======================================================================
  function veq(r)
    implicit none
    real(8),intent(in):: r
    real(8):: veq
!.....TODO: code...
    veq = 0d0
    return
  end function veq
!=======================================================================
  function xi(x)
    implicit none
    real(8),intent(in):: x
    real(8):: xi

    xi= 0.1818d0*exp(-3.2d0*x) &
         +0.5099d0*exp(-0.9423d0*x) &
         +0.2802d0*exp(-0.4029d0*x) &
         +0.02817d0*exp(-0.2016d0*x)
    return
  end function xi
!=======================================================================
  function dxi(x)
    implicit none
    real(8),intent(in):: x
    real(8):: dxi

    dxi= -0.58176d0*exp(-3.2d0*x) &
         -0.48047877d0*exp(-0.9423d0*x) &
         -0.11289258d0*exp(-0.4029d0*x) &
         -0.005679072d0*exp(-0.2016d0*x)
    return
  end function dxi
!=======================================================================
  function zeta(x)
    implicit none
    real(8),intent(in):: x
    real(8):: zeta

    zeta = (3d0*x**5 -10d0*x**3 +15d0*x +8d0)/16d0
    return
  end function zeta
!=======================================================================
  function dzeta(x)
    implicit none
    real(8),intent(in):: x
    real(8):: dzeta

    dzeta = (15d0*x**4 -30d0*x**2 +15d0)/16d0
    return
  end function dzeta
end module EAM
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd lib"
!     End:
