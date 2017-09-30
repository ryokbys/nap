module EAM
!-----------------------------------------------------------------------
!                     Last modified: <2017-09-30 19:19:34 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
!  Parallel implementation of the EAM pontential.
!-----------------------------------------------------------------------
  implicit none

  character(len=128):: paramsdir = '.'
  character(len=128),parameter:: paramsfname = 'in.params.EAM'

  integer,parameter:: ioprms = 20
  
!.....Max num of species
  integer,parameter:: msp = 9
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

  real(8),allocatable,dimension(:,:):: ea_a,ea_b,ea_c,ea_re, &
       ea_alp,ea_beta,ea_xi,ea_rc
  logical,allocatable:: interact(:,:)

  logical:: lprmset = .false.

!.....Parameters from fitpot
  integer:: nprms
  real(8),allocatable:: prms(:)

contains
  subroutine init_EAM()
!
!  Allocate and initialize parameters to be used.
!
    if( .not.allocated(ea_a) ) then
      allocate(ea_a(msp,msp), ea_b(msp,msp), ea_c(msp,msp), &
           ea_re(msp,msp), ea_alp(msp,msp), ea_beta(msp,msp),&
           ea_xi(msp,msp), interact(msp,msp), ea_rc(msp,msp))
    endif
    
  end subroutine init_EAM
!=======================================================================
  subroutine read_params_EAM(myid_md,mpi_md_world,iprint)
!
!  Read parameters from file
!
    include 'mpif.h'
    integer,intent(in):: myid_md,mpi_md_world,iprint

    integer:: isp,jsp,ierr
    real(8):: a,b,c,re,alp,beta,xi,rc
    character(len=128):: cline,fname

    if( myid_md.eq.0 ) then
      fname = trim(paramsdir)//'/'//trim(paramsfname)
      open(ioprms,file=trim(fname),status='old')
      ea_a(1:msp,1:msp) = 0d0
      ea_b(1:msp,1:msp) = 0d0
      ea_c(1:msp,1:msp) = 0d0
      ea_re(1:msp,1:msp) = 0d0
      ea_alp(1:msp,1:msp) = 0d0
      ea_beta(1:msp,1:msp) = 0d0
      ea_xi(1:msp,1:msp) = 0d0
      ea_rc(1:msp,1:msp) = 0d0
      interact(1:msp,1:msp) = .false.

      do while(.true.)
        read(ioprms,*,end=10) cline
        if( cline(1:1).eq.'#' .or. cline(1:1).eq.'1' ) cycle
        backspace(ioprms)
        read(ioprms,*) isp,jsp,a,b,c,re,alp,beta,xi,rc
        if( isp.gt.msp .or. jsp.gt.msp ) then
          print *,'Warning @read_params_EAM: isp/jsp is greater than msp,'&
               //' skip reading the line.'
          cycle
        endif

        ea_a(isp,jsp) = a
        ea_b(isp,jsp) = b
        ea_c(isp,jsp) = c
        ea_re(isp,jsp) = re
        ea_alp(isp,jsp) = alp
        ea_beta(isp,jsp) = beta
        ea_xi(isp,jsp) = xi
        ea_rc(isp,jsp) = rc
        interact(isp,jsp) = .true.
!.....Symmetrize parameters
        ea_a(jsp,isp) = a
        ea_b(jsp,isp) = b
        ea_c(jsp,isp) = c
        ea_re(jsp,isp) = re
        ea_alp(jsp,isp) = alp
        ea_beta(jsp,isp) = beta
        ea_xi(jsp,isp) = xi
        ea_rc(jsp,isp) = rc
        interact(jsp,isp) = .true.
      enddo
10    continue
      close(ioprms)
      if( iprint.ne.0 ) then
        print *,'Finished reading '//trim(fname)
        print *,''
      endif
    endif
    call mpi_bcast(ea_a,msp*msp,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(ea_b,msp*msp,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(ea_c,msp*msp,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(ea_re,msp*msp,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(ea_alp,msp*msp,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(ea_beta,msp*msp,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(ea_xi,msp*msp,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(ea_rc,msp*msp,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(interact,msp*msp,mpi_integer,0,mpi_md_world,ierr)
    
  end subroutine read_params_EAM
!=======================================================================
  subroutine set_paramsdir_EAM(dname)
!
!  Accessor routine for setting paramsdir
!
    implicit none
    character(len=*),intent(in):: dname

    paramsdir = trim(dname)
    return
  end subroutine set_paramsdir_EAM
!=======================================================================
  subroutine set_params_EAM(nprms_in,prms_in)
!
!  Accessor routine to set EAM parameters from outside.
!  Curretnly this routine is supposed to be called only on serial run.
!
    integer,intent(in):: nprms_in
    real(8),intent(in):: prms_in(nprms_in)

    nprms = nprms_in
    if( .not.allocated(prms) ) allocate(prms(nprms))
    prms(1:nprms) = prms_in(1:nprms_in)
    lprmset = .true.
    return
  end subroutine set_params_EAM
!=======================================================================
  subroutine update_params_EAM()
!
!  Update EAM parameters by taking parameter values from params array.
!  This routine would be called only from externally within fitpot.
!
    integer:: i,inc,is,js
    real(8),parameter:: tiny = 1d-10

    if( .not.lprmset ) then
      print *,'ERROR: params have not been set yet.'
      stop
    endif

    ea_a(1:msp,1:msp) = 0d0
    ea_b(1:msp,1:msp) = 0d0
    ea_c(1:msp,1:msp) = 0d0
    ea_re(1:msp,1:msp) = 0d0
    ea_alp(1:msp,1:msp) = 0d0
    ea_beta(1:msp,1:msp) = 0d0
    ea_xi(1:msp,1:msp) = 0d0
    interact(1:msp,1:msp) = .false.

    inc = 0
    do is=1,nsp
      do js=is,nsp
        inc = inc + 1
        ea_a(is,js) = prms(inc)
        ea_a(js,is) = ea_a(is,js)
        inc = inc + 1
        ea_b(is,js) = prms(inc)
        ea_b(js,is) = ea_b(is,js)
        inc = inc + 1
        ea_c(is,js) = prms(inc)
        ea_c(js,is) = ea_c(is,js)
        inc = inc + 1
        ea_re(is,js) = prms(inc)
        ea_re(js,is) = ea_re(is,js)
        inc = inc + 1
        ea_alp(is,js) = prms(inc)
        ea_alp(js,is) = ea_alp(is,js)
        inc = inc + 1
        ea_beta(is,js) = prms(inc)
        ea_beta(js,is) = ea_beta(is,js)
        inc = inc + 1
        ea_xi(is,js) = prms(inc)
        ea_xi(js,is) = ea_xi(is,js)
        if( abs(ea_a(is,js)).gt.tiny .or. abs(ea_b(is,js)).gt.tiny &
             .or. abs(ea_c(is,js)).gt.tiny ) then
          interact(is,js) = .true.
          interact(js,is) = .true.
        endif
      enddo
    enddo
    return
  end subroutine update_params_EAM
!=======================================================================
  subroutine force_EAM_SM(namax,natm,tag,ra,nnmax,aa,strs,h,hi,tcom &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
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
!!$    include "params_SM_Al.h"
    integer,intent(in):: namax,natm,nnmax,nismax,lspr(0:nnmax,namax)&
         ,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_md_world,myid_md,nex(3)
    real(8),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),sv(3,6) &
         ,acon(nismax),rc,tag(namax)
    real(8),intent(inout):: tcom
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical:: lstrs

    integer:: i,j,k,l,m,n,ierr,is,js,ixyz,jxyz
    real(8):: xij(3),rij,dfi,dfj,drhoij,drdxi(3),drdxj(3),r,dphi,at(3)
    real(8):: x,y,z,xi(3),epotl,epott,tmp
    logical,save:: l1st=.true.
    real(8),allocatable,save:: sqrho(:)
    real(8),allocatable,save:: strsl(:,:,:)
!!$    real(8),save:: exrc,phic,dphic
    real(8),allocatable,save:: exrc(:,:),phic(:,:),dphic(:,:)

    if( l1st ) then
!.....Check cutoff radius
      if( allocated(sqrho) ) deallocate(sqrho)
      allocate(sqrho(namax+nbmax))
      if( allocated(strsl) ) deallocate(strsl)
      allocate(strsl(3,3,namax))
      allocate(exrc(msp,msp),phic(msp,msp),dphic(msp,msp))
      exrc(:,:) = 0d0
      phic(:,:) = 0d0
      dphic(:,:)= 0d0
      do is=1,msp
        do js=is,msp
          if( .not.interact(is,js) ) cycle
          if( rc.lt.ea_rc(is,js) ) then
            if( myid_md.eq.0 ) then
              write(6,*) ' Error: rc is smaller than one of EAM rcs'
            endif
            call mpi_finalize(ierr)
            stop
          endif
!-----smoothing embeded term
          exrc(is,js)= exp(-ea_beta(is,js)*(ea_rc(is,js)-ea_re(is,js)))
!-----smoothing 2-body term
          r= ea_rc(is,js) -ea_re(is,js)
          phic(is,js) = 2d0*ea_b(is,js) &
               *exp(-0.5d0*ea_beta(is,js)*r) &
               -ea_c(is,js)*(1d0+ea_alp(is,js)*r) &
               *exp(-ea_alp(is,js)*r)
          dphic(is,js)= -ea_beta(is,js)*ea_b(is,js) &
               *exp(-0.5d0*ea_beta(is,js)*r)  &
               +ea_c(is,js)*ea_alp(is,js)**2 &
               *r*exp(-ea_alp(is,js)*r)
!.....Symmetrize
          exrc(js,is) = exrc(is,js)
          phic(js,is) = phic(is,js)
          dphic(js,is)= dphic(is,js)
        enddo
      enddo
      l1st=.false.
    endif

    if( size(strsl).lt.3*3*namax ) then
      deallocate(strsl)
      allocate(strsl(3,3,namax))
    endif
    if( size(sqrho).lt.namax+nbmax ) then
      deallocate(sqrho)
      allocate(sqrho(namax+nbmax))
    endif
    
    epotl= 0d0
    sqrho(1:natm)= 0d0
    strsl(1:3,1:3,1:natm) = 0d0

!-----rho(i)
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is = int(tag(i))
      do k=1,lspr(0,i)
        j=lspr(k,i)
        if(j.eq.0) exit
        js = int(tag(j))
        x= ra(1,j) -xi(1)
        y= ra(2,j) -xi(2)
        z= ra(3,j) -xi(3)
        xij(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
        rij=sqrt(xij(1)*xij(1)+ xij(2)*xij(2) +xij(3)*xij(3))
        if( rij.gt.ea_rc(is,js) ) cycle
        sqrho(i)=sqrho(i) +exp(-ea_beta(is,js)*(rij-ea_re(is,js))) &
             -exrc(is,js) -(rij-ea_rc(is,js))*(-ea_beta(is,js))*exrc(is,js)
      enddo
      sqrho(i)= dsqrt(sqrho(i))
    enddo

!-----copy rho of boundary atoms
    call copy_dba_fwd(tcom,namax,natm,nb,nbmax,lsb,nex,&
         lsrc,myparity,nn,sv,mpi_md_world,sqrho,1)
!!$    call copy_rho_ba(tcom,namax,natm,nb,nbmax,lsb,lsrc,myparity,nn,sv &
!!$         ,mpi_md_world,sqrho)

!-----dE/dr_i
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is = int(tag(i))
      dfi= -0.5d0*ea_a(is,js)/sqrho(i)
      do k=1,lspr(0,i)
        j=lspr(k,i)
        if(j.eq.0) exit
        if(j.le.i) cycle
        js = int(tag(j))
        x= ra(1,j) -xi(1)
        y= ra(2,j) -xi(2)
        z= ra(3,j) -xi(3)
        xij(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
        rij=sqrt(xij(1)**2+ xij(2)**2 +xij(3)**2)
        if( rij.gt.ea_rc(is,js) ) cycle
        drdxi(1:3)= -xij(1:3)/rij
        r= rij -ea_re(is,js)
!.....2-body term
        tmp= 0.5d0 *( 2d0*ea_b(is,js)*exp(-0.5d0*ea_beta(is,js)*r) &
             -ea_c(is,js)*(1d0+ea_alp(is,js)*r)*exp(-ea_alp(is,js)*r) &
             -phic(is,js) -(rij-ea_rc(is,js))*dphic(is,js) )
        epi(i)= epi(i) +tmp
        epi(j)= epi(j) +tmp
        if(j.le.natm) then
          epotl=epotl +tmp +tmp
        else
          epotl=epotl +tmp
        endif
        dphi= -ea_beta(is,js)*ea_b(is,js)*exp(-0.5d0*ea_beta(is,js)*r)  &
             +ea_c(is,js)*ea_alp(is,js)*ea_alp(is,js)*r &
             *exp(-ea_alp(is,js)*r) &
             -dphic(is,js)
        aa(1:3,i)=aa(1:3,i) -dphi*drdxi(1:3)
        aa(1:3,j)=aa(1:3,j) +dphi*drdxi(1:3)
!.....Atomic stress for 2-body terms
        if( lstrs ) then
          do ixyz=1,3
            do jxyz=1,3
              strsl(jxyz,ixyz,i)=strsl(jxyz,ixyz,i) &
                   -0.5d0*dphi*xij(ixyz)*(-drdxi(jxyz))
              strsl(jxyz,ixyz,j)=strsl(jxyz,ixyz,j) &
                   -0.5d0*dphi*xij(ixyz)*(-drdxi(jxyz))
            enddo
          enddo
        endif
!.....Embedded term
        drhoij= -ea_beta(is,js)*exp(-ea_beta(is,js)*r) &
             +ea_beta(is,js)*exrc(is,js)
        dfj= -0.5d0*ea_a(is,js)/sqrho(j)
        tmp = (dfi+dfj)*drhoij
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
      epi(i)=epi(i) -ea_a(is,js)*sqrho(i)
      epotl=epotl -ea_a(is,js)*sqrho(i)
    enddo

    if( lstrs ) then
      call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
           ,nn,mpi_md_world,strsl,9)
      strs(1:3,1:3,1:natm)= strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)
    endif

!-----gather epot
    call mpi_allreduce(epotl,epott,1,mpi_real8 &
         ,mpi_sum,mpi_md_world,ierr)
    epot= epot +epott

!      deallocate(sqrho)
  end subroutine force_EAM_SM
end module EAM
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
