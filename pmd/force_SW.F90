module SW
!-----------------------------------------------------------------------
!                     Last modified: <2021-03-02 13:35:58 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
  use pmdvars,only: nspmax
  include "./const.h"
  
  integer,parameter:: ioprms = 50
  character(len=128):: paramsdir = '.'
  character(len=128),parameter:: paramsfname = 'in.params.SW'
!-----Si mass (to be multiplied by umass)
  real(8),parameter:: am_si = 28.0855d0
!.....length scaling factor for matching this potential to VASP
!  real(8),parameter:: sfac  = 1.0062662d0
  real(8),parameter:: sfac  = 1d0
!.....number of parameters
  integer,parameter:: nprms = 10
!.....Small enough value for some criterion
  real(8),parameter:: eps = 1d-10

!-----SW unit energy in eV
  real(8):: swe   = 2.1678d0
!-----SW unit length in Ang
  real(8):: swl   = 2.0951d0*sfac
!-----si-si
  real(8):: swa   = 7.049556277d0
  real(8):: swb   = 0.6022245584d0
  real(8):: swp   = 4.d0
  real(8):: swq   = 0.d0
  real(8):: swc   = 1.d0
  real(8):: swrc  = 1.8d0
!-----si-si-si
  real(8):: sws   = 21.d0
  real(8):: swt   = 1.2d0

  integer,parameter:: msp = nspmax
  integer:: nsp

  logical:: interact(msp,msp)
  logical:: interact3(msp,msp,msp)
  real(8):: aswe,aswl
  real(8):: aswa(msp,msp)
  real(8):: aswb(msp,msp)
  real(8):: aswp(msp,msp)
  real(8):: aswq(msp,msp)
  real(8):: aswc(msp,msp)
  real(8):: aswrc(msp,msp)
  real(8):: asws(msp,msp,msp)
  real(8):: aswt(msp,msp,msp)

contains
  subroutine force_SW(namax,natm,tag,ra,nnmax,aa,strs,h,hi,tcom &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_world,myid,epi,epot,nismax,specorder,lstrs,iprint)
!-----------------------------------------------------------------------
!  Parallel implementation of SW(Si) force calculation for pmd
!    - 2014.04.07 by R.K.
!      Parameters are loaded at the first call.
!    - 2010.03.29 by R.K.
!      1st version.
!-----------------------------------------------------------------------
    implicit none
    include "mpif.h"
    include "./params_unit.h"
!    include "params_SW_Si.h"
    integer,intent(in):: namax,natm,nnmax,nismax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_world,myid,lspr(0:nnmax,namax),nex(3)
    real(8),intent(in):: ra(3,namax),tag(namax) &
         ,h(3,3),hi(3,3),sv(3,6),rc
    real(8),intent(inout):: tcom
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    character(len=3),intent(in):: specorder(msp)
    logical,intent(in):: lstrs

!-----local
    integer:: i,j,k,l,m,n,ixyz,jxyz,is,js,ks,ierr,nbl
    real(8):: rij,rik,riji,riki,rij2,rik2,rc2,src,src2,srcij,srcij2,srcik,srcik2
    real(8):: tmp,tmpj(3),tmpk(3),vexp,df2,csn,tcsn,tcsn2,dhrij,dhrik &
         ,dhcsn,vol,voli,volj,volk,drij(3),rcmax
    real(8):: drik(3),dcsni(3),dcsnj(3),dcsnk(3),drijc,drikc,x,y,z,bl &
         ,xi(3),xj(3),xk(3),xij(3),xik(3),at(3)
    real(8):: epotl,epotl2,epotl3,epott
    real(8),save:: swli,a8d3r3,rcmax2
    real(8),save,allocatable:: aa2(:,:),aa3(:,:)
    real(8),allocatable,save:: strsl(:,:,:)
!-----1st call
    logical,save:: l1st=.true.

!-----only at 1st call
    if( l1st ) then
      call read_params_SW(myid,mpi_world,iprint,specorder)
      allocate(aa2(3,namax),aa3(3,namax),strsl(3,3,namax))
!-------check rc
      rcmax = 0d0
      do is=1,msp
        do js=1,msp
          rcmax = max(rcmax,aswrc(is,js)*aswl)
        enddo
      enddo
      rcmax2= rcmax*rcmax
      if( myid.eq.0 .and. iprint.ge.ipl_basic ) then
        write(6,'(a,es12.4)') ' rc of input         =',rc
        write(6,'(a,es12.4)') ' rc of this potential=',rcmax
      endif
      if( rc .lt. rcmax ) then
!!$      if( int(rc*100d0) &
!!$           .ne.int(swrc*swl*100d0) ) then
        if( myid.eq.0 ) then
          write(6,'(1x,a)') "ERROR: Cutoff radius is not appropriate !!!"
          write(6,'(1x,a,es12.4)') "  rc should be longer than ", rcmax
        endif
        call mpi_finalize(ierr)
        stop
      endif
      swli= 1d0/swl
!!$      a8d3r3= 8d0/(3d0*sqrt(3d0))
!!$      avol= 5.427d0**3/8
!-------finally set l1st
      l1st=.false.
    endif

    if( size(aa2).ne.3*namax ) then
      deallocate(aa2,aa3,strsl)
      allocate(aa2(3,namax),aa3(3,namax),strsl(3,3,namax))
    endif

    epotl= 0d0
    epi(1:natm+nb)= 0d0
    strsl(1:3,1:3,1:natm+nb)= 0d0

!-----2 body term
    epotl2= 0d0
    aa2(1:3,1:natm+nb)=0d0
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is= int(tag(i))
      do k=1,lspr(0,i)
        j=lspr(k,i)
        if(j.eq.0) exit
        if(j.le.i) cycle
        js= int(tag(j))
        if( .not. interact(is,js) ) cycle
        src= aswrc(is,js)
        xj(1:3)= ra(1:3,j)
        x = xj(1) -xi(1)
        y = xj(2) -xi(2)
        z = xj(3) -xi(3)
        xij(1:3)= (h(1:3,1)*x +h(1:3,2)*y +h(1:3,3)*z)/aswl
        rij2 = xij(1)*xij(1) +xij(2)*xij(2) +xij(3)*xij(3)
        rij = dsqrt(rij2)
        if( rij.ge.src ) cycle
!!$        rij = dlspr(0,k,i) /aswl
!!$        if( rij.ge.src ) cycle
!!$        xij(1:3) = dlspr(1:3,k,i) /aswl
        riji= 1d0/rij
        drijc= 1d0/(rij-src)
        vexp=exp(aswc(is,js)*drijc)
!---------potential
        tmp= 0.5d0 *aswe *aswa(is,js) *vexp &
             *(aswb(is,js)*riji**aswp(is,js) -riji**aswq(is,js))
        epi(i)= epi(i) +tmp
        epotl2= epotl2 +tmp
        if( j.le.natm ) then
          epi(j)= epi(j) +tmp
          epotl2= epotl2 +tmp
        endif
!---------force
        df2= -aswe*aswa(is,js)*vexp*(aswp(is,js)*aswb(is,js) &
             *(riji**(aswp(is,js)+1d0)) &
             -aswq(is,js)*(riji**(aswq(is,js)+1d0)) &
             +(aswb(is,js)*(riji**aswp(is,js))&
             -riji**aswq(is,js))*aswc(is,js)*drijc*drijc)
        drij(1:3) = -xij(1:3)*riji /aswl
        aa2(1:3,i)= aa2(1:3,i) -df2*drij(1:3)
        aa2(1:3,j)= aa2(1:3,j) +df2*drij(1:3)
!-----------Stress
        if( .not. lstrs ) cycle
        if( j.le.natm ) then
          do jxyz=1,3
            strsl(1:3,jxyz,i)= strsl(1:3,jxyz,i) &
                 -0.5d0*xij(jxyz)*aswl*(-df2*drij(1:3))
            strsl(1:3,jxyz,j)= strsl(1:3,jxyz,j) &
                 -0.5d0*xij(jxyz)*aswl*(-df2*drij(1:3))
          enddo
        else
          do jxyz=1,3
            strsl(1:3,jxyz,i)= strsl(1:3,jxyz,i) &
                 -0.5d0*xij(jxyz)*aswl*(-df2*drij(1:3))
          enddo
        endif

      enddo
    enddo

!-----3 body term
    epotl3= 0d0
    aa3(1:3,1:natm+nb)=0d0
!-----atom (i)
    do i=1,natm
      xi(1:3)=ra(1:3,i)
      is= int(tag(i))
      do n=1,lspr(0,i)
!---------atom (j)
        j=lspr(n,i)
        if(j.eq.0) exit
        if( j.eq.i ) cycle
        js= int(tag(j))
        srcij= aswrc(is,js)
        xj(1:3)= ra(1:3,j)
        x = xj(1) -xi(1)
        y = xj(2) -xi(2)
        z = xj(3) -xi(3)
        xij(1:3)= (h(1:3,1)*x +h(1:3,2)*y +h(1:3,3)*z)/aswl
        rij2 = xij(1)*xij(1) +xij(2)*xij(2) +xij(3)*xij(3)
        rij = dsqrt(rij2)
        if( rij.ge.srcij ) cycle
!!$        rij= dsqrt(rij2)
!!$        rij = dlspr(0,n,i) /aswl
!!$        xij(1:3) = dlspr(1:3,n,i) /aswl
        riji= 1d0/rij
        drijc= 1d0/(rij-srcij)
!---------atom (k)
        do m=1,lspr(0,i)
          k=lspr(m,i)
          if(k.eq.0) exit
          if( k.le.j .or. k.eq.i ) cycle
          ks= int(tag(k))
          if( .not. interact3(is,js,ks) ) cycle
          srcik= aswrc(is,ks)
          xk(1:3)= ra(1:3,k)
          x = xk(1) -xi(1)
          y = xk(2) -xi(2)
          z = xk(3) -xi(3)
          xik(1:3)= (h(1:3,1)*x +h(1:3,2)*y +h(1:3,3)*z)/aswl
          rik2 = xik(1)*xik(1) +xik(2)*xik(2) +xik(3)*xik(3)
          rik = dsqrt(rik2)
          if( rik.ge.srcik ) cycle
!!$          rik = dlspr(0,m,i) /aswl
!!$          if( rik.ge.srcik ) cycle
!!$          xik(1:3) = dlspr(1:3,m,i) /aswl
          riki= 1d0/rik
          drikc= 1d0/(rik-srcik)
!-----------common term
          csn=(xij(1)*xik(1) +xij(2)*xik(2) +xij(3)*xik(3)) &
               * (riji*riki)
          tcsn = csn +1d0/3d0
          tcsn2= tcsn*tcsn
          vexp= dexp(aswt(is,js,ks)*drijc +aswt(is,js,ks)*drikc)
!-----------potential
          tmp= aswe *asws(is,js,ks) *vexp *tcsn2
          epi(i)= epi(i) +tmp
          epotl3= epotl3 +tmp
!-----------force
          dhrij= -asws(is,js,ks) *aswt(is,js,ks) *vexp *tcsn2 *drijc*drijc
          dhrik= -asws(is,js,ks) *aswt(is,js,ks) *vexp *tcsn2 *drikc*drikc
          dhcsn= 2d0 *asws(is,js,ks) *vexp *tcsn 
          drij(1:3)= -xij(1:3)*riji /aswl
          drik(1:3)= -xik(1:3)*riki /aswl
          dcsnj(1:3)= (-xij(1:3)*csn*(riji*riji) +xik(1:3)*(riji*riki)) /aswl
          dcsnk(1:3)= (-xik(1:3)*csn*(riki*riki) +xij(1:3)*(riji*riki)) /aswl
          dcsni(1:3)= -dcsnj(1:3) -dcsnk(1:3)
          tmpj(1:3)= aswe*(dhcsn*dcsnj(1:3) +dhrij*(-drij(1:3)))
          tmpk(1:3)= aswe*(dhcsn*dcsnk(1:3) +dhrik*(-drik(1:3)))
          aa3(1:3,i)= aa3(1:3,i) +(tmpj(1:3)+tmpk(1:3))
          aa3(1:3,j)= aa3(1:3,j) -tmpj(1:3)
          aa3(1:3,k)= aa3(1:3,k) -tmpk(1:3)
!-------------Stress
          if( .not. lstrs ) cycle
          do jxyz=1,3
            strsl(1:3,jxyz,i)=strsl(1:3,jxyz,i) &
                 -0.5d0*xij(jxyz)*aswl*tmpj(1:3) & !*volj &
                 -0.5d0*xik(jxyz)*aswl*tmpk(1:3) !*volk
            strsl(1:3,jxyz,j)=strsl(1:3,jxyz,j) &
                 -0.5d0*xij(jxyz)*aswl*tmpj(1:3) !*volj
            strsl(1:3,jxyz,k)=strsl(1:3,jxyz,k) &
                 -0.5d0*xik(jxyz)*aswl*tmpk(1:3) !*volk
          enddo

        enddo
      enddo
    enddo

    call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,aa3,3)
    call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,epi,1)
    aa(1:3,1:natm)= aa(1:3,1:natm) +aa2(1:3,1:natm) +aa3(1:3,1:natm)
    
    if( lstrs ) then
      call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
           ,nn,mpi_world,strsl,9)
      strs(1:3,1:3,1:natm) = strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)
    endif

!-----gather epot
    epotl= epotl2 +epotl3
!!$    epotl= epotl3
    call mpi_allreduce(epotl,epott,1,mpi_real8,mpi_sum,mpi_world,ierr)
    epot= epot +epott
    if( iprint.ge.ipl_info ) print *,'SW epot = ',epott
    return
  end subroutine force_SW
!=======================================================================
  subroutine read_params_SW(myid,mpi_world,iprint,specorder)
    use util, only: num_data
    implicit none
    include 'mpif.h'
    integer,intent(in):: myid,mpi_world,iprint
    character(len=3),intent(in):: specorder(msp)

    integer:: itmp,ierr,isp,jsp,ksp,nd
    real(8):: rctmp,tswa,tswb,tswp,tswq,tswc,tswrc,tsws,tswt
    logical:: lexist
    character(len=128):: cfname,ctmp,cline
    character(len=3):: cspi,cspj

!!$    integer,external:: num_data

!.....read parameters at the 1st call
    if( myid.eq.0 ) then
!.....Interact between only Si
      interact(1:msp,1:msp) = .false.
      interact3(1:msp,1:msp,1:msp) = .false.
      do isp=1,msp
        cspi = specorder(isp)
        if( trim(cspi).ne.'Si' ) cycle
        interact(isp,isp) = .true.
        interact3(isp,isp,isp) = .true.
      enddo
!.....Initialize parameters
      aswe = swe
      aswl = swl
      aswa(1,1) = swa
      aswb(1,1) = swb
      aswp(1,1) = swp
      aswq(1,1) = swq
      aswc(1,1) = swc
      aswrc(1,1) = swrc
      asws(1,1,1) = sws
      aswt(1,1,1) = swt
!.....Check whether the file exists      
      cfname = trim(paramsdir)//'/'//trim(paramsfname)
      inquire(file=cfname,exist=lexist)
      if( .not. lexist ) then
        if( iprint.ge.ipl_warn ) then
          write(6,'(a)') ' WARNING: in.params.SW does not exist !!!.'
          write(6,'(a)') '           Default parameters will be used.'
        endif
        goto 20
      endif
!.....Read file if exists
      if( iprint.ne.0 ) write(6,'(/,a)') ' SW parameters read from file:'
      open(ioprms,file=cfname,status='old')
      do while(.true.)
        read(ioprms,'(a)',end=10) cline
        nd = num_data(cline,' ')
        if( nd.eq.0 ) cycle
        if( cline(1:1).eq.'!' .or. cline(1:1).eq.'#' ) cycle
        isp = 0
        jsp = 0
        ksp = 0
        if( index(cline,'unit').ne.0 .and. nd.eq.3 ) then
          backspace(ioprms)
          read(ioprms,*) ctmp, aswe, aswl
          if( iprint.ne.0 ) &
               write(6,'(a,2es14.4)') '   unit (energy,length) = ',aswe,aswl
        else if( nd.eq.8 ) then  ! two body
          backspace(ioprms)
          read(ioprms,*) isp,jsp,tswa,tswb,tswp,tswq,tswc,tswrc
          if( isp.gt.msp .or. jsp.gt.msp ) then
            if( iprint.ne.0 ) then
              write(6,*) 'WARNING@read_params_SW: isp/jsp greater than msp, ' &
                   //'so skip the line.'
            endif
            cycle
          endif
          if( iprint.ne.0 ) &
               write(6,'(2i4,6es14.4)') isp,jsp,tswa,tswb,tswp,tswq,tswc,tswrc
          interact(isp,jsp) = .true.
          aswa(isp,jsp) = tswa
          aswb(isp,jsp) = tswb
          aswp(isp,jsp) = tswp
          aswq(isp,jsp) = tswq
          aswc(isp,jsp) = tswc
          aswrc(isp,jsp) = tswrc
!.....Symmetrize if needed
          if( isp.ne.jsp ) then
            interact(isp,jsp) = interact(jsp,isp)
            aswa(isp,jsp) = aswa(jsp,isp)
            aswb(isp,jsp) = aswb(jsp,isp)
            aswp(isp,jsp) = aswp(jsp,isp)
            aswq(isp,jsp) = aswq(jsp,isp)
            aswc(isp,jsp) = aswc(jsp,isp)
            aswrc(isp,jsp) = aswrc(jsp,isp)
          endif
        else if( nd.eq.5 ) then  ! three body
          backspace(ioprms)
          read(ioprms,*) isp,jsp,ksp,tsws,tswt
          if( isp.gt.msp .or. jsp.gt.msp .or. ksp.gt.msp ) then
            if( iprint.ne.0 ) &
                 write(6,*) 'WARNING@read_params_SW: isp/jsp/ksp greater than msp, ' &
                 //'so skip the line.'
            cycle
          endif
          if( iprint.ne.0 ) write(6,'(3i4,6es14.4)') isp,jsp,ksp,tsws,tswt
          interact3(isp,jsp,ksp) = .true.
          asws(isp,jsp,ksp) = tsws
          aswt(isp,jsp,ksp) = tswt
!.....Symmetrize if needed
          if( jsp.ne.ksp ) then
            interact3(isp,ksp,jsp) = interact3(isp,jsp,ksp)
            asws(isp,ksp,jsp) = asws(isp,jsp,ksp)
            aswt(isp,ksp,jsp) = aswt(isp,jsp,ksp)
          endif
        else
          if( iprint.ne.0 ) then
            write(6,*) 'WARNING@read_params_SW: number of entry wrong, ' &
                 //'so skip the line.'
          endif
          cycle
        endif
      enddo
10    close(ioprms)
      
!!$      read(50,*) itmp,rctmp
!!$      if( itmp.ne.nprms ) then
!!$        write(6,'(a)') ' [Error] itmp.ne.nprms'
!!$        write(6,'(a,i3)') '  itmp =',itmp
!!$        stop
!!$      endif
!!$      read(50,*) swe
!!$      read(50,*) swl
!!$      read(50,*) swa
!!$      read(50,*) swb
!!$      read(50,*) swp
!!$      read(50,*) swq
!!$      read(50,*) swc
!!$      read(50,*) swrc
!!$      read(50,*) sws
!!$      read(50,*) swt
!!$      close(50)
    endif

20  continue
    call mpi_bcast(interact,msp*msp,mpi_logical,0,mpi_world,ierr)
    call mpi_bcast(aswe,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(aswl,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(aswa,msp*msp,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(aswb,msp*msp,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(aswp,msp*msp,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(aswq,msp*msp,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(aswc,msp*msp,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(aswrc,msp*msp,mpi_real8,0,mpi_world,ierr)

    call mpi_bcast(interact3,msp*msp*msp,mpi_logical,0,mpi_world,ierr)
    call mpi_bcast(asws,msp*msp*msp,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(aswt,msp*msp*msp,mpi_real8,0,mpi_world,ierr)

!!$    call mpi_bcast(swe,1,mpi_double_precision,0,mpi_world,ierr)
!!$    call mpi_bcast(swl,1,mpi_double_precision,0,mpi_world,ierr)
!!$    call mpi_bcast(swa,1,mpi_double_precision,0,mpi_world,ierr)
!!$    call mpi_bcast(swb,1,mpi_double_precision,0,mpi_world,ierr)
!!$    call mpi_bcast(swp,1,mpi_double_precision,0,mpi_world,ierr)
!!$    call mpi_bcast(swq,1,mpi_double_precision,0,mpi_world,ierr)
!!$    call mpi_bcast(swc,1,mpi_double_precision,0,mpi_world,ierr)
!!$    call mpi_bcast(swrc,1,mpi_double_precision,0,mpi_world,ierr)
!!$    call mpi_bcast(sws,1,mpi_double_precision,0,mpi_world,ierr)
!!$    call mpi_bcast(swt,1,mpi_double_precision,0,mpi_world,ierr)
    return
  end subroutine read_params_SW
end module SW
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
