module BMH
!-----------------------------------------------------------------------
!                     Last modified: <2021-07-07 14:56:57 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
!  Parallel implementation of fitpot Born-Mayer-Huggins (BMH) potential.
!  This potential should be used with Coulomb and dipole potentials.
!-----------------------------------------------------------------------
  use pmdvars, only: nspmax
  use util,only: csp2isp
  implicit none
  save
  include "./const.h"
  character(len=128):: paramsdir = '.'
  character(len=128),parameter:: cprmfname = 'in.params.BMH'

  integer,parameter:: ioprms = 20

!.....Max number of species available in the potential
  integer:: nspcs
  real(8):: bmh_aij(nspmax,nspmax),bmh_alpij(nspmax,nspmax)
  logical:: interact(nspmax,nspmax)

  logical:: lprmset_BMH
  
!.....params
  integer:: nprms
  real(8),allocatable:: params(:)

contains
  subroutine force_BMH(namax,natm,tag,ra,nnmax,aa,strs,h,hi,tcom &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid,epi,epot,nismax,specorder,lstrs,iprint,l1st)
    use util,only: itotOf
    implicit none
    include "mpif.h"
    include "./params_unit.h"
    integer,intent(in):: namax,natm,nnmax,nismax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),lspr(0:nnmax,namax),nex(3)
    integer,intent(in):: mpi_md_world,myid
    real(8),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),rc &
         ,tag(namax),sv(3,6)
    real(8),intent(inout):: tcom
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: l1st
    character(len=3),intent(in):: specorder(nspmax)
    logical:: lstrs

    integer:: i,j,k,l,m,n,ierr,is,js,ixyz,jxyz
    real(8):: xi(3),xij(3),rij(3),dij,diji,dvdr,dij2 &
         ,drdi(3),drdj(3),x,y,z,epotl,epott,at(3),tmp,tmp2 &
         ,aij,alpij,vs2b &
         ,vrc,dvdrc,dvs2b
    real(8):: vs2bc,dvs2bc
    real(8),save:: vrcs(nspmax,nspmax),dvdrcs(nspmax,nspmax)
    real(8),allocatable,save:: strsl(:,:,:)
    

    real(8),save:: rcmax2

    if( l1st ) then
      if( allocated(strsl) ) deallocate(strsl)
      allocate(strsl(3,3,namax))
      rcmax2 = rc*rc
!.....Initialize smooth cutoff
      vrcs(:,:) = 0d0
      dvdrcs(:,:)= 0d0
      do is=1,nspmax
        do js=is,nspmax
          if( .not. interact(is,js) ) cycle
          aij= bmh_aij(is,js)
          alpij= bmh_alpij(is,js)
          vs2bc = aij *exp(-alpij*rc)
          dvs2bc= -alpij*aij*exp(-alpij*rc)
          vrcs(is,js) = vs2bc
          vrcs(js,is) = vrcs(is,js)
          dvdrcs(is,js) = dvs2bc
          dvdrcs(js,is) = dvdrcs(is,js)
        enddo
      enddo
    endif

    if( size(strsl).lt.3*3*namax ) then
      deallocate(strsl)
      allocate(strsl(3,3,namax))
    endif

    epotl= 0d0
    strsl(1:3,1:3,1:natm+nb) = 0d0

!-----loop over resident atoms
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is = int(tag(i))
      do k=1,lspr(0,i)
        j=lspr(k,i)
        if(j.eq.0) exit
        if(j.le.i) cycle
        js = int(tag(j))
        if( .not.interact(is,js) ) cycle
        xij(1:3)= ra(1:3,j) -xi(1:3)
        rij(1:3)= h(1:3,1,0)*xij(1) +h(1:3,2,0)*xij(2) +h(1:3,3,0)*xij(3)
        dij2= rij(1)**2+ rij(2)**2 +rij(3)**2
        if( dij2.gt.rcmax2 ) cycle
        dij = dsqrt(dij2)
        diji= 1d0/dij
        drdi(1:3)= -rij(1:3)*diji
        drdj(1:3)=  rij(1:3)*diji
        vrc = vrcs(is,js)
        dvdrc = dvdrcs(is,js)
        aij = bmh_aij(is,js)
        alpij = bmh_alpij(is,js)
!.....Potential
        vs2b = aij*dexp(-alpij*dij)
        tmp = vs2b
        tmp2 = 0.5d0 *(tmp -vrc -dvdrc*(dij-rc))
        if(j.le.natm) then
          epi(i)=epi(i) +tmp2
          epi(j)=epi(j) +tmp2
          epotl= epotl +tmp2 +tmp2
        else
          epi(i)=epi(i) +tmp2
          epotl= epotl +tmp2
        endif
!.....Force
        dvs2b = -alpij*aij*dexp(-alpij*dij)
        dvdr = dvs2b -dvdrc
        aa(1:3,i)= aa(1:3,i) -drdi(1:3)*dvdr
        aa(1:3,j)= aa(1:3,j) -drdj(1:3)*dvdr
!!$        print '(a,4i5,10es11.3)','i,j,is,js,dij,vs2b,dvs2b,aai,aaj=' &
!!$             ,i,j,is,js,dij,vs2b,dvs2b,aa(1:3,i),aa(1:3,j)
!.....Stress
        if( .not.lstrs ) cycle
        do ixyz=1,3
          do jxyz=1,3
            strsl(jxyz,ixyz,i)=strsl(jxyz,ixyz,i) &
                 -0.5d0*dvdr*rij(ixyz)*(-drdi(jxyz))
            strsl(jxyz,ixyz,j)=strsl(jxyz,ixyz,j) &
                 -0.5d0*dvdr*rij(ixyz)*(-drdi(jxyz))
          enddo
        enddo
      enddo
    enddo

    if( lstrs ) then
!!$      call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
!!$           ,nn,mpi_md_world,strsl,9)
      strs(1:3,1:3,1:natm)= strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)
    endif

    epott = 0d0
    call mpi_allreduce(epotl,epott,1,mpi_real8 &
         ,mpi_sum,mpi_md_world,ierr)
    epot= epot +epott
    if( iprint.ge.ipl_info ) print '(a,es15.7)',' epot BMH = ',epott
    
  end subroutine force_BMH
!=======================================================================
  subroutine set_paramsdir_BMH(dname)
!
!  Accessor routine to set paramsdir.
!
    implicit none
    character(len=*),intent(in):: dname

    paramsdir = trim(dname)
    return
  end subroutine set_paramsdir_BMH
!=======================================================================
  subroutine set_params_BMH(ndimp,params_in,ctype,interact_in)
!
!  Accessor routine to set BMH parameters from outside.
!  Curretnly this routine is supposed to be called only on serial run.
!
    integer,intent(in):: ndimp
    real(8),intent(in):: params_in(ndimp)
    character(len=*),intent(in):: ctype
    logical,intent(in):: interact_in(nspmax,nspmax)

    integer:: i,j,inc,itmp,nspt,nint

    nprms = ndimp
    if( .not.allocated(params) ) allocate(params(nprms))
    params(1:nprms) = params_in(1:ndimp)
    lprmset_BMH = .true.

    interact(:,:) = interact_in(:,:)
    nint = 0
    do i=1,nspmax
      do j=i,nspmax
        if( .not.interact(i,j) ) cycle
        nint = nint +1
      enddo
    enddo
    if( nint*2.ne.nprms ) then
      write(6,*) 'ERROR @set_params_BMH: nint*2.ne.nprms !!!'
      write(6,*) '  nint,nprms= ',nint,nprms
      write(6,*) 'Probably you need to set interactions correctly...'
      write(6,*) 'interact:'
      do i=1,nspmax
        do j=i,nspmax
          write(6,'(a,2i5,l4)') '  isp,jsp,interact(isp,jsp)=',i,j,interact(i,j)
        enddo
      enddo
      stop
    endif

    bmh_aij(:,:)= 0d0
    bmh_alpij(:,:)= 0d0

    inc = 0
    do i=1,nspmax
      do j=i,nspmax
        if( .not.interact(i,j) ) cycle
        inc= inc +1
        bmh_aij(i,j) = params(inc)
        bmh_aij(j,i) = bmh_aij(i,j)
        inc= inc +1
        bmh_alpij(i,j) = params(inc)
        bmh_alpij(j,i) = bmh_alpij(i,j)
      enddo
    enddo

    return
  end subroutine set_params_BMH
!=======================================================================
  subroutine read_params_BMH(myid_md,mpi_md_world,iprint,specorder)
    use util, only: num_data
    implicit none
    include 'mpif.h'
    include './params_unit.h'
    integer,intent(in):: myid_md,mpi_md_world,iprint
    character(len=3),intent(in):: specorder(nspmax)

    integer:: isp,jsp,ierr,ni,nj
    real(8):: aij,alpij
    character(len=3):: csp,cspi,cspj
    character(len=128):: cline,fname

    if( myid_md.eq.0 ) then
      bmh_aij(:,:) = 0d0
      bmh_alpij(:,:)= 0d0
      interact(:,:) = .false.
      fname = trim(paramsdir)//'/'//trim(cprmfname)
      open(ioprms,file=trim(fname),status='old')
      if( iprint.ge.ipl_basic ) write(6,'(/,a)') ' BMH parameters:'
      do while(.true.)
        read(ioprms,'(a)',end=10) cline
        if( num_data(cline,' ').eq.0 ) cycle
        if( cline(1:1).eq.'#' .or. cline(1:1).eq.'!' ) cycle
        if( num_data(cline,' ').eq.4 ) then
          backspace(ioprms)
          read(ioprms,*) cspi,cspj, aij,alpij
          isp = csp2isp(cspi)
          jsp = csp2isp(cspj)
          if( isp.gt.nspmax .or. jsp.gt.nspmax ) then
            write(6,*) ' Warning @read_params: since isp/jsp is greater than nspmax,'&
                 //' skip reading the line.'
            cycle
          endif
          bmh_aij(isp,jsp)= aij
          bmh_alpij(isp,jsp)= alpij
          interact(isp,jsp)= .true.
!.....Symmetrize
          bmh_aij(jsp,isp)= aij
          bmh_alpij(jsp,isp)= alpij
          interact(jsp,isp)= .true.
        endif
      enddo  ! while reading...
10    close(ioprms)
!.....Set parameters
      if( iprint.ge.ipl_basic ) print *,'  cspi,cspj,aij,alpij:'
      do isp=1,nspmax
        do jsp=isp,nspmax
          if( .not. interact(isp,jsp) ) cycle
          if( iprint.ge.ipl_basic ) then
            print '(4x,2(1x,a3),6es12.4)' &
               ,specorder(isp),specorder(jsp) &
               ,bmh_aij(isp,jsp),bmh_alpij(isp,jsp)
          endif
        enddo
      enddo
    endif  ! myid_md.eq.0

    call mpi_bcast(bmh_aij,nspmax*nspmax,mpi_integer,0,mpi_md_world,ierr)
    call mpi_bcast(bmh_alpij,nspmax*nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(interact,nspmax*nspmax,mpi_logical,0,mpi_md_world,ierr)
    
  end subroutine read_params_BMH
end module BMH
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
