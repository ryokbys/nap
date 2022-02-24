module BMH
!-----------------------------------------------------------------------
!                     Last modified: <2022-02-11 09:41:52 KOBAYASHI Ryo>
!-----------------------------------------------------------------------
!  Parallel implementation of fitpot Born-Mayer-Huggins (BMH) potential.
!  This potential should be used with Coulomb potential.
!  
!  Potential form (param_type==pair) is defined as:
!    V(rij) = fij*bij *exp((aij-rij)/bij) -c6ij/rij**6 -c8ij/rij**8
!  and for (param_type==species) as:
!    V(rij) = f0*(bi+bj) *exp((ai+aj-rij)/(bi+bj)) -c6i*c6j/rij**6 -c8i*c8j/rij**8
!-----------------------------------------------------------------------
  use pmdvars, only: nspmax
  use util,only: csp2isp, num_data, itotOf
  implicit none
  save
  include "./const.h"
  character(len=128):: paramsdir = '.'
  character(len=128),parameter:: cprmfname = 'in.params.BMH'

  integer,parameter:: ioprms = 20

!.....Max number of species available in the potential
  integer:: nspcs
  real(8):: bmh_fij(nspmax,nspmax), bmh_rc
  real(8):: bmh_aij(nspmax,nspmax),bmh_bij(nspmax,nspmax), &
       bmh_c6ij(nspmax,nspmax),bmh_c8ij(nspmax,nspmax)
  logical:: interact(nspmax,nspmax)
  character(len=12):: param_type = 'pair'

  logical:: lprmset_BMH
  
!.....params
  integer:: nprms
  real(8),allocatable:: params(:)

contains
  subroutine force_BMH(namax,natm,tag,ra,nnmax,aa,strs,h,hi &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rcin,lspr &
       ,mpi_md_world,myid,epi,epot,nismax,specorder,lstrs,iprint,l1st)
    implicit none
    include "mpif.h"
    include "./params_unit.h"
    integer,intent(in):: namax,natm,nnmax,nismax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),lspr(0:nnmax,namax),nex(3)
    integer,intent(in):: mpi_md_world,myid
    real(8),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),rcin &
         ,tag(namax),sv(3,6)
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: l1st
    character(len=3),intent(in):: specorder(nspmax)
    logical:: lstrs

    integer:: i,j,k,l,m,n,ierr,is,js,ixyz,jxyz
    real(8):: xi(3),xij(3),rij(3),dij,diji,diji2,diji6,diji8,dvdr,dij2 &
         ,drdi(3),drdj(3),x,y,z,epotl,epott,at(3),tmp,tmp2 &
         ,fij,aij,bij,c6ij,c8ij,vs2b &
         ,vrc,dvdrc,dvs2b,rc
    real(8):: vs2bc,dvs2bc,rc6,rc8
    real(8),save:: vrcs(nspmax,nspmax),dvdrcs(nspmax,nspmax)
    real(8),allocatable,save:: strsl(:,:,:)
    
    real(8),save:: rcmax2

    rc = rcin
    if( bmh_rc.gt.0d0 ) rc = bmh_rc
    
    if( l1st ) then
      if( allocated(strsl) ) deallocate(strsl)
      allocate(strsl(3,3,namax))
      rcmax2 = rc*rc
!.....Initialize smooth cutoff
      vrcs(:,:) = 0d0
      dvdrcs(:,:)= 0d0
      rc6 = rc**6
      rc8 = rc6 *rc**2
      do is=1,nspmax
        do js=is,nspmax
          if( .not. interact(is,js) ) cycle
          fij= bmh_fij(is,js)
          aij= bmh_aij(is,js)
          bij= bmh_bij(is,js)
          c6ij= bmh_c6ij(is,js)
          c8ij= bmh_c8ij(is,js)
          vs2bc = fij*bij *exp((aij-rc)/bij) -c6ij/rc6 -c8ij/rc8
          dvs2bc= -fij *exp((aij-rc)/bij) +6d0*c6ij/rc6/rc +8d0*c8ij/rc8/rc
          vrcs(is,js) = vs2bc
          vrcs(js,is) = vs2bc
          dvdrcs(is,js) = dvs2bc
          dvdrcs(js,is) = dvs2bc
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
        fij = bmh_fij(is,js)
        aij = bmh_aij(is,js)
        bij = bmh_bij(is,js)
        c6ij= bmh_c6ij(is,js)
        c8ij= bmh_c8ij(is,js)
!.....Potential
        diji2 = diji*diji
        diji6 = diji2*diji2*diji2
        diji8 = diji6*diji2
        vs2b = fij*bij*dexp((aij-dij)/bij) -c6ij*diji6 -c8ij*diji8
        tmp2 = 0.5d0 *(vs2b -vrc -dvdrc*(dij-rc))
        if(j.le.natm) then
          epi(i)=epi(i) +tmp2
          epi(j)=epi(j) +tmp2
          epotl= epotl +tmp2 +tmp2
        else
          epi(i)=epi(i) +tmp2
          epotl= epotl +tmp2
        endif
!.....Force
        dvs2b = -fij *dexp((aij-dij)/bij) +6d0*c6ij*diji6*diji +8d0*c8ij*diji8*diji
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
!!$      call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
!!$           ,nn,mpi_md_world,strsl,9)
      strs(1:3,1:3,1:natm)= strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)
    endif

!!$    print *,'strs BMH:'
!!$    print *,' 1:  ',strsl(1,1,1),strsl(2,2,1),strsl(3,3,1)
!!$    print *,'65:  ',strsl(1,1,65),strsl(2,2,65),strsl(3,3,65)

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
    bmh_bij(:,:)= 0d0

    inc = 0
    do i=1,nspmax
      do j=i,nspmax
        if( .not.interact(i,j) ) cycle
        inc= inc +1
        bmh_aij(i,j) = params(inc)
        bmh_aij(j,i) = bmh_aij(i,j)
        inc= inc +1
        bmh_bij(i,j) = params(inc)
        bmh_bij(j,i) = bmh_bij(i,j)
      enddo
    enddo

    return
  end subroutine set_params_BMH
!=======================================================================
  subroutine read_params_BMH(myid_md,mpi_md_world,iprint,specorder)
    implicit none
    include 'mpif.h'
    include './params_unit.h'
    integer,intent(in):: myid_md,mpi_md_world,iprint
    character(len=3),intent(in):: specorder(nspmax)

    integer:: isp,jsp,ierr,ni,nj,nd
    real(8):: aij,bij,c6ij,c8ij,fij,rc
    character(len=3):: cspi,cspj
    character(len=128):: cline,cfname,mode,ctmp

    if( myid_md.eq.0 ) then
!.....Initialize params
      bmh_rc = -1d0    ! apply global rcut if this is negative
      bmh_aij(:,:)= 0d0
      bmh_bij(:,:)= 0d0
      bmh_c6ij(:,:)= 0d0
      bmh_c8ij(:,:)= 0d0
      bmh_fij(:,:)= 0d0
      interact(:,:) = .false.
      cfname = trim(paramsdir)//'/'//trim(cprmfname)

!.....Detect parameter type: pair or species
      open(ioprms,file=trim(cfname),status='old')
      do while(.true.)
        read(ioprms,'(a)',end=1) cline
        nd = num_data(cline,' ')
        if( nd.eq.0 ) cycle
        if( cline(1:1).eq.'!' .or. cline(1:1).eq.'#' ) cycle
        if( index(cline,'param_type').eq.0 ) cycle
        if( nd.lt.2 ) then
          print *,'ERROR: entry for param_type should have at least 1 arugment: ' &
               //' pair or species.'
          print *,'  e.g.)  param_type  pair'
          stop 1
        endif
        backspace(ioprms)
        read(ioprms,*) ctmp, param_type
        exit
      enddo
1     close(ioprms)

      if( iprint.ge.ipl_basic ) then
        write(6,'(/,a)') ' BMH parameters:'
        write(6,'(a,a)') '   param_type = ',trim(param_type)
      endif

      open(ioprms,file=trim(cfname),status='old')
      if( param_type(1:4).eq.'pair' ) then
        do while(.true.)
          read(ioprms,'(a)',end=10) cline
          nd = num_data(cline,' ')
          if( nd.eq.0 ) cycle
          if( cline(1:1).eq.'#' .or. cline(1:1).eq.'!' ) cycle
          if( nd.eq.7 ) then
            read(cline,*) cspi, cspj, aij,bij,c6ij,c8ij,fij
            isp = csp2isp(cspi)
            jsp = csp2isp(cspj)
            if( isp.le.0 .or. jsp.le.0 ) then
              if( iprint.ge.ipl_warn ) then
                write(6,*) ' Warning @read_params: since isp/jsp <= 0,'&
                     //' skip reading the line.'
              endif
              cycle
            endif
            bmh_fij(isp,jsp)= fij
            bmh_aij(isp,jsp)= aij
            bmh_bij(isp,jsp)= bij
            bmh_c6ij(isp,jsp)= c6ij
            bmh_c8ij(isp,jsp)= c8ij
            interact(isp,jsp)= .true.
!.....Symmetrize
            bmh_fij(jsp,isp)= fij
            bmh_aij(jsp,isp)= aij
            bmh_bij(jsp,isp)= bij
            bmh_c6ij(jsp,isp)= c6ij
            bmh_c8ij(jsp,isp)= c8ij
            interact(jsp,isp)= .true.
            if( iprint.ge.ipl_basic ) then
              print '(a,2a4,2f8.4,3es12.4)','   cspi,cspj,aij,bij,c6ij,c8ij,fij=',&
                   trim(cspi),trim(cspj),aij,bij,c6ij,c8ij,fij
            endif
          else if( index(cline,'cutoff').ne.0 ) then
            read(cline,*) ctmp, bmh_rc
          endif
        enddo  ! while reading...

      else if( param_type(1:4).eq.'spec' ) then
        mode = 'none'
        do while(.true.)
          read(ioprms,'(a)',end=10) cline
          nd = num_data(cline,' ')
          if( nd.eq.0 ) cycle
          if( cline(1:1).eq.'#' .or. cline(1:1).eq.'!' ) cycle
          if( nd.eq.5 ) then
            mode = 'none'
            read(cline,*) cspi, aij,bij,c6ij,c8ij
            isp = csp2isp(cspi)
            if( isp.le.0 ) then
              if( iprint.ge.ipl_warn ) then
                write(6,*) ' Warning @read_params: since isp <= 0,'&
                     //' skip reading the line.'
              endif
              cycle
            endif
            bmh_aij(isp,isp)= aij
            bmh_bij(isp,isp)= bij
            bmh_c6ij(isp,isp)= c6ij
            bmh_c8ij(isp,isp)= c8ij
            if( iprint.ge.ipl_basic ) then
              print '(a,a4,2f8.5,2es12.4)','   cspi,ai,bi,c6i,c8i=',&
                   trim(cspi),aij,bij,c6ij,c8ij
            endif
          else if( index(cline,'f0').ne.0 .and. nd.ge.2 ) then
            mode = 'none'
            read(cline,*) ctmp, fij
            if( iprint.ge.ipl_basic ) then
              print '(a,es12.4)', '   f0= ',fij
            endif
            bmh_fij(:,:) = fij
          else if( index(cline,'interact').ne.0 ) then
            mode = 'interact'
          else if( mode.eq.'interact' .and. nd.eq.2 ) then
            read(cline,*) cspi,cspj
            isp = csp2isp(cspi)
            jsp = csp2isp(cspj)
            if( isp.le.0 .or. isp.gt.nspmax .or. jsp.le.0 .or. jsp.gt.nspmax ) then
              print *,' WARNING: cspi,cspj,isp,jsp not available:', &
                   cspi,cspj,isp,jsp
            endif
            interact(isp,jsp) = .true.
            interact(jsp,isp) = .true.
          else
            mode = 'none'
          endif
        enddo  ! while reading...
      else
        print *,'ERROR: no such param_type available '//trim(param_type)
        stop 1
      endif
10    close(ioprms)
      if( iprint.ge.ipl_basic ) then
        if( bmh_rc.gt.0d0 ) print '(a,f7.3)','  cutoff radius = ',bmh_rc
      endif
    endif  ! myid_md.eq.0

    call mpi_bcast(bmh_rc,1,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(bmh_fij,nspmax*nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(bmh_aij,nspmax*nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(bmh_bij,nspmax*nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(bmh_c6ij,nspmax*nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(bmh_c8ij,nspmax*nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(interact,nspmax*nspmax,mpi_logical,0,mpi_md_world,ierr)

  end subroutine read_params_BMH
end module BMH
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
