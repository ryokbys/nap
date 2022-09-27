module Pellenq
!-----------------------------------------------------------------------
!                     Last modified: <2022-09-27 13:42:11 KOBAYASHI Ryo>
!-----------------------------------------------------------------------
!  Parallel implementation of the  pontential by Pellenq and Nicholson.
!    - Pellenq & Nicholson, J. Phys. Chem. 98, 13339–13349 (1994)
!-----------------------------------------------------------------------
  use pmdvars,only: nspmax
  use util,only: csp2isp
  use memory,only: accum_mem
  use vector,only: dot
  implicit none
  include "./const.h"
  save

  character(len=128):: paramsdir = '.'
  character(len=128),parameter:: paramsfname = 'in.params.Pellenq'

  integer,parameter:: ioprms = 20

!.....Max number of species available in this potential
!!$  integer,parameter:: nspmax = 9
  integer:: nsp
  logical:: interact(nspmax,nspmax)

  character(len=10):: cprmtype = 'element'
  real(8):: Aij(nspmax,nspmax),rhoij(nspmax,nspmax)
  real(8):: c6ij(nspmax,nspmax),c8ij(nspmax,nspmax),c10ij(nspmax,nspmax)
  real(8):: rcij(nspmax,nspmax)

!.....Smooth cutoff
  real(8):: vrcs(nspmax,nspmax), dvdrcs(nspmax,nspmax)

  real(8),allocatable:: strsl(:,:,:)

  logical:: lprmset_Pellenq = .false.

contains
  subroutine read_params_Pellenq(myid_md,mpi_md_world,iprint,specorder)
!
!  Read pair parameters for Pellenq potential from file
!
    include 'mpif.h'
    integer,intent(in):: myid_md,mpi_md_world,iprint
    character(len=3),intent(in):: specorder(nspmax)
    integer:: i,j,isp,jsp,id,ierr
    character(len=128):: cline,fname
    character(len=3):: cspi,cspj
    real(8):: Ai,rhoi,c6i,c8i,c10i,rci

    if( myid_md.eq.0 ) then
      fname = trim(paramsdir)//'/'//trim(paramsfname)
      open(ioprms,file=trim(fname),status='old')
      interact(1:nspmax,1:nspmax) = .false.
      Aij(:,:)= 0d0
      rhoij(:,:)= 0d0
      c6ij(:,:)= 0d0
      c8ij(:,:)= 0d0
      c10ij(:,:)= 0d0
      rcij(:,:)= 0d0
      cprmtype = 'element'
      if( iprint.ne.0 ) write(6,'(/,a)') ' Pellenq parameters:'
      do while(.true.)
        read(ioprms,'(a)',end=10) cline
        if( cline(1:1).eq.'#' .or. cline(1:1).eq.'!' ) then
          call parse_option(cline)
          cycle
        endif
        backspace(ioprms)
        if( trim(cprmtype) == 'pair' ) then
          read(ioprms,*,end=10) cspi,cspj,Ai,rhoi,c6i,c8i,c10i,rci
          isp = csp2isp(cspi)
          jsp = csp2isp(cspj)
          if( isp.gt.0 .and. jsp.gt.0 ) then
            Aij(isp,jsp) = Ai
            rhoij(isp,jsp) = rhoi
            c6ij(isp,jsp) = c6i
            c8ij(isp,jsp) = c8i
            c10ij(isp,jsp) = c10i
            rcij(isp,jsp) = rci
            interact(isp,jsp) = .true.
            if( iprint.ge.ipl_basic ) then
              write(6,'(a,2a4,6f10.3)') '   cspi,cspj,Aij,rhoij,c6ij,c8ij,c10ij,rcij = ', &
                   trim(cspi),trim(cspj),Ai,rhoi,c6i,c8i,c10i,rci
            endif
          else
            if( iprint.ge.ipl_info ) then
              print *,' Pellenq parameter read but not used: cspi,cspj=', &
                   trim(cspi),trim(cspj)
            endif
          endif
        else  ! cprmtype not 'pair', then it must be 'element' (default)
          read(ioprms,*,end=10) cspi,Ai,rhoi,c6i,c8i,c10i,rci
          isp = csp2isp(cspi)
          if( isp.gt.0 ) then
            Aij(isp,isp) = Ai
            rhoij(isp,isp) = rhoi
            c6ij(isp,isp) = c6i
            c8ij(isp,isp) = c8i
            c10ij(isp,isp) = c10i
            rcij(isp,isp) = rci
            interact(isp,isp) = .true.
            if( iprint.ge.ipl_basic ) then
              write(6,'(a,a4,6f10.3)') '   cspi,Ai,rhoi,c6i,c8i,c10i,rci = ', &
                   trim(cspi),Ai,rhoi,c6i,c8i,c10i,rci
            endif
          else
            if( iprint.ge.ipl_info ) then
              print *,' Pellenq parameter read but not used: cspi=',cspi
            endif
          endif
        endif
      enddo
10    close(ioprms)

!.....Mixing parameters from elemental ones
      if( trim(cprmtype).eq.'element' ) then
        do isp=1,nspmax-1
          if( .not.interact(isp,isp) ) cycle
          do jsp=isp+1,nspmax
            if( .not.interact(jsp,jsp) ) cycle
            interact(isp,jsp) = .true.
            Aij(isp,jsp) = sqrt(Aij(isp,isp) *Aij(jsp,jsp))
!.....Mixing lengths by averaging
            rhoij(isp,jsp)= (rhoij(isp,isp) +rhoij(jsp,jsp))/2
            rcij(isp,jsp) = (rcij(isp,isp) +rcij(jsp,jsp))/2
!.....Mixiing c6,c8,c10 by just sqrt, which may not be appropriate...
            c6ij(isp,jsp) = sqrt(c6ij(isp,isp) *c6ij(jsp,jsp))
            c8ij(isp,jsp) = sqrt(c8ij(isp,isp) *c8ij(jsp,jsp))
            c10ij(isp,jsp)= sqrt(c10ij(isp,isp) *c10ij(jsp,jsp))
            if( iprint.ge.ipl_basic ) then
              cspi = specorder(isp)
              cspj = specorder(jsp)
              write(6,'(a,2a4,6f10.3)') &
                   '   cspi,cspj,Ai,rhoi,c6i,c8i,c10i,rci = ', &
                   trim(cspi),trim(cspj),Aij(isp,jsp),rhoij(isp,jsp),&
                   c6ij(isp,jsp),c8ij(isp,jsp),c10ij(isp,jsp), &
                   rcij(isp,jsp)
            endif
!.....Symmetrize
            interact(jsp,isp) = .true.
            Aij(jsp,isp) = Aij(isp,jsp)
            rhoij(jsp,isp)= rhoij(isp,jsp)
            rcij(jsp,isp)= rcij(isp,jsp)
            c6ij(jsp,isp) = c6ij(isp,jsp)
            c8ij(jsp,isp) = c8ij(isp,jsp)
            c10ij(jsp,isp)= c10ij(isp,jsp)
          enddo
        enddo
      else  ! if cprmtype == 'pair', only symmetrize
        do isp=1,nspmax-1
          do jsp=isp+1,nspmax
            if( .not. (interact(isp,jsp).or.interact(jsp,isp)) ) cycle
            if( interact(isp,jsp) ) then
              interact(jsp,isp) = .true.
              Aij(jsp,isp) = Aij(isp,jsp)
              rhoij(jsp,isp)= rhoij(isp,jsp)
              rcij(jsp,isp)= rcij(isp,jsp)
              c6ij(jsp,isp) = c6ij(isp,jsp)
              c8ij(jsp,isp) = c8ij(isp,jsp)
              c10ij(jsp,isp)= c10ij(isp,jsp)
            else
              interact(isp,jsp) = .true.
              Aij(isp,jsp) = Aij(jsp,isp)
              rhoij(isp,jsp)= rhoij(jsp,isp)
              rcij(isp,jsp)= rcij(jsp,isp)
              c6ij(isp,jsp) = c6ij(jsp,isp)
              c8ij(isp,jsp) = c8ij(jsp,isp)
              c10ij(isp,jsp)= c10ij(jsp,isp)
            endif
          enddo
        enddo
      endif
    endif  ! myid_md.eq.0

    call mpi_bcast(Aij,nspmax*nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(rhoij,nspmax*nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(rcij,nspmax*nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(c6ij,nspmax*nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(c8ij,nspmax*nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(c10ij,nspmax*nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(interact,nspmax*nspmax,mpi_logical,0,mpi_md_world,ierr)

    if( iprint.ge.ipl_debug .and. myid_md.eq.0 ) then
      write(6,'(a)') ' Finished reading '//trim(fname)
      write(6,*) ''
    endif
  end subroutine read_params_Pellenq
!=======================================================================
  subroutine parse_option(cline)
!
!  Parse options from a comment line.
!  Lines starting from ! or # are treated as comment lines,
!  and options can be given at the comment lines.
!  The option words should be put after these comment characters with
!  one or more spaces between them for example,
!
!  param_type: element
!
!  Currently available options are:
!    - "param_type:", Parameter type, "element" (default) or "pair"
!
    use util, only: num_data
    include "./const.h"
    character(len=*),intent(in):: cline

    integer:: num
    character(len=10):: c1,copt

    if( index(cline,'param_type:').ne.0 ) then
      num = num_data(trim(cline),' ')
      if( num.lt.3 ) stop 'ERROR: num of entry wrong...'
      read(cline,*) c1, copt, cprmtype
      print *,'  param_type: ',trim(cprmtype)
    endif
    
  end subroutine parse_option
!=======================================================================
  subroutine force_Pellenq(namax,natm,tag,ra,nnmax,aa,strs,h,hi &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc_global,lspr,d2lspr &
       ,mpi_md_world,myid,epi,epot,nismax,lstrs,iprint,l1st)
    use util,only: itotOf
    implicit none
    include "mpif.h"
    include "./params_unit.h"
    integer,intent(in):: namax,natm,nnmax,nismax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),lspr(0:nnmax,namax),nex(3)
    integer,intent(in):: mpi_md_world,myid
    real(8),intent(in):: ra(3,namax),h(3,3),hi(3,3),rc_global &
         ,tag(namax),sv(3,6),d2lspr(nnmax,namax)
    real(8),intent(inout):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: l1st
    logical:: lstrs

    integer:: i,j,k,l,m,n,ierr,is,js,ixyz,jxyz
    real(8):: xi(3),xj(3),xij(3),rij(3),dij,diji,dedr,epott, &
         dxdi(3),dxdj(3),x,y,z,epotl,at(3),tmp,tmp2, &
         dij2,vrc,dvdrc,expbrc,expbr,A,rho,c6,c8,c10, &
         r2,r10,r10i,r8i,r6i,f6,f8,f10,df6,df8,df10,rc
    real(8),save:: rc2
    real(8),external:: fcut1,dfcut1

    if( l1st ) then
      if( allocated(strsl) ) then
        call accum_mem('force_Pellenq',-8*size(strsl))
        deallocate(strsl)
      endif
      allocate(strsl(3,3,namax))
      call accum_mem('force_Pellenq',8*size(strsl))
      rc2 = -1d0
!.....Initialize smooth cutoff
      vrcs(:,:) = 0d0
      dvdrcs(:,:) = 0d0
      do is=1,nspmax
        do js=is,nspmax
          if( .not.(interact(is,js).or.interact(js,is)) ) cycle
          rc = rcij(is,js)
          rc2 = max(rc2,rc**2)
          call compute_f2n(rc,is,js,f6,f8,f10,df6,df8,df10)
          expbrc = exp(-rc/rhoij(is,js))
          vrc = Aij(is,js)*expbrc &
               -f6*c6ij(is,js)/rc**6 &
               -f8*c8ij(is,js)/rc**8 &
               -f10*c10ij(is,js)/rc**10
!!$          vrc = Aij(is,js)*expbrc
          vrcs(is,js) = vrc
          vrcs(js,is) = vrc
          dvdrc= -Aij(is,js)*expbrc/rhoij(is,js) &
               -c6ij(is,js)/rc**6 *(df6 -6d0*f6/rc) &
               -c8ij(is,js)/rc**8 *(df8 -8d0*f8/rc) &
               -c10ij(is,js)/rc**10 *(df10 -10d0*f10/rc)
!!$          dvdrc= -Aij(is,js)*expbrc/rhoij(is,js)
          dvdrcs(is,js) = dvdrc
          dvdrcs(js,is) = dvdrc
        enddo
      enddo
      if( rc2.gt.rc_global**2 ) then
        if( myid.eq.0 ) print *, &
             'ERROR: rcmax of Pellenq pot is greater than global rc'
        stop 1
      endif
    endif

    if( .not.allocated(strsl) ) then
      allocate(strsl(3,3,namax))
      call accum_mem('force_Pellenq',8*size(strsl))
    else if( size(strsl).lt.3*3*namax ) then
      call accum_mem('force_Pellenq',-8*size(strsl))
      deallocate(strsl)
      allocate(strsl(3,3,namax))
      call accum_mem('force_Pellenq',8*size(strsl))
    endif

    epotl= 0d0
    strsl(1:3,1:3,1:namax) = 0d0

!.....Loop over resident atoms
!$omp parallel
!$omp do private(i,xi,is,k,j,js,xj,xij,rij,dij2,dij,diji,dxdi, &
!$omp     vrc,dvdrc,tmp,tmp2,dedr,ixyz,jxyz,A,rho,c6,c8,c10, &
!$omp     r2,r10,r10i,r8i,r6i,expbr,f6,f8,f10,df6,df8,df10) &
!$omp     reduction(+:epotl)
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is=int(tag(i))
      do k=1,lspr(0,i)
        if( d2lspr(k,i).ge.rc2 ) cycle
        j=lspr(k,i)
!!$        if( j.le.i ) cycle
        js= int(tag(j))
!.....Check if these two species interact
        if( .not. interact(is,js) ) cycle
        xj(1:3)= ra(1:3,j)
        xij(1:3)= xj(1:3)-xi(1:3)
        rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij2 = rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3)
!!$        if( dij2.gt.rc2 ) cycle
        dij= sqrt(dij2)
        rc = rcij(is,js)
        if( dij.ge.rc ) cycle
        diji= 1d0/dij
        dxdi(1:3)= -rij(1:3)*diji
        A = Aij(is,js)
        rho= rhoij(is,js)
        c6 = c6ij(is,js)
        c8 = c8ij(is,js)
        c10= c10ij(is,js)
        vrc = vrcs(is,js)
        dvdrc = dvdrcs(is,js)
        call compute_f2n(dij,is,js,f6,f8,f10,df6,df8,df10)
        r2 = dij**2
        r10 = r2**5
        r10i = 1d0/r10
        r8i = r10i *r2
        r6i = r8i *r2
        expbr = exp(-dij/rho)
!.....potential
        tmp= A*expbr -f6*c6*r6i -f8*c8*r8i -f10*c10*r10i
        tmp2 = 0.5d0 *(tmp -vrc -dvdrc*(dij-rc))
        epi(i)= epi(i) +tmp2
        epotl= epotl +tmp2
!!$        if( j.le.natm ) then
!!$          epi(i)= epi(i) +tmp2
!!$          epotl= epotl +tmp2 +tmp2
!!$!$omp atomic
!!$          epi(j)= epi(j) +tmp2
!!$        else
!!$          epi(i)= epi(i) +tmp2
!!$          epotl= epotl +tmp2
!!$        endif
!.....force
        dedr= -A*expbr/rho &
             -c6*r6i*(df6 -6d0*f6/dij) &
             -c8*r8i*(df8 -8d0*f8/dij) &
             -c10*r10i*(df10 -10d0*f10/dij) &
             -dvdrc
        aa(1:3,i)= aa(1:3,i) -dxdi(1:3)*dedr
!!$        do ixyz=1,3
!!$!$omp atomic
!!$          aa(ixyz,i)= aa(ixyz,i) +dxdi(ixyz)*dedr
!!$        enddo
!.....stress
        do ixyz=1,3
          do jxyz=1,3
            strsl(jxyz,ixyz,i)= strsl(jxyz,ixyz,i) &
                 -0.5d0 *dedr*rij(ixyz)*(-dxdi(jxyz))
!!$!$omp atomic
!!$            strsl(jxyz,ixyz,j)= strsl(jxyz,ixyz,j) &
!!$                 -0.5d0 *dedr*rij(ixyz)*(-dxdi(jxyz))
          enddo
        enddo
      enddo
    enddo
!$omp end do
!$omp end parallel

    strs(1:3,1:3,1:natm)= strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)

!-----gather epot
    epott= 0d0
    call mpi_allreduce(epotl,epott,1,mpi_real8 &
         ,mpi_sum,mpi_md_world,ierr)
    epot= epot +epott
    if( iprint.ge.ipl_info ) print *,'epot Pellenq = ',epott
    return
  end subroutine force_Pellenq
!=======================================================================
  subroutine compute_f2n(r,is,js,f6,f8,f10,df6,df8,df10)
    real(8),intent(in):: r
    integer,intent(in):: is,js
    real(8),intent(out):: f6,f8,f10,df6,df8,df10

    integer:: k
    real(8):: br,expbr,fack,brk,ftmp,dftmp,rhoi,tmp

    rhoi = 1d0/rhoij(is,js)
    br = r *rhoi
    expbr = exp( -br)
    fack = 1d0
    brk = 1d0
    ftmp = 1d0 -1d0*expbr
    dftmp = +expbr *rhoi
    do k=1,6
      fack = fack*k
      brk = brk*br
      tmp = brk/fack *expbr
      ftmp = ftmp -tmp
      dftmp = dftmp -tmp *(dble(k)/r -rhoi)
    enddo
    f6 = ftmp
    df6 = dftmp
    do k=7,8
      fack = fack*k
      brk = brk*br
      tmp = brk/fack *expbr
      ftmp = ftmp -tmp
      dftmp = dftmp -tmp *(dble(k)/r -rhoi)
    enddo
    f8 = ftmp
    df8= dftmp
    do k=9,10
      fack = fack*k
      brk = brk*br
      tmp = brk/fack *expbr
      ftmp = ftmp -tmp
      dftmp = dftmp -tmp *(dble(k)/r -rhoi)
    enddo
    f10 = ftmp
    df10= dftmp
    return
  end subroutine compute_f2n
  
end module Pellenq
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd lib"
!     End:
