module repel
!-----------------------------------------------------------------------
!                     Last modified: <2022-09-28 14:15:03 KOBAYASHI Ryo>
!-----------------------------------------------------------------------
!  Parallel implementation of repulsion pontential
!-----------------------------------------------------------------------
  use pmdvars,only: nspmax
  use util,only: csp2isp, is_numeric
  use memory,only: accum_mem
  implicit none
  include "./const.h"
  save

  character(len=128):: paramsdir = '.'
  character(len=128),parameter:: paramsfname = 'in.params.repel'

  integer,parameter:: ioprms = 20

!.....Max number of species available in this potential
  integer:: nsp
  logical:: interact(nspmax,nspmax)

  character(len=10):: ctype = 'exp'
  real(8):: Aij(nspmax,nspmax), rhoij(nspmax,nspmax), &
       sgmij(nspmax,nspmax), rcij(nspmax,nspmax)

!.....Smooth cutoff
  real(8):: vrcs(nspmax,nspmax), dvdrcs(nspmax,nspmax)

  real(8),allocatable:: strsl(:,:,:)

  logical:: lprmset_repel = .false.

contains
  subroutine read_params_repel(myid_md,mpi_md_world,iprint,specorder)
!
!  Read pair parameters for repel potential from file
!
    include 'mpif.h'
    integer,intent(in):: myid_md,mpi_md_world,iprint
    character(len=3),intent(in):: specorder(nspmax)
    integer:: i,j,isp,jsp,id,ierr
    character(len=128):: cline,fname
    character(len=3):: cspi,cspj
    real(8):: Ai,rhoi,c6i,c8i,c10i,rci,sgmi

    if( myid_md.eq.0 ) then
      fname = trim(paramsdir)//'/'//trim(paramsfname)
      open(ioprms,file=trim(fname),status='old')
      interact(1:nspmax,1:nspmax) = .false.
      Aij(:,:)= 0d0
      rhoij(:,:)= 0d0
      sgmij(:,:)= 0d0
      rcij(:,:)= 0d0
      ctype = 'exp'
      if( iprint.ne.0 ) write(6,'(/,a)') ' repel parameters:'
      do while(.true.)
        read(ioprms,'(a)',end=10) cline
        if( cline(1:1).eq.'#' .or. cline(1:1).eq.'!' ) then
          call parse_option(cline)
          cycle
        endif
        backspace(ioprms)
        if( trim(ctype) == 'exp' ) then
          read(ioprms,*,end=10) cspi,cspj,Ai,rhoi,sgmi,rci
          isp = csp2isp(cspi)
          jsp = csp2isp(cspj)
          if( isp.gt.0 .and. jsp.gt.0 ) then
            Aij(isp,jsp) = Ai
            rhoij(isp,jsp) = rhoi
            sgmij(isp,jsp) = sgmi
            rcij(isp,jsp) = rci
            interact(isp,jsp) = .true.
            if( iprint.ge.ipl_basic ) then
              write(6,'(a,2a4,6f10.3)') '   cspi,cspj,Aij,rhoij,sgmij,rcij = ', &
                   trim(cspi),trim(cspj),Ai,rhoi,sgmi,rci
            endif
          else
            if( iprint.ge.ipl_info ) then
              print *,' repel parameter read but not used: cspi,cspj=', &
                   trim(cspi),trim(cspj)
            endif
          endif  ! isp.gt.0 ...
        else
          print *,' repel type '//trim(ctype)//' not available...'
        endif  ! ctype
      enddo
10    close(ioprms)

      do isp=1,nspmax-1
        do jsp=isp+1,nspmax
          if( .not. (interact(isp,jsp).or.interact(jsp,isp)) ) cycle
          if( interact(isp,jsp) ) then
            interact(jsp,isp) = .true.
            Aij(jsp,isp) = Aij(isp,jsp)
            rhoij(jsp,isp)= rhoij(isp,jsp)
            sgmij(jsp,isp)= sgmij(isp,jsp)
            rcij(jsp,isp)= rcij(isp,jsp)
          else
            interact(isp,jsp) = .true.
            Aij(isp,jsp) = Aij(jsp,isp)
            rhoij(isp,jsp)= rhoij(jsp,isp)
            sgmij(isp,jsp)= sgmij(jsp,isp)
            rcij(isp,jsp)= rcij(jsp,isp)
          endif
        enddo
      enddo
    endif  ! myid_md.eq.0

    call mpi_bcast(Aij,nspmax*nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(rhoij,nspmax*nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(sgmij,nspmax*nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(rcij,nspmax*nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(interact,nspmax*nspmax,mpi_logical,0,mpi_md_world,ierr)

    if( iprint.ge.ipl_debug .and. myid_md.eq.0 ) then
      write(6,'(a)') ' Finished reading '//trim(fname)
      write(6,*) ''
    endif
  end subroutine read_params_repel
!=======================================================================
  subroutine parse_option(cline)
!
!  Parse options from a comment line.
!  Lines starting from ! or # are treated as comment lines,
!  and options can be given at the comment lines.
!  The option words should be put after these comment characters with
!  one or more spaces between them for example,
!
!  repel_type: exp
!
!  Currently available options are:
!    - "repel_type:", repel potential type, "exp" (default)
!
    use util, only: num_data
    include "./const.h"
    character(len=*),intent(in):: cline

    integer:: num
    character(len=10):: c1,copt

    if( index(cline,'repel_type:').ne.0 ) then
      num = num_data(trim(cline),' ')
      if( num.lt.3 ) stop 'ERROR: num of entry wrong...'
      read(cline,*) c1, copt, ctype
      print *,'  repel_type: ',trim(ctype)
    endif
    
  end subroutine parse_option
!=======================================================================
  subroutine force_repel(namax,natm,tag,ra,nnmax,aa,strs,h,hi &
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
         dij2,vrc,dvdrc,expbrc,expbr,A,rho,rc,sgm
    real(8),save:: rc2
    real(8),external:: fcut1,dfcut1

    if( l1st ) then
      if( allocated(strsl) ) then
        call accum_mem('force_repel',-8*size(strsl))
        deallocate(strsl)
      endif
      allocate(strsl(3,3,namax))
      call accum_mem('force_repel',8*size(strsl))
      rc2 = -1d0
!.....Initialize smooth cutoff
      vrcs(:,:) = 0d0
      dvdrcs(:,:) = 0d0
      do is=1,nspmax
        do js=is,nspmax
          if( .not.(interact(is,js).or.interact(js,is)) ) cycle
          rc = rcij(is,js)
          rc2 = max(rc2,rc**2)
          expbrc = exp(-(rc-sgmij(is,js))/rhoij(is,js))
          vrc = Aij(is,js)*expbrc
          vrcs(is,js) = vrc
          vrcs(js,is) = vrc
          dvdrc= -Aij(is,js)*expbrc/rhoij(is,js)
          dvdrcs(is,js) = dvdrc
          dvdrcs(js,is) = dvdrc
        enddo
      enddo
      if( rc2.gt.rc_global**2 ) then
        if( myid.eq.0 ) print *, &
             'ERROR: rcmax of repel pot is greater than global rc'
        stop 1
      endif
    endif

    if( .not.allocated(strsl) ) then
      allocate(strsl(3,3,namax))
      call accum_mem('force_repel',8*size(strsl))
    else if( size(strsl).lt.3*3*namax ) then
      call accum_mem('force_repel',-8*size(strsl))
      deallocate(strsl)
      allocate(strsl(3,3,namax))
      call accum_mem('force_repel',8*size(strsl))
    endif

    epotl= 0d0
    strsl(1:3,1:3,1:namax) = 0d0

!.....Loop over resident atoms
!$omp parallel
!$omp do private(i,xi,is,k,j,js,xj,xij,rij,dij2,dij,diji,dxdi, &
!$omp     vrc,dvdrc,tmp,tmp2,dedr,ixyz,jxyz,A,rho,expbr) &
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
        sgm= sgmij(is,js)
        vrc = vrcs(is,js)
        dvdrc = dvdrcs(is,js)
        expbr = exp(-(dij-sgm)/rho)
!.....potential
        tmp= A*expbr
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
        dedr= -A*expbr/rho -dvdrc
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
    if( iprint.ge.ipl_info ) print *,'epot repel = ',epott
    return
  end subroutine force_repel
!=======================================================================
end module repel
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd lib"
!     End:
