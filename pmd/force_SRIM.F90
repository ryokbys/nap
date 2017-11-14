module SRIM
!-----------------------------------------------------------------------
!                     Last modified: <2017-11-14 12:57:38 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
!  Parallel implementation of SRIM repulsive potential.
!  See Ziegler, Biersack, Littmark (ZBL) The Stopping and Range of Ions
!  in Matter (SRIM), 1985.
!-----------------------------------------------------------------------
  implicit none

  character(len=128):: paramsdir = '.'
  character(len=128),parameter:: paramsfname = 'in.params.SRIM'
  integer,parameter:: ioprms = 20
  
!.....Max num of species
  integer,parameter:: msp = 9

  logical:: interact(msp,msp)
  real(8):: qnucl(msp)
  real(8):: srim_rc

contains
!=======================================================================
  subroutine read_params_SRIM(myid,mpi_world,iprint)
    implicit none
    include "mpif.h"
    integer,intent(in):: myid,mpi_world,iprint

    integer:: isp,jsp,ierr,ndata
    real(8):: qnucli
    character(len=128):: cline,fname
    real(8),parameter:: qtiny = 1d-10

    ndata = 0
    if( myid.eq.0 ) then
      fname = trim(paramsdir)//'/'//trim(paramsfname)
      open(ioprms,file=trim(fname),status='old')
      interact(1:msp,1:msp) = .false.
      qnucl(1:msp) = 0d0

      if( iprint.ne.0 ) write(6,'(/,a)') ' SRIM parameters:'
      do while(.true.)
        read(ioprms,*,end=10) cline
        if( cline(1:1).eq.'#' .or. cline(1:1).eq.'!' ) cycle
        backspace(ioprms)
        if( ndata.eq.0 ) then
          read(ioprms,*) ndata, srim_rc
          if( iprint.ne.0 ) write(6,'(a,f0.3)') '   Cutoff radius = ',srim_rc
        else
          read(ioprms,*) isp, qnucli
          if( isp.gt.msp ) then
            write(6,*) ' Warning @read_params: since isp is greater than msp,'&
                 //' skip reading the line.'
          endif
          qnucl(isp) = qnucli
          if( iprint.ne.0 ) write(6,'(a,i4,f0.3)') '   isp,qnucl = ',isp,qnucli
        endif
      enddo
10    close(ioprms)
      do isp=1,msp
        if( abs(qnucl(isp)).lt.qtiny ) cycle
        do jsp=isp,msp
          if( abs(qnucl(jsp)).lt.qtiny ) cycle
          interact(isp,jsp) = .true.
          interact(jsp,isp) = .true.
        enddo
      enddo
    endif
    
    call mpi_bcast(srim_rc,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(qnucl,msp,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(interact,msp*msp,mpi_logical,0,mpi_world,ierr)
    return
  end subroutine read_params_SRIM
!=======================================================================
  subroutine force_SRIM(namax,natm,tag,ra,nnmax,aa,strs,h,hi,tcom &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint,l1st)
    implicit none
    include "mpif.h"
    include "./params_unit.h"
    integer,intent(in):: namax,natm,nnmax,nismax,lspr(0:nnmax,namax)&
         ,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_md_world,myid_md,nex(3)
    real(8),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),sv(3,6) &
         ,acon(nismax),rc,tag(namax)
    real(8),intent(inout):: tcom
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: l1st
    logical:: lstrs

    integer:: i,j,k,l,m,n,ierr,is,js,ixyz,jxyz
    real(8):: xij(3),rij,rcij,dfi,dfj,drdxi(3),drdxj(3),r,at(3)
    real(8):: x,y,z,xi(3),epotl,epott,tmp,dtmp,drhoi,drhoj
    real(8),allocatable,save:: rho(:)
    real(8),allocatable,save:: strsl(:,:,:)

    real(8),save:: srim_rc2

    if( l1st ) then
      if( rc.lt.srim_rc ) then
        if( myid_md.eq.0 .and. iprint.gt.0 ) then
          print '(/,a)',' Input cutoff radius is smaller than rc of SRIM potential.'
          print '(a,f0.3)', '   Input rc     = ',rc
          print '(a,f0.3)', '   Potential rc = ',srim_rc
        endif
        call mpi_finalize(ierr)
        stop
      endif
      if( allocated(strsl) ) deallocate(strsl)
      allocate(strsl(3,3,namax))
      srim_rc2 = srim_rc*srim_rc
    endif

    if( size(strsl).lt.3*3*namax ) then
      deallocate(strsl)
      allocate(strsl(3,3,namax))
    endif
    
    epotl= 0d0
    strsl(1:3,1:3,1:namax) = 0d0

!-----dE/dr_i
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is = int(tag(i))
      do k=1,lspr(0,i)
        j=lspr(k,i)
        if(j.eq.0) exit
        if(j.le.i) cycle
        js = int(tag(j))
        if( .not. interact(is,js) ) cycle
        x= ra(1,j) -xi(1)
        y= ra(2,j) -xi(2)
        z= ra(3,j) -xi(3)
        xij(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
        rij= xij(1)**2+ xij(2)**2 +xij(3)**2
        if( rij.gt.srim_rc2 ) cycle
        rij = sqrt(rij)
        drdxi(1:3)= -xij(1:3)/rij
!.....2-body term
        tmp = 0.5d0 *vnucl(is,js,rij)
        epi(i)= epi(i) +tmp
        epi(j)= epi(j) +tmp
        if(j.le.natm) then
          epotl=epotl +tmp +tmp
        else
          epotl=epotl +tmp
        endif
        dtmp = dvnucl(is,js,rij)
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
      enddo
    enddo

    if( lstrs ) then
      strs(1:3,1:3,1:natm)= strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)
    endif

!-----gather epot
    call mpi_allreduce(epotl,epott,1,mpi_real8 &
         ,mpi_sum,mpi_md_world,ierr)
    epot= epot +epott

  end subroutine force_SRIM
!=======================================================================
  function vnucl(is,js,rij)
!
!  Repulsive potential between nuclei
!
    implicit none
    integer,intent(in):: is,js
    real(8),intent(in):: rij
    real(8):: vnucl
    real(8):: rs,qi,qj

    qi = qnucl(is)
    qj = qnucl(js)
    rs = 0.4683766d0  /(qi**(2d0/3) +qj**(2d0/3))
    vnucl = qi*qj/rij *xi(rij/rs)
    return
  end function vnucl
!=======================================================================
  function dvnucl(is,js,rij)
!
!  Derivative of the repulsive potential between nuclei
!
    implicit none
    integer,intent(in):: is,js
    real(8),intent(in):: rij
    real(8):: dvnucl
    real(8):: rs,qi,qj

    qi = qnucl(is)
    qj = qnucl(js)
    rs = 0.4683766d0  /(qi**(2d0/3) +qj**(2d0/3))
    dvnucl = qi*qj/rij* ( -1d0/rij*xi(rij/rs) &
         +dxi(rij/rs)/rs )
    return
  end function dvnucl
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
  
end module SRIM
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
