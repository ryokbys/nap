module SW_Si
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

contains
  subroutine force_SW_Si(namax,natm,tag,ra,nnmax,aa,strs,h,hi,tcom &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_world,myid,epi,epot,nismax,acon,lstrs,iprint)
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
    real(8),intent(in):: ra(3,namax),tag(namax),acon(nismax) &
         ,h(3,3),hi(3,3),sv(3,6),rc
    real(8),intent(inout):: tcom
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical:: lstrs

!-----local
    integer:: i,j,k,l,m,n,ixyz,jxyz,is,js,ks,ierr,nbl
    real(8):: rij,rik,riji,riki,rij2,rik2,rc2,src,src2
    real(8):: tmp,tmp1(3),tmp2(3),vexp,df2,csn,tcsn,tcsn2,dhrij,dhrik &
         ,dhcsn,vol,voli,volj,volk,drij(3)
    real(8):: drik(3),dcsni(3),dcsnj(3),dcsnk(3),drijc,drikc,x,y,z,bl
    real(8):: epotl,epotl2,epotl3,epott
    real(8),save:: swli,a8d3r3
    real(8),save,allocatable:: aa2(:,:),aa3(:,:)
    real(8),save,allocatable,dimension(:):: xi,xj,xk,xij,xik,at,bli
    real(8),allocatable,save:: strsl(:,:,:)
!-----1st call
    logical,save:: l1st=.true.

!-----only at 1st call
    if( l1st ) then
      call read_params(myid,mpi_world,iprint)
      allocate(aa2(3,namax),aa3(3,namax),strsl(3,3,namax))
      allocate(xi(3),xj(3),xk(3),xij(3),xik(3),at(3),bli(namax))
!-------check rc
      if( myid.eq.0 .and. iprint.gt.0 ) then
        write(6,'(a,es12.4)') ' rc of input         =',rc
        write(6,'(a,es12.4)') ' rc of this potential=',swrc*swl
      endif
      if( rc .lt. swrc*swl ) then
!!$      if( int(rc*100d0) &
!!$           .ne.int(swrc*swl*100d0) ) then
        if( myid.eq.0 ) then
          write(6,'(1x,a)') "!!! Cutoff radius is not appropriate !!!"
          write(6,'(1x,a,es12.4)') "rc should be longer than ", swrc*swl
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

    epotl= 0d0
    epi(1:natm+nb)= 0d0
    strsl(1:3,1:3,1:natm+nb)= 0d0

!-----2 body term
    epotl2= 0d0
    aa2(1:3,1:natm+nb)=0d0
    src2 = swrc*swrc
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is= int(tag(i))
!!$      bl= 0d0
!!$      nbl= 0
      do k=1,lspr(0,i)
        j=lspr(k,i)
        if(j.eq.0) exit
        if(j.le.i) cycle
        js= int(tag(j))
        xj= ra(1:3,j)-xi(1:3)
        xij(1:3)= ( h(1:3,1)*xj(1) +h(1:3,2)*xj(2) &
             +h(1:3,3)*xj(3) )*swli
        rij2= xij(1)*xij(1) +xij(2)*xij(2) +xij(3)*xij(3)
        src= swrc
        if( rij2.ge.src2 ) cycle
        rij= dsqrt(rij2)
!!$        nbl=nbl +1
!!$        bl=bl +rij
        riji= 1d0/rij
        drijc= 1d0/(rij-src)
        vexp=exp(swc*drijc)
!---------potential
        tmp= 0.5d0*swe *swa*vexp*(swb*riji**swp -riji**swq)
        epi(i)= epi(i) +tmp
!!$        if( i.eq.1 ) print *,'i,j,rij,tmp,epi=',i,j,rij,tmp,epi(i)
        epotl2= epotl2 +tmp
        if( j.le.natm ) then
          epi(j)= epi(j) +tmp
          epotl2= epotl2 +tmp
        endif
!---------force
        df2= -swe*swa*vexp*(swp*swb*(riji**(swp+1d0)) &
             -swq*(riji**(swq+1d0)) &
             +(swb*(riji**swp) -riji**swq)*swc*drijc*drijc)
        drij(1:3) = -xij(1:3)*riji*swli
        aa2(1:3,i)= aa2(1:3,i) -df2*drij(1:3)
        aa2(1:3,j)= aa2(1:3,j) +df2*drij(1:3)
!-----------Stress
        if( j.le.natm ) then
          do jxyz=1,3
            strsl(1:3,jxyz,i)= strsl(1:3,jxyz,i) &
                 -0.5d0*xij(jxyz)*swl*(-df2*drij(1:3))
            strsl(1:3,jxyz,j)= strsl(1:3,jxyz,j) &
                 -0.5d0*xij(jxyz)*swl*(-df2*drij(1:3))
          enddo
        else
          do jxyz=1,3
            strsl(1:3,jxyz,i)= strsl(1:3,jxyz,i) &
                 -0.5d0*xij(jxyz)*swl*(-df2*drij(1:3))
          enddo
        endif

!!$        if( i.eq.1 ) then
!!$          print '(a,2i5,7es12.4)','i,j,rij,aa2,strs= ',i,j,rij,aa2(1:3,i)&
!!$               ,strsl(1,1,i),strsl(2,2,i),strsl(3,3,i)
!!$        endif
      enddo
!!$      write(6,'(i6,9f10.3)') i,strs(1:3,1:3,i)
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
        xj(1:3)= ra(1:3,j) -xi(1:3)
        xij(1:3)= ( h(1:3,1)*xj(1) +h(1:3,2)*xj(2) &
             +h(1:3,3)*xj(3) )*swli
        rij2= xij(1)*xij(1) +xij(2)*xij(2) +xij(3)*xij(3)
        src= swrc
!!$        src2= src*src
        if( rij2.ge.src2 ) cycle
        rij= dsqrt(rij2)
        riji= 1d0/rij
        drijc= 1d0/(rij-src)
!!$        vol= a8d3r3*rij**3
!!$        volj= 1d0/vol
!---------atom (k)
        do m=1,lspr(0,i)
          k=lspr(m,i)
          if(k.eq.0) exit
          if( k.le.j .or. k.eq.i ) cycle
          ks= int(tag(k))
          xk(1:3)= ra(1:3,k) -xi(1:3)
          xik(1:3)= ( h(1:3,1)*xk(1) +h(1:3,2)*xk(2) &
               +h(1:3,3)*xk(3) )*swli
          rik2= xik(1)*xik(1)+xik(2)*xik(2)+xik(3)*xik(3)
          src= swrc
!!$          src2= src*src
          if( rik2.ge.src2 ) cycle
          rik=dsqrt(rik2)
          riki= 1d0/rik
          drikc= 1d0/(rik-src)
!!$          vol= a8d3r3*rik**3
!!$          volk= 1d0/vol
!-----------common term
          csn=(xij(1)*xik(1) +xij(2)*xik(2) +xij(3)*xik(3)) &
               * (riji*riki)
          tcsn = csn +1d0/3d0
          tcsn2= tcsn*tcsn
          vexp= dexp(swt*drijc +swt*drikc)
!-----------potential
          tmp= swe *sws *vexp *tcsn2
          epi(i)= epi(i) +tmp
          epotl3= epotl3 +tmp
!-----------force
          dhrij= -sws *swt *vexp *tcsn2 *drijc*drijc
          dhrik= -sws *swt *vexp *tcsn2 *drikc*drikc
          dhcsn= 2d0 *sws *vexp *tcsn 
          drij(1:3)= -xij(1:3)*riji*swli
          drik(1:3)= -xik(1:3)*riki*swli
          dcsnj(1:3)= (-xij(1:3)*csn*(riji*riji) +xik(1:3)*(riji*riki))*swli
          dcsnk(1:3)= (-xik(1:3)*csn*(riki*riki) +xij(1:3)*(riji*riki))*swli
          dcsni(1:3)= -dcsnj(1:3) -dcsnk(1:3)
!!$          aa3(1:3,i)=aa3(1:3,i) -swe*(dhcsn*dcsni(1:3) +dhrij*drij(1:3) &
!!$                 +dhrik*drik(1:3))

          tmp1(1:3)= swe*(dhcsn*dcsnj(1:3) +dhrij*(-drij(1:3)))
          tmp2(1:3)= swe*(dhcsn*dcsnk(1:3) +dhrik*(-drik(1:3)))
          aa3(1:3,i)= aa3(1:3,i) +(tmp1(1:3)+tmp2(1:3))
          aa3(1:3,j)= aa3(1:3,j) -tmp1(1:3)
          aa3(1:3,k)= aa3(1:3,k) -tmp2(1:3)
!-------------Stress
          do jxyz=1,3
            strsl(1:3,jxyz,i)=strsl(1:3,jxyz,i) &
                 -0.5d0*xij(jxyz)*swl*tmp1(1:3) & !*volj &
                 -0.5d0*xik(jxyz)*swl*tmp2(1:3) !*volk
            strsl(1:3,jxyz,j)=strsl(1:3,jxyz,j) &
                 -0.5d0*xij(jxyz)*swl*tmp1(1:3) !*volj
            strsl(1:3,jxyz,k)=strsl(1:3,jxyz,k) &
                 -0.5d0*xik(jxyz)*swl*tmp2(1:3) !*volk
          enddo

        enddo
      enddo
    enddo

    call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,strsl,9)
    call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,aa3,3)
    call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,epi,1)

!-----sum
    aa(1:3,1:natm)= aa2(1:3,1:natm) +aa3(1:3,1:natm)
!!$    aa(1:3,1:natm)= aa3(1:3,1:natm)
    strs(1:3,1:3,1:natm) = strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)

!-----gather epot
    epotl= epotl2 +epotl3
!!$    epotl= epotl3
    call mpi_allreduce(epotl,epott,1,MPI_DOUBLE_PRECISION &
         ,MPI_SUM,mpi_world,ierr)
    epot= epot +epott
    return
  end subroutine force_SW_Si
!=======================================================================
  subroutine read_params(myid,mpi_world,iprint)
    implicit none
    include 'mpif.h'

    integer,intent(in):: myid,mpi_world,iprint
    integer:: itmp,ierr
    real(8):: rctmp
    logical:: lexist

!.....read parameters at the 1st call
    inquire(file='in.params.SW_Si',exist=lexist)
    if( .not. lexist ) then
      if( myid.eq.0 .and. iprint.gt.0 ) then
        write(6,'(a)') ' WARNING: in.params.SW_Si does not exist !!!.'
        write(6,'(a)') '           Default parameters will be used.'
      endif
      return
    endif
    if( myid.eq.0 ) then
      open(50,file='in.params.SW_Si',status='old')
      read(50,*) itmp,rctmp
      if( itmp.ne.nprms ) then
        write(6,'(a)') ' [Error] itmp.ne.nprms'
        write(6,'(a,i3)') '  itmp =',itmp
        stop
      endif
      read(50,*) swe
      read(50,*) swl
      read(50,*) swa
      read(50,*) swb
      read(50,*) swp
      read(50,*) swq
      read(50,*) swc
      read(50,*) swrc
      read(50,*) sws
      read(50,*) swt
      close(50)
    endif

    call mpi_bcast(swe,1,mpi_double_precision,0,mpi_world,ierr)
    call mpi_bcast(swl,1,mpi_double_precision,0,mpi_world,ierr)
    call mpi_bcast(swa,1,mpi_double_precision,0,mpi_world,ierr)
    call mpi_bcast(swb,1,mpi_double_precision,0,mpi_world,ierr)
    call mpi_bcast(swp,1,mpi_double_precision,0,mpi_world,ierr)
    call mpi_bcast(swq,1,mpi_double_precision,0,mpi_world,ierr)
    call mpi_bcast(swc,1,mpi_double_precision,0,mpi_world,ierr)
    call mpi_bcast(swrc,1,mpi_double_precision,0,mpi_world,ierr)
    call mpi_bcast(sws,1,mpi_double_precision,0,mpi_world,ierr)
    call mpi_bcast(swt,1,mpi_double_precision,0,mpi_world,ierr)
    return
  end subroutine read_params
end module SW_Si
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
