module dipole
!-----------------------------------------------------------------------
!                     Last modified: <2021-11-24 11:48:24 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
!  Parallel implementation of dipole-dipole and quadrupole-dipole interactions.
!-----------------------------------------------------------------------
  use pmdmpi
  use mod_precision
  use pmdvars, only: nspmax
  use util,only: csp2isp
  implicit none
  save
  character(len=128):: paramsdir = '.'
  character(len=128),parameter:: cprmfname = 'in.params.dipole'

  integer,parameter:: ioprms = 20
!.....Coulomb's constant, acc = 1.0/(4*pi*epsilon0)
  real(rp),parameter:: acc  = 14.3998554737_rp
  real(rp):: acc2 = acc*acc
!.....permittivity of vacuum
  real(rp),parameter:: eps0 = 0.00552634939836_rp  ! e^2 /Ang /eV

!.....Max number of species available in the potential
  integer:: nspcs
  integer:: dip_ni(nspmax)
  real(rp):: dip_ei(nspmax),dip_cij(nspmax,nspmax),dip_dij(nspmax,nspmax)
  logical:: interact(nspmax,nspmax)

  logical:: lprmset_dipole
  
!.....params
  integer:: nprms
  real(rp),allocatable:: params(:)

contains
  subroutine force_dipole(namax,natm,tag,ra,nnmax,aa,strs,h,hi &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid,epi,epot,nismax,specorder,lstrs,iprint,l1st)
    use util,only: itotOf
    implicit none
    include "./params_unit.h"
    integer,intent(in):: namax,natm,nnmax,nismax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),lspr(0:nnmax,namax),nex(3)
    integer,intent(in):: mpi_md_world,myid
    real(rp),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),rc &
         ,tag(namax),sv(3,6)
    real(rp),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: l1st
    character(len=3),intent(in):: specorder(nspmax)
    logical:: lstrs

    integer:: i,j,k,l,m,n,ierr,is,js,ixyz,jxyz
    real(rp):: xi(3),xij(3),rij(3),dij,diji,dvdr,dij2,dij6,dij8 &
         ,drdi(3),drdj(3),x,y,z,epotl,epott,at(3),tmp,tmp2 &
         ,c6ij,c8ij,dv6,dv8,v6,v8 &
         ,vrc,dvdrc
    real(rp):: vrc6,vrc8,dvdrc6,dvdrc8
    real(rp),save:: vrcs(nspmax,nspmax),dvdrcs(nspmax,nspmax)
    real(rp),allocatable,save:: strsl(:,:,:)
    

    real(rp),save:: rcmax2

    if( l1st ) then
      if( allocated(strsl) ) deallocate(strsl)
      allocate(strsl(3,3,namax))
      rcmax2 = rc*rc
!.....Initialize smooth cutoff
      vrcs(:,:) = 0.0_rp
      dvdrcs(:,:)= 0.0_rp
      do is=1,nspmax
        do js=is,nspmax
          if( .not. interact(is,js) ) cycle
          c6ij= dip_cij(is,js)
          c8ij= dip_dij(is,js)
          vrc6 = -acc2*c6ij/rc**6
          vrc8 = -acc2*c8ij/rc**8
          dvdrc6 = 6.0_rp*acc2*c6ij/rc**7
          dvdrc8 = 8.0_rp*acc2*c8ij/rc**9
          vrcs(is,js) = vrc6 +vrc8
          vrcs(js,is) = vrcs(is,js)
          dvdrcs(is,js) = dvdrc6 +dvdrc8
          dvdrcs(js,is) = dvdrcs(is,js)
        enddo
      enddo
    endif

    if( size(strsl).lt.3*3*namax ) then
      deallocate(strsl)
      allocate(strsl(3,3,namax))
    endif

    epotl= 0.0_rp
    strsl(1:3,1:3,1:natm+nb) = 0.0_rp

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
        dij = sqrt(dij2)
        diji= 1.0_rp/dij
        drdi(1:3)= -rij(1:3)*diji
        drdj(1:3)=  rij(1:3)*diji
        c6ij = dip_cij(is,js)
        c8ij = dip_dij(is,js)
        dij6 = dij2*dij2*dij2
        dij8 = dij6*dij2
        vrc = vrcs(is,js)
        dvdrc = dvdrcs(is,js)
!.....Potential
        v6 = -acc2*c6ij/dij6
        v8 = -acc2*c8ij/dij8
        tmp = v6 +v8
        tmp2 = 0.5_rp *(tmp -vrc -dvdrc*(dij-rc))
!!$        tmp2 = 0.5d0 *tmp
        if(j.le.natm) then
          epi(i)=epi(i) +tmp2
          epi(j)=epi(j) +tmp2
          epotl= epotl +tmp2 +tmp2
        else
          epi(i)=epi(i) +tmp2
          epotl= epotl +tmp2
        endif
!.....Force
        dv6 = 6.0_rp*acc2*c6ij/dij6*diji
        dv8 = 8.0_rp*acc2*c8ij/dij8*diji
        dvdr = dv6 +dv8 -dvdrc
        aa(1:3,i)= aa(1:3,i) -drdi(1:3)*dvdr
        aa(1:3,j)= aa(1:3,j) -drdj(1:3)*dvdr
!.....Stress
        if( .not.lstrs ) cycle
        do ixyz=1,3
          do jxyz=1,3
            strsl(jxyz,ixyz,i)=strsl(jxyz,ixyz,i) &
                 -0.5_rp*dvdr*rij(ixyz)*(-drdi(jxyz))
            strsl(jxyz,ixyz,j)=strsl(jxyz,ixyz,j) &
                 -0.5_rp*dvdr*rij(ixyz)*(-drdi(jxyz))
          enddo
        enddo
      enddo
    enddo

    if( lstrs ) then
!!$      call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
!!$           ,nn,mpi_md_world,strsl,9)
      strs(1:3,1:3,1:natm)= strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)
    endif

    epott = 0.0_rp
    call mpi_allreduce(epotl,epott,1,mpi_real_rp &
         ,mpi_sum,mpi_md_world,ierr)
    epot= epot +epott
    if( iprint.gt.2 ) print '(a,es15.7)',' epot dipole = ',epott
    
  end subroutine force_dipole
!=======================================================================
  subroutine set_paramsdir_dipole(dname)
!
!  Accessor routine to set paramsdir.
!
    implicit none
    character(len=*),intent(in):: dname

    paramsdir = trim(dname)
    return
  end subroutine set_paramsdir_dipole
!=======================================================================
  subroutine read_params_dipole(myid_md,mpi_md_world,iprint,specorder,amass)
    use util, only: num_data
    implicit none
    include './params_unit.h'
    integer,intent(in):: myid_md,mpi_md_world,iprint
    character(len=3),intent(in):: specorder(nspmax)
    real(rp),intent(in):: amass(nspmax)

    integer:: isp,jsp,ierr,ni,nj
    real(rp):: ei,ej,ai,aj,cij,dij,aij,alpij,bij,betij
    character(len=3):: csp,cspi,cspj
    character(len=128):: cline,fname

    if( myid_md.eq.0 ) then
      dip_ni(:) = -1
      dip_ei(:) = -1.0_rp
      dip_cij(:,:) = 0.0_rp
      dip_dij(:,:) = 0.0_rp
      interact(:,:) = .false.
      fname = trim(paramsdir)//'/'//trim(cprmfname)
      open(ioprms,file=trim(fname),status='old')
      if( iprint.gt.0 ) write(6,'(/,a)') ' dipole parameters:'
      do while(.true.)
        read(ioprms,'(a)',end=10) cline
        if( num_data(cline,' ').eq.0 ) cycle
        if( cline(1:1).eq.'#' .or. cline(1:1).eq.'!' ) cycle
        if( num_data(cline,' ').eq.3 ) then
          backspace(ioprms)
!.....species, # of total electrons, 1st ionization energy
          read(ioprms,*) csp, ni, ei
          isp = csp2isp(csp)
          if( isp.gt.nspmax ) then
            write(6,*) ' Warning @read_params: since isp is greater than nspmax,'&
                 //' skip reading the line.'
            cycle
          endif
          dip_ni(isp) = ni
          dip_ei(isp) = ei
          if( iprint.gt.0 ) print '(a,a4,i5,es12.4)','   csp,ni,ei=',csp,ni,ei
        endif
      enddo  ! while reading...
10    close(ioprms)
!.....Set parameters
      if( iprint.gt.0 ) print *,'  cspi,cspj,cij,dij:'
      do isp=1,nspmax
        if( dip_ni(isp).le.0 ) cycle
        ni = dip_ni(isp)
        ei = dip_ei(isp)
        ai = ni*plankh**2 /amass(isp) /ei**2
        do jsp=isp,nspmax
          if( dip_ni(jsp).le.0 ) cycle
          nj = dip_ni(jsp)
          ej = dip_ei(jsp)
          aj = nj*plankh**2 /amass(jsp) /ei**2
          cij = 3.0_rp*ai*aj*ei*ej /2.0_rp/(ei+ej)
          dij = 9.0_rp*cij /4.0_rp *(ai*ei/ni +aj*ej/nj)
          dip_cij(isp,jsp) = cij
          dip_dij(isp,jsp) = dij
          interact(isp,jsp) = .true.
!.....Symmetrize
          dip_cij(jsp,isp) = cij
          dip_dij(jsp,isp) = dij
          interact(jsp,isp) = .true.
          if( iprint.gt.0 ) then
            print '(4x,2(1x,a3),6es12.4)' &
               ,specorder(isp),specorder(jsp),cij,dij
          endif
        enddo
      enddo
    endif  ! myid_md.eq.0

    call mpi_bcast(dip_ni,nspmax*nspmax,mpi_integer,0,mpi_md_world,ierr)
    call mpi_bcast(dip_ei,nspmax*nspmax,mpi_real_rp,0,mpi_md_world,ierr)
    call mpi_bcast(dip_cij,nspmax*nspmax,mpi_real_rp,0,mpi_md_world,ierr)
    call mpi_bcast(dip_dij,nspmax*nspmax,mpi_real_rp,0,mpi_md_world,ierr)
    call mpi_bcast(interact,nspmax*nspmax,mpi_logical,0,mpi_md_world,ierr)
    
  end subroutine read_params_dipole
end module dipole
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
