subroutine run(ntot0,rtot,vtot,atot,stot,ekitot,epitot, &
     naux,auxtot,hmat,ispcs,ifmvs,ekin,epot,stnsr,linit)
  use pmdvars,only: nx,ny,nz,iprint,nstp,ifpmd
  use util,only: calc_nfmv
  implicit none
  integer,intent(in):: ntot0,ispcs(ntot0),ifmvs(ntot0),naux
  real(8),intent(inout):: rtot(3,ntot0),vtot(3,ntot0),hmat(3,3,0:1)
!f2py integer,intent(hide),depend(ispcs):: ntot0=shape(ispcs,0)
!f2py intent(in,out):: rtot,vtot,hmat
  real(8),intent(out):: atot(3,ntot0),stot(3,3,ntot0), &
       ekitot(3,3,ntot0),epitot(ntot0),auxtot(naux,ntot0),ekin, &
       epot,stnsr(3,3)
  logical,intent(in):: linit

  integer:: tagtot_isp(ntot0),tagtot_ifmv(ntot0),tagtot_igrp(4,ntot0)
  integer:: tagtot_itot(ntot0)
  real(8):: hunit
  
  if( nx.lt.0 .or. ny.lt.0 .or. nz.lt.0 ) then
    print *,'Some pmdvars should be set before calling run().'
    stop
  endif

  call get_tagtot(ntot0,ispcs,ifmvs,tagtot_isp,tagtot_ifmv, &
       tagtot_igrp,tagtot_itot)
  call calc_nfmv(ntot0,tagtot_ifmv,tagtot_igrp)
  
  hunit = 1d0
!!$  print *,'nstp,iprint=',nstp,iprint
!!$  print *,'iprint,ntot0,rtot(:,ntot0)=',iprint,ntot0,rtot(:,ntot0)
  call pmd_core(hunit,hmat,ntot0,tagtot_isp,tagtot_ifmv, &
       tagtot_igrp,tagtot_itot,rtot,vtot,atot,stot,ekitot,epitot, &
       auxtot,epot,ekin,stnsr)
  return
end subroutine run
!=======================================================================
subroutine get_tagtot(ntot,ispcs,ifmvs,tagtot_isp,tagtot_ifmv, &
     tagtot_igrp,tagtot_itot)
  implicit none 
  integer,intent(in):: ntot,ispcs(ntot),ifmvs(ntot)
  integer,intent(out):: tagtot_isp(ntot),tagtot_ifmv(ntot)
  integer,intent(out):: tagtot_igrp(4,ntot),tagtot_itot(ntot)

  integer:: i

  do i=1,ntot
    tagtot_isp(i) = ispcs(i)
    tagtot_ifmv(i) = ifmvs(i)
    tagtot_igrp(:,i) = 0
    tagtot_itot(i) = i
  enddo
  return
end subroutine get_tagtot
!=======================================================================
subroutine get_naux(naux0)
  use pmdvars,only: naux
  implicit none 
  integer,intent(out):: naux0

  naux0 = naux
  return
end subroutine get_naux
!=======================================================================
subroutine set_pmdvars(nsp0,ns,ls,cspcs,nf,lf,cfrcs,rc0,rbuf0, &
     iprint0,nstp0,dt0, & !,naux0,laux,cauxarr0
     ifdmp0,dmp0,eps_conv0,n_conv0, &
     lctctl,ctctl0,tinit0,tfin0,ttgt0,trlx0,nrmtrans0, &
     lcpctl,cpctl0,ptgt0,stgt0,srlx0,lcellfix0, &
     ifpmd0,npmd0,nerg0,nnmax0,lrealloc0)
!
!  Set variables to be stored in pmdvars module that are required 
!  to call pmd_core.
!
  use pmdvars,only: specorder,has_specorder,nspmax,nsp,dt,rbuf,rc1nn,rc,iprint, &
       am,nstp,naux,ifpmd,npmd,nerg,cauxarr,ifdmp,dmp,eps_conv,n_conv, &
       cpctl,ptgt,stgt,srlx,lcellfix,nnmax,lrealloc, &
       ctctl,tinit,tfin,ttgt,trlx,nrmtrans,nx,ny,nz
  use force
  use element
  implicit none
  integer,intent(in):: nsp0
  integer,intent(in):: ns,ls
  character(len=1),intent(in):: cspcs(ns,ls)
!f2py integer,intent(hide),depend(cspcs):: ns=shape(cspcs,0),ls=shape(cspcs,1)
  integer,intent(in):: nf,lf
  character(len=1),intent(in):: cfrcs(nf,lf)
!f2py integer,intent(hide),depend(cfrcs):: nf=shape(cfrcs,0),lf=shape(cfrcs,1)
  real(8),intent(in):: rc0,rbuf0,dt0,dmp0,eps_conv0
  integer,intent(in):: iprint0,nstp0,ifdmp0,n_conv0
!!$  integer,intent(in):: naux0,laux
!!$  character(len=1),intent(in):: cauxarr0(naux0,laux)
!!$!f2py integer,intent(hide),depend(cauxarr0):: laux=shape(cauxarr0,1)
  integer,intent(in):: lctctl
  character(len=1),intent(in):: ctctl0(lctctl)
!f2py integer,intent(hide),depend(ctctl0):: lctctl=shape(ctctl0,0)
  real(8),intent(in):: tinit0, tfin0, trlx0
  real(8),intent(in):: ttgt0(9)
  integer,intent(in):: nrmtrans0
  real(8),intent(in):: ptgt0,srlx0,stgt0(3,3)
  logical,intent(in):: lcellfix0(3,3)
  integer,intent(in):: lcpctl
  character(len=1),intent(in):: cpctl0(lcpctl)
!f2py integer,intent(hide),depend(cpctl0):: lcpctl=shape(cpctl0,0)
  integer,intent(in):: ifpmd0,npmd0,nerg0,nnmax0
  logical:: lrealloc0

  integer:: i,j
  character:: c3*3, c128*128, c6*6, c20*20
  type(atom):: elem
  logical:: lcoulomb = .false.

  call init_element()
  iprint = iprint0

!.....Set specorder
  nsp = nsp0
  if( ls.ne.3 ) then
    print *,' The length of specorder char should be 3, ls = ',ls
    stop
  else if( ns.ne.size(specorder) ) then
    print *,' The length of specorder array should be 9, ns = ',ns
    stop
  endif
  do i=1,ns
    write(c3,'(3a1)') cspcs(i,1:ls)
    specorder(i) = trim(c3)
!!$    specorder(i) = trim(csp)
!!$    write(specorder(i),'(3a1)') cspcs(i,1:ls)
  enddo
  has_specorder = .true.

!.....Set force_list
  if( lf.ne.128 ) then
    print *,' The length of force_list char should be 128, lf = ',lf
    stop
  endif
  num_forces = nf
  do i=1,num_forces
    write(c128,'(128a1)') cfrcs(i,1:lf)
    force_list(i) = trim(c128)
  end do

  rc1nn = 3d0
  rbuf = rbuf0
!.....RC should be set before calling init_force
  rc = rc0

!.....It is required to call init_force and read some in.params.XXX to define aux array
  call init_force(.true.)
!.....Before allocating auxiliary array, set naux (num of auxiliary data)
  call set_cauxarr()

  write(c20,'(20a1)') cpctl0(1:lcpctl)
  cpctl = trim(c20)
  write(c20,'(20a1)') ctctl0(1:lctctl)
  ctctl = trim(c20)

!!$  naux = naux0
  nstp = nstp0
  do i=1,nspmax
    c3 = specorder(i)
    if( trim(c3).ne.'x' ) then
      elem = get_element(trim(c3))
      am(i) = elem%mass
    endif
  enddo
  dt = dt0
  nx = 1
  ny = 1
  nz = 1
  nerg = nerg0
  ifpmd = ifpmd0
  npmd = npmd0
  ifdmp = ifdmp0
  dmp = dmp0
  eps_conv = eps_conv0
  n_conv = n_conv0
  ptgt = ptgt0
  stgt(:,:) = stgt0(:,:)
  srlx = srlx0
  lcellfix(:,:) = lcellfix0(:,:)
  nnmax = nnmax0
  lrealloc = lrealloc0
  tinit = tinit0
  tfin = tfin0
  trlx = trlx0
  ttgt(:) = ttgt0(:)
  nrmtrans = nrmtrans0

end subroutine set_pmdvars
!=======================================================================
subroutine set_mpivars(mpi_comm,nodes,myid)
  use pmdvars,only: mpi_md_world,myid_md,nodes_md
  integer,intent(in):: mpi_comm,nodes,myid
  
  nodes_md = nodes
  myid_md = myid
  mpi_md_world = mpi_comm
end subroutine set_mpivars
!=======================================================================
subroutine wrap_calc_rdf(natm,ra,tag,h,hi,rmax,rmin,l1st, &
     lpairwise,nbin,msp,dists,rdfs)
!
! Wrapper for calling distfunc.calc_rdf from python 
!
  use distfunc,only: calc_rdf
  use pairlist,only: mk_lspr_sngl
  use util,only: get_vol
  implicit none
  integer,intent(in):: natm,nbin,msp
  real(8),intent(in):: ra(3,natm),tag(natm),h(3,3),hi(3,3)
  real(8),intent(in):: rmax,rmin
  logical,intent(in):: l1st,lpairwise
  real(8),intent(out):: dists(nbin)
  real(8),intent(out):: rdfs(nbin,0:msp,0:msp)

  real(8),parameter:: pi = 3.14159265358979d0
  integer:: i,iprint,nnmax
  real(8):: vol,rho
  integer,allocatable:: lspr(:,:)
  integer,allocatable:: tag_isp(:)

!.....Estimate nnmax
  vol = get_vol(h)
  rho = max(dble(natm)/vol, 0.2d0)
  nnmax = int(1.2d0 *rho *4d0*pi*rmax**3 /3)  ! margin 20 %
  allocate(lspr(0:nnmax,natm))
  allocate(tag_isp(natm))
  do i=1,natm
    tag_isp(i) = int(tag(i))
  enddo
  
  iprint = 1
  call mk_lspr_sngl(natm,natm,nnmax,tag_isp,ra,rmax,h,hi,lspr,iprint,l1st)

  call calc_rdf(natm,nnmax,tag,h,ra,rmax,rmin,lspr,iprint,l1st, &
       lpairwise,msp,nbin,dists,rdfs)
  deallocate(tag_isp)
  deallocate(lspr)
end subroutine wrap_calc_rdf
!=======================================================================
subroutine wrap_calc_adf(natm,ra,tag,h,hi,rmax,ntrpl,itriples, &
     nbin,angs,adfs,l1st)
!
! Wrapper for calling distfunc.calc_adf from python 
!
  use distfunc,only: calc_adf
  use pairlist,only: mk_lspr_sngl
  use util,only: get_vol
  implicit none
  integer,intent(in):: natm,nbin,ntrpl
  real(8),intent(in):: ra(3,natm),tag(natm),h(3,3),hi(3,3)
  real(8),intent(in):: rmax
  integer,intent(in):: itriples(3,ntrpl)
  logical,intent(in):: l1st
  real(8),intent(out):: angs(nbin)
  real(8),intent(out):: adfs(nbin,ntrpl)

  real(8),parameter:: pi = 3.14159265358979d0
  integer:: i,iprint,nnmax
  real(8):: vol,rho,dang
  integer,allocatable:: lspr(:,:)
  integer,allocatable:: tag_isp(:)

!.....Estimate nnmax
  vol = get_vol(h)
  rho = max(dble(natm)/vol, 0.2d0)
  nnmax = int(1.2d0 *rho *4d0*pi*rmax**3 /3)  ! margin 20 %
  allocate(lspr(0:nnmax,natm))
  allocate(tag_isp(natm))
  do i=1,natm
    tag_isp(i) = int(tag(i))
  enddo
  
  iprint = 1
  call mk_lspr_sngl(natm,natm,nnmax,tag_isp,ra,rmax,h,hi,lspr,iprint,l1st)

  dang = 180d0 /nbin
  call calc_adf(natm,nnmax,tag,h,ra,rmax,lspr,ntrpl,itriples, &
       dang,nbin,angs,adfs)
  deallocate(tag_isp)
  deallocate(lspr)
end subroutine wrap_calc_adf
!=======================================================================
subroutine wrap_lspr_sngl(natm,ra,tag,h,hi,rcut,iprint,l1st,nnmax,lspr)
!
! Wrapper for calling distfunc.calc_adf from python 
!
  use pairlist,only: mk_lspr_sngl
  implicit none
  integer,intent(in):: natm,nnmax,iprint
  integer:: i
  real(8),intent(in):: ra(3,natm),tag(natm),h(3,3),hi(3,3)
  real(8),intent(in):: rcut
  logical,intent(in):: l1st
  integer,intent(out):: lspr(0:nnmax,natm)
  integer:: tag_isp(natm)

  do i=1,natm
    tag_isp(i) = int(tag(i))
  enddo
  call mk_lspr_sngl(natm,natm,nnmax,tag_isp,ra,rcut,h,hi,lspr,iprint,l1st)
end subroutine wrap_lspr_sngl
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make"
!     End:
