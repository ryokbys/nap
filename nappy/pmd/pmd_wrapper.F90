subroutine run(ntot0,rtot,vtot,atot,stot,ekitot,epitot, &
     naux,auxtot,hmat,ispcs,ekin,epot,stnsr,linit)
  use pmdvars,only: nx,ny,nz,iprint
  implicit none
  integer,intent(in):: ntot0,naux,ispcs(ntot0)
  real(8),intent(inout):: rtot(3,ntot0),vtot(3,ntot0),hmat(3,3,0:1)
!f2py intent(in,out):: rtot,vtot,hmat
  real(8),intent(out):: atot(3,ntot0),stot(3,3,ntot0),ekitot(3,3,ntot0), &
       epitot(ntot0),auxtot(naux,ntot0),ekin,epot,stnsr(3,3)
  logical,intent(in):: linit

  integer:: ntot
  real(8):: hunit,tagtot(ntot0)
  
  if( nx.lt.0 .or. ny.lt.0 .or. nz.lt.0 ) then
    print *,'Some pmdvars should be set before calling run().'
    stop
  endif

  call get_tagtot(ntot0,ispcs,tagtot)
  
  ntot = ntot0
  hunit = 1d0
!!$  print *,'iprint,ntot0,rtot(:,ntot0)=',iprint,ntot0,rtot(:,ntot0)
!!$  call pmd_core(hunit,hmat,ntot0,ntot,tagtot,rtot,vtot,atot,stot, &
!!$       ekitot,epitot,auxtot,epot,ekin,stnsr)
  call oneshot(hunit,hmat,ntot0,tagtot,rtot,vtot,atot,stot, &
       ekitot,epitot,auxtot,ekin,epot,stnsr,linit)
  return
end subroutine run
!=======================================================================
subroutine get_tagtot(ntot,ispcs,tagtot)
  integer,intent(in):: ntot,ispcs(ntot)
  real(8),intent(out):: tagtot(ntot)

  integer:: i,ifmv,isp

  do i=1,ntot
    isp = ispcs(i)
    ifmv = 1
    tagtot(i) = isp*1d0 +ifmv*1d-1 +i*1d-14
  enddo
  return
end subroutine get_tagtot
!=======================================================================
subroutine set_pmdvars(ns,ls,cspcs,nf,lf,cfrcs,rc0, &
     mpi_comm,myid,nodes,iprint0,nstp0,naux0,laux,cauxarr0)
!
!  Set variables to be stored in pmdvars module that are required 
!  to call pmd_core.
!
  use pmdvars,only: specorder,nspmax,dt,rbuf,rc1nn,rc,nx,ny,nz,iprint, &
       mpi_md_world,myid_md,nodes_md,am,nstp,naux,ifpmd,nerg,cauxarr
  use force
  use element
  implicit none 
  integer,intent(in):: ns,ls
  character(len=1),intent(in):: cspcs(ns,ls)
!f2py integer,intent(hide),depend(cspcs):: ns=shape(cspcs,0),ls=shape(cspcs,1)
  integer,intent(in):: nf,lf
  character(len=1),intent(in):: cfrcs(nf,lf)
!f2py integer,intent(hide),depend(cfrcs):: nf=shape(cfrcs,0),lf=shape(cfrcs,1)
  real(8),intent(in):: rc0
  integer,intent(in):: mpi_comm,myid,nodes,naux0,iprint0,nstp0
  integer,intent(in):: laux
  character(len=1),intent(in):: cauxarr0(naux0,laux)
!f2py integer,intent(hide),depend(cauxarr0):: laux=shape(cauxarr0,1)

  integer:: i,j
  character:: c3*3, c128*128, c6*6
  type(atom):: elem
  logical:: lcoulomb = .false.

!.....Set specorder
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

!.....Set cauxarr0
  if( laux.ne.6 ) then
    print *,' The length of cauxarr char should be 6, laux = ',laux
    stop
  endif
  if( lcoulomb .and. naux0.lt.2 ) then
    print *,' naux should be greater than 1 when Coulomb potential is used.'
    stop
  endif
  if( allocated(cauxarr) .and. size(cauxarr).ne.naux0 ) deallocate(cauxarr)
  if( .not.allocated(cauxarr) ) allocate(cauxarr(naux0))
  
  do i=1,naux0
    write(c6,'(6a1)') cauxarr0(i,1:laux)
    cauxarr(i) = trim(c6)
  end do

  naux = naux0
  nstp = nstp0
  do i=1,nspmax
    c3 = specorder(i)
    if( trim(c3).ne.'x' ) then
      elem = get_element(trim(c3))
      am(i) = elem%mass
    endif
  enddo
  dt = 1d0
  rbuf = 0d0
  rc1nn = 3d0
  rc = rc0
  nx = 1
  ny = 1
  nz = 1
  iprint = iprint0
  ifpmd = 0
  nerg = 0

  nodes_md = nodes
  myid_md = myid
  mpi_md_world = mpi_comm

end subroutine set_pmdvars
!=======================================================================
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make"
!     End:
