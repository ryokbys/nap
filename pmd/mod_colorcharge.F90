module clrchg
!
!  Module for color charge NEMD
!
  use pmdvars,only: nspmax
  use util,only: csp2isp
  implicit none
  include "./const.h"
  save

  character(len=128),parameter:: cfclrini = 'clrini'
  integer,parameter:: ioini = 68

  logical:: lclrchg = .false.  ! flag to apply external force
  logical:: initialized = .false.
  character(len=3):: cspc_clrchg = 'non'
  integer:: ispc_clrchg
  character(len=20):: clr_set = 'random'
  real(8):: clrfield(3) = (/ 0d0, 0d0, 0d0 /)   ! in [eV/Ang]
  real(8):: clraccel(3) = (/ 0d0, 0d0, 0d0 /)   ! scaled acceleration
  real(8):: vacc(3) = (/ 0d0, 0d0, 0d0 /)
  real(8):: clrregion(3,2)  ! regional condition
  data clrregion(1:3,1) / -1d0, -1d0, -1d0 / ! sufficiently smaller than 0
  data clrregion(1:3,2) /  2d0,  2d0,  2d0 / ! sufficiently greater than 1
  
contains
!=======================================================================
  subroutine init_clrchg(specorder,ntot,clrtot,tagtot,myid,iprint)
!
!  Initialization for color charge.
!  This routine should be called after broadcasting data from node-0,
!  but before calling space_decomp.
!
    use force,only: luse_charge
    integer,intent(in):: myid,iprint,ntot
    character(len=3),intent(in):: specorder(nspmax)
    real(8),intent(in):: tagtot(ntot)
    real(8),intent(inout):: clrtot(ntot)

!!$    real(8),external:: urnd

    if( trim(cspc_clrchg).eq.'non' ) then
      stop 'ERROR: spcs_clrchg must be specified.'
    else
      ispc_clrchg = csp2isp(trim(cspc_clrchg))
    endif

    clrtot(:) = 0d0
    if( trim(clr_set).eq.'read' ) then  ! read from clrini
      call read_clr(ntot,clrtot,myid)
      initialized = .true.
    else if( trim(clr_set).eq.'all_one' ) then  ! set all the clr == 1.0
      call set_clr_one(ntot,tagtot,clrtot,myid,iprint)
      initialized = .true.
    else if( trim(clr_set).eq.'random' ) then ! random
      call set_clr_random(ntot,tagtot,clrtot,myid,iprint)
      initialized = .true.
    else if( index(clr_set,'chg').ne.0 ) then  ! charge related
      if( .not. luse_charge ) stop 'ERROR: clr_set is set as chg related, ' &
           //'but the potential does not use charges.'
    endif
! else do nothing for the moment

    return
  end subroutine init_clrchg
!=======================================================================
  subroutine read_clr(ntot,clrtot,myid)
!
!  Read clr of atoms from file.
!
    integer,intent(in):: ntot,myid
    real(8),intent(inout):: clrtot(ntot)
    integer:: n,i,j

    if( myid.eq.0 ) then
      open(ioini,file=trim(cfclrini),status='old')
      read(ioini,*) n
      if( n.ne.ntot ) then
        print *,'ERROR: n in clrini is not same with ntot in pmdini.'
        stop
      endif
      do i=1,ntot
        read(ioini,*) j,clrtot(j)
      enddo
      close(ioini)
    endif
    return
  end subroutine read_clr
!=======================================================================
  subroutine set_clr_one(ntot,tagtot,clrtot,myid,iprint)
!
!  Set all the clr one
!
    integer,intent(in):: ntot,myid,iprint
    real(8),intent(in):: tagtot(ntot)
    real(8),intent(inout):: clrtot(ntot)

    integer:: i,is
    
    if( myid.eq.0 ) then
      clrtot(:) = 0d0
      do i=1,ntot
        is = int(tagtot(i))
        if( is.eq.ispc_clrchg ) then
          clrtot(i) = 1d0
        endif
      enddo
    endif

  end subroutine set_clr_one
!=======================================================================
  subroutine set_clr_round_chg(natm,tag,clr,chg,myid,iprint)
!
!  Set clr as round(chg)
!
    integer,intent(in):: natm,myid,iprint
    real(8),intent(in):: tag(natm),chg(natm)
    real(8),intent(inout):: clr(natm)

    integer:: i,is

    do i=1,natm
      is = int(tag(i))
      clr(i) = anint(chg(i))
    enddo
    return
  end subroutine set_clr_round_chg
!=======================================================================
  subroutine set_clr_random(ntot,tagtot,clrtot,myid,iprint)
!
!  Set clr of each atom randomly.
!
    integer,intent(in):: ntot,myid,iprint
    real(8),intent(in):: tagtot(ntot)
    real(8),intent(inout):: clrtot(ntot)

    integer:: i,is,nclr,np,nm
    real(8):: r
    
    if( myid.eq.0 ) then
!.....Count num of ions which color charges are to be set.
      nclr = 0
      do i=1,ntot
        is = int(tagtot(i))
        if( is.eq.ispc_clrchg ) then
          nclr = nclr +1
        endif
      enddo
      if( mod(nclr,2) .ne. 0 ) stop 'ERROR: num of ions for color charge' &
           //' must be even.'
!.....Set color charges to those ions
      np = nclr/2
      nm = nclr/2
      clrtot(:) = 0d0
      do i=1,ntot
        is = int(tagtot(i))
        if( is.ne.ispc_clrchg ) cycle
!!$        r = urnd()  !<=== something wrong with urnd() func...
        call random_number(r)
        if( r.lt.0.5d0 ) then  ! +1
          if( np.le.0 ) then
            clrtot(i) = -1d0
            nm = nm -1
          else
            clrtot(i) = 1d0
            np = np -1
          endif
        else  ! -1
          if( nm.le.0 ) then
            clrtot(i) = 1d0
            np = np -1
          else
            clrtot(i) = -1d0
            nm = nm -1
          endif
        endif
      enddo

    endif

  end subroutine set_clr_random
!=======================================================================
  subroutine clrchg_force(namax,natm,tag,aa,aux,hi,specorder,myid,iprint)
!
!  Add external forces on atoms of specified species
!
    use pmdvars,only: naux,iaux_chg,iaux_clr,ra,sorg
    integer,intent(in):: namax,natm,myid,iprint
    character(len=3),intent(in):: specorder(nspmax)
    real(8),intent(in):: tag(namax),hi(3,3)
    real(8),intent(inout):: aa(3,namax),aux(naux,namax)
    logical,save:: l1st = .true.

    integer:: i,is
    real(8):: ri(3)

    if( l1st ) then
      if( myid.eq.0 ) then
!.....Write some settings
        if( iprint.ge.ipl_basic ) then
          print *,''
          print '(a)', ' Color charge NEMD:'
          print '(a,a)', '   color charge type:', trim(clr_set)
          print '(a,i3,a5)', '   Specified species = ',ispc_clrchg,trim(cspc_clrchg)
          print '(a,3f10.5)', '   Applied field [eV/A] = ',clrfield(1:3)
          if( trim(clr_set).eq.'region' ) then
            print '(a)', '   regional condition:'
            do i=1,3
              print '(i6,2f6.3)', i, max(clrregion(i,1),0d0), min(clrregion(i,2),1d0)
            enddo
          endif
        endif
      endif
      l1st = .false.
    endif

!!$!.....It is not necessary to do this every step,
!!$!     but for simplicity it is done every time
!!$    clraccel(1:3)= hi(1:3,1)*clrfield(1) +hi(1:3,2)*clrfield(2) &
!!$         +hi(1:3,3)*clrfield(3)
!.....Now the aa(:,:) is in real unit, not normalized unit,
!     no need to multiply hi(:,:) matrix

    if( trim(clr_set).eq.'round_chg' ) then
      do i=1,natm
        aux(iaux_clr,i) = anint(aux(iaux_chg,i))
      enddo
    else if( trim(clr_set).eq.'chg_cation' ) then
      !...Only cations with greater 0.1 e charges
      do i=1,natm
        if( aux(iaux_chg,i).gt.0.1d0 ) then
          aux(iaux_clr,i) = 1d0
        else
          aux(iaux_clr,i) = 0d0
        endif
      enddo
!.....Set color charges on the atoms inside specified region
    else if( trim(clr_set).eq.'region') then
      aux(iaux_clr,1:natm) = 0d0
      do i=1,natm
        if( int(tag(i)).ne.ispc_clrchg ) cycle
        ri(1:3) = ra(1:3,i) +sorg(1:3)
        if( clrregion(1,1).lt.ri(1) .and. ri(1).lt.clrregion(1,2) .and. &
            clrregion(2,1).lt.ri(2) .and. ri(2).lt.clrregion(2,2) .and. &
            clrregion(3,1).lt.ri(3) .and. ri(3).lt.clrregion(3,2) ) then
          aux(iaux_clr,i) = 1d0
        endif
      enddo
    endif

    do i=1,natm
      is= int(tag(i))
      if( is.ne.ispc_clrchg ) cycle
      aa(1:3,i)= aa(1:3,i) +aux(iaux_clr,i)*clrfield(1:3)
    enddo
    return
  end subroutine clrchg_force
!=======================================================================
  subroutine accum_vels(namax,natm,hmat,va,clr,istp,nouterg,dt)
!
!  Accumurate velocities with mutiplying color charge.
!
    integer,intent(in):: namax,natm,istp,nouterg
    real(8),intent(in):: hmat(3,3),va(3,namax),clr(namax),dt

    integer:: i
    real(8):: vi(3)

    do i=1,natm
      vi(1:3) = hmat(1:3,1)*va(1,i) +hmat(1:3,2)*va(2,i) &
           +hmat(1:3,3)*va(3,i)
      vacc(1:3) = vacc(1:3) +vi(1:3)*clr(i)
    enddo
    if( mod(istp,nouterg).eq.0 ) then
      write(90,'(i10,4es12.4)') istp, dt*istp, vacc(1:3)
    endif
    
  end subroutine accum_vels
!=======================================================================
  subroutine rm_trans_clrchg(natm,tag,va,am,mpi_world,myid,iprint)
!
!  Remove translational momentum of all the atoms except specified species.
!
    use pmdvars,only: nrmtrans
    include 'mpif.h'
    integer,intent(in):: natm,mpi_world,myid,iprint
    real(8),intent(in):: tag(natm),am(nspmax)
    real(8),intent(inout):: va(3,natm)
    
    integer:: i,is,ierr
    real(8):: sump(3),tmps(3),tmp,amtot,ami

!.....If nrmtrans < 0, not to remove translation,
!.....because the user must set that intentionally.
    if( nrmtrans.lt.0 ) return

    sump(:)= 0d0
    amtot = 0d0
    do i=1,natm
      is = int(tag(i))
      if( is.eq.ispc_clrchg ) cycle
      ami= am(is)
      sump(1:3)= sump(1:3) +ami*va(1:3,i)
      amtot= amtot +ami
    enddo
    tmps(:)= sump(:)
    call mpi_allreduce(tmps,sump,3,mpi_real8,mpi_sum,mpi_world,ierr)
    tmp= amtot
    call mpi_allreduce(tmp,amtot,1,mpi_real8,mpi_sum,mpi_world,ierr)
    if( amtot.lt.1d-1 ) then
      print *,'Error: amtot.le.0.1 !, myid=',myid
      stop
    endif
    do i=1,natm
      is = int(tag(i))
      if( is.eq.ispc_clrchg ) cycle
      va(1:3,i)= va(1:3,i) -sump(1:3)/amtot
    enddo
  end subroutine rm_trans_clrchg
end module clrchg
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd lib"
!     End:
