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
  character(len=20):: clr_init = 'random'
  real(8):: clrfield(3) = (/ 0d0, 0d0, 0d0 /)   ! in [eV/Ang]
  real(8):: clraccel(3) = (/ 0d0, 0d0, 0d0 /)   ! scaled acceleration
  real(8):: vacc(3) = (/ 0d0, 0d0, 0d0 /)
  
contains
!=======================================================================
  subroutine init_clrchg(specorder,ntot,clrtot,tagtot,myid,iprint)
!
!  Initialization for color charge.
!  This routine should be called after broadcasting data from node-0,
!  but before calling space_decomp.
!
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
    if( trim(clr_init).eq.'read' ) then  ! read from clrini
      call read_clr(ntot,clrtot,myid)
    else if(  trim(clr_init).eq.'all_one' ) then ! set all the clr == 1.0
      call set_clr_one(ntot,tagtot,clrtot,myid,iprint)
    else  ! random
      call set_clr_random(ntot,tagtot,clrtot,myid,iprint)
    endif

99  initialized = .true.
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
  subroutine clrchg_force(namax,natm,tag,aa,clr,hi,specorder,myid,iprint)
!
!  Add external forces on atoms of specified species
!
    integer,intent(in):: namax,natm,myid,iprint
    character(len=3),intent(in):: specorder(nspmax)
    real(8),intent(in):: tag(namax),hi(3,3),clr(namax)
    real(8),intent(inout):: aa(3,namax)
    logical,save:: l1st = .true.

    integer:: i,is

    if( l1st ) then
      if( myid.eq.0 ) then
!.....Write some settings
        if( iprint.ge.ipl_basic ) then
          print *,''
          print '(a)', ' Color charge NEMD:'
          print '(a,i3,a5)', '   Specified species = ',ispc_clrchg,trim(cspc_clrchg)
          print '(a,3f10.5)', '   Applied field [eV/A] = ',clrfield(1:3)
        endif
      endif
      l1st = .false.
    endif

!.....It is not necessary to do this every step,
!     but for simplicity it is done every time
    clraccel(1:3)= hi(1:3,1)*clrfield(1) +hi(1:3,2)*clrfield(2) &
         +hi(1:3,3)*clrfield(3)

    do i=1,natm
      is= int(tag(i))
      if( is.ne.ispc_clrchg ) cycle
      aa(1:3,i)= aa(1:3,i) +clr(i)*clraccel(1:3)
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
    include 'mpif.h'
    integer,intent(in):: natm,mpi_world,myid,iprint
    real(8),intent(in):: tag(natm),am(nspmax)
    real(8),intent(inout):: va(3,natm)
    
    integer:: i,is,ierr
    real(8):: sump(3),tmps(3),tmp,amtot,ami

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
!     compile-command: "make pmd"
!     End:
