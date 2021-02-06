module extforce
  use pmdio,only: csp2isp,nspmax
  implicit none
  include "./const.h"
  save

  logical:: lextfrc = .false.  ! flag to apply external force
  logical:: initialized = .false.
  character(len=3):: cspc_extfrc = 'all' ! [default: 'all']
  integer:: ispc_extfrc
  real(8):: extfrc(3) = (/ 0d0, 0d0, 0d0 /)   ! in [eV/Ang]
  real(8):: extaa(3) = (/ 0d0, 0d0, 0d0 /)
  
contains
!=======================================================================
  subroutine init_extfrc(specorder,myid,iprint)
    integer,intent(in):: myid,iprint
    character(len=3),intent(in):: specorder(nspmax)

    if( trim(cspc_extfrc).eq.'all' ) then
      ispc_extfrc = 0
    else
      ispc_extfrc = csp2isp(trim(cspc_extfrc),specorder)
    endif

!.....Write some settings
    if( myid.eq.0 .and. iprint.ge.ipl_basic ) then
      print *,''
      print '(a)', ' External force:'
      print '(a,a5)', '   Specified species = ',trim(cspc_extfrc)
      print '(a,3f10.5)', '   Applied force [eV/A] = ',extfrc(1:3)
    endif
    
    initialized = .true.
  end subroutine init_extfrc
!=======================================================================
  subroutine add_extfrc(natm,tag,aa,hi,specorder,myid,iprint)
!
!  Add external forces on atoms of specified species
!
    integer,intent(in):: natm,myid,iprint
    character(len=3),intent(in):: specorder(nspmax)
    real(8),intent(in):: tag(natm),hi(3,3)
    real(8),intent(inout):: aa(3,natm)

    integer:: i,is

    if( .not.initialized ) call init_extfrc(specorder,myid,iprint)

!.....It is not necessary to do this every step,
!     but for simplicity it is done every time
    extaa(1:3)= hi(1:3,1)*extfrc(1) +hi(1:3,2)*extfrc(2) &
         +hi(1:3,3)*extfrc(3)

    do i=1,natm
      is= int(tag(i))
      if( is.ne.ispc_extfrc ) cycle
      aa(1:3,i)= aa(1:3,i) +extaa(1:3)
    enddo
    return
  end subroutine add_extfrc
!=======================================================================
  subroutine rm_trans_extfrc(natm,tag,va,am,mpi_world,myid,iprint)
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
      if( is.eq.ispc_extfrc ) cycle
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
      if( is.eq.ispc_extfrc ) cycle
      va(1:3,i)= va(1:3,i) -sump(1:3)/amtot
    enddo
  end subroutine rm_trans_extfrc
end module extforce
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
