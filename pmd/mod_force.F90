module force
!-----------------------------------------------------------------------
!                     Last-modified: <2019-05-24 16:55:10 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
  implicit none
  save
  integer:: num_forces = -1
  character(len=128),allocatable:: force_list(:)

!.....Force index list
  integer,parameter:: N_FORCES = 29
  character(len=128):: force_index_list(N_FORCES)

  logical:: luse_force(N_FORCES)

  logical:: luse_charge

!.....Overlay main potential with a nuclear repulsive one, usually ZBL.
  logical:: loverlay = .false.
!.....Overlay type: pair or atom
  character(len=128):: overlay_type = 'pair'
  character(len=128):: overlay_force = 'ZBL'
  type overlay
    character(len=3):: csp1,csp2
    real(8):: rin,rout
  end type overlay
  type(overlay),allocatable:: overlays(:,:)

contains
!=======================================================================
  subroutine init_mod_force()

    force_index_list(1)  =  "LJ"
    force_index_list(2)  =  "Ito3_WHe"
    force_index_list(3)  =  "RK_WHe"
    force_index_list(4)  =  "RK_FeH"
    force_index_list(5)  =  "Ramas_FeH"
    force_index_list(6)  =  "Ackland_Fe"
    force_index_list(7)  =  "SW_Si"
    force_index_list(8)  =  "EDIP_Si"
    force_index_list(9)  =  "Brenner"
    force_index_list(10) =  "Brenner_vdW"
    force_index_list(11) =  "Lu_WHe"
    force_index_list(12) =  "Branicio_AlN"
    force_index_list(13) =  "Mishin_Al"
    force_index_list(14) =  "AFS_W"
    force_index_list(15) =  "SC_Fe"
    force_index_list(16) =  "SM_Al"
    force_index_list(17) =  "EAM"
    force_index_list(18) =  "linreg"
    force_index_list(19) =  "NN"
    force_index_list(20) =  "Morse"
    force_index_list(21) =  "Morse_repul"
    force_index_list(22) =  "vcMorse"
    force_index_list(23) =  "Buckingham"
    force_index_list(24) =  "Bonny_WRe"
    force_index_list(25) =  "ZBL"
    force_index_list(26) =  "screened_Coulomb"
    force_index_list(27) =  "Ewald"
    force_index_list(28) =  "Ewald_long"
    force_index_list(29) =  "Coulomb"

    luse_force(:) = .false.

!!$    do i=1,num_forces
!!$      
!!$    enddo

  end subroutine init_mod_force
!=======================================================================
  function use_force(force_name)
    implicit none
    character(len=*),intent(in):: force_name
    logical:: use_force
    integer:: i

    use_force = .false.
    do i=1,num_forces
      if( trim(force_name).eq.trim(force_list(i)) ) then
        use_force = .true.
        return
      endif
    enddo
    return
  end function use_force
!=======================================================================
  subroutine set_use_charge()

    luse_charge = .false.
    if(  use_force('screened_Coulomb') .or. &
         use_force('Ewald') .or. &
         use_force('Ewald_long') .or. &
         use_force('Coulomb') ) luse_charge = .true.

  end subroutine set_use_charge
!=======================================================================
  subroutine write_forces(myid)
    implicit none
    integer,intent(in):: myid
    integer:: i

    if( myid.eq.0 ) then
      write(6,'(a)',advance='no') ' Use the following force-fields:'
      do i=1,num_forces
        write(6,'(2x,a)',advance='no') trim(force_list(i))
      enddo
      print *,''
    endif
  end subroutine write_forces
!=======================================================================
  subroutine bcast_force()
!
!   Broadcast variables related to mod_force.
!
    use pmdmpi
    include 'mpif.h'
    integer:: ierr

    call mpi_bcast(num_forces,1,mpi_integer,0,mpicomm,ierr)
    if( .not.allocated(force_list) ) then
      allocate(force_list(num_forces))
    endif
    call mpi_bcast(force_list,128*num_forces,mpi_character &
         ,0,mpicomm,ierr)
    
  end subroutine bcast_force
end module force
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
