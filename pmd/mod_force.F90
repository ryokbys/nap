module force
!-----------------------------------------------------------------------
!                     Last-modified: <2024-07-28 11:47:09 KOBAYASHI Ryo>
!-----------------------------------------------------------------------
  use pmdvars,only: nspmax
  implicit none
  private
  save

  public:: init_mod_force,use_force,set_use_charge,set_use_elec_temp, &
       write_forces,bcast_force,calc_overlay
  public:: luse_charge, luse_elec_temp, force_list, num_forces, loverlay, &
       ol_type, ol_force, ol_ranges, ol_alphas, ol_dalphas, ol_pair
  
!.....Force index list
  integer,parameter:: N_FORCES = 35
  character(len=128):: force_index_list(N_FORCES)

  integer:: num_forces = -1
  character(len=128):: force_list(N_FORCES)

  logical:: luse_force(N_FORCES)
  logical:: luse_charge
  logical:: luse_elec_temp

!.....Overlay main potential with a nuclear repulsive one, usually ZBL.
  logical:: loverlay = .false.

!.....Overlay type: pot or force
  character(len=128):: ol_type = 'pot'
  character(len=128):: ol_force = 'ZBL'
  real(8):: ol_ranges(2,nspmax)
  real(8),allocatable:: ol_alphas(:,:),ol_dalphas(:,:)

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
    force_index_list(30) =  "FPC"
    force_index_list(31) =  "cspline"
    force_index_list(32) =  "Tersoff"
    force_index_list(33) =  "angular"
    force_index_list(34) =  "DNN"
    force_index_list(35) =  "fdesc"

    luse_force(:) = .false.

    ol_ranges(:,:) = -1d0
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
  subroutine set_use_elec_temp()

    luse_elec_temp = .false.
    if(  use_force('Tersoff')  ) luse_elec_temp= .true.

  end subroutine set_use_elec_temp
!=======================================================================
  subroutine write_forces(myid)
    implicit none
    integer,intent(in):: myid
    integer:: i

    if( myid.eq.0 ) then
      print *,''
      write(6,'(a)',advance='no') ' Use the following force-fields:'
      do i=1,num_forces
        write(6,'(2x,a)',advance='no') trim(force_list(i))
      enddo
      print *,''
    endif
  end subroutine write_forces
!=======================================================================
  subroutine bcast_force(mpicomm)
!
!   Broadcast variables related to mod_force.
!
    integer,intent(in):: mpicomm
    include 'mpif.h'
    integer:: ierr

    call mpi_bcast(num_forces,1,mpi_integer,0,mpicomm,ierr)
    call mpi_bcast(force_list,128*N_FORCES,mpi_character &
         ,0,mpicomm,ierr)

!.....Overlay
    call mpi_bcast(loverlay,1,mpi_logical,0,mpicomm,ierr)
    call mpi_bcast(ol_type,128,mpi_character,0,mpicomm,ierr)
    call mpi_bcast(ol_force,128,mpi_character,0,mpicomm,ierr)
    call mpi_bcast(ol_ranges,2*nspmax,mpi_real8,0,mpicomm,ierr)
      
  end subroutine bcast_force
!=======================================================================
  function ol_pair(ir,isp,jsp)
!
!  Return rin(ir==1), rout(ir==2) for pair (isp,jsp).
!
    integer:: ir,isp,jsp
    real(8):: ol_pair

    ol_pair = (ol_ranges(ir,isp)+ol_ranges(ir,jsp)) /2
    return
  end function ol_pair
!=======================================================================
  subroutine ol_allocate(namax,nnmax)
    integer,intent(in):: namax,nnmax
    
    if( .not.allocated(ol_alphas) ) then
      allocate(ol_alphas(0:nnmax,namax),ol_dalphas(nnmax,namax))
    else if( size(ol_alphas).ne.(nnmax+1)*namax ) then
      deallocate(ol_alphas,ol_dalphas)
      allocate(ol_alphas(0:nnmax,namax),ol_dalphas(nnmax,namax))
    endif
    return
  end subroutine ol_allocate
!=======================================================================
  subroutine get_fol_dfol(r,isp,jsp,fol,dfol)
    include "params_unit.h"
!!$    real(8),parameter:: pi = 3.14159265358979d0
    real(8),intent(in):: r
    integer,intent(in):: isp,jsp
    real(8),intent(out):: fol,dfol

    real(8):: ri,ro,x

    ri = ol_pair(1,isp,jsp)
    ro = ol_pair(2,isp,jsp)

    if( r.ge.ro ) then
      fol = 1d0
      dfol = 0d0
    else if( r.ge.ri .and. r.lt.ro ) then
      x = (r-ro)/(ri-ro)*pi
      fol = 0.5d0 *(1d0 +cos(x))
      dfol = -0.5d0*pi /(ri-ro) *sin(x)
    else
      fol = 0d0
      dfol = 0d0
    endif
    
    return
  end subroutine get_fol_dfol
!=======================================================================
  subroutine calc_overlay(namax,natm,nb,nnmax,h,tag,ra,lspr &
       ,l1st,iprint)
!
!  Compute overlay coefficients of each pair and atom.
!
    integer,intent(in):: namax,natm,nb,nnmax,lspr(0:nnmax,namax)
    integer,intent(in):: iprint
    real(8),intent(in):: h(3,3),tag(namax),ra(3,namax)
    logical,intent(in):: l1st

    integer:: ia,ja,jj,is,js
    real(8):: xi(3),xj(3),xij(3),rij(3),dij2,dij,ri,ro,fol,dfol

!.....Check rcut of lspr, which should be larger than 2*rout.

    call ol_allocate(namax,nnmax)
    ol_alphas(1:nnmax,:) = 0d0
    ol_alphas(0,:) = 1d0
    ol_dalphas(:,:) = 1d0

    do ia=1,natm+nb
      xi(1:3) = ra(1:3,ia)
      is = int(tag(ia))
      do jj=1,lspr(0,ia)
        ja = lspr(jj,ia)
        js = int(tag(ja))
        xj(1:3)= ra(1:3,ja)
        xij(1:3)= xj(1:3) -xi(1:3)
        rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij2= rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3)
        dij = dsqrt(dij2)
        call get_fol_dfol(dij,is,js,fol,dfol)
        ol_alphas(jj,ia)= fol
        ol_alphas(0,ia) = ol_alphas(0,ia) *fol
        ol_dalphas(jj,ia)= dfol
      enddo  ! jj=...
    enddo  ! ia=...

    return
  end subroutine calc_overlay
  
end module force
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
