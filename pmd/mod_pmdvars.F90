module pmdvars
!-----------------------------------------------------------------------
!                    Last modified: <2017-06-19 17:53:07 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
  implicit none
!=======================================================================
! PARAMETERS
!=======================================================================
!.....max. num. of atoms in a node
  integer:: namax = 10000
!.....max. num. of boundary-particles
  integer:: nbmax = 10000
!.....max. num. of neighbors
  integer:: nnmax = 1000

!.....output #
  integer,parameter:: ioerg = 11
  integer,parameter:: iotemp= 12
  integer,parameter:: iostrs= 13
  integer,parameter:: iotdst= 14

!=======================================================================
! VARIABLES
!=======================================================================
  integer:: nouterg,noutpmd,istp &
       ,iocntpmd,iocnterg
  integer:: natm,nb,nsp
  real(8):: tcpu,tcpu1,tcpu2,tcom,tspdcmp
  real(8):: epot0,vmaxold,vmax
  real(8):: tgmm
!!$  real(8):: tgmm,tfac(9),ediff(9),ediff0(9),temp(9),ekl(9)
  real(8),allocatable:: tfac(:),ediff(:),ediff0(:),temp(:),ekl(:)
  integer,allocatable:: ndof(:)
!!$  integer:: ndof(9)
  integer:: nxmlt
!.....Search time and expiration time
  real(8):: ts,te
  integer:: istps,istpe
!.....simulation box
  real(8):: hi(3,3),vol,sgm(3,3),al(3),avol
  real(8):: ht(3,3,0:1),hti(3,3),dh
!.....positions, velocities, and accelerations
  real(8),allocatable:: ra(:,:),va(:,:),aa(:,:),ra0(:,:) &
       ,strs(:,:,:),stt(:,:,:)
!!$    real(8):: ra(3,namax),va(3,namax),aa(3,namax),ra0(3,namax) &
!!$         ,strs(3,3,namax),stt(3,3,namax)
!.....real*8 identifier which includes species, index of FMV, total id
  real(8),allocatable:: tag(:)
!!$    real(8):: tag(namax)
  integer,allocatable:: lspr(:,:)
!!$    integer:: lspr(0:nnmax,namax)
!.....potential and kinetic energy per atoms
  real(8),allocatable:: epi(:),eki(:,:,:),stp(:,:,:)
!!$    real(8):: epi(namax),eki(3,3,namax),stp(3,3,namax)
!.....mass, prefactors
  real(8),allocatable:: acon(:),fack(:)
!.....atomic strain
  real(8),allocatable:: stn(:,:,:)
!.....atomic charge and electronegativity
  real(8),allocatable:: chg(:),chi(:)
!!$    real(8):: stn(3,3,namax)

!.....Shear stress
  real(8):: shrfx

!.....Barostat
  real(8):: phyd,ah(3,3),aht(3,3) &
       ,g(3,3,0:1),gt(3,3,0:1),gi(3,3),gg(3,3)

!.....FIRE parameters
  integer:: nmin_fire = 5
  real(8):: finc_fire = 1.1
  real(8):: fdec_fire = 0.5
  real(8):: alp0_fire = 0.1
  real(8):: falp_fire = 0.99
  real(8):: dtmax_fire = 10.0
! factor to be multiplied to dt to get dtmax_fire
  real(8):: dtmfctr_fire = 1000.0

!.....temperature distribution along x
  real(8),allocatable:: tdst(:)
  integer,allocatable:: nadst(:)

!.....space decomposition
  integer,allocatable:: lsb(:,:),lsex(:,:)
  integer:: nn(6),myparity(3),lsrc(6),nex(3),myx,myy,myz
  real(8):: sv(3,6),sorg(3),anxi,anyi,anzi

contains
  subroutine initialize_pmdvars(nspmax)
    integer,intent(in):: nspmax
    
    if( .not. allocated(tfac) ) then
      allocate(tfac(nspmax),ediff(nspmax),ediff0(nspmax),temp(nspmax)&
           ,ekl(nspmax),ndof(nspmax))
    endif
  end subroutine initialize_pmdvars
end module pmdvars
