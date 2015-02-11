  module variables
!-----------------------------------------------------------------------
!                        Time-stamp: <2015-02-10 17:26:53 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
    implicit none
    save

!=======================================================================
! PARAMETERS
!=======================================================================
!.....max. num. of atoms
    integer,parameter:: namax = 2000
!.....max. num. of species
    integer,parameter:: nismax= 9
!.....max. num. of boundary-particles
    integer,parameter:: nbmax = 2000
!.....max. num. of neighbors
    integer,parameter:: nnmax = 200

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
    integer:: natm,nb,ntot,nis
    real(8):: tcpu,tcpu1,tcpu2,tcom
    real(8):: epot,ekin,epot0,vmaxold,vmax
    integer:: nstp = 0
    integer:: nerg = 1000
    integer:: ifpmd= 1
    integer:: npmd = 10
    integer:: ifdmp= 0
    real(8):: dmp  = 0.9d0
    real(8):: dt = 1d0
    real(8):: rc = 5.0d0
    real(8):: rbuf= 0d0
    character(len=6):: ciofmt='ascii '
    character(len=20):: cforce='LJ_Ar'
!.....temperature
    character(len=20):: ctctl='none'
    integer:: iftctl= 0
    real(8):: tinit= -1d0
    real(8):: ttgt(9)
    data ttgt / 300d0, 300d0, 300d0, 300d0, 300d0, 300d0, &
         300d0, 300d0, 300d0 /
    real(8):: trlx = 100d0
    real(8):: tgmm,tfac(9),temp(9),ekl(9),ediff(9),ediff0(9)
    integer:: ndof(9)
!.....temperature distribution on x-direction
    logical:: ltdst= .false.
    integer:: ntdst= 1
    integer:: nxmlt
    real(8),allocatable:: tdst(:)
    integer,allocatable:: nadst(:)
!.....Search time and expiration time
    real(8):: ts,te
    integer:: istps,istpe
!.....parallel-related variables
    integer:: nxyz
    integer:: nn(6),myparity(3),lsrc(6),lsb(0:nbmax,6) &
         ,myid_md,nodes_md,mpi_md_world,myx,myy,myz,ierr
    real(8):: sv(3,6),sorg(3),anxi,anyi,anzi
    integer:: nx = 1
    integer:: ny = 1
    integer:: nz = 1
!.....simulation box
    real(8):: hunit,h(3,3,0:1),hi(3,3),vol,sgm(3,3),al(3),avol
    real(8):: ht(3,3,0:1),hti(3,3),dh
!.....factors on each moving direction
    real(8):: fmv(3,0:9)
    data fmv &
         / 0d0, 0d0, 0d0, &
           1d0, 1d0, 1d0, &
           1d0, 1d0, 1d0, &
           1d0, 1d0, 1d0, &
           1d0, 1d0, 1d0, &
           1d0, 1d0, 1d0, &
           1d0, 1d0, 1d0, &
           1d0, 1d0, 1d0, &
           1d0, 1d0, 1d0, &
           1d0, 1d0, 1d0 &
         /
!.....positions, velocities, and accelerations
    real(8):: ra(3,namax),va(3,namax),aa(3,namax),ra0(3,namax) &
         ,strs(3,3,namax),stt(3,3,namax)
!.....real*8 identifier which includes species, index of FMV, total id
    real(8):: tag(namax)
    integer:: lspr(0:nnmax,namax)
!.....potential and kinetic energy per atoms
    real(8):: epi(namax),eki(3,3,namax),stp(3,3,namax)
!.....mass, prefactors
    real(8):: acon(nismax),fack(nismax)
    real(8):: am(1:nismax)= 12.0d0
!.....strain
    real(8):: stn(3,3,namax)

!.....zload type
    character(len=5):: czload_type= 'no'
    real(8):: strfin = 0.0d0
!.....Shear stress
    real(8):: shrst = 0.0d0
    real(8):: shrfx


!.....Isobaric
    integer:: ifpctl= 0 ! 0:no  1:Parrinello-Rahman  2:Andersen
    real(8):: ptgt   = 0d0
    real(8):: vmcoeff= 1d0
    real(8):: voldmp = 1d0
    real(8):: stgt(1:3,1:3)= 0d0
    real(8):: phyd,vm,ah(3,3),aht(3,3),ptnsr(3,3) &
         ,g(3,3,0:1),gt(3,3,0:1),gi(3,3),gg(3,3)
    
  end module variables
