  module variables
    implicit none
    save

!.....max. num. of atoms
    integer,parameter:: namax = 100000
!.....max. num. of species
    integer,parameter:: nismax= 9
!.....max. num. of boundary-particles
    integer,parameter:: nbmax = 100000
!.....max. num. of neighbors
    integer,parameter:: nnmax = 200

    integer:: nstp,nouterg,noutpmd,istp,nerg &
         ,ifpmd,npmd,iftctl,ifdmp,iocntpmd,iocnterg
    integer:: natm,nb,ntot,nis
    real(8):: tcpu,tcpu1,tcpu2,tcom
    real(8):: dt,rc,dmp,treq,trlx,temp,epot,ekin,epot0,vmaxold,vmax
    real(8):: tinit= -1d0
    real(8):: rbuf= 0d0
    character(len=6):: ciofmt='ascii '
!.....Search time and expiration time
    real(8):: ts,te
    integer:: istps,istpe
!.....parallel-related variables
    integer:: nx,ny,nz,nxyz
    integer:: nn(6),myparity(3),lsrc(6),lsb(0:nbmax,6) &
         ,myid_md,nodes_md,mpi_md_world,myx,myy,myz,ierr
    real(8):: sv(3,6),sorg(3),anxi,anyi,anzi
!.....simulation box
    real(8):: h(3,3,0:1),hi(3,3),vol,sgm(3,3),al(3),avol
    real(8):: ht(3,3,0:1),hti(3,3),dh
!.....factors on each moving direction
    real(8):: fmv(3,0:9)
!.....positions, velocities, and accelerations
    real(8):: ra(3,namax),va(3,namax),aa(3,namax),ra0(3,namax) &
         ,strs(3,3,namax),stt(3,3,namax)
!.....real*8 identifier which includes species, index of FMV, total id
    real(8):: tag(namax)
    integer:: lspr(0:nnmax,namax)
!.....potential and kinetic energy per atoms
    real(8):: epi(namax),eki(3,3,namax),stp(3,3,namax)
!.....mass, prefactors
    real(8):: am(nismax),acon(nismax),fack(nismax)
!.....strain
    real(8):: stn(3,3,namax)

!.....Final strain value
    real(8):: strfin
!.....Shear stress
    real(8):: shrst,shrfx

!.....Isobaric
    integer:: ifpctl= 0 ! 0:no  1:Parrinello-Rahman  2:Andersen
    real(8):: ptgt   = 0d0
    real(8):: vmcoeff= 1d0
    real(8):: voldmp = 1d0
    real(8):: phyd,vm,ah(3,3),aht(3,3),ptnsr(3,3) &
         ,g(3,3,0:1),gt(3,3,0:1),gi(3,3),gg(3,3)
    
  end module variables
