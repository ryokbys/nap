!-----------------------------------------------------------------------
!                     Last-modified: <2016-11-04 14:07:42 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
module pmc
! 
! Module includes variables commonly used in pmc.
!
! Species IDs are fixed: 0) Vac, 1) Al, 2) Mg, 3) Si
!
!.....input file name
  character(len=128):: cinpfname = 'in.pmc'
  integer:: ionum_inp = 10
!.....kinetic or not
  logical:: lkinetic = .true.
!.....number of MC steps
  integer:: num_steps = 10
!.....size of the system in multiple of FCC cell
  integer:: ncx = 6
  integer:: ncy = 6
  integer:: ncz = 6
  integer:: natm
  real(8):: hmat(3,3),hmati(3,3)
  real(8),allocatable:: pos0(:,:),pos(:,:),tagmc(:),epimc(:)
!.....symbols and SIDs array
  character,allocatable:: csymbols(:)
  integer(1),allocatable:: i1sids(:)
!.....lattice constant of unit cell
  real(8):: alat = 4.0448d0
!.....number of solute atoms
  integer:: num_Mg = 10
  integer:: num_Si = 5
  integer:: num_Vac= 1
!.....initial structure: 1) random, 2) clustered,
!      3) read from file (restart)
  integer:: init_strct = 2
!.....maximum steps for relaxation in pmd
  integer:: nstps_relax = 10
!.....atoms to be moved
  character:: species(0:3) = (/ 'V','A','M','S' /)
  logical:: lmove(0:3) = (/ &
       .true., &
       .false., &
       .false., &
       .false. /)
!.....frequency prefactors in 1/sec
!.....values taken from Mantina et al., Acta Mater. 57 (2009)
  real(8):: prefreq(1:3) = (/ &
       16.6d+12, &
       18.6d+12, &
       15.7d+12 /)
!.....average migration barriers in eV
!.....also taken from Mantina et al.
  real(8):: demig(1:3) = (/ &
       0.58d0, &
       0.42d0, &
       0.55d0 /)
!.....pair-list for MC only, not to conflict with pmd
  integer,parameter:: nnmaxmc = 20
  integer,allocatable:: lsprmc(:,:)
!.....parallel setting for pmd
  integer:: nx = 1
  integer:: ny = 1
  integer:: nz = 1
  
contains
!=======================================================================
  subroutine symbols2sids(natm,csymbols,i1sids)
!
!  Change symbol array to SID array
!
    implicit none 
    integer,intent(in):: natm
    character,intent(in):: csymbols(natm)
    integer(1),intent(out):: i1sids(natm)
!.....local
    integer:: i
    character(len=1):: c

    do i = 1,natm
      c = csymbols(i)
      i1sids(i) = symbol2sid(c)
    enddo
    return

  end subroutine symbols2sids
!=======================================================================
  function symbol2sid(csymbol) result(sid)
!
!  Symbol to SID converter for the system of AlMgSi.
!  This function should be changed once you change the components.
!
    implicit none 
    character(len=1),intent(in):: csymbol
    integer(1):: sid
    if( csymbol == 'V' ) then
      sid = 0
    elseif( csymbol == 'A' ) then
      sid = 1
    elseif( csymbol == 'M' ) then
      sid = 2
    elseif( csymbol == 'S' ) then
      sid = 3
    else
      stop 'There is no symbol specified by sid.'
    endif
    return
  end function symbol2sid
end module pmc
!=======================================================================
program prec_mc
!-----------------------------------------------------------------------
! (Kinetic) MC simulation using parallel computation of force
! calculation used in pmd.
!-----------------------------------------------------------------------
! INPUT FILES:
! ------------
!   in.pmc:  Input parameters for this MC
!
! OUTPUT FILES:
! -------------
!   dat.erg:    MC step, real time, epot
!   dat.symbols: MC step, symbol array 
!   POSCAR_#####: Cell info and atom coordinations of a certain steps.
!-----------------------------------------------------------------------
  use pmc
  implicit none 
  include "mpif.h"
  include "./params_unit.h"

  integer:: i,j,k,l,m,n,ierr
  integer:: mpi_md_world,nodes_md,myid_md,myx,myy,myz
  real(8):: rc,anxi,anyi,anzi,sorg(3)

!.....initialize parallel
  call init_parallel(mpi_md_world,nodes_md,myid_md)
  if( myid_md.eq.0 ) write(6,'(a)') ' program pmc starts...'
!.....read input parameters
  call read_in_pmc(ionum_inp, cinpfname)
!.....parallel setting for pmd
  call parallel_setting(nx,ny,nz,myid_md,myx,myy,myz,anxi,anyi,anzi,sorg)

!.....create 1st Al fcc crystal
  natm = ncx*ncy*ncz*4
  allocate(csymbols(natm),i1sids(natm),pos0(3,natm),pos(3,natm) &
       ,tagmc(natm),epimc(natm))
  call create_Al_fcc(ncx,ncy,ncz,natm,alat,pos0,csymbols,hmat)
  call symbols2sids(natm,csymbols,i1sids)
  hmati(1:3,1:3) = 0d0
  hmati(1,1) = 1d0/hmat(1,1)
  hmati(2,2) = 1d0/hmat(2,2)
  hmati(3,3) = 1d0/hmat(3,3)

!.....create neighbor list only once here, and no longer needed after
  rc = 3.0
  allocate(lsprmc(0:nnmaxmc,natm))
  call make_tag(natm,csymbols,tagmc)
  call mk_lspr_sngl(natm,natm,nnmaxmc,tagmc,pos0,rc,hmat,hmati,lsprmc)
!.....check restart,
!     if not restart, initialize system eigher random or clustered
  if( init_strct.eq.1 ) then  ! random
    call random_symbols(natm,csymbols,num_Mg,num_Si,num_Vac)
  elseif( init_strct.eq.2 ) then  ! clustered
    call clustered_symbols(natm,pos0,csymbols,num_Mg,num_Si,num_Vac &
         ,nnmaxmc,lsprmc)
  elseif( init_strct.eq.3 ) then  ! read from file (restart)
    call read_symbols(11,'dat.symbols',natm,csymbols)
  endif

  call write_POSCAR('POSCAR_000000',natm,csymbols,pos0,hmat,species)
  
  if( myid_md.eq.0 ) write(6,'(a)') ' program pmc ends'
  call mpi_finalize(ierr)

end program prec_mc
!=======================================================================
subroutine kinetic_mc()
!
! Kinetic MC simulation using
!
  use pmc
  
!.....initialize some values here
  
!.....compute chemical potentials of solutes

!.....compute energy of the initial configuration

!.....main loop starts..................................................

!.....store previous symbols

!.....pick one solutes to be moved (usually vacancy)

!.....look for neighbor site and check if they are different species

!.....loop for possible neighbor sites and compute migration barriers.
!.....to reduce computation cost, look at history of symbol array
!.....if the symbol matches one of the history, no need to calc energy
!.....of the system and just take from the history.

!.....compute the every and total probabilities

!.....pick one event from the event list

!.....proceed real-time clock

!.....output if needed
  
  return
end subroutine kinetic_mc
!=======================================================================
subroutine create_Al_fcc(nx,ny,nz,natm,alat,pos,csymbols,hmat)
!
! Create Al fcc crystalline structure as an initial template.
!
  implicit none 
!.....arguments
  integer,intent(in):: nx,ny,nz,natm
  real(8),intent(in):: alat
  real(8),intent(out):: pos(3,natm),hmat(3,3)
  character,intent(out):: csymbols(natm)
!.....local variables
  integer:: ix,iy,iz,m,inc
  real(8):: upos(3,4)

!.....set h-matrix
  hmat(1:3,1:3) = 0d0
  hmat(1,1) = alat*nx
  hmat(2,2) = alat*ny
  hmat(3,3) = alat*nz

!.....positions of unit cell of FCC
  upos(1:3,1) = (/ 0.0d0, 0.0d0, 0.0d0 /)
  upos(1:3,2) = (/ 0.0d0, 0.5d0, 0.5d0 /)
  upos(1:3,3) = (/ 0.5d0, 0.0d0, 0.5d0 /)
  upos(1:3,4) = (/ 0.5d0, 0.5d0, 0.0d0 /)

!.....extend unit cell with nx,ny,nz
  inc = 0
  do ix = 0,nx-1
    do iy = 0,ny-1
      do iz = 0,nz-1
        do m = 1,4
          inc = inc + 1
          if( inc.gt.natm ) then
            print *, 'Error: inc.gt.natm'
            stop
          endif
          pos(1,inc) = dble(upos(1,m) + ix)/nx
          pos(2,inc) = dble(upos(2,m) + iy)/ny
          pos(3,inc) = dble(upos(3,m) + iz)/nz
          csymbols(inc) = 'A'
        enddo
      enddo
    enddo
  enddo
  return
  
end subroutine create_Al_fcc
!=======================================================================
subroutine read_in_pmc(ionum,cfname)
!
! Read frexible input format
!
  implicit none
  integer,intent(in):: ionum
  character(len=*),intent(in):: cfname
  character(len=128):: c1st

  write(6,'(a)') ' reading '//trim(cfname)//'...'

  open(ionum,file=trim(cfname))
  do
!.....Read 1st word in each line
    read(ionum,*,end=10) c1st
!.....Skip comment line
    if( c1st(1:1).eq.'!' .or. &
         c1st(1:1).eq.'#' .or. &
!.....Skip lines starting from digits or sign
         c1st(1:1).eq.'0' .or. &
         c1st(1:1).eq.'1' .or. &
         c1st(1:1).eq.'2' .or. &
         c1st(1:1).eq.'3' .or. &
         c1st(1:1).eq.'4' .or. &
         c1st(1:1).eq.'5' .or. &
         c1st(1:1).eq.'6' .or. &
         c1st(1:1).eq.'7' .or. &
         c1st(1:1).eq.'8' .or. &
         c1st(1:1).eq.'9' .or. &
         c1st(1:1).eq.'+' .or. &
         c1st(1:1).eq.'-' ) cycle
!        write(6,'(a)') c1st
    call read_in_pmc_core(ionum,c1st)
  enddo
  close(ionum)
10 write(6,'(a)') ' read '//trim(cfname)//' done'
end subroutine read_in_pmc
!=======================================================================
subroutine read_in_pmc_core(ionum,cname)
  use pmc
  implicit none
  integer,intent(in):: ionum
  character(len=*) ,intent(in):: cname
  
  character(len=128):: ctmp
  integer:: ndata,nrow,is,itmp
  
  if( trim(cname).eq.'kinetic' ) then
    call read_l1(ionum,lkinetic)
    return
  elseif( trim(cname).eq.'num_steps' ) then
    call read_i1(ionum,num_steps)
    return
  elseif( trim(cname).eq.'ncopy_x' ) then
    call read_i1(ionum,ncx)
    return
  elseif( trim(cname).eq.'ncopy_y' ) then
    call read_i1(ionum,ncy)
    return
  elseif( trim(cname).eq.'ncopy_z' ) then
    call read_i1(ionum,ncz)
    return
  elseif( trim(cname).eq.'lattice_constant' ) then
    call read_r1(ionum,alat)
    return
  elseif( trim(cname).eq.'num_Mg' ) then
    call read_i1(ionum,num_Mg)
    return
  elseif( trim(cname).eq.'num_Si' ) then
    call read_i1(ionum,num_Si)
    return
  elseif( trim(cname).eq.'num_Vac' ) then
    call read_i1(ionum,num_Vac)
    return
  elseif( trim(cname).eq.'initial_structure' ) then
    call read_i1(ionum,init_strct)
    return
  elseif( trim(cname).eq.'parallel_x' ) then
    call read_i1(ionum,nx)
    return
  elseif( trim(cname).eq.'parallel_y' ) then
    call read_i1(ionum,ny)
    return
  elseif( trim(cname).eq.'parallel_z' ) then
    call read_i1(ionum,nz)
    return
  elseif( trim(cname).eq.'num_relax_steps' ) then
    call read_i1(ionum,nstps_relax)
    return
  elseif( trim(cname).eq.'atoms_to_be_moved' ) then
    backspace(ionum)
    read(ionum,*) ctmp, lmove(0:3)
    return
  endif
end subroutine read_in_pmc_core
!=======================================================================
subroutine init_parallel(mpi_world,nodes,myid)
  implicit none
  include "mpif.h"
  integer,intent(out):: mpi_world,nodes,myid
  
  integer:: ierr
  
!.....initialize the MPI environment
  call mpi_init(ierr)
!.....total number of MD-nodes
  call mpi_comm_size(mpi_comm_world, nodes, ierr)
!.....my rank in MD-nodes
  call mpi_comm_rank(mpi_comm_world, myid, ierr)
  mpi_world= mpi_comm_world

  return
end subroutine init_parallel
!=======================================================================
subroutine parallel_setting(nx,ny,nz,myid_md,myx,myy,myz &
     ,anxi,anyi,anzi,sorg)
  implicit none
  integer,intent(in):: nx,ny,nz,myid_md
  integer,intent(out):: myx,myy,myz
  real(8),intent(out):: anxi,anyi,anzi,sorg(3)

  anxi= 1d0/nx
  anyi= 1d0/ny
  anzi= 1d0/nz
!.....vector node indices: range [0:nx-1]
  myx=myid_md/(ny*nz)
  myy=mod(myid_md/nz,ny)
  myz=mod(myid_md,nz)
!.....reduced node origin
  sorg(1)= anxi*myx
  sorg(2)= anyi*myy
  sorg(3)= anzi*myz
  return

end subroutine parallel_setting
!=======================================================================
subroutine make_tag(natm,csymbols,tag)
  use pmc, only: symbol2sid
  implicit none
  integer,intent(in):: natm
  character,intent(in):: csymbols(natm)
  real(8),intent(out):: tag(natm)

  integer:: i
  character(len=1):: c

  do i = 1,natm
    c = csymbols(i)
    tag(i) = dble(symbol2sid(c)) +0.1d0 +1d-14*i
  end do
  return
  
end subroutine make_tag
!=======================================================================
subroutine random_symbols(natm,csymbols,num_Mg,num_Si,num_Vac)
  implicit none
  integer,intent(in):: natm,num_Mg,num_Si,num_Vac
  character,intent(inout):: csymbols(natm)

  integer:: i,img,isi,ivac,irnd
  real(8),external:: urnd

  img = num_Mg
  isi = num_Si
  ivac= num_Vac

!.....Mg
  do while(.true.)
    irnd = int(natm*urnd())+1
    if( csymbols(irnd).eq.'A' ) then
      csymbols(irnd) = 'M'
      img = img - 1
    endif
    if( img.eq.0 ) exit
  enddo
!.....Si
  do while(.true.)
    irnd = int(natm*urnd())+1
    if( csymbols(irnd).eq.'A' ) then
      csymbols(irnd) = 'S'
      isi = isi - 1
    endif
    if( isi.eq.0 ) exit
  enddo
!.....Vac
  do while(.true.)
    irnd = int(natm*urnd())+1
    if( csymbols(irnd).eq.'A' ) then
      csymbols(irnd) = 'V'
      ivac = ivac - 1
    endif
    if( ivac.eq.0 ) exit
  enddo
  return
  
end subroutine random_symbols
!=======================================================================
subroutine clustered_symbols(natm,pos0,csymbols &
     ,num_Mg,num_Si,num_Vac,nnmax,lspr)
  implicit none
  integer,intent(in):: natm,num_Mg,num_Si,num_Vac &
       ,nnmax,lspr(0:nnmax,natm)
  character,intent(inout):: csymbols(natm)
  real(8),intent(in):: pos0(3,natm)

  integer:: i,jj,j,irnd,nsol,isol,icntr
  real(8):: cntr(3),dmin,d
  real(8),external:: urnd
  character,allocatable:: carr(:) 

!.....1st, pick one site close to the center
  cntr(1:3) = (/ 0.5d0, 0.5d0, 0.5d0 /)
  icntr = -1
  dmin = 1d+30
  do i = 1,natm
    d = (cntr(1)-pos0(1,i))**2 &
         +(cntr(2)-pos0(2,i))**2 &
         +(cntr(3)-pos0(3,i))**2
    if( d < dmin ) then
      dmin = d
      icntr = i
    endif
  enddo
  
!.....make random array of symbols to be replaced
  nsol = num_Mg +num_Si +num_Vac
  allocate(carr(nsol))
  do i = 1,nsol
    if( i.le.num_Mg ) then
      carr(i) = 'M'
    else if( i.le.num_Mg+num_Si) then
      carr(i) = 'S'
    else
      carr(i) = 'V'
    endif
  enddo

  isol = 1
  csymbols(icntr) = carr(isol)

  do while(.true.)
    do i = 1,natm
      if( csymbols(i).eq.'A' ) cycle
      do jj = 1,lspr(0,i)
        j = lspr(jj,i)
        if( csymbols(j).eq.'A' ) then
          isol = isol + 1
          csymbols(j) = carr(isol)
        endif
        if( isol.eq.nsol) return
      enddo
    enddo
  enddo
  stop 'Error: something wrong @clustered_symbols'
  
end subroutine clustered_symbols
!=======================================================================
subroutine loads_symbols(natm,txt,csymbols)
!
! Load symbols from a character string
!
  implicit none
  integer,intent(in):: natm
  character,intent(in):: txt(natm)
  character,intent(out):: csymbols(natm)

  integer:: i

  do i = 1,natm
    csymbols(i) = txt(i)
  enddo
  return
  
end subroutine loads_symbols
!=======================================================================
subroutine read_symbols(ionum,fname,natm,csymbols)
  implicit none
  integer,intent(in):: ionum,natm
  character(len=*),intent(in):: fname
  character,intent(out):: csymbols(natm)

  character:: txt(natm)
  integer:: itmp,ios

!.....read the last symbols from file
  open(ionum,file='fname',status='old')
  do while(.true.)
    read(ionum,iostat=ios) itmp, txt
    if( ios.lt.0 ) exit
  end do
  close(ionum)

  call loads_symbols(natm,txt,csymbols)
  return
  
end subroutine read_symbols
!=======================================================================
subroutine write_POSCAR(cfname,natm,csymbols,pos,hmat,species)
!
!  Write system config in POSCAR format.
!  Vacancies are written as an atom of species Vanadium in order to
!  make them visible on purpose.
!
  implicit none
  character(len=*),intent(in):: cfname
  integer,intent(in):: natm
  character,intent(in):: csymbols(natm),species(0:3)
  real(8),intent(in):: pos(3,natm),hmat(3,3)

  integer:: i,m,id,ns(0:3),idorder(natm)
  integer:: date_time(8)
  character(len=10):: sys_time(3)
  character(len=1):: cs
  character(len=2):: sname(0:3) = (/ 'V ','Al','Mg','Si' /)
  
  call date_and_time(sys_time(1), sys_time(2), sys_time(3), date_time)

  id = 0
  do m = 0,3
    cs = species(m)
    ns(m) = 0
    do i = 1,natm
      if( csymbols(i).eq.cs ) then
        id = id + 1
        idorder(id) = i
        ns(m) = ns(m) + 1
      endif
    end do
  enddo
  
  open(90,file=trim(cfname))
  write(90,'(3a)') 'written by PMC at ', &
       trim(sys_time(1)),trim(sys_time(3))
  write(90,'(a)') '   1.00000000'
  write(90,'(3es15.7)') hmat(1:3,1)
  write(90,'(3es15.7)') hmat(1:3,2)
  write(90,'(3es15.7)') hmat(1:3,3)
  write(90,'(4a8)') sname(0:3)
  write(90,'(4i8)') ns(0:3)
  write(90,'(a)') 'Direct'
  do m = 1,natm
    i = idorder(m)
    write(90,'(3f15.8)') pos(1:3,i)
  enddo
  close(90)
  
end subroutine write_POSCAR
!=======================================================================
subroutine run_pmd(hmat,natm,pos0,csymbols,epimc,epotmc &
     ,nstps_pmd,nx,ny,nz,mpi_md_world,nodes_md,myid_md)
  implicit none
  integer,intent(in):: natm,nstps_pmd,nx,ny,nz&
       ,mpi_md_world,nodes_md,myid_md
  real(8),intent(in):: hmat(3,3),pos0(3,natm)
  real(8),intent(out):: epimc(natm),epotmc
  character,intent(in):: csymbols(natm)

  integer:: i,inc
  integer,parameter:: nismax = 9
  integer:: nstp,nerg,npmd,ifpmd,minstp,ntdst,n_conv,ifsort,iprint &
       ,ifdmp
  real(8):: hunit,h(3,3,0:1),am(nismax),dt,rc,dmp,tinit,tfin,ttgt(9)&
       ,trlx,stgt(3,3),ptgt,srlx,stbeta,strfin,fmv(3,0:9),ptnsr(3,3) &
       ,epot,ekin,eps_conv,rbuf
  character:: ciofmt*6,cforce*20,ctctl*20,cpctl*20,czload_type*5,csi*1
  logical:: ltdst
  
  logical,save:: l1st = .true.
  integer,save:: ntot
  real(8),save,allocatable:: tagtot(:),rtot(:,:),vtot(:,:),atot(:,:) &
       ,epitot(:),ekitot(:,:,:),stot(:,:,:)

!.....at the 1st call, evaluate number of total atoms to be used in pmd
!     and allocate total system arrays
  if( l1st ) then
    inc = 0
    do i = 1,natm
      csi = csymbols(i)
      if( csi.eq.'V' ) cycle
      inc = inc + 1
    enddo
    ntot = inc
    allocate(tagtot(ntot),rtot(3,ntot),vtot(3,ntot),atot(3,ntot) &
         ,epitot(ntot),ekitot(3,3,ntot),stot(3,3,ntot))
  endif

  call make_tag(natm,csymbols,tagtot)
  hunit = 1d0
  h(1:3,1:3,0) = hmat(1:3,1:3)
  
  inc = 0
  do i = 1,natm
    csi = csymbols(i)
    if( csi.eq.'V' ) cycle
    inc = inc + 1
    rtot(1:3,inc) = pos0(1:3,i)
    vtot(1:3,inc) = 0d0
    atot(1:3,inc) = 0d0
    stot(1:3,1:3,inc) = 0d0
    ekitot(1:3,1:3,inc) = 0d0
    epitot(inc) = 0d0
  enddo
  nstp = ntot
  nerg = nstp
  npmd = 1
  am(1:9) = 1d0
  am(1) = 26.982  ! Al
  am(2) = 24.305  ! Mg
  am(3) = 28.085  ! Si
  dt = 5d0
  ciofmt = 'ascii'
  ifpmd = 1
  cforce = 'NN'
  rc = 5.8d0
  rbuf = 0.2d0
  ifdmp = 2  ! FIRE
  dmp = 0.99d0
  minstp = 3
  tinit = 100d0
  tfin = 1d0
  ctctl = 'none'
  ttgt(1:9) = 300d0
  trlx = 100d0
  ltdst = .false.
  ntdst = 1
  cpctl = 'none'
  stgt(1:3,1:3) = 0d0
  ptgt = 0d0
  srlx = 100d0
  stbeta = 1d-1
  strfin = 0d0
  fmv(1:3,0) = (/ 0d0, 0d0, 0d0 /)
  fmv(1:3,1:9) = 1d0
  ptnsr(1:3,1:3) = 0d0
  epot = 0d0
  ekin = 0d0
  n_conv = 1
  czload_type = 'no'
  eps_conv = 1d-3
  ifsort = 1
  iprint = 1

!.....call pmd_core to perfom MD
  call pmd_core(hunit,h,ntot,tagtot,rtot,vtot,atot,stot &
       ,ekitot,epitot,nstp,nerg,npmd,myid_md,mpi_md_world,nodes_md &
       ,nx,ny,nz &
       ,nismax,am,dt,ciofmt,ifpmd,cforce,rc,rbuf,ifdmp,dmp,minstp &
       ,tinit,tfin,ctctl,ttgt,trlx,ltdst,ntdst,cpctl,stgt,ptgt &
       ,srlx,stbeta,strfin &
       ,fmv,ptnsr,epot,ekin,n_conv &
       ,czload_type,eps_conv,ifsort,iprint)
  
  
end subroutine run_pmd
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmc"
!     End:
