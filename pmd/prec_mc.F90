module pmc
!-----------------------------------------------------------------------
!                     Last-modified: <2018-10-08 12:51:54 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
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
  integer:: nstps_mc = 10
  integer:: noutint  = 10
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
  integer:: num_Al_clst= 0
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
!!$  real(8):: demig(1:3) = (/ &
!!$       0.58d0, &
!!$       0.42d0, &
!!$       0.55d0 /)
!.....obtained by own DFT calculations
  real(8):: demig(1:3) = (/ &
       0.569d0, &
       0.450d0, &
       0.479d0 /)
!.....pair-list for MC only, not to conflict with pmd
  integer,parameter:: nnmaxmc = 20
  integer,allocatable:: lsprmc(:,:)
!.....parallel setting for pmd
  integer:: nx = 1
  integer:: ny = 1
  integer:: nz = 1
!.....temperature
  real(8):: temp = 300d0
!.....random seed
  real(8):: dseed0 = 12345d0
  
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
!   out.mc.erg:     MC step, real time, epot
!   out.mc.symbols: MC step, symbol array 
!   poscars/POSCAR_######: Cell info and atom coordinations
!-----------------------------------------------------------------------
  use pmc
  implicit none 
  include "mpif.h"
  include "./params_unit.h"

  interface
    function urnd(dseed0)
      real(8),intent(in),optional:: dseed0
      real(8):: urnd
    end function urnd
  end interface

  integer:: i,j,k,l,m,n,ierr
  integer:: ihour,imin,isec,iday
  integer:: mpi_md_world,nodes_md,myid_md,myx,myy,myz
  real(8):: rc,anxi,anyi,anzi,sorg(3),t0,t1
  character:: cnum*128

!.....initialize parallel
  call init_parallel(mpi_md_world,nodes_md,myid_md)
  t0 = mpi_wtime()
  if( myid_md.eq.0 ) then
    write(6,'(a)') ' program pmc starts...'
!.....read input parameters
    call read_in_pmc(ionum_inp, cinpfname)
  endif
!.....set random seed
  t1 = urnd(dseed0)
  
  call bcast_params(myid_md,mpi_md_world,lkinetic, &
       nstps_mc,ncx,ncy,ncz,alat,num_Mg,num_Si,num_Vac, &
       init_strct,nx,ny,nz,nstps_relax,lmove,temp,num_Al_clst)
!.....parallel setting for pmd
  call parallel_setting(nx,ny,nz,myid_md,myx,myy,myz,anxi,anyi,anzi,sorg)

  if( myid_md.eq.0 ) then
    call write_init_params(lkinetic,nstps_mc, ncx,ncy,ncz,alat, &
         num_Mg,num_Si,num_Vac,num_Al_clst, &
         init_strct,nx,ny,nz,nstps_relax,lmove,temp)
  endif
  
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

  if( myid_md.eq.0 ) then
    write(6,'(a,i8)') ' Number of total atoms = ',natm
  endif

!.....create neighbor list only once here, and no longer needed after
  rc = 3.0
  allocate(lsprmc(0:nnmaxmc,natm))
  call make_tag(natm,csymbols,tagmc)
  if( myid_md.eq.0 ) then
    write(cnum,'(i0)') 0
    call write_POSCAR('poscars/POSCAR_'//trim(cnum),natm,csymbols,pos0,hmat,species)
  endif
  call mk_lspr_sngl(natm,natm,nnmaxmc,tagmc,pos0,rc,hmat,hmati,lsprmc, &
       0,.true.)

!.....check restart,
!     if not restart, initialize system eigher random or clustered
  if( myid_md.eq.0 ) then
    if( init_strct.eq.0 ) then  ! read from file (restart)
      call read_symbols(11,'dat.symbols',natm,csymbols)
    elseif( init_strct.eq.1 ) then  ! random
      call random_symbols(natm,csymbols,num_Mg,num_Si,num_Vac)
    elseif( init_strct.eq.2 ) then  ! clustered
      call clustered_symbols(natm,pos0,csymbols,num_Mg,num_Si,num_Vac, &
           nnmaxmc,lsprmc,num_Al_clst)
    elseif( init_strct.eq.3 ) then  ! pairs between same species
      call paired_symbols(natm,pos0,csymbols,num_Mg,num_Si,num_Vac,&
           nnmaxmc,lsprmc,'XX')
    elseif( init_strct.eq.4 ) then  ! pairs between different species
      call paired_symbols(natm,pos0,csymbols,num_Mg,num_Si,num_Vac,&
           nnmaxmc,lsprmc,'XY')
    endif
  endif

  call kinetic_mc(mpi_md_world,nodes_md,myid_md,myx,myy,myz &
       ,nx,ny,nz,anxi,anyi,anzi,sorg, hmat,natm,pos0,csymbols&
       ,nstps_mc,nstps_relax,noutint,nnmaxmc,lsprmc,species,temp &
       ,demig,prefreq)

  t1 = mpi_wtime() - t0
  if( myid_md.eq.0 ) then
    call write_POSCAR('poscars/POSCAR_final',natm,csymbols,pos0,hmat,species)
    iday  = int(t1/86400)
    ihour = int(t1/3600)
    imin  = int((t1-ihour*3600)/60)
    isec  = int(t1 -ihour*3600 -imin*60)
    write(6,*) ''
    if( iday.gt.0 ) then
      write(6,'(a,f12.1,a,i3,"D",i2,"h",i2.2,"m",i2.2,"s")') &
           " time =",t1, &
           " sec  = ",iday,ihour,imin,isec
    else
      write(6,'(a,f14.1,a,i5,"h",i2.2,"m",i2.2,"s")') &
           " time =",t1, &
           " sec  = ",ihour,imin,isec
    endif
  endif

  call mpi_finalize(ierr)

end program prec_mc
!=======================================================================
subroutine bcast_params(myid_md,mpi_md_world,lkinetic, &
     nstps_mc, ncx,ncy,ncz,alat,num_Mg,num_Si,num_Vac, &
     init_strct,nx,ny,nz,nstps_relax,lmove,temp,num_Al_clst)
  implicit none
  include 'mpif.h'
  integer,intent(in):: myid_md,mpi_md_world
  integer,intent(inout):: nx,ny,nz,ncx,ncy,ncz,nstps_mc, &
       num_Mg,num_Si,num_Vac,init_strct,nstps_relax,num_Al_clst
  real(8),intent(inout):: alat,temp
  logical,intent(inout):: lmove(0:3),lkinetic
  
  integer:: ierr

  call mpi_bcast(nstps_mc,1,mpi_integer,0,mpi_md_world,ierr)
  call mpi_bcast(nstps_relax,1,mpi_integer,0,mpi_md_world,ierr)

  call mpi_bcast(nx,1,mpi_integer,0,mpi_md_world,ierr)
  call mpi_bcast(ny,1,mpi_integer,0,mpi_md_world,ierr)
  call mpi_bcast(nz,1,mpi_integer,0,mpi_md_world,ierr)

  call mpi_bcast(ncx,1,mpi_integer,0,mpi_md_world,ierr)
  call mpi_bcast(ncy,1,mpi_integer,0,mpi_md_world,ierr)
  call mpi_bcast(ncz,1,mpi_integer,0,mpi_md_world,ierr)

  call mpi_bcast(num_Mg,1,mpi_integer,0,mpi_md_world,ierr)
  call mpi_bcast(num_Si,1,mpi_integer,0,mpi_md_world,ierr)
  call mpi_bcast(num_Vac,1,mpi_integer,0,mpi_md_world,ierr)
  call mpi_bcast(num_Al_clst,1,mpi_integer,0,mpi_md_world,ierr)

  call mpi_bcast(init_strct,1,mpi_integer,0,mpi_md_world,ierr)

  call mpi_bcast(lkinetic,1,mpi_logical,0,mpi_md_world,ierr)
  call mpi_bcast(lmove,4,mpi_logical,0,mpi_md_world,ierr)

  call mpi_bcast(alat,1,mpi_double_precision,0,mpi_md_world,ierr)
  call mpi_bcast(temp,1,mpi_double_precision,0,mpi_md_world,ierr)

end subroutine bcast_params
!=======================================================================
subroutine kinetic_mc(mpi_md_world,nodes_md,myid_md,myx,myy,myz &
     ,nx,ny,nz,anxi,anyi,anzi,sorg, hmat,natm,pos0,csymbols &
     ,nstps_mc,nstps_relax,noutint,nnmaxmc,lsprmc,species,temp &
     ,demig,prefreq)
!
! Kinetic MC simulation using
!
  implicit none
  include 'mpif.h'
  integer,intent(in):: mpi_md_world,nodes_md,myid_md,myx,myy,myz &
       ,nx,ny,nz,natm,nstps_mc,nstps_relax,noutint &
       ,nnmaxmc,lsprmc(0:nnmaxmc,natm)
  real(8),intent(in):: anxi,anyi,anzi,sorg(3),hmat(3,3),pos0(3,natm), &
       temp,demig(3),prefreq(3)
  character,intent(in):: species(0:3)
  character,intent(inout):: csymbols(natm) 

  integer:: i,ic,ievent,ihist,iorder,istp,jc,jj,js,ncalc,ncandidate, &
       nhist,nspcs,ierr,isc
  integer:: nstps_pmd,maxhist,mem,nstps_done
  real(8):: epotmc,epotmc0,de,dt,epot,ergp,ptmp,ptot,rand, &
       tclck,p,efrm,efrm0
  real(8),allocatable:: epimc(:),ecpot(:),erghist(:),ergtmp(:), &
       probtmp(:)
  integer,allocatable:: nstptmp(:)
  character:: ci*1,cj*1,cfmt*10,cergtxt*1024,cnum*128,csi*1
  character,allocatable:: csymprev(:),csymhist(:,:),csymtmp(:,:) &
       ,cjtmp(:)
  integer,external:: cs2is,check_history
  real(8),external:: epot2efrm
  interface
    function urnd(dseed0)
      real(8),intent(in),optional:: dseed0
      real(8):: urnd
    end function urnd
  end interface


  real(8),parameter:: fkb = 8.61733035d-5  ! eV/K
  integer,parameter:: ioerg = 30
  integer,parameter:: iosym = 31

  maxhist = nstps_mc * 12 + 1
  if( myid_md.eq.0 ) then
    allocate(epimc(natm),ecpot(0:3),csymprev(natm), &
         csymhist(natm,maxhist),erghist(maxhist), &
         csymtmp(natm,nnmaxmc),ergtmp(nnmaxmc), &
         cjtmp(nnmaxmc),probtmp(nnmaxmc),nstptmp(nnmaxmc))
    mem = 8*(natm +4 +maxhist +nnmaxmc +nnmaxmc) &
         +1*(natm +natm*maxhist +natm*nnmaxmc +nnmaxmc)
    write(6,'(a,f10.3,a)') ' allocated array in kinetic_mc = ', &
         dble(mem)/1000/1000, ' MB'
  else
    allocate(epimc(natm),ecpot(0:3),csymprev(natm), &
         csymhist(natm,1),erghist(1), &
         csymtmp(natm,1),ergtmp(1), &
         cjtmp(1),probtmp(1),nstptmp(1))
  endif
  

!.....test run pmd

!!$  print *,'natm,nstps_relax, =',natm,nstps_relax
!!$  print *,'nx,ny,nz = ',nx,ny,nz
!!$  call run_pmd(hmat,natm,pos0,csymbols,epimc,epotmc &
!!$       ,nstps_relax,nx,ny,nz,mpi_md_world,nodes_md,myid_md)
!!$  print *,'epotmc = ',epotmc
!!$  print *,'epimc(1) =',epimc(1)
!!$  print *,'epimc(N) =',epimc(natm)

!.....initialize some values here
  nstps_pmd = nstps_relax
  if( nstps_pmd.lt.0 ) then
    nstps_pmd = natm
  endif

!.....Open output files
  if( myid_md.eq.0 ) then
    open(ioerg,file='out.mc.erg',status='replace')
    open(iosym,file='out.mc.symbols',status='replace')
  endif

  nspcs = 3
!.....compute chemical potentials of solutes
  call calc_chem_pot(nspcs,species,ecpot,hmat,natm,pos0 &
       ,nstps_relax,nx,ny,nz,mpi_md_world,nodes_md,myid_md)  

!.....compute energy of the initial configuration
  call run_pmd(hmat,natm,pos0,csymbols,epimc,epotmc &
       ,nstps_pmd,nx,ny,nz,mpi_md_world,nodes_md,myid_md,nstps_done)
  efrm = epot2efrm(natm,ecpot,csymbols,epotmc)
  efrm0 = efrm
  ergp = efrm0

  if( myid_md.eq.0 ) then
    write(cfmt,'(i10)') natm
    write(iosym,'(i8,3x,'//trim(cfmt)//'a)') 0,csymbols(1:natm)
    flush(iosym)
    write(6,*) ''
    write(6,'(a)') ' chemical potentials:'
    do i=0,3
      write(6,'(1x,i4,a2,f12.3)') i, species(i), ecpot(i)
    enddo
    write(6,*) ''
    write(6,'(a)') ' base migration barriers (eV) and prefactors (Hz):'
    do i=1,3
      write(6,'(1x,i4,a2,f12.3,es12.3)') i, species(i),demig(i),prefreq(i)
    enddo
    write(6,*) ''
    write(6,'(a,es15.7)') ' initial formation energy = ',efrm0
    write(6,*) ''
!.....Register initial structure to history list
    nhist = 1
    csymhist(1:natm,nhist) = csymbols(1:natm)
    erghist(nhist) = efrm

!.....main loop starts..................................................
    do istp = 1, nstps_mc

!.....store previous symbols
      csymprev(1:natm) = csymbols(1:natm)

!.....pick solutes/vacancies to be moved
      do i=1,natm
        if( csymbols(i).eq.'V' ) then  ! only vacancies
          ic = i
          exit
        endif
      enddo
      ci = csymbols(ic)

!.....Look for neighbor site and check if they are different species
      ncandidate = 0
      do jj=1,lsprmc(0,ic)
        jc = lsprmc(jj,ic)
        cj = csymbols(jc)
        if( cj.ne.ci ) then
          ncandidate = ncandidate + 1
        endif
      enddo

!.....Loop for possible neighbor sites and compute migration barriers.
!.....To reduce computation cost, look at history of symbol array.
!.....If the symbol matches one of the history, no need to calc energy
!.....of the system and just take from the history.
      ncalc = 0
      do jj=1,lsprmc(0,ic)
        jc = lsprmc(jj,ic)
        cj = csymbols(jc)
        if( cj.eq.ci ) cycle
!.....Exchange positions or migrate vacancy
!!$        print *,' jj,ic,jc = ',jj,ic,jc
        csymtmp(1:natm,jj)= csymbols(1:natm)
        csymtmp(ic,jj) = cj
        csymtmp(jc,jj) = ci
!.....Look for the history to check if we need to compute energy
        ihist = check_history(natm,nhist,csymtmp(1,jj),csymhist)
        if( ihist.eq.-1 ) then  ! no same symbols
!.....Send the order to run_pmd
          iorder = 1
          call mpi_bcast(iorder,1,mpi_integer,0,mpi_md_world,ierr)
          call run_pmd(hmat,natm,pos0,csymtmp(1,jj),epimc,epotmc &
               ,nstps_pmd,nx,ny,nz,mpi_md_world,nodes_md,myid_md&
               ,nstps_done)
          nstptmp(jj) = nstps_done
          ncalc = ncalc + 1
          nhist = nhist + 1
          if( nhist.gt.maxhist ) then
            stop ' Error: nhist.gt.maxhist'
          endif
          csymhist(1:natm,nhist) = csymtmp(1:natm,jj)
          efrm = epot2efrm(natm,ecpot,csymtmp(1,jj),epotmc)
          erghist(nhist) = efrm
        else  ! there is one same symbols
          efrm = erghist(ihist)
        endif
        ergtmp(jj) = efrm
!        csymtmp(1:natm,jj) = csymbols(1:natm)
        cjtmp(jj) = cj
!.....Since ci should be vacancy, cj is the migrating atom.
        js = cs2is(cj)
        de = demig(js) + (efrm-ergp)/2
        p = prefreq(js) *exp(-de/(temp*fkb))
!!$        print *,'jj,cj,ergp,erg,de,p = ',jj,cj,ergp,efrm,de,p
        probtmp(jj) = p
      enddo  ! loop over nearest neighbors, jj

!.....Compute total probability
      ptot = 0d0
      do jj=1,lsprmc(0,ic)
        ptot = ptot +probtmp(jj)
      enddo
!.....pick one event from the event list
      rand = urnd()*ptot
      ievent = lsprmc(0,ic)
      ptmp = 0d0
      do jj=1,lsprmc(0,ic)
        ptmp = ptmp +probtmp(jj)
        if( rand.lt.ptmp ) then
          ievent = jj
          exit
        endif
      enddo
      ergp = ergtmp(ievent)
      nstps_done = nstptmp(ievent)
!!$      do jj=1,lsprmc(0,ic)
!!$        print *,' jj,csymtmp = ',jj,csymtmp(1:natm,jj)
!!$      enddo
!!$      print *,' ievent = ',ievent
      csymbols(1:natm) = csymtmp(1:natm,ievent)
!.....Proceed real-time clock
      rand = urnd()
      dt = -log(rand)/ptot
      tclck = tclck +dt

!.....Output if needed
!.....Write symbols
      write(cfmt,'(i10)') natm
      write(iosym,'(i8,3x,'//trim(cfmt)//'a)') istp,csymbols(1:natm)
!.....Write POSCAR if migrating atom is not Al
      if( (noutint.le.0 .and. cjtmp(ievent).ne.'A') .or. &
           (noutint.gt.0 .and. mod(istp,noutint).eq.0) ) then
        write(cnum,'(i0)') istp
        call write_POSCAR('poscars/POSCAR_'//trim(cnum),natm,csymbols,pos0, &
             hmat,species)
      endif
!.....Write energy
      write(cergtxt,'(i8,2es15.7,a,i3,a,i2,a,a,a,i3)') &
           istp,tclck&
           ,ergp,', ncalc=',ncalc,', ievent=',ievent,', ',cjtmp(ievent) &
           ,', nstp=',nstps_done
      write(ioerg,'(a)') trim(cergtxt)
      write(6,'(a)') trim(cergtxt)
      flush(iosym)
      flush(ioerg)

    enddo  ! end of loop: istp
    iorder = -1
    call mpi_bcast(iorder,1,mpi_integer,0,mpi_md_world,ierr)

  else  ! myid_md.ne.0
    do while(.true.)
!.....Recieve order from the master
      call mpi_bcast(iorder,1,mpi_integer,0,mpi_md_world,ierr)
!.....Perform something depending on the order
!     iorder ==  1 : run_pmd
!     iorder == -1 : exit
      if( iorder.eq.1 ) then
        call run_pmd(hmat,natm,pos0,csymtmp,epimc,epot &
             ,nstps_pmd,nx,ny,nz,mpi_md_world,nodes_md,myid_md&
             ,nstps_done)
      else if( iorder.eq.-1 ) then
        exit
      endif
    end do
  endif

  if( myid_md.eq.0 ) then
    close(ioerg)
    close(iosym)
  endif

  deallocate(epimc,ecpot,csymprev,csymhist,erghist)

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
10 write(6,'(a)') ' reading '//trim(cfname)//' done'
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
    call read_i1(ionum,nstps_mc)
    return
  elseif( trim(cname).eq.'ncopy_x' ) then
    call read_i1(ionum,ncx)
    return
  elseif( trim(cname).eq.'num_out_interval' ) then
    call read_i1(ionum,noutint)
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
  elseif( trim(cname).eq.'num_Al_in_cluster' ) then
    call read_i1(ionum,num_Al_clst)
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
  elseif( trim(cname).eq.'atoms_to_move' ) then
    backspace(ionum)
    read(ionum,*) ctmp, lmove(0:3)
    return
  elseif( trim(cname).eq.'temperature' ) then
    call read_r1(ionum,temp)
    return
  elseif( trim(cname).eq.'random_seed' ) then
    call read_r1(ionum,dseed0)
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
subroutine write_init_params(lkinetic,nstps_mc, ncx,ncy,ncz,alat, &
     num_Mg,num_Si,num_Vac,num_Al_clst, &
     init_strct,nx,ny,nz,nstps_relax,lmove,temp)
  implicit none
  integer,intent(in):: nstps_mc,ncx,ncy,ncz,num_Mg,num_Si,num_Vac &
       ,init_strct,nx,ny,nz,nstps_relax,num_Al_clst
  real(8),intent(in):: alat,temp
  logical,intent(in):: lmove(0:3),lkinetic

  write(6,*) '=============== initial parameters =================='
  write(6,*) ''
  write(6,'(1x,a20,1x,i8)') 'num_steps',nstps_mc
  write(6,'(1x,a20,1x,i8)') 'num_relax_steps',nstps_relax
  write(6,*) ''
  write(6,'(1x,a20,1x,i2)') 'ncopy_x',ncx
  write(6,'(1x,a20,1x,i2)') 'ncopy_y',ncy
  write(6,'(1x,a20,1x,i2)') 'ncopy_z',ncz
  write(6,*) ''
  write(6,'(1x,a20,1x,i4)') 'num_Mg',num_Mg
  write(6,'(1x,a20,1x,i4)') 'num_Si',num_Si
  write(6,'(1x,a20,1x,i4)') 'num_Vac',num_Vac
  write(6,'(1x,a20,1x,i4)') 'num_Al_in_cluster',num_Al_clst
  write(6,*) ''
  write(6,'(1x,a20,1x,i2)') 'parallel_x',nx
  write(6,'(1x,a20,1x,i2)') 'parallel_y',ny
  write(6,'(1x,a20,1x,i2)') 'parallel_z',nz
  write(6,'(1x,a20,1x,i2)') 'initial_structure',init_strct
  write(6,*) ''
  write(6,'(1x,a20,1x,f8.3)') 'lattice_constant',alat
  write(6,'(1x,a20,1x,f8.3)') 'temperature',temp
  write(6,*) ''
  write(6,'(1x,a20,1x,l1)') 'kinetic',lkinetic
  write(6,'(1x,a20,4(1x,l1))') 'atoms_to_be_moved',lmove(0:3)
  write(6,*) ''
  write(6,*) '====================================================='
  
end subroutine write_init_params
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
  interface
    function urnd(dseed0)
      real(8),intent(in),optional:: dseed0
      real(8):: urnd
    end function urnd
  end interface

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
     ,num_Mg,num_Si,num_Vac,nnmaxmc,lsprmc,num_Al_clst)
  implicit none
  integer,intent(in):: natm,num_Mg,num_Si,num_Vac, &
       nnmaxmc,lsprmc(0:nnmaxmc,natm),num_Al_clst
  character,intent(inout):: csymbols(natm)
  real(8),intent(in):: pos0(3,natm)

  interface
    function urnd(dseed0)
      real(8),intent(in),optional:: dseed0
      real(8):: urnd
    end function urnd
  end interface
  integer:: i,jj,j,k,irnd,nsol,isol,icntr,nmg,nsi,nvac,inc,nal
  real(8):: cntr(3),dmin,d,r,rMg,rSi,rAl,rc
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
  nsol = num_Mg +num_Si +num_Vac +num_Al_clst
  nmg= 0
  nsi= 0
  nal= 0
  nvac= 0
  allocate(carr(nsol))
  inc = 0
  do while(.true.)
    if( inc.eq.nsol ) exit
    rMg = dble(num_Mg-nmg)/(nsol-inc)
    rSi = dble(num_Mg-nmg +num_Si-nsi)/(nsol-inc)
    rAl = dble(num_Mg-nmg +num_Si-nsi +num_Al_clst-nal)/(nsol-inc)
    r = urnd()
    if( r.lt.rMg ) then
      if( nmg.ge.num_Mg ) cycle
      inc= inc +1
      carr(inc) = 'M'
      nmg= nmg +1
    else if( r.le.rSi ) then
      if( nsi.ge.num_Si ) cycle
      inc= inc +1
      carr(inc) = 'S'
      nsi= nsi +1
    else if( r.le.rAl ) then
      if( nal.ge.num_Al_clst ) cycle
      inc= inc +1
      carr(inc) = 'X'
      nal= nal +1
    else
      if( nvac.ge.num_Vac ) cycle
      inc= inc +1
      carr(inc) = 'V'
      nvac= nvac +1
    endif
  enddo

  isol = 1
  csymbols(icntr) = carr(isol)

  do while(.true.)
    do i = 1,natm
      if( csymbols(i).eq.'A' ) cycle
      do jj = 1,lsprmc(0,i)
        j = lsprmc(jj,i)
        if( csymbols(j).eq.'A' ) then
          isol = isol + 1
          csymbols(j) = carr(isol)
        endif
        if( isol.eq.nsol ) then
!.....Replace X with A to restore Al atoms
          do k=1,natm
            if( csymbols(k).eq.'X' ) csymbols(k) = 'A'
          enddo
          return
        endif
      enddo
    enddo
  enddo
  stop 'Error: something wrong @clustered_symbols'
  
end subroutine clustered_symbols
!=======================================================================
subroutine paired_symbols(natm,pos0,csymbols,&
     num_Mg,num_Si,num_Vac,nnmax,lspr,cpair_type)
  implicit none
  integer,intent(in):: natm,num_Mg,num_Si,num_Vac &
       ,nnmax,lspr(0:nnmax,natm)
  character,intent(in):: cpair_type*2
  character,intent(inout):: csymbols(natm)
  real(8),intent(in):: pos0(3,natm)

  interface
    function urnd(dseed0)
      real(8),intent(in),optional:: dseed0
      real(8):: urnd
    end function urnd
  end interface

  integer:: npair,ipair,irnd,n,ichosen,jchosen,jj,i,j
  character:: cpairs(2,natm),c1,c2
  logical,external:: solute_in_neighbors
  integer:: isite_avail(natm),navail

  call create_pairs(natm,csymbols,num_Mg,num_Si,num_Vac, &
       npair,cpairs,cpair_type)

  write(6,'(a,i5)') ' npair = ',npair
  do i=1,npair
    write(6,'(a,i4,2a4)') ' i,cpairs(1:2,i) = ',i,cpairs(1:2,i)
  enddo

!.....Choose sites where Mg or Si or Vac will be put
!     which should be separated far enough not to connect each other
  ipair = 0
  do while(.true.)
    if( ipair.eq.npair ) exit
    isite_avail(1:natm) = 1
    do i=1,natm
      if( csymbols(i).ne.'A' ) then
        isite_avail(i) = 0
        do jj=1,lspr(0,i)
          j= lspr(jj,i)
          isite_avail(j)= 0
        enddo
      endif
    enddo
    navail= 0
    do i=1,natm
      if( isite_avail(i).eq.1 ) navail= navail +1
    enddo
    if( navail.eq.0 ) then
      write(6,'(a)') ' No more available sites for non-overlapping distribution.'
      stop
    endif

    irnd = int(navail*urnd()) +1
    n = 0
!.....Pick a random site that is available
    do i=1,natm
      if( isite_avail(i).eq.1 ) then
        n=n+1
        if( n.ge.irnd ) then
          ichosen = i
          exit
        endif
      endif
    enddo
!.....Check whether there is no Mg, Si, or Vac among neighbors
    if( solute_in_neighbors(ichosen,natm,csymbols,nnmax,lspr) ) cycle
!.....Put solute or vacancy on the chosen site
    ipair = ipair +1
    c1 = cpairs(1,ipair)
    c2 = cpairs(2,ipair)
    csymbols(ichosen)= c1
    if( c2.eq.'0' ) cycle
!.....Pick one one of neighbor sites and put c2 into it
    jj = lspr(0,ichosen)*urnd() +1
    jchosen = lspr(jj,ichosen)
    csymbols(jchosen)= c2
  enddo

end subroutine paired_symbols
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

!!$  character:: txt(natm)
  character(len=natm):: txt
  integer:: itmp,ios

!.....read the last symbols from file
  open(ionum,file=trim(fname),status='old')
  read(ionum,*) itmp, txt
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
  character(len=8):: date,cdate
  character(len=10):: time
  character(len=5):: zone
  character(len=1):: cs
  character(len=2):: sname(0:3) = (/ 'V ','Al','Mg','Si' /)
  
  call date_and_time(date, time, zone, date_time)

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
  write(cdate,'(i4.4,i2.2,i2.2)') date_time(1),date_time(2),date_time(3)
  write(90,'(3a)') 'written by PMC at ',trim(cdate)
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
     ,nstps_pmd,nx,ny,nz,mpi_md_world,nodes_md,myid_md,nstps_done)
  use pmc, only: symbol2sid
  implicit none
  integer,intent(in):: natm,nstps_pmd,nx,ny,nz&
       ,mpi_md_world,nodes_md,myid_md
  integer,intent(out):: nstps_done
  real(8),intent(in):: hmat(3,3),pos0(3,natm)
  real(8),intent(out):: epimc(natm),epotmc
  character,intent(in):: csymbols(natm)

  integer:: i,inc
  integer,parameter:: nismax = 9
  integer:: nstp,nerg,npmd,ifpmd,minstp,ntdst,n_conv,ifsort,iprint &
       ,ifdmp,ifcoulomb,numff,nrmtrans
  real(8):: hunit,h(3,3,0:1),am(nismax),dt,rc,dmp,tinit,tfin,ttgt(9)&
       ,trlx,stgt(3,3),ptgt,srlx,stbeta,strfin,fmv(3,0:9),ptnsr(3,3) &
       ,epot,ekin,eps_conv,rbuf,pini,pfin
  character:: ciofmt*6,cforce*20,ctctl*20,cpctl*20,czload_type*5,csi*1&
       ,boundary*3
  character(len=20):: cffs(1)
  logical:: ltdst,lstrs,lcellfix(3,3),lvc

  logical,save:: l1st = .true.
  integer,save:: ntot = 0
  real(8),save,allocatable:: tagtot(:),rtot(:,:),vtot(:,:),atot(:,:) &
       ,epitot(:),ekitot(:,:,:),stot(:,:,:),chgtot(:),chitot(:)

!.....at the 1st call, evaluate number of total atoms to be used in pmd
!     and allocate total system arrays
  if( myid_md.eq.0 ) then
    inc = 0
    do i = 1,natm
      csi = csymbols(i)
      if( csi.eq.'V' ) cycle
      inc = inc + 1
    enddo
    if( ntot.ne.inc ) then
      ntot = inc
      if( allocated(tagtot) ) then
        deallocate(tagtot,rtot,vtot,atot,epitot,ekitot,stot,chgtot,chitot)
      endif
      allocate(tagtot(ntot),rtot(3,ntot),vtot(3,ntot),atot(3,ntot) &
           ,epitot(ntot),ekitot(3,3,ntot),stot(3,3,ntot) &
           ,chgtot(ntot),chitot(ntot))
    endif
  else
    if( ntot.ne.1 ) then
      ntot = 1
      if( allocated(tagtot) ) then
        deallocate(tagtot,rtot,vtot,atot,epitot,ekitot,stot,chgtot,chitot)
      endif
      allocate(tagtot(ntot),rtot(3,ntot),vtot(3,ntot),atot(3,ntot) &
           ,epitot(ntot),ekitot(3,3,ntot),stot(3,3,ntot) &
           ,chgtot(ntot),chitot(ntot))
    endif
  endif

  hunit = 1d0
  h(1:3,1:3,0) = hmat(1:3,1:3)
  lcellfix(1:3,1:3) = .false.

  if( myid_md.eq.0 ) then
    inc = 0
    do i = 1,natm
      csi = csymbols(i)
      if( csi.eq.'V' ) cycle
      inc = inc + 1
      tagtot(inc) = dble(symbol2sid(csi)) +0.1d0 +1d-14*inc
      rtot(1:3,inc) = pos0(1:3,i)
      vtot(1:3,inc) = 0d0
      atot(1:3,inc) = 0d0
      stot(1:3,1:3,inc) = 0d0
      ekitot(1:3,1:3,inc) = 0d0
      epitot(inc) = 0d0
    enddo
  end if
  nerg = nstps_pmd
  npmd = 1
  am(1:9) = 1d0
  am(1) = 26.982  ! Al
  am(2) = 24.305  ! Mg
  am(3) = 28.085  ! Si
  dt = 5d0
  ciofmt = 'ascii'
  ifpmd = 1
  numff = 1
  cforce = 'NN'
  cffs(1) = 'NN'
  rc = 5.8d0
  rbuf = 0.2d0
  ifdmp = 2  ! FIRE
  dmp = 0.99d0
  minstp = 3
  tinit = 0d0
  tfin = 1d0
  ctctl = 'none'
  ttgt(1:9) = 300d0
  trlx = 100d0
  ltdst = .false.
  ntdst = 1
  nrmtrans = 1
  lstrs = .false.
  cpctl = 'none'
  stgt(1:3,1:3) = 0d0
  ptgt = 0d0
  pini = 0d0
  pfin = 0d0
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
  iprint = 0
  ifcoulomb = 0
  lvc = .false.
  boundary = 'ppp'

!.....call pmd_core to perfom MD
!!$  print *,'nstps_pmd = ',nstps_pmd
!!$  print *,'minstp = ',minstp
  call pmd_core(hunit,h,ntot,tagtot,rtot,vtot,atot,stot &
       ,ekitot,epitot,chgtot,chitot,nstps_pmd,nerg,npmd &
       ,myid_md,mpi_md_world,nodes_md,nx,ny,nz &
       ,nismax,am,dt,ciofmt,ifpmd,numff,cffs,rc,rbuf,ifdmp,dmp,minstp &
       ,tinit,tfin,ctctl,ttgt,trlx,ltdst,ntdst,nrmtrans,cpctl,stgt,ptgt &
       ,pini,pfin,srlx,stbeta,strfin,lstrs,lcellfix &
       ,fmv,ptnsr,epot,ekin,n_conv,ifcoulomb &
       ,czload_type,eps_conv,ifsort,iprint,nstps_done,lvc&
       ,boundary)
!!$  print *,'nstps_pmd,minstp,nstps_done = ' &
!!$       ,nstps_pmd,minstp,nstps_done
  if( myid_md.eq.0 ) then
    inc = 0
    epimc(1:natm) = 0d0
    do i = 1,natm
      csi = csymbols(i)
      if( csi.eq.'V' ) cycle
      inc = inc + 1
      epimc(i) = epitot(inc)
    enddo
    epotmc = epot
  endif
!!$  print *,'epot after pmd =',epot

  l1st = .false.

end subroutine run_pmd
!=======================================================================
subroutine calc_chem_pot(nspcs,species,ecpot,hmat,natm,pos0 &
     ,nstps_pmd,nx,ny,nz,mpi_md_world,nodes_md,myid_md)
!
! Calculate chemical potentials
!
  implicit none
  integer,intent(in):: nspcs,natm,nx,ny,nz,nstps_pmd&
       ,mpi_md_world,nodes_md,myid_md
  real(8),intent(in):: hmat(3,3),pos0(3,natm)
  character,intent(in):: species(0:nspcs)
  real(8),intent(out):: ecpot(0:nspcs)

  integer:: ispcs,i,nstps_done
  character:: cspcs*1
  character,allocatable:: csymtmp(:)
  real(8):: epot
  real(8),allocatable:: epi(:)

  allocate(csymtmp(natm),epi(natm))

  do ispcs = 0,nspcs
    do i=1,natm
      csymtmp(i) = 'A'
    enddo
    cspcs = species(ispcs)
    csymtmp(1) = cspcs
    call run_pmd(hmat,natm,pos0,csymtmp,epi,epot &
         ,nstps_pmd,nx,ny,nz,mpi_md_world,nodes_md,myid_md&
         ,nstps_done)
    ecpot(ispcs) = epot
  enddo
  
!.....chemical potential of Al
  ecpot(1) = ecpot(1) / natm

!.....subtract epot of Al bulk
  do ispcs = 0,nspcs
    if( ispcs.eq.1 ) cycle
    do i=1,natm
      csymtmp(i) = 'A'
    enddo
    cspcs = species(ispcs)
    csymtmp(1) = cspcs
    do i = 1,natm
      if( csymtmp(i).eq.'A' ) then
        ecpot(ispcs) = ecpot(ispcs) -ecpot(1)
      endif
    enddo
  enddo
  
end subroutine calc_chem_pot
!=======================================================================
function check_history(natm,nhist,csymbols,csymhist) result(ihist)
  implicit none
  integer,intent(in):: natm,nhist
  character,intent(in):: csymbols(natm),csymhist(natm,nhist)
  integer:: ihist

  integer:: i
  logical,external:: symbols_same

!!$  print *, 'check_history:'
!!$  print *, 'csymbols = ',csymbols(1:natm)
  ihist = -1
  do i=1,nhist
!!$    print *, 'csymhist(:,i) = ',csymhist(1:natm,i)
    if( symbols_same(natm,csymbols,csymhist(1,i)) ) then
      ihist=i
      return
    endif
  end do
!.....If there is no same structure in the history, return -1
  return

end function check_history
!=======================================================================
function symbols_same(ndim,csym1,csym2) result(lsame)
  implicit none
  integer,intent(in):: ndim
  character,intent(in):: csym1(ndim),csym2(ndim)

  logical:: lsame
  integer:: i

  do i=1,ndim
    if( csym1(i).ne.csym2(i) ) then
      lsame = .false.
      return
    endif
  enddo
  lsame = .true.
  return
  
end function symbols_same
!=======================================================================
function cs2is(cspcs) result(ispcs)
  use pmc
  implicit none
  character,intent(in):: cspcs*1

  integer:: ispcs
  integer:: i

  do i=0,3
    if( species(i).eq.cspcs ) then
      ispcs = i
    endif
  enddo
  return
end function cs2is
!=======================================================================
function epot2efrm(natm,ecpot,csym,epot)
  implicit none
  integer,intent(in):: natm
  real(8),intent(in):: ecpot(0:4),epot
  character,intent(in):: csym(natm)
  real(8):: epot2efrm

  integer:: i,isc
  character:: csi
  integer,external:: cs2is

  epot2efrm = epot
  do i=1,natm
    csi = csym(i)
    if( csi.eq.'V' ) cycle
    isc = cs2is(csi)
    epot2efrm = epot2efrm - ecpot(isc)
  enddo
  return
end function epot2efrm
!=======================================================================
function solute_in_neighbors(isite,natm,csymbols,nnmax,lspr)
  implicit none
  integer,intent(in):: isite,natm,nnmax,lspr(0:nnmax,natm)
  character,intent(in):: csymbols(natm)
  logical:: solute_in_neighbors

  integer:: j,jj

  solute_in_neighbors = .false.

  do jj=1,lspr(0,isite)
    j=lspr(jj,isite)
    if( csymbols(j).ne.'A' ) then
      solute_in_neighbors = .true.
      return
    endif
  enddo
  return
  
end function solute_in_neighbors
!=======================================================================
subroutine create_pairs(natm,csymbols,num_Mg,num_Si,num_Vac, &
     npair,cpairs,cpair_type)
!
!  Create pairs of X-X or X-Y if possible.
!  If it is not possible, put '0' in cpair(2,i).
!
  implicit none
  integer,intent(in):: natm,num_Mg,num_Si,num_Vac
  character,intent(in):: csymbols(natm),cpair_type*2
  integer,intent(out):: npair
  character,intent(out):: cpairs(2,natm)

  integer:: img,isi,ivac
  
  img = num_Mg
  isi = num_Si
  ivac= num_Vac

  npair = 0
  cpairs(1:2,1:natm) = '0'

  if( cpair_type.eq.'XX' ) then  ! pairs of same species
!.....Mg
    do while(.true.)
      if( img.le.0 ) exit
      npair = npair +1
      cpairs(1,npair) = 'M'
      img = img -1
      if( img.le.0 ) exit
      cpairs(2,npair) = 'M'
      img = img -1
    enddo
!.....Si
    do while(.true.)
      if( isi.le.0 ) exit
      npair = npair +1
      cpairs(1,npair) = 'S'
      isi = isi -1
      if( isi.le.0 ) exit
      cpairs(2,npair) = 'S'
      isi = isi -1
    enddo
!.....Vac
    do while(.true.)
      if( ivac.le.0 ) exit
      npair = npair +1
      cpairs(1,npair) = 'V'
      ivac = ivac -1
      if( ivac.le.0 ) exit
      cpairs(2,npair) = 'V'
      ivac = ivac -1
    enddo
  else if( cpair_type.eq.'XY' ) then  ! pairs of different species
!.....Mg
    do while(.true.)
      if( img.le.0 ) exit
      npair = npair +1
      cpairs(1,npair) = 'M'
      img = img -1
      if( isi.le.0 ) exit
      cpairs(2,npair) = 'S'
      isi = isi -1
    enddo
!.....Si
    do while(.true.)
      if( isi.le.0 ) exit
      npair = npair +1
      cpairs(1,npair) = 'S'
      isi = isi -1
      if( img.le.0 ) exit
      cpairs(2,npair) = 'M'
      img = img -1
    enddo
!.....Vac
    do while(.true.)
      if( ivac.le.0 ) exit
      npair = npair +1
      cpairs(1,npair) = 'V'
      ivac = ivac -1
      if( ivac.le.0 ) exit
      cpairs(2,npair) = 'V'
      ivac = ivac -1
    enddo
    
  endif
  
end subroutine create_pairs
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmc"
!     End:
