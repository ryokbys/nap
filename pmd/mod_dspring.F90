module dspring
!-----------------------------------------------------------------------
!                     Last-modified: <2023-01-18 14:14:36 KOBAYASHI Ryo>
!-----------------------------------------------------------------------
!  Spring force in descriptor space.
!  Originally for the purpose of restricting structure, at 2021-05-17, by R.K.
!-----------------------------------------------------------------------
  use pmdvars,only: nspmax, nsp
  use util,only: csp2isp
  implicit none
  include 'mpif.h'
  include "./const.h"
  save

  character(len=128):: paramsdir = '.'
  character(len=128),parameter:: cfname = 'in.dspring'
  integer,parameter:: ionum = 68

  logical:: ldspring = .false.  ! Flag to use dspring.
  logical:: initialized = .false.

  integer:: ndim_desc = -1
  integer:: ifmvdsp = -1
  real(8):: scnst = 1.0d0
  real(8),allocatable:: desctgt(:,:),dspfrc(:,:),dspstrs(:,:,:)
  real(8):: descstd(nspmax),epotdsp
  logical:: ldspc(nspmax)

contains
!=======================================================================
  subroutine init_dspring(myid,mpi_world,iprint)
!
!  Initialize dspring module.
!
    integer,intent(in):: myid,mpi_world,iprint
    
    integer:: inc,ix,iy,iz,mx,my,mz,ixyz,nxyz
    real(8):: anxi,anyi,anzi,fext
    integer:: istat(mpi_status_size),itag,ierr

    if( myid.eq.0 .and. iprint.ge.ipl_basic ) then
      print *,''
      print '(a)',' Desc-spring ON:'
    endif

    call read_dspring_params(myid,mpi_world,iprint)

!...Write out some information for users
    if( myid.eq.0 ) then
      if( iprint.ge.ipl_basic ) then
        print '(a,i6)','   Descriptor dimension = ',ndim_desc
        print '(a,i6)','   Group of ifmv        = ',ifmvdsp
        print '(a,es11.3)','   Spring constant      = ',scnst
      endif
    endif

    initialized = .true.
    return
  end subroutine init_dspring
!=======================================================================
  subroutine read_dspring_params(myid,mpi_world,iprint)
!
!  Read some info from in.dspring:
!    - dimension of the descriptor that must be consistent with in.params.desc
!    - spring constant
!    - target descriptor values of the reference structure
!    - normalization factors to divide the descriptors before
!      they are used in the spring force
!    - to which group of atoms the forces exert on (specified by ifmv)
!
!  The format of in.dspring is like the following:
!-----------------------------------------------------------------------
!  !  Comment line if begins with ! or #
!  !
!  group_ifmv     2
!  spring_constant    1.0   (common for all the descriptors)
!  desc_std   0.0405  0.0222  0.0171
!  desc_dimension     55
!  !  isf  tgt(1)    tgt(2)   tgt(3)
!     1    0.1245    0.1234   0.2234
!     2    0.1245    0.1234   0.2234
!     3    0.1245    0.1234   0.2234
!     ....
!     55   0.1245    0.1234   0.2234
!-----------------------------------------------------------------------
!
    use util, only: num_data
    integer,intent(in):: myid,mpi_world,iprint
    
    integer:: nentry,isp,ierr,isf,idesc
    real(8):: tmp
    character:: cline*128, c1st*128, csp*3, cmode*128
    real(8),allocatable:: desctmp(:)

    if( myid.eq.0 ) then ! only at master, node-0.
      ldspc(:) = .true.
      open(ionum,file=trim(cfname),status='old')
      cmode = 'none'
      do while(.true.)
        read(ionum,'(a)',end=10) cline
        nentry = num_data(cline,' ')
        if( nentry.eq.0 ) cycle  ! blank line
        if( cline(1:1).eq.'!' .or. cline(1:1).eq.'#' ) cycle  ! comment line
!!$        print *,trim(cline)
        if( trim(cmode).eq.'none' ) then
          read(cline,*) c1st
          if( trim(c1st).eq.'spring_constant' ) then
            backspace(ionum)
            read(ionum,*) c1st, scnst
          else if( trim(c1st).eq.'group_ifmv' ) then
            backspace(ionum)
            read(ionum,*) c1st, ifmvdsp
          else if( trim(c1st).eq.'desc_std' ) then
            backspace(ionum)
            read(ionum,*) c1st, (descstd(isp),isp=1,nsp)
          else if( trim(c1st).eq.'desc_dimension' ) then
            backspace(ionum)
            read(ionum,*) c1st, ndim_desc
            allocate(desctmp(ndim_desc),desctgt(ndim_desc,nspmax))
            cmode = 'read_descs'
          else
            print *,'There is no such dspring keyword: ', trim(c1st)
          endif
        else if( trim(cmode).eq.'read_descs' ) then
          backspace(ionum)
          read(ionum,*) idesc, (desctgt(idesc,isp),isp=1,nsp)
          if( idesc.gt.ndim_desc ) then
            print *,'ERROR: desc-ID exceeds the descriptor dimension.'
            stop 1
          endif
        endif
      enddo
!!$      read(ionum,*) tmp
10    close(ionum)
    endif

!.....Normalize desctgt by dividing by descnrm
    do isp=1,3
      do isf=1,ndim_desc
        desctgt(isf,isp) = desctgt(isf,isp) / descstd(isp)
      enddo
    enddo

!.....Broadcast some parameters
    call mpi_barrier(mpi_world,ierr)
    call mpi_bcast(ndim_desc,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(scnst,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(ifmvdsp,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(ldspc,nspmax,mpi_logical,0,mpi_world,ierr)
    if( myid.ne.0 ) allocate(desctgt(ndim_desc,nspmax))
    call mpi_bcast(desctgt,ndim_desc*nspmax,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(descstd,nspmax,mpi_real8,0,mpi_world,ierr)

    if( allocated(desctmp) ) deallocate(desctmp)
    return
  end subroutine read_dspring_params
!=======================================================================
  subroutine force_dspring(namax,natm,nnmax,lspr,rcin,h,hi,tag,ra, &
       aa,edsp,nb,nbmax,lsb,nex,lsrc, &
       myparity,nn,myid,mpi_world,iprint,l1st)
!
!  Write the code of spring force in the descriptor space,
!  and add them to the variable, aa.
!  Note that aa should be normalized by h-matrix.
!
    use descriptor,only: gsfi,dgsfi,calc_desci,pre_desci,make_gsf_arrays,nsf
    use util,only: ifmvOf
    integer,intent(in):: namax,natm,nnmax,lspr(0:nnmax,namax), &
         myid,mpi_world,iprint,nn(6),nex(3), &
         nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3)
    real(8),intent(in):: rcin,tag(namax),h(3,3),hi(3,3),ra(3,namax)
    real(8),intent(inout):: aa(3,namax),edsp(namax)
    logical,intent(in):: l1st

    integer:: ia,ifmv,isf,isp,jsp,jj,ja,i,k,ixyz,jxyz,ierr
    real(8):: tmp,at(3),ave_aa,ave_dsp,dsq, &
         xi(3),xj(3),xij(3),rij(3),dij,sij

    if( .not.allocated(dspfrc) ) then
      allocate(dspfrc(3,namax),dspstrs(3,3,namax))
    else if( size(dspfrc).lt.3*namax ) then
      deallocate(dspfrc,dspstrs)
      allocate(dspfrc(3,namax),dspstrs(3,3,namax))
    endif
    
!.....Compute descriptor values of atoms
    call pre_desci(namax,natm,nnmax,lspr,iprint,rcin)
    call make_gsf_arrays(l1st,namax,natm,tag,nnmax,lspr,myid,mpi_world,iprint)

    epotdsp = 0d0
    edsp(:) = 0d0
    dspfrc(:,:) = 0d0
    dspstrs(:,:,:) = 0d0
    do ia=1,natm
      ifmv = ifmvOf(tag(ia))
      if( ifmv.ne.ifmvdsp ) cycle  ! only specified atoms pass here
      isp = int(tag(ia))
      if( .not. ldspc(isp) ) cycle  ! only specified species pass here
      call calc_desci(ia,namax,natm,nnmax,h,tag,ra,lspr,rcin,iprint)
      dsq = 0d0
      xi(1:3) = ra(1:3,ia)
      do isf=1,nsf
!.....Normalize gsfi with stddev of desc vector of species.
        gsfi(isf) = gsfi(isf) / descstd(isp)
        dsq = dsq +(gsfi(isf)-desctgt(isf,isp))**2
!.....Since dgsfi is not normalized, normalize the prefactor here instead.
        tmp = scnst*(gsfi(isf)-desctgt(isf,isp)) /descstd(isp)
!.....Compute spring forces
        do jj=1,lspr(0,ia)
          ja = lspr(jj,ia)
!.....dspring force
          dspfrc(1:3,ja) = dspfrc(1:3,ja) -dgsfi(1:3,isf,jj)*tmp
!.....dspring stress
          xj(1:3) = ra(1:3,ja)
          xij(1:3) = xj(1:3) -xi(1:3)
          rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
          dij = sqrt( rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3) )
          do ixyz=1,3
            do jxyz=1,3
              sij = -dgsfi(jxyz,isf,jj)*tmp*rij(ixyz)
              dspstrs(ixyz,jxyz,ja) = dspstrs(ixyz,jxyz,ja) +sij
              dspstrs(ixyz,jxyz,ia) = dspstrs(ixyz,jxyz,ia) +sij
            enddo
          enddo
        enddo
        dspfrc(1:3,ia) = dspfrc(1:3,ia) -dgsfi(1:3,isf,0)*tmp
      enddo
      edsp(ia) = edsp(ia) +0.5d0 *scnst *dsq
      epotdsp = epotdsp +0.5d0 *scnst *dsq
    enddo

!.....Send back forces on immigrants to the neighboring nodes
    call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,dspfrc,3)
!.....Normalize forces by h-matrix
    do i=1,natm
      at(1:3)= dspfrc(1:3,i)
      dspfrc(1:3,i)= hi(1:3,1)*at(1) +hi(1:3,2)*at(2) +hi(1:3,3)*at(3)
    enddo
!.....Add dspring forces to the original forces    
    aa(1:3,1:natm) = aa(1:3,1:natm) +dspfrc(1:3,1:natm)
    return
  end subroutine force_dspring
!=======================================================================
  subroutine add_dspring_strs(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
       ,nn,mpi_world,strs)
!
!  Add dspring stress per atom to the total stress per atom.
!  It requires that force_dspring is called beforehand and stress are computed already.
!
    integer,intent(in):: namax,natm,nbmax,nb,myparity(3),nn(6),nex(3), &
         lsb(0:nbmax,6),lsrc(6),mpi_world
    real(8),intent(inout):: strs(3,3,namax)
    
    call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,dspstrs,9)
    strs(1:3,1:3,1:natm) = strs(1:3,1:3,1:natm) +dspstrs(1:3,1:3,1:natm)*0.5d0
    return
  end subroutine add_dspring_strs
!=======================================================================
  subroutine add_dspring_epot(epot,mpi_world)
!
!  Add dspring epot to that of the system.
!  It requires that force_dspring is called beforehand.
!
    real(8),intent(inout):: epot
    integer,intent(in):: mpi_world

    integer:: ierr
    real(8):: tmp
    
!.....Gather energy
    tmp = epotdsp
    call mpi_allreduce(tmp,epotdsp,1,mpi_real8,mpi_sum,mpi_world,ierr)
    epot = epot +epotdsp
    return
  end subroutine add_dspring_epot
!=======================================================================
  subroutine final_dspring(myid)
!
!  Finalize dspring if needed.
!
    integer,intent(in):: myid

    deallocate(desctgt)
  end subroutine final_dspring
!=======================================================================
end module dspring
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd lib"
!     End:
