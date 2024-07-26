module fdesc
!-----------------------------------------------------------------------
!                     Last-modified: <2024-07-26 10:36:37 KOBAYASHI Ryo>
!-----------------------------------------------------------------------
!  Potential in descriptor space.
!  Originally for the purpose of restricting structure, at 2021-05-17, by R.K.
!-----------------------------------------------------------------------
  use pmdvars,only: nspmax, nsp
  use util,only: csp2isp
  implicit none
  include 'mpif.h'
  include "./const.h"
  save

  character(len=128):: paramsdir = '.'
  character(len=128),parameter:: cfname = 'in.fdesc'
  integer,parameter:: ionum = 68

  logical:: lfdesc = .false.  ! Flag to use fdesc.
  logical:: initialized = .false.
  logical:: pc1_as_edesc = .false.  ! pc1 value at edesci

  integer:: ndim_desc = -1
  integer:: giddesc = -1
  real(8):: scnst = 1.0d0
  real(8):: gcoef = 0.1d0
  real(8):: gsgm  = 1.0d0
  real(8),allocatable:: desctgt(:),descpca(:), descfrc(:,:),descstrs(:,:,:)
  real(8):: edesc
!  logical:: ldspc(nspmax)

contains
!=======================================================================
  subroutine init_fdesc(myid,mpi_world,iprint)
!
!  Initialize fdesc module.
!
    integer,intent(in):: myid,mpi_world,iprint
    
    integer:: inc,ix,iy,iz,mx,my,mz,ixyz,nxyz
    real(8):: anxi,anyi,anzi,fext
    integer:: istat(mpi_status_size),itag,ierr

    if( myid.eq.0 .and. iprint.ge.ipl_basic ) then
      print *,''
      print '(a)',' Desc-spring ON:'
    endif

    call read_fdesc_params(myid,mpi_world,iprint)

!...Write out some information for users
    if( myid.eq.0 ) then
      if( iprint.ge.ipl_basic ) then
        print '(a,i6)','   Descriptor dimension = ',ndim_desc
        print '(a,i6)','   Group ID for fdesc   = ',giddesc
!        print '(a,es11.3)','   Spring constant      = ',scnst
        print '(a,f7.3)','   Gaussian coeff.      = ',gcoef
        print '(a,f7.3)','   Gaussian sigma       = ',gsgm
      endif
    endif

    initialized = .true.
    return
  end subroutine init_fdesc
!=======================================================================
  subroutine read_fdesc_params(myid,mpi_world,iprint)
!
!  Read some info from in.fdesc:
!    - dimension of the descriptor that must be consistent with in.params.desc
!    - spring constant
!    - target descriptor values of the reference structure
!    - PCA 1st coefficients of descriptor components
!    - to which group of atoms the forces exert on (specified by gid)
!
!  The format of in.fdesc is like the following:
!-----------------------------------------------------------------------
!  !  Comment line if begins with ! or #
!  !
!  group_id  1
!  spring_constant    1.0   (common for all the descriptors)
!  ! gauss_coef         1.0   (common for all the descriptors)
!  ! gauss_sigma        1.0   (common for all the descriptors)
!  desc_dimension     55
!  !  isf  desctgt  descpca
!     1    0.1245   0.345  
!     2    0.1452   0.456  
!     3    0.1542   0.567  
!     ....
!     55   0.1245   1.234  
!-----------------------------------------------------------------------
!
    use util, only: num_data
    integer,intent(in):: myid,mpi_world,iprint
    
    integer:: nentry,isp,ierr,isf,idesc
    real(8):: tmp
    character:: cline*128, c1st*128, csp*3, cmode*128
    real(8),allocatable:: desctmp(:)

    if( myid.eq.0 ) then ! only at master, node-0.
!      ldspc(:) = .true.
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
          else if( trim(c1st).eq.'group_id' ) then
            backspace(ionum)
            read(ionum,*) c1st, giddesc
          else if( trim(c1st).eq.'gauss_coef' ) then
            backspace(ionum)
            read(ionum,*) c1st, gcoef
          else if( trim(c1st).eq.'gauss_sigma' ) then
            backspace(ionum)
            read(ionum,*) c1st, gsgm
          else if( trim(c1st).eq.'pc1_as_edesc' ) then
            backspace(ionum)
            read(ionum,*) c1st, pc1_as_edesc
          else if( trim(c1st).eq.'desc_dimension' ) then
            backspace(ionum)
            read(ionum,*) c1st, ndim_desc
            if( allocated(desctgt) ) deallocate(desctgt,descpca)
            allocate(desctmp(ndim_desc),desctgt(ndim_desc), &
                 descpca(ndim_desc))
            cmode = 'read_descs'
          else
            print *,'There is no such fdesc keyword: ', trim(c1st)
          endif
        else if( trim(cmode).eq.'read_descs' ) then
          backspace(ionum)
          read(ionum,*) idesc, desctgt(idesc), descpca(idesc)
          if( idesc.gt.ndim_desc ) then
            print *,'ERROR: desc-ID exceeds the descriptor dimension.'
            stop 1
          endif
        endif
      enddo
!!$      read(ionum,*) tmp
10    close(ionum)
    endif

!.....Broadcast some parameters
    call mpi_barrier(mpi_world,ierr)
    call mpi_bcast(ndim_desc,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(scnst,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(giddesc,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(pc1_as_edesc,1,mpi_logical,0,mpi_world,ierr)
!    call mpi_bcast(ldspc,nspmax,mpi_logical,0,mpi_world,ierr)
    if( myid.ne.0) then
      if( allocated(desctgt) ) deallocate(desctgt, descpca)
      allocate(desctgt(ndim_desc), descpca(ndim_desc))
    endif
    call mpi_bcast(desctgt,ndim_desc,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(descpca,ndim_desc,mpi_real8,0,mpi_world,ierr)

    if( allocated(desctmp) ) deallocate(desctmp)
    return
  end subroutine read_fdesc_params
!=======================================================================
  subroutine force_desc(namax,natm,nnmax,lspr,rcin,h,hi,tag,ra, &
       aa,epot,edesci,strs,nb,nbmax,lsb,nex,lsrc, &
       myparity,nn,myid,mpi_world,iprint,l1st)
!
!  Write the code of spring force in the descriptor space,
!  and add them to the variable, aa.
!  Note that aa should be normalized by h-matrix.
!
    use descriptor,only: gsfi,dgsfi,calc_desci,pre_desci,make_gsf_arrays,nsf
    use util,only: igvarOf
    integer,intent(in):: namax,natm,nnmax,lspr(0:nnmax,namax), &
         myid,mpi_world,iprint,nn(6),nex(3), &
         nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3)
    real(8),intent(in):: rcin,tag(namax),h(3,3),hi(3,3),ra(3,namax)
    real(8),intent(inout):: aa(3,namax),edesci(namax),epot,strs(3,3,namax)
    logical,intent(in):: l1st

    integer:: ia,igv,isf,isp,jsp,jj,ja,i,k,ixyz,jxyz,ierr
    real(8):: tmp,at(3),ave_aa,ave_dsp,dsq,dpca1,pca1,aexp, &
         xi(3),xj(3),xij(3),rij(3),dij,sij

    if( .not.allocated(descfrc) ) then
      allocate(descfrc(3,namax),descstrs(3,3,namax))
    else if( size(descfrc).lt.3*namax ) then
      deallocate(descfrc,descstrs)
      allocate(descfrc(3,namax),descstrs(3,3,namax))
    endif
    
!.....Compute descriptor values of atoms
    call pre_desci(namax,natm,nnmax,lspr,iprint,rcin)
    call make_gsf_arrays(l1st,namax,natm,tag,nnmax,lspr,myid,mpi_world,iprint)

    edesc = 0d0
    edesci(:) = 0d0
    descfrc(:,:) = 0d0
    descstrs(:,:,:) = 0d0
    do ia=1,natm
      igv = igvarOf(tag(ia),giddesc)
      if( igv.eq.0 ) cycle  ! fdesc works only on atoms of igv > 0 
      isp = int(tag(ia))
!      if( .not. ldspc(isp) ) cycle  ! only specified species pass here
      call calc_desci(ia,namax,natm,nnmax,h,tag,ra,lspr,rcin,iprint)
      xi(1:3) = ra(1:3,ia)
      dpca1 = 0d0
      pca1 = 0d0
      do isf=1,nsf
        dpca1 = dpca1 + (gsfi(isf)-desctgt(isf)) * descpca(isf)
        pca1 = pca1 + gsfi(isf) * descpca(isf)
      enddo
      dsq = dpca1**2
      aexp = -gcoef * exp(-dsq /2 /gsgm**2)
      do isf=1,nsf
!!$        tmp = scnst*dpca1 * descpca(isf)
        tmp = -aexp *dpca1 /gsgm**2 *descpca(isf)
!.....Compute spring forces
        do jj=1,lspr(0,ia)
          ja = lspr(jj,ia)
!.....desc force
          descfrc(1:3,ja) = descfrc(1:3,ja) -dgsfi(1:3,isf,jj)*tmp
!.....desc stress
          xj(1:3) = ra(1:3,ja)
          xij(1:3) = xj(1:3) -xi(1:3)
          rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
          dij = sqrt( rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3) )
          do ixyz=1,3
            do jxyz=1,3
              sij = -dgsfi(jxyz,isf,jj)*tmp*rij(ixyz)
              descstrs(ixyz,jxyz,ja) = descstrs(ixyz,jxyz,ja) +sij
              descstrs(ixyz,jxyz,ia) = descstrs(ixyz,jxyz,ia) +sij
            enddo
          enddo ! ixyz
        enddo ! jj
        descfrc(1:3,ia) = descfrc(1:3,ia) -dgsfi(1:3,isf,0)*tmp
      enddo ! isf
      if( pc1_as_edesc ) then
        edesci(ia) = edesci(ia) + dpca1
      else
        edesci(ia) = edesci(ia) + aexp
      endif
!!$      edesc = edesc +0.5d0 *scnst *dsq
      edesc = edesc +aexp
    enddo ! ia

!.....Send back forces on immigrants to the neighboring nodes
    call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,descfrc,3)
!!$!.....Normalize forces by h-matrix
!!$    do i=1,natm
!!$      at(1:3)= descfrc(1:3,i)
!!$      descfrc(1:3,i)= hi(1:3,1)*at(1) +hi(1:3,2)*at(2) +hi(1:3,3)*at(3)
!!$    enddo
    if( iprint.ge.ipl_debug .and. myid.eq.0 ) then
      print '(a,es12.3)','edesc = ',edesc
      do i=1,10
        print '(i6, 6(2x, es12.3))', i, descfrc(1:3,i), aa(1:3,i)
      enddo
    endif
!.....Add desc forces to the original forces
    aa(1:3,1:natm) = aa(1:3,1:natm) +descfrc(1:3,1:natm)
!.....Gather epot
    tmp = edesc
    call mpi_allreduce(tmp,edesc,1,mpi_real8,mpi_sum,mpi_world,ierr)
    epot = epot + edesc
!.....Send back stresses
    call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,descstrs,9)
    strs(1:3,1:3,1:natm) = strs(1:3,1:3,1:natm) +descstrs(1:3,1:3,1:natm)*0.5d0

    return
  end subroutine force_desc
!=======================================================================
  subroutine add_fdesc_strs(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
       ,nn,mpi_world,strs)
!
!  Add desc stress per atom to the total stress per atom.
!  It requires that force_desc is called beforehand and stress are computed already.
!
    integer,intent(in):: namax,natm,nbmax,nb,myparity(3),nn(6),nex(3), &
         lsb(0:nbmax,6),lsrc(6),mpi_world
    real(8),intent(inout):: strs(3,3,namax)
    
    call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,descstrs,9)
    strs(1:3,1:3,1:natm) = strs(1:3,1:3,1:natm) +descstrs(1:3,1:3,1:natm)*0.5d0
    return
  end subroutine add_fdesc_strs
!=======================================================================
  subroutine add_fdesc_epot(epot,mpi_world)
!
!  Add desc epot to that of the system.
!  It requires that force_desc is called beforehand.
!
    real(8),intent(inout):: epot
    integer,intent(in):: mpi_world

    integer:: ierr
    real(8):: tmp
    
!.....Gather energy
    tmp = edesc
    call mpi_allreduce(tmp,edesc,1,mpi_real8,mpi_sum,mpi_world,ierr)
    epot = epot +edesc
    return
  end subroutine add_fdesc_epot
!=======================================================================
  subroutine final_fdesc(myid)
!
!  Finalize fdesc if needed.
!
    integer,intent(in):: myid

    deallocate(desctgt,descpca)
  end subroutine final_fdesc
!=======================================================================
end module fdesc
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd lib"
!     End:
