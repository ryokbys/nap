module group
!-----------------------------------------------------------------------
!                     Last modified: <2024-07-12 19:45:30 KOBAYASHI Ryo>
!-----------------------------------------------------------------------
!  Module for grouping atoms.
!-----------------------------------------------------------------------
  use pmdvars, only: max_group, nspmax
  use util, only: ispOf, igvarOf, replace_igvar
  implicit none
  private
  include 'const.h'
  save

  public:: grouping, bcast_group, init_group, register_group

!.....Group types:
!       1) species
!       2) sphere
  integer:: gtype(max_group)
!.....Grouping timing:
!       0) never,
!       1) every step,
!       n) per each n steps,
!       -1) only at the 1st time
  integer:: gtiming(max_group) = 0

!.....species grouping
  integer:: ispc(max_group,nspmax)
!.....sphere grouping
  real(8):: sph_org(max_group,3) ! origin vector in h-mat unit
  real(8):: sph_rad2(max_group) ! radius**2
  integer:: isph(max_group,2)  ! 1) inside, 2) outside

contains
!=======================================================================
  subroutine init_group()

!.....Initialize variables
    gtype(:) = 0
    ispc(:,:) = 0
    isph(:,:) = 0
    sph_org(:,:) = 0d0
    sph_rad2(:) = 0d0

  end subroutine init_group
!=======================================================================
  function gtypeOf(gname)
    character(len=*),intent(in):: gname
    integer:: gtypeOf

    if( trim(gname).eq.'species' ) gtypeOf = 1
    if( trim(gname).eq.'sphere' )  gtypeOf = 2

  end function gtypeOf
!=======================================================================
  subroutine register_group(cline)
!
!  Register a group by a given line in in.pmd file.
!  This routine is supposed to be called only at the MPI node-0.
!
    use util,only: num_data
    character(len=*),intent(in):: cline

    character(len=20):: c1st,cgname
    integer:: ndat,i,igrp,kin,kout,npair
    integer,allocatable:: itmp(:,:)
    real(8):: rad,orgx,orgy,orgz

    ndat = num_data(cline, ' ')
    read(cline,*) c1st
    if( trim(c1st).ne.'group' ) then
      print *,'WARNING: '//trim(c1st)//' is not group.'
      return
    endif
    if( ndat.lt.2 ) then
      print *,'WARNING: wrong number of input entry, ndat=',ndat
      return
    endif

    read(cline,*) c1st, igrp, cgname
    if( igrp.gt.4 ) then
      print *,'ERROR: group-ID for '//trim(cgname)//' exceeds the limit, 4.'
      print *,'       group-ID= ',igrp
      return
    endif
    if( trim(cgname).eq.'species' ) then
! e.g.) group  1  species  1 0  2 1  3 2
!              ^                         : group-ID (<=4)
!                 ^^^^^^^                : grouping name
!                          ^    ^    ^   : species
!                            ^    ^    ^ : group variable (0-9)
      if( mod(ndat-3, 2).ne.0 ) then
        print '(a,i0,a)',' ERROR: group ', igrp,  &
             ' '//trim(cgname)//'; number of entry items is wrong.'
        return
      endif
      npair = int((ndat-3)/2)
      allocate(itmp(npair,2))
      itmp(:,:) = 0
      read(cline,*) c1st, igrp, cgname, (itmp(i,1),itmp(i,2),i=1,npair)
      gtype(igrp) = gtypeOf(cgname)
      do i=1,npair
        ispc(igrp,itmp(i,1)) = itmp(i,2)
      enddo
      gtiming(igrp) = -1  ! only at the 1st time
      deallocate(itmp)
    else if( trim(cgname).eq.'sphere' ) then
! e.g.) group  1  sphere  20.0  0.5  0.5  0.5  1  2
!              ^                                    : group-ID (<=4)
!                 ^^^^^^                            : grouping name
!                         ^^^^                      : radius
!                               ^^^  ^^^  ^^^       : origian
!                                              ^  ^ : group variables for inside & outside
      if( ndat.lt.9 ) then
        print '(a,i0,a)',' ERROR: group ', igrp,  &
             ' '//trim(cgname)//'; number of entry items is wrong.'
        return
      endif
      read(cline,*) c1st, igrp, cgname, rad, orgx, orgy, orgz, kin, kout
      gtype(igrp) = gtypeOf(cgname)
      sph_org(igrp,1:3) = (/ orgx, orgy, orgz /)
      sph_rad2(igrp) = rad**2
      isph(igrp,1:2) = (/ kin, kout /)
      gtiming(igrp) = 1
    else
      print *,'WARNING: no such group name: '//trim(cgname)
      return
    endif
    
    return
  end subroutine register_group
!=======================================================================
  subroutine bcast_group(mpicomm)
    include 'mpif.h'
    integer,intent(in):: mpicomm
    integer:: ierr

    call mpi_bcast(gtype,max_group,mpi_integer,0,mpicomm,ierr)
    call mpi_bcast(gtiming,max_group,mpi_integer,0,mpicomm,ierr)
    call mpi_bcast(ispc,max_group*nspmax,mpi_integer,0,mpicomm,ierr)
    call mpi_bcast(isph,max_group*2,mpi_integer,0,mpicomm,ierr)
    call mpi_bcast(sph_org,max_group*3,mpi_real8,0,mpicomm,ierr)
    call mpi_bcast(sph_rad2,max_group*3,mpi_real8,0,mpicomm,ierr)
    
  end subroutine bcast_group
!=======================================================================
  subroutine grouping(namax,natm,h,tag,ra,sorg,istp,myid,mpicomm,iprint)
!
!  Grouping according to the given grouping method name.
!
    integer,intent(in):: namax,natm,myid,mpicomm,iprint,istp
    real(8),intent(in):: h(3,3),ra(3,namax),sorg(3)
    real(8),intent(inout):: tag(namax)
    integer:: gid,isp,ia,gti,igvar
    real(8):: ti,ri(3),dri(3),rad2

    do gid=1,4
      if( gtiming(gid).eq.0 ) cycle
      if( gtiming(gid).lt.0 .and. abs(gtiming(gid)).ne.istp ) cycle
      if( gtiming(gid).gt.1 .and. mod(istp,gtiming(gid)).ne.1 ) cycle
      gti = gtype(gid)
      if( gti.eq.1 ) then  ! species grouping
        do ia=1,natm
          ti = tag(ia)
          isp = ispOf(ti)
          igvar = ispc(gid,isp)
          call replace_igvar(ti,gid,igvar)
!.....Store the modified tag
          tag(ia) = ti
        enddo  ! ia loop
      else if( gti.eq.2 ) then  ! sphere grouping
        do ia=1,natm
          ti = tag(ia)
          ri(1:3) = ra(1:3,ia) +sorg(1:3) -sph_org(gid,1:3)
          dri(1:3) = h(1:3,1)*ri(1) +h(1:3,2)*ri(2) +h(1:3,3)*ri(3)
          rad2 = dri(1)**2 +dri(2)**2 +dri(3)**2
          if( rad2.le.sph_rad2(gid) ) then
            igvar = isph(gid,1)
          else
            igvar = isph(gid,2)
          endif
          call replace_igvar(ti,gid,igvar)
          tag(ia) = ti
        enddo  ! ia loop
      endif
    enddo  ! gid loop
    return
  end subroutine grouping
end module group
!-----------------------------------------------------------------------
!  Local Variables:
!  compile-command: "make pmd lib"
!  End:
