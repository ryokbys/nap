module zload
  implicit none
  save
!.....initial z top and bottom positions
  real(8):: ztop0, zbot0, zlen0, dzl
  integer:: nztop, nzbot
!.....variables that change during z-loading simulation
  real(8):: zlen, ztop, zbot
  
contains
  subroutine set_zload_atoms(natm,ra,tag,h,fmv,sorg,strfin,nstp &
       ,myid_md,mpi_md_world)
!
!     Choose atoms to be applied z-loading.
!     The criterion is 5% of top and bottom atoms in z-direction.
!     ifmv's of these atoms are set to 9.
!
    implicit none
    include 'mpif.h'
    integer,intent(in):: natm,myid_md,mpi_md_world,nstp
    real(8),intent(in):: ra(3,natm),sorg(3),strfin,h(3,3)
    real(8),intent(inout):: tag(natm),fmv(3,0:9)
    integer:: ierr,i,nztopl,nzbotl
    real(8):: ztopl,zbotl,dzlfin

!.....Set 9th fmv to fix y and z directions, not x
    fmv(1:3,9) = (/ 1d0, 0d0, 0d0 /)
    
!.....Detect initial top and bottom positions
    ztopl = 0d0
    zbotl = 1d0
    do i=1,natm
      zbotl= min(zbotl,ra(3,i)+sorg(3))
      ztopl= max(ztopl,ra(3,i)+sorg(3))
    enddo
    ztop0= 0d0
    zbot0= 0d0
    call mpi_allreduce(ztopl,ztop0,1,mpi_double_precision,mpi_max &
         ,mpi_md_world,ierr)
    call mpi_allreduce(zbotl,zbot0,1,mpi_double_precision,mpi_min &
         ,mpi_md_world,ierr)
    zlen0 = ztop0 -zbot0
    dzlfin = zlen0 *strfin
    dzl= dzlfin /nstp /2
    zlen = zlen0
    ztop = ztop0
    zbot = zbot0
    if( myid_md.eq.0 ) then
      write(6,'(a,2es12.4)') ' zlen0,dzl = ',zlen0,dzl
    endif

!.....Define top and bottom atoms z-motions of which are controlled
    nztopl = 0
    nzbotl = 0
    do i=1,natm
      if( ra(3,i)+sorg(3) .gt. ztop-1.5d0/h(3,3) ) then
        call replaceTag('ifmv',9,tag(i))
        nztopl = nztopl +1
      else if( ra(3,i)+sorg(3) .lt. zbot+1.5d0/h(3,3) ) then
        call replaceTag('ifmv',9,tag(i))
        nzbotl = nzbotl +1
      endif
    enddo

    call mpi_allreduce(nztopl,nztop,1,mpi_integer,mpi_max &
         ,mpi_md_world,ierr)
    call mpi_allreduce(nzbotl,nzbot,1,mpi_integer,mpi_max &
         ,mpi_md_world,ierr)
    if( myid_md.eq.0 ) then
      write(6,'(a,2i8)') ' nztop, nzbot = ',nztop,nzbot
    endif

  end subroutine set_zload_atoms
!=======================================================================
  subroutine zload_atoms(natm,ra,tag,nstp,strfin,strnow &
       ,sorg,myid_md,mpi_md_world)
!
!  Apply z-loading by moving atoms of top and
!  bottom layers towards opposite directions.
!
!  ifmv of top and bottom atoms should be 9.
!  And the z-motions of those atoms are controled.
!     
    implicit none
    include 'mpif.h'
    integer,intent(in):: natm,nstp,myid_md,mpi_md_world
    real(8),intent(in):: strfin,sorg(3)
    real(8),intent(inout):: ra(3,natm),strnow,tag(natm)

    integer:: i,l,ierr,ifmv

    do i=1,natm
      ifmv= int(mod(tag(i)*10,10d0))
      if( ifmv.eq.9 .and. ra(3,i)+sorg(3).gt.0.5d0 ) then  !top
        ra(3,i)= ra(3,i) +dzl
      else if( ifmv.eq.9 .and. ra(3,i)+sorg(3).le.0.5d0 ) then  !bottom
        ra(3,i)= ra(3,i) -dzl
      endif
    enddo
    zlen= zlen +2d0*dzl
    strnow= (zlen-zlen0)/zlen0
!      write(6,'(a,2es12.3e3)') ' zl,strnow=',zl,strnow

  end subroutine zload_atoms
!=======================================================================
  subroutine get_forces_on_base(natm,ra,aa,tag,h,fa2v,ftop,fbot &
       ,sorg,myid,mpi_md_world)
!
!  Compute forces on atoms at top and bottom bases,
!  and get stress on bases dividing by area.
!
    implicit none
    include 'mpif.h'
    integer,intent(in):: natm,myid,mpi_md_world
    real(8),intent(in):: ra(3,natm),aa(3,natm),tag(natm),h(3,3,0:1) &
         ,fa2v,sorg(3)
    real(8),intent(out):: ftop,fbot
    integer:: i,ifmv,ierr
    real(8):: ftopl,fbotl

    ftopl= 0d0
    fbotl= 0d0
    do i=1,natm
      ifmv= int(mod(tag(i)*10,10d0))
      if( ifmv.eq.9 .and. ra(3,i)+sorg(3).gt.0.5d0 ) then !top layer
        ftopl=ftopl +aa(3,i)
      else if( ifmv.eq.9 .and. ra(3,i)+sorg(3).le.0.5d0 ) then !bottom
        fbotl=fbotl +aa(3,i)
      endif
    enddo

    ftop= 0d0
    fbot= 0d0
    call mpi_allreduce(ftopl,ftop,1,mpi_double_precision &
         ,mpi_sum,mpi_md_world,ierr)
    call mpi_allreduce(fbotl,fbot,1,mpi_double_precision &
         ,mpi_sum,mpi_md_world,ierr)

!.....force on top/bottom in eV/A
    ftopl= ftop /fa2v
    fbotl= fbot /fa2v
    ftop = h(3,1,0)*ftopl +h(3,2,0)*ftopl +h(3,3,0)*ftopl
    fbot = h(3,1,0)*fbotl +h(3,2,0)*fbotl +h(3,3,0)*fbotl

  end subroutine get_forces_on_base
!=======================================================================
  subroutine zload_box(natm,nstp,istp,dt,strfin,strnow,h,myid)
!
!  Apply tensile loading by rescaling z-length of h-matrix.
!
    implicit none
    integer,intent(in):: natm,nstp,istp,myid
    real(8),intent(in):: strfin,dt
    real(8),intent(inout):: strnow,h(1:3,1:3,0:1)

    logical,save:: l1st= .true.
    real(8),save:: zl0,dzlfin,dzl,zl,zv0(3)

!.....at first call, determine velocity of top and bottom layers
    if( l1st ) then
!.....get z-length
      zv0(1:3)= h(1:3,3,0)
      zl0= sqrt(zv0(1)**2 +zv0(2)**2 +zv0(3)**2)
      if( myid.eq.0 ) then
        write(6,'(a)') ' z-loading parameters:'
        write(6,'(a,i5,3es12.4)') '   zv0= ',zv0(1:3)
        write(6,'(a,es12.4,a)') '   strain rate=' &
             ,strfin/100/(nstp*dt*1d-15),' /s'
      endif
      l1st=.false.
    endif

    h(1:3,3,0)= zv0(1:3) *(1d0 +strfin/100/nstp*istp)
    strnow= strfin /nstp *istp
    return
  end subroutine zload_box
!=======================================================================
end module zload
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
