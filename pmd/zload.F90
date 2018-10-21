module zload
!-----------------------------------------------------------------------
!  Module for loading on plane perpendicular to z-axis
!                     Last-modified: <2018-10-21 20:05:57 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
  implicit none
  save
!.....initial z top and bottom positions
  real(8):: ztop0, zbot0, zlen0, dlen0, dzl, dl
  integer:: nztop, nzbot
!.....variables that change during z-loading simulation
  real(8):: zlen, ztop, zbot
  real(8):: dlen
!.....unit vector of shear direction (ux,uy)
  real(8):: uvx, uvy
!.....deviation along a1,a2 vectors in case of shear
  real(8):: d1,d2
  
contains
  subroutine set_zload_atoms(natm,ra,tag,h,fmv,sorg,strfin,nstp &
       ,zskin_width,myid_md,mpi_md_world,iprint)
!
!     Choose atoms to be applied z-loading.
!     ifmv's of top and bottom atoms within zskin_width in z-direction
!     are set to 9.
!
    implicit none
    include 'mpif.h'
    integer,intent(in):: natm,myid_md,mpi_md_world,nstp,iprint
    real(8),intent(in):: ra(3,natm),sorg(3),strfin,h(3,3),zskin_width
    real(8),intent(inout):: tag(natm),fmv(3,0:9)
    integer:: ierr,i,nztopl,nzbotl
    real(8):: ztopl,zbotl,dzlfin

    if( myid_md.eq.0 .and. iprint.gt.0 ) then
      print *,''
      print *,'Setup zload_atoms'
    endif

!.....Set z-fmv of 9th to 0d0, which means no motion along z
    fmv(3,9) = 0d0
    
!.....Detect initial top and bottom positions
    ztopl = 0d0
    zbotl = 1d0
    do i=1,natm
      zbotl= min(zbotl,ra(3,i)+sorg(3))
      ztopl= max(ztopl,ra(3,i)+sorg(3))
    enddo
    ztop0= 0d0
    zbot0= 0d0
    call mpi_allreduce(ztopl,ztop0,1,mpi_real8,mpi_max &
         ,mpi_md_world,ierr)
    call mpi_allreduce(zbotl,zbot0,1,mpi_real8,mpi_min &
         ,mpi_md_world,ierr)
    zlen0 = ztop0 -zbot0
    dzlfin = zlen0 *strfin
    dzl= dzlfin /nstp /2
    zlen = zlen0
    ztop = ztop0
    zbot = zbot0
    if( myid_md.eq.0 .and. iprint.gt.0 ) then
      write(6,'(a,2es12.4)') ' zlen0,dzl = ',zlen0,dzl
    endif

!.....Define top and bottom atoms z-motions of which are controlled
    nztopl = 0
    nzbotl = 0
    do i=1,natm
      if( ra(3,i)+sorg(3) .gt. ztop -zskin_width/h(3,3) ) then
        call replaceTag('ifmv',9,tag(i))
        nztopl = nztopl +1
      else if( ra(3,i)+sorg(3) .lt. zbot +zskin_width/h(3,3) ) then
        call replaceTag('ifmv',9,tag(i))
        nzbotl = nzbotl +1
      endif
    enddo

    call mpi_allreduce(nztopl,nztop,1,mpi_integer,mpi_max &
         ,mpi_md_world,ierr)
    call mpi_allreduce(nzbotl,nzbot,1,mpi_integer,mpi_max &
         ,mpi_md_world,ierr)
    if( myid_md.eq.0 .and. iprint.gt.0 ) then
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
  subroutine set_shear(natm,ra,tag,h,fmv,sorg,strfin,nstp &
       ,zskin_width,zshear_angle,myid_md,mpi_md_world,iprint)
!
!  Choose atoms to be applied shear loading.
!  ifmv's of top and bottom atoms within zskin_width in z-direction
!  are set to 9.
!
    implicit none
    include 'mpif.h'
    integer,intent(in):: natm,myid_md,mpi_md_world,nstp,iprint
    real(8),intent(in):: ra(3,natm),sorg(3),strfin,h(3,3),zskin_width &
         ,zshear_angle
    real(8),intent(inout):: tag(natm),fmv(3,0:9)
    integer:: ierr,i,nztopl,nzbotl
    real(8):: ztopl,zbotl,dlfin,zskin,angle
    real(8):: x1,y1,x2,y2,amati(2,2),det
    real(8),external:: absv
    real(8),parameter:: pi = 3.14159265358979d0

    if( myid_md.eq.0 .and. iprint.gt.0 ) then
      print *,''
      print *,'Setup shear'
    endif

!.....Set fmv of 9th to zero, which means no motion except z
    fmv(1:2,9) = 0d0

!.....Unit vector along shear from zshear_angle
    angle = mod(zshear_angle,180.0)
    angle = angle/180.0 *pi
    uvx = cos(angle)
    uvy = sin(angle)

!.....Detect initial top and bottom positions
    ztopl = 0d0
    zbotl = 1d0
    do i=1,natm
      zbotl= min(zbotl,ra(3,i)+sorg(3))
      ztopl= max(ztopl,ra(3,i)+sorg(3))
    enddo
    ztop0= 0d0
    zbot0= 0d0
    call mpi_allreduce(ztopl,ztop0,1,mpi_real8,mpi_max &
         ,mpi_md_world,ierr)
    call mpi_allreduce(zbotl,zbot0,1,mpi_real8,mpi_min &
         ,mpi_md_world,ierr)
    zlen0 = (ztop0 -zbot0)*absv(3,h(1:3,3))
    dlfin = zlen0 *strfin
    dl= (dlfin /nstp /2)
!.....Conversion of shear from x,y to a1,a2
    x1 = h(1,1)
    y1 = h(2,1)
    x2 = h(1,2)
    y2 = h(2,2)
    det = x1*y2 -x2*y1
    amati(1,1) =  y2 /det
    amati(1,2) = -x2 /det
    amati(2,1) = -y1 /det
    amati(2,2) =  x1 /det
    d1 = amati(1,1) *uvx*dl +amati(1,2)*uvy*dl
    d2 = amati(2,1) *uvx*dl +amati(2,2)*uvy*dl
   
    dlen0 = 0d0
    dlen = dlen0
    ztop = ztop0
    zbot = zbot0
    if( myid_md.eq.0 .and. iprint.gt.0 ) then
      write(6,'(a,3es12.4)') ' zlen0,dlen0,dl = ',zlen0,dlen0,dl
      write(6,'(a,5es12.3)') ' angle,uvx,uvy,d1,d2 = ',angle,uvx,uvy,d1,d2
    endif

!.....Define top and bottom atoms z-motions of which are controlled
    nztopl = 0
    nzbotl = 0
    zskin = zskin_width /h(3,3)
    do i=1,natm
      if( ra(3,i)+sorg(3) .gt. ztop -zskin ) then
        call replaceTag('ifmv',9,tag(i))
        nztopl = nztopl +1
      else if( ra(3,i)+sorg(3) .lt. zbot +zskin ) then
        call replaceTag('ifmv',9,tag(i))
        nzbotl = nzbotl +1
      endif
    enddo

    call mpi_allreduce(nztopl,nztop,1,mpi_integer,mpi_max &
         ,mpi_md_world,ierr)
    call mpi_allreduce(nzbotl,nzbot,1,mpi_integer,mpi_max &
         ,mpi_md_world,ierr)
    if( myid_md.eq.0 .and. iprint.gt.0 ) then
      write(6,'(a,2i8)') ' nztop, nzbot = ',nztop,nzbot
    endif

  end subroutine set_shear
!=======================================================================
  subroutine shear_atoms(natm,ra,tag,nstp,strfin,strnow &
       ,sorg,myid_md,mpi_md_world)
!
!  Apply shear by moving atoms of top-z layers within the plane normal to z.
!
!  ifmv of top and bottom atoms should be 9.
!  And the xy-motions of those atoms are controled.
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
        ra(1,i)= ra(1,i) +d1
        ra(2,i)= ra(2,i) +d2
      else if( ifmv.eq.9 .and. ra(3,i)+sorg(3).le.0.5d0 ) then  !bottom
        ra(1,i)= ra(1,i) -d1
        ra(2,i)= ra(2,i) -d2
      endif
    enddo
    dlen= dlen +2d0*dl
    strnow= (dlen-dlen0)/zlen0
!!$    if( myid_md.eq.0 ) write(6,'(a,2es12.3e3)') ' xlen,strnow=',xlen,strnow

  end subroutine shear_atoms
!=======================================================================
  subroutine get_forces_on_base(natm,ra,aa,tag,h,ftop,fbot &
       ,sorg,myid,mpi_md_world,iprint,czload_type)
!
!  Compute forces on atoms at top and bottom bases,
!  and get stress on bases dividing by area.
!  Assumes that the z-axis (3) is normal to xy-plane.
!
    implicit none
    include 'mpif.h'
    include './params_unit.h'
    integer,intent(in):: natm,myid,mpi_md_world,iprint
    real(8),intent(in):: ra(3,natm),aa(3,natm),tag(natm),h(3,3) &
         ,sorg(3)
    real(8),intent(out):: ftop,fbot
    character(len=*),intent(in):: czload_type
    integer:: i,is,ifmv,ierr,ixyz
    real(8):: ftopl,fbotl,xyarea,a(3),b(3),axb(3),fx,fy
    real(8),external:: absv

    ftopl= 0d0
    fbotl= 0d0
    if( trim(czload_type).eq.'atoms' ) then
      do i=1,natm
        ifmv= int(mod(tag(i)*10,10d0))
!.....Scaled force to real force to get eV/Ang unit
        if( ifmv.eq.9 .and. ra(3,i)+sorg(3).gt.0.5d0 ) then !top layer
          ftopl=ftopl +(h(3,1)*aa(1,i) +h(3,2)*aa(2,i) &
               +h(3,3)*aa(3,i) )
        else if( ifmv.eq.9 .and. ra(3,i)+sorg(3).le.0.5d0 ) then !bottom
          fbotl=fbotl +(h(3,1)*aa(1,i) +h(3,2)*aa(2,i) &
               +h(3,3)*aa(3,i) )
        endif
      enddo
    else if( trim(czload_type).eq.'shear' ) then
      do i=1,natm
        ifmv= int(mod(tag(i)*10,10d0))
!.....Scaled force to real force to get eV/Ang unit
        if( ifmv.eq.9 .and. ra(3,i)+sorg(3).gt.0.5d0 ) then !top layer
          fx=h(1,1)*aa(1,i) +h(1,2)*aa(2,i) &
               +h(1,3)*aa(3,i)
          fy=h(1,1)*aa(1,i) +h(1,2)*aa(2,i) &
               +h(1,3)*aa(3,i)
          ftopl = ftopl +(fx*uvx +fy*uvy)
        else if( ifmv.eq.9 .and. ra(3,i)+sorg(3).le.0.5d0 ) then !bottom
          fx=h(1,1)*aa(1,i) +h(1,2)*aa(2,i) &
               +h(1,3)*aa(3,i)
          fy=h(1,1)*aa(1,i) +h(1,2)*aa(2,i) &
               +h(1,3)*aa(3,i)
          fbotl = fbotl +(fx*uvx +fy*uvy)
        endif
      enddo
    endif

    ftop= 0d0
    fbot= 0d0
    call mpi_allreduce(ftopl,ftop,1,mpi_real8 &
         ,mpi_sum,mpi_md_world,ierr)
    call mpi_allreduce(fbotl,fbot,1,mpi_real8 &
         ,mpi_sum,mpi_md_world,ierr)

!.....Divide forces on top/bottom by xy-area to get stress
    a(1:3)= h(1:3,1)
    b(1:3)= h(1:3,2)
    call vprod(a,b,axb)
    xyarea = absv(3,axb)
    ftop = ftop /xyarea *up2gpa
    fbot = fbot /xyarea *up2gpa

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
