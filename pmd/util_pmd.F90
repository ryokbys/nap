subroutine write_pmd0_bin(ionum,cfname,cstat,natm,tag,ra,va,h &
     ,hunit,eki,epi,strs)
  integer,intent(in):: ionum,natm
  real(8),intent(in):: tag(natm),ra(3,natm),va(3,natm) &
       ,eki(natm),epi(natm),strs(3,3,natm),h(3,3,0:1),hunit
  character(len=*),intent(in):: cfname,cstat

  integer:: ia,ib,l,i

  open(ionum,file=trim(cfname),form='unformatted' &
       ,status=trim(cstat))
  write(ionum) hunit
  write(ionum) (((h(ia,ib,l)/hunit,ia=1,3),ib=1,3),l=0,1)
  write(ionum) natm
  write(ionum) (tag(i),ra(1:3,i),va(1:3,i) &
       ,eki(i),epi(i) &
       ,strs(1,1,i),strs(2,2,i),strs(3,3,i) &
       ,strs(2,3,i),strs(3,1,i),strs(1,2,i) &
       ,i=1,natm)
  close(ionum)

end subroutine write_pmd0_bin
!=======================================================================
subroutine write_pmd0_ascii(ionum,cfname,cstat,natm,tag,ra,va,h &
     ,hunit,eki,epi,strs)
  integer,intent(in):: ionum,natm
  real(8),intent(in):: tag(natm),ra(3,natm),va(3,natm) &
       ,eki(natm),epi(natm),strs(3,3,natm),h(3,3,0:1),hunit
  character(len=*),intent(in):: cfname,cstat

  integer:: ia,ib,l,i

  open(ionum,file=trim(cfname),status=trim(cstat))
  write(ionum,'(es23.14e3)') hunit
  write(ionum,'(3es23.14e3)') (((h(ia,ib,l)/hunit,ia=1,3) &
       ,ib=1,3),l=0,1)
  write(ionum,'(i10)') natm
  do i=1,natm
    write(ionum,'(7es23.14e3,11es12.4)') tag(i),ra(1:3,i),va(1:3,i) &
         ,eki(i),epi(i) &
         ,strs(1,1,i),strs(2,2,i),strs(3,3,i) &
         ,strs(2,3,i),strs(3,1,i),strs(1,2,i)
  enddo
  close(ionum)

end subroutine write_pmd0_ascii
!=======================================================================
subroutine write_pmd_bin(ionum,cfname &
     ,natm,h,hunit,tag,ra,va,eki,epi,strs,sorg,dt)
  implicit none
  include './params_unit.h'
  integer,intent(in):: ionum
  character(len=*),intent(in) :: cfname
  integer,intent(in):: natm
  real(8),intent(in):: h(3,3,0:1),tag(natm),ra(3,natm),va(3,natm) &
       ,eki(3,3,natm),epi(natm),strs(3,3,natm),sorg(3),dt,hunit

  integer:: ia,ib,l,i

  open(ionum,file=cfname,form='unformatted' &
       ,status='replace')
  write(ionum) hunit
  write(ionum) (((h(ia,ib,l)/hunit,ia=1,3),ib=1,3),l=0,1)
  write(ionum) natm
  write(ionum) (tag(i),ra(1:3,i)+sorg(1:3),va(1:3,i) & !/dt
       ,eki(1,1,i)+eki(2,2,i)+eki(3,3,i) &
       ,epi(i) &
       ,strs(1,1,i),strs(2,2,i),strs(3,3,i) &
       ,strs(2,3,i),strs(3,1,i),strs(1,2,i) &
       ,i=1,natm)
  close(ionum)

end subroutine write_pmd_bin
!=======================================================================
subroutine read_pmd_bin(ionum,cfname &
     ,namax,natm,h,hunit,tag,ra,va,eki,epi,strs)
  implicit none
  integer,intent(in):: ionum,namax
  character(len=*),intent(in):: cfname
  integer,intent(out):: natm
  real(8),intent(out):: hunit,h(3,3,0:1),tag(namax),ra(3,namax) &
       ,va(3,namax),eki(3,3,namax),epi(namax),strs(3,3,namax)

  integer:: ia,ib,l,i

  open(ionum,file=trim(cfname),form='unformatted',status='old')
!-----natm: num. of particles in this node
  read(ionum) hunit
  read(ionum) (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
  h(1:3,1:3,0:1)= h(1:3,1:3,0:1)*hunit
  read(ionum) natm
  read(ionum) (tag(i),ra(1:3,i),va(1:3,i) &
       ,eki(1,1,i),epi(i) &
       ,strs(1,1,i),strs(2,2,i),strs(3,3,i) &
       ,strs(2,3,i),strs(3,1,i),strs(1,2,i) &
       ,i=1,natm)
  close(ionum)

end subroutine read_pmd_bin
!=======================================================================
subroutine write_pmd_ascii(ionum,cfname &
     ,natm,h,hunit,tag,ra,va,eki,epi,strs,sorg,dt)
  implicit none
  include './params_unit.h'
  integer,intent(in):: ionum
  character(len=*),intent(in) :: cfname
  integer,intent(in):: natm
  real(8),intent(in):: h(3,3,0:1),tag(natm),ra(3,natm),va(3,natm) &
       ,eki(3,3,natm),epi(natm),strs(3,3,natm),sorg(3),dt,hunit

  integer:: ia,ib,l,i

  open(ionum,file=cfname,status='replace')
  write(ionum,'(es23.14e3)') hunit
  write(ionum,'(3es23.14e3)') (((h(ia,ib,l)/hunit,ia=1,3) &
       ,ib=1,3),l=0,1)
  write(ionum,'(i10)') natm
  do i=1,natm
    write(ionum,'(7es23.14e3,11es22.14)') tag(i),ra(1:3,i) &
         +sorg(1:3),va(1:3,i) & !/dt
         ,eki(1,1,i)+eki(2,2,i)+eki(3,3,i) &
         ,epi(i) &
         ,strs(1,1,i),strs(2,2,i),strs(3,3,i) &
         ,strs(2,3,i),strs(3,1,i),strs(1,2,i)
  enddo
  close(ionum)

end subroutine write_pmd_ascii
!=======================================================================
subroutine read_pmd_ascii(ionum,cfname &
     ,namax,natm,h,hunit,tag,ra,va,eki,epi,strs)
  implicit none
  integer,intent(in):: ionum,namax
  character(len=*),intent(in):: cfname
  integer,intent(out):: natm
  real(8),intent(out):: hunit,h(3,3,0:1),tag(namax),ra(3,namax) &
       ,va(3,namax),eki(3,3,namax),epi(namax),strs(3,3,namax)

  integer:: ia,ib,l,i

  open(ionum,file=trim(cfname),status='old')
!-----natm: num. of particles in this node
  read(ionum,*) hunit
  read(ionum,*) (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
  h(1:3,1:3,0:1)= h(1:3,1:3,0:1)*hunit
  read(ionum,*) natm
  read(ionum,*) (tag(i),ra(1:3,i),va(1:3,i) &
       ,eki(1,1,i),epi(i) &
       ,strs(1,1,i),strs(2,2,i),strs(3,3,i) &
       ,strs(2,3,i),strs(3,1,i),strs(1,2,i) &
       ,i=1,natm)
  close(ionum)

end subroutine read_pmd_ascii
!=======================================================================
subroutine read_pmdtot_ascii(ionum,cfname)
  use pmdio
  implicit none
  integer,intent(in):: ionum
  character(len=*),intent(in):: cfname

  integer:: ia,ib,l,i

  open(ionum,file=trim(cfname),status='old')
  read(ionum,*) hunit
  read(ionum,*) (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
  h(1:3,1:3,0:1)= h(1:3,1:3,0:1)*hunit
  read(ionum,*) ntot0
  ntot = ntot0
  allocate(tagtot(ntot0),rtot(3,ntot0),vtot(3,ntot0),epitot(ntot0) &
       ,ekitot(3,3,ntot0),stot(3,3,ntot0),atot(3,ntot0))
  do i=1,ntot0
    read(ionum,*) tagtot(i),rtot(1:3,i),vtot(1:3,i) &
         ,ekitot(1,1,i),epitot(i) &
         ,stot(1,1,i),stot(2,2,i),stot(3,3,i) &
         ,stot(2,3,i),stot(3,1,i),stot(1,2,i)
  enddo
  close(ionum)

end subroutine read_pmdtot_ascii
!=======================================================================
subroutine write_pmdtot_ascii(ionum,cfname)
  use pmdio
  implicit none
  include './params_unit.h'
  integer,intent(in):: ionum
  character(len=*),intent(in) :: cfname

  integer:: ia,ib,l,i

  open(ionum,file=cfname,status='replace')
  write(ionum,'(es23.14e3)') hunit
  write(ionum,'(3es23.14e3)') (((h(ia,ib,l)/hunit,ia=1,3) &
       ,ib=1,3),l=0,1)
  write(ionum,'(i10)') ntot
  do i=1,ntot
    write(ionum,'(7es23.14e3,11es13.4e3)') tagtot(i) &
         ,rtot(1:3,i) &
         ,vtot(1:3,i) & !/dt
         ,ekitot(1,1,i)+ekitot(2,2,i)+ekitot(3,3,i) &
         ,epitot(i) &
         ,stot(1,1,i),stot(2,2,i),stot(3,3,i) &
         ,stot(2,3,i),stot(3,1,i),stot(1,2,i)

  enddo
  close(ionum)

end subroutine write_pmdtot_ascii
!=======================================================================
subroutine read_pmdtot_bin(ionum,cfname)
  use pmdio
  implicit none
  integer,intent(in):: ionum
  character(len=*),intent(in):: cfname

  integer:: ia,ib,l,i

  open(ionum,file=trim(cfname),form='unformatted',status='old')
!-----natm: num. of particles in this node
  read(ionum) hunit
  read(ionum) (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
  h(1:3,1:3,0:1)= h(1:3,1:3,0:1)*hunit
  read(ionum) ntot0
  ntot = ntot0
  allocate(tagtot(ntot0),rtot(3,ntot0),atot(3,ntot0) &
       ,vtot(3,ntot0),epitot(ntot0) &
       ,ekitot(3,3,ntot0),stot(3,3,ntot0))
  read(ionum) (tagtot(i),rtot(1:3,i),vtot(1:3,i) &
       ,ekitot(1,1,i),epitot(i) &
       ,stot(1,1,i),stot(2,2,i),stot(3,3,i) &
       ,stot(2,3,i),stot(3,1,i),stot(1,2,i) &
       ,i=1,ntot0)
  close(ionum)

end subroutine read_pmdtot_bin
!=======================================================================
subroutine write_pmdtot_bin(ionum,cfname)
  use pmdio
  implicit none
  include './params_unit.h'
  integer,intent(in):: ionum
  character(len=*),intent(in) :: cfname

  integer:: ia,ib,l,i

  open(ionum,file=cfname,form='unformatted' &
       ,status='replace')
  write(ionum) hunit
  write(ionum) (((h(ia,ib,l)/hunit,ia=1,3),ib=1,3),l=0,1)
  write(ionum) ntot
  do i=1,ntot
    write(ionum) tagtot(i),rtot(1:3,i),vtot(1:3,i) & !/dt
         ,ekitot(1,1,i)+ekitot(2,2,i)+ekitot(3,3,i) &
         ,epitot(i) &
         ,stot(1,1,i),stot(2,2,i),stot(3,3,i) &
         ,stot(2,3,i),stot(3,1,i),stot(1,2,i)
  enddo
  close(ionum)

end subroutine write_pmdtot_bin
!=======================================================================
subroutine read_chgtot_ascii(ionum,cfname)
  use pmdio
  implicit none
  integer,intent(in):: ionum
  character(len=*),intent(in):: cfname

  integer:: ia,ib,l,i,ntott,itmp

  open(ionum,file=trim(cfname),status='old')
  read(ionum,*) ntott
  allocate(chgtot(ntott),chitot(ntott))
  do i=1,ntott
    read(ionum,*) itmp,chgtot(i),chitot(i)
  enddo
  close(ionum)

end subroutine read_chgtot_ascii
!=======================================================================
subroutine read_chgtot_bin(ionum,cfname)
  use pmdio
  implicit none
  integer,intent(in):: ionum
  character(len=*),intent(in):: cfname

  integer:: ia,ib,l,i,ntott,itmp

  open(ionum,file=trim(cfname),form='unformatted',status='old')
  read(ionum) ntott
  allocate(chgtot(ntott),chitot(ntott))
  do i=1,ntott
    read(ionum) itmp,chgtot(i),chitot(i)
  enddo
  close(ionum)

end subroutine read_chgtot_bin
!=======================================================================
subroutine write_chgtot_ascii(ionum,cfname)
  use pmdio
  use util,only: itotOf
  implicit none
  integer,intent(in):: ionum
  character(len=*),intent(in) :: cfname
!!$  integer,external:: itotOf

  integer:: ia,ib,l,i

  open(ionum,file=cfname,status='replace')
  write(ionum,'(i10)') ntot
  do i=1,ntot
    write(ionum,'(es23.14e3,2es13.4e3)') itotOf(tagtot(i)) &
         ,chgtot(i),chitot(i)
  enddo
  close(ionum)

end subroutine write_chgtot_ascii
!=======================================================================
subroutine write_chgtot_bin(ionum,cfname)
  use pmdio
  use util,only: itotOf
  implicit none
  integer,intent(in):: ionum
  character(len=*),intent(in) :: cfname
!!$  integer,external:: itotOf

  integer:: ia,ib,l,i

  open(ionum,file=cfname,form='unformatted',status='replace')
  write(ionum) ntot
  do i=1,ntot
    write(ionum) itotOf(tagtot(i)) &
         ,chgtot(i),chitot(i)
  enddo
  close(ionum)

end subroutine write_chgtot_bin
!=======================================================================
subroutine write_dump(ionum,cfname)
!
!     Write atomic configuration in LAMMPS-dump format file.
!
  use pmdio
  use util,only: itotOf
  implicit none 
  integer,intent(in):: ionum
  character(len=*),intent(in) :: cfname

  integer:: i,k,l
  real(8):: xi(3),ri(3),eki,epi,xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz, &
       xlo_bound,xhi_bound,ylo_bound,yhi_bound, &
       zlo_bound,zhi_bound,st(3,3)
!!$  integer,external:: itotOf
  real(8),allocatable,save:: rlmp(:,:)

  real(8),parameter:: tiny = 1d-14

  if( .not. allocated(rlmp) ) allocate(rlmp(3,ntot))
  if( size(rlmp).ne.3*ntot ) then
    deallocate(rlmp)
    allocate(rlmp(3,ntot))
  endif

  open(ionum,file=trim(cfname),status='replace')
  write(ionum,'(a)') 'ITEM: TIMESTEP'
  write(ionum,'(i10)') 0
  write(ionum,'(a)') 'ITEM: NUMBER OF ATOMS'
  write(ionum,'(i10)') ntot
  write(ionum,'(a)') 'ITEM: BOX BOUNDS xy xz yz'
  call pmd2lammps(h,ntot,rtot,rlmp,xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz)
  xlo_bound = xlo +min(0d0, xy, xz, xy+xz)
  xhi_bound = xhi +max(0d0, xy, xz, xy+xz)
  ylo_bound = ylo +min(0d0, yz)
  yhi_bound = yhi +max(0d0, yz)
  zlo_bound = zlo
  zhi_bound = zhi
  write(ionum,'(3f15.4)') xlo_bound, xhi_bound, xy
  write(ionum,'(3f15.4)') ylo_bound, yhi_bound, xz
  write(ionum,'(3f15.4)') zlo_bound, zhi_bound, yz
  write(ionum,'(a)') 'ITEM: ATOMS id type x y z ekin epot' &
       //' sxx syy szz syz sxz sxy chg chi'
  do i=1,ntot
!        xi(1:3)= rtot(1:3,i)
!        ri(1:3)= h(1:3,1,0)*xi(1) +h(1:3,2,0)*xi(2) +h(1:3,3,0)*xi(3)
    eki = ekitot(1,1,i) +ekitot(2,2,i) +ekitot(3,3,i)
    epi = epitot(i)
    st(1:3,1:3) = stot(1:3,1:3,i)
    if( eki.lt.tiny ) eki = 0d0
    if( abs(epi).lt.tiny ) epi = 0d0
    do l=1,3
      do k=l,3
        if( abs(st(l,k)).lt.tiny ) st(l,k) = 0d0
      enddo
    enddo
    write(ionum,'(i8,i3,3f12.5,8es11.3,f9.4,f9.2)') &
         itotOf(tagtot(i)),int(tagtot(i)),rlmp(1:3,i),eki, &
         epi, &
         st(1,1),st(2,2),st(3,3), &
         st(2,3),st(1,3),st(1,2), &
         chgtot(i),chitot(i)
  enddo

  close(ionum)
end subroutine write_dump
!=======================================================================
subroutine pmd2lammps(h,ntot,rtot,rlmp, &
     xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz)
!
!     Convert pmd data format to LAMMPS format not only cell
!     but also atomic positions.
!     Parameters to be output:
!       xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz
!     LAMMPS cell should be defined as,
!       a = ( xhi-xlo,       0,       0 )
!       b = (      xy, yhi-hlo,       0 )
!       c = (      xz,      yz, zhi-zlo )
!     See, http://lammps.sandia.gov/doc/Section_howto.html, for detail.
!
  implicit none
  integer,intent(in):: ntot
  real(8),intent(in):: h(3,3),rtot(3,ntot)
  real(8),intent(out):: xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz &
       ,rlmp(3,ntot)

  integer:: i,lxy,lxz,lyz
  real(8):: a0(3),b0(3),c0(3),a1(3),a2(3),a3(3) &
       ,b1(3),b2(3),b3(3),rt(3),amat(3,3),bmat(3,3) &
       ,x,y,z,a23(3),a31(3),a12(3),vol
  real(8):: a,b,c,alpha,beta,gamma
  real(8),external:: absv,sprod

  xlo = 0d0
  ylo = 0d0
  zlo = 0d0
  a0(1:3) = h(1:3,1)
  b0(1:3) = h(1:3,2)
  c0(1:3) = h(1:3,3)
  a = absv(3,a0)
  b = absv(3,b0)
  c = absv(3,c0)
  alpha = acos(sprod(3,b0,c0)/b/c)
  beta  = acos(sprod(3,a0,c0)/a/c)
  gamma = acos(sprod(3,a0,b0)/a/b)
  xhi = a
  xy = b*cos(gamma)
  xz = c*cos(beta)
  yhi= sqrt(b*b -xy*xy)
  yz = (b*c*cos(alpha) -xy*xz)/yhi
  zhi= sqrt(c*c -xz*xz -yz*yz)
  x = xhi -xlo
  y = yhi -ylo
  z = zhi -zlo
  lxy = 0
  if( xy.gt.xhi/2 ) then
    xy = xy -xhi
    lxy = -1
  else if( xy.lt.-xhi/2 ) then
    xy = xy +xhi
    lxy = 1
  endif
  lxz = 0
  if( xz.gt.xhi/2 ) then
    xz = xz -xhi
    lxz = -1
  else if( xz.lt.-xhi/2 ) then
    xz = xz +xhi
    lxz = 1
  endif
  lyz = 0
  if( yz.gt.yhi/2 ) then
    yz = yz -yhi
    lyz = -1
  else if( yz.lt.-yhi/2 ) then
    yz = yz +yhi
    lyz = 1
  endif

  a1(1:3) = h(1:3,1)
  a2(1:3) = h(1:3,2)
  a3(1:3) = h(1:3,3)
  call vprod(a2,a3,a23)
  call vprod(a3,a1,a31)
  call vprod(a1,a2,a12)
  vol = abs( sprod(3,a1,a23) )
  amat(1:3,1:3) = 0d0
  amat(1,1:3) = a23(1:3)
  amat(2,1:3) = a31(1:3)
  amat(3,1:3) = a12(1:3)
  b1(1:3) = (/ x, 0d0, 0d0 /)
  b2(1:3) = (/ xy,  y, 0d0 /)
  b3(1:3) = (/ xz, yz,   z /)
  bmat(1:3,1:3) = 0d0
  bmat(1:3,1) = b1(1:3)
  bmat(1:3,2) = b2(1:3)
  bmat(1:3,3) = b3(1:3)
  do i=1,ntot
    rlmp(1:3,i) = 0d0
    call shift_pos_for_lammps(rtot(1,i),rlmp(1,i),lxy,lxz,lyz &
         ,x,y,z,yz,xz,xy)
    rt = matmul(h,rlmp(1:3,i))
    rlmp(1:3,i) = matmul(bmat,matmul(amat,rt))/vol
  enddo

  return
end subroutine pmd2lammps
!=======================================================================
subroutine shift_pos_for_lammps(r,rn,lxy,lxz,lyz,x,y,z,yz,xz,xy)
  use pmdio,only: boundary
  implicit none
  real(8),intent(in):: r(3),x,y,z,yz,xz,xy
  integer,intent(in):: lxy,lxz,lyz
  real(8),intent(out):: rn(3)

  integer:: i
  real(8):: xyp
  real(8),external:: pbc

  xyp = xy -lxy*x
  rn(1:3) = r(1:3)
  rn(2) = rn(2) -lyz*r(3)
  rn(1) = rn(1) -lxz*r(3) +(r(2)*xyp -rn(2)*xy)/x
  do i=1,3
    if( boundary(i:i).eq.'p' ) then
      rn(i) = pbc(rn(i))
    endif
  enddo
  return
end subroutine shift_pos_for_lammps
!=======================================================================
function pbc(x)
  implicit none
  real(8),intent(in):: x
  real(8):: pbc

  if( x.lt.0d0 ) then
    pbc = x -int(x) +1d0
  else if( x.ge.1d0 ) then
    pbc = x -int(x)
  else
    pbc = x
  endif
  return
end function pbc
!=======================================================================
subroutine hmat2lammps(h,xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz)
!
!     Convert h-matrix to LAMMPS cell vectors.
!     Parameters to be output:
!       xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz
!     LAMMPS cell should be defined as,
!       a = ( xhi-xlo,       0,       0 )
!       b = (      xy, yhi-hlo,       0 )
!       c = (      xz,      yz, zhi-zlo )
!     See, http://lammps.sandia.gov/doc/Section_howto.html, for detail.
!
  implicit none
  real(8),intent(in):: h(3,3)
  real(8),intent(out):: xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz

  real(8):: a0(3),b0(3),c0(3)
  real(8):: a,b,c,alpha,beta,gamma
  real(8),external:: absv,sprod

  xlo = 0d0
  ylo = 0d0
  zlo = 0d0
  a0(1:3) = h(1:3,1)
  b0(1:3) = h(1:3,2)
  c0(1:3) = h(1:3,3)
  a = absv(3,a0)
  b = absv(3,b0)
  c = absv(3,c0)
  alpha = acos(sprod(3,b0,c0)/b/c)
  beta  = acos(sprod(3,a0,c0)/a/c)
  gamma = acos(sprod(3,a0,b0)/a/b)
  xhi = a
  xy = b*cos(gamma)
  xz = c*cos(beta)
  yhi= sqrt(b*b -xy*xy)
  yz = (b*c*cos(alpha) -xy*xz)/yhi
  zhi= sqrt(c*c -xz*xz -yz*yz)
  return
end subroutine hmat2lammps
!=======================================================================
function num_data(str,delim)
  implicit none
  character(len=*),intent(in):: str
  character(len=1),intent(in):: delim
  integer:: num_data

  integer:: i

  i=1
  num_data = 0
  do
    if( i.gt.len(str) ) exit
    if( str(i:i).ne.delim ) then
      num_data = num_data + 1
      do
        i = i + 1
        if( i.gt.len(str) ) exit
        if( str(i:i).eq.delim ) exit
      end do
    end if
    i = i + 1
  end do
  return
end function num_data
!=======================================================================
function is_numeric(str)
!
!  Check if the given string is numeric or not.
!
  character(len=*),intent(in):: str
  logical:: is_numeric
  integer:: i
  i = verify(trim(str),'0123456789')
  is_numeric = i==0
  return
end function is_numeric
!=======================================================================
subroutine time_stamp(prefix)
  implicit none
  character(len=*),intent(in):: prefix
  character(len=10):: c1,c2,c3
  integer:: date(8)
  character(len=128):: cdate,ctime

  call date_and_time(c1,c2,c3,date)

  write(ctime,'(i0.2,a,i0.2,a,i0.2)') date(5),':',date(6) &
       ,':',date(7)
  write(cdate,'(i4,a,i0.2,a,i0.2)') date(1),'-',date(2),'-',date(3)
  write(6,'(a,1x,a,1x,a)') prefix, 'at '//trim(ctime), &
       'on '//trim(cdate)
  return
end subroutine time_stamp
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
