      program put_interstitials
c-----------------------------------------------------------------------
c Put interstitial atoms into the given crystal/poly-crystal system
c read from pmd-style file.
c-----------------------------------------------------------------------
c Usage:
c   $ ./intrstl num-of-interstitials pmd00000 outfile
c-----------------------------------------------------------------------
      implicit none
      integer,parameter:: namax= 100000
      integer,parameter:: nnmax= 20
      real(8),parameter:: rcut = 3.0d0
      character,parameter:: ciofmt*6 = 'ascii '
      integer,parameter:: is_int= 2

      integer:: nargc,natm,ntet,ntetmax
      real(8):: tag(namax),ra(3,namax),va(3,namax),eki(3,3,namax)
     &     ,epi(namax),strs(3,3,namax),sorg(3),dt,hunit,h(3,3,0:1)
     &     ,hi(3,3),tagi,ri(3)
      integer:: lspr(0:nnmax,namax)
      real(8),allocatable:: ptet(:,:),vtet(:)

      integer:: i,nint,num
      character(len=128):: ciname,coname,cnum

      integer,external:: iargc

c.....Read an argument, ipmd
      nargc=iargc()
      if( nargc.ne.3 ) then
        stop 'Usage: $ ./intrstl num-of-intestials pmd00000'//
     &       ' pmd00000put'
      endif
      call getarg(1,cnum)
      read(cnum,*) nint
      call getarg(2,ciname)
      call getarg(3,coname)

      write(6,'(a,i6)') ' num of interstitials=',nint

      if( ciofmt .eq. 'ascii') then
        call read_pmd_ascii(10,ciname,namax,natm,h,hunit,tag,ra,va
     &       ,eki,epi,strs)
      else if( ciofmt .eq. 'binary' ) then
        call read_pmd_bin(10,ciname,namax,natm,h,hunit,tag,ra,va
     &       ,eki,epi,strs)
      else
        write(6,'(a)') ' [Error] ciofmt must be ascii or binary !!!'
        stop
      endif

c.....make hi by inverting h
      call invert(h,hi)

c.....make pair list
      call mk_lspr_sngl(namax,natm,nnmax,tag,ra,rcut,h,hi,lspr)

c.....search tetrahedra around each atom and compute its volume
      ntetmax= natm *4
      allocate(ptet(3,ntetmax),vtet(ntetmax))
      call make_tetrahedra(natm,nnmax,lspr,ra,h,ntetmax,ntet,ptet,vtet)

c      write(6,'(a)') ' before sorting'
c      do i=1,ntet
c        write(6,'(4f10.5)') ptet(1:3,i),vtet(i)
c      enddo
c.....sort tetra according to its volume
      call heapsort(ntet,ntetmax,vtet,3,ptet)
c      write(6,'(a)') ' after sorting'
c      do i=1,ntet
c        write(6,'(4f10.5)') ptet(1:3,i),vtet(i)
c      enddo

      if( nint.gt.ntet ) then
        write(6,'(a)') ' [Error] nint.gt.ntet !!!'
        stop
      endif

c.....choose 'num-of-interstitial' sites to be occupied
      do i=1,nint
        natm= natm+1
        ri(1:3)= ptet(1:3,ntet-(i-1))
        tagi= 1d0*is_int +1d-1 +1d-14*natm
        ra(1:3,natm)= ri(1:3)
        tag(natm)= tagi
      enddo

      write(6,'(a)') ' here2'
      if( ciofmt .eq. 'ascii') then
        call write_pmd0_ascii(20,coname,'replace',natm,tag
     &       ,ra,va,h,hunit,eki,epi,strs)
      else if( ciofmt .eq. 'binary' ) then
        call write_pmd0_bin(20,coname,'replace',natm,tag
     &       ,ra,va,h,hunit,eki,epi,strs)
      endif

      end program put_interstitials
c=======================================================================
      subroutine make_tetrahedra(natm,nnmax,lspr,ra,h
     &     ,ntetmax,ntet,ptet,vtet)
      implicit none
      integer,intent(in):: natm,nnmax,lspr(0:nnmax,natm),ntetmax
      real(8),intent(in):: ra(3,natm),h(3,3)
      integer,intent(out):: ntet
      real(8),intent(out):: ptet(3,ntetmax) ! position
      real(8),intent(out):: vtet(ntetmax)   ! volume

      integer:: i,j,k,l,jj,kk,ll,itet
      integer,allocatable:: iqrtt(:,:) ! quartet of indices
      real(8):: a(3),b(3),c(3),ha(3),hb(3),hc(3),hbc(3)
     &     ,ri(3),rj(3),rk(3),rl(3)

      logical,external:: lexists
      real(8),external:: dot

      allocate(iqrtt(4,ntetmax))

      ntet= 0
c.....every atom-i
      do i=1,natm-3
        ri(1:3)= ra(1:3,i)
c.....j: neighbor of atom-i
        do jj=1,lspr(0,i)
          j= lspr(jj,i)
c.....only i<j<k<l
          if( j.le.i ) cycle
          rj(1:3)= ra(1:3,j)
c.....loop for neighbors of atom-i
          do kk=1,lspr(0,i)
            k= lspr(kk,i)
            if( .not. lexists(k,nnmax,natm,lspr,j) ) cycle
c.....only i<j<k<l
            if( k.le.j ) cycle
            do ll=1,lspr(0,i)
              l= lspr(ll,i)
c.....only i<j<k<l
              if( l.le.k ) cycle
              if( .not. lexists(l,nnmax,natm,lspr,j) ) cycle
              ntet= ntet +1
              iqrtt(1:4,ntet)= (/i,j,k,l/)
            enddo
          enddo
        enddo
      enddo
      write(6,'(a,i8)') ' num of tetrahedron=',ntet
      
c.....get tetrahedron positions and volumes
      do itet=1,ntet
c        write(6,'(4i5)') iqrtt(1:4,itet)
        ri(1:3)= ra(1:3,iqrtt(1,itet))
        rj(1:3)= ra(1:3,iqrtt(2,itet))
        rk(1:3)= ra(1:3,iqrtt(3,itet))
        rl(1:3)= ra(1:3,iqrtt(4,itet))
c.....center of the tetrahedron
        ptet(1:3,itet)= (ri(1:3)+rj(1:3)+rk(1:3)+rl(1:3))/4
c.....volume= (1/6)*a*(b x c)
        a(1:3)=rj(1:3)-ri(1:3)-anint(rj(1:3)-ri(1:3))
        b(1:3)=rk(1:3)-ri(1:3)-anint(rk(1:3)-ri(1:3))
        c(1:3)=rl(1:3)-ri(1:3)-anint(rl(1:3)-ri(1:3))
        ha(1:3)= h(1:3,1)*a(1) +h(1:3,2)*a(2) +h(1:3,3)*a(3)
        hb(1:3)= h(1:3,1)*b(1) +h(1:3,2)*b(2) +h(1:3,3)*b(3)
        hc(1:3)= h(1:3,1)*c(1) +h(1:3,2)*c(2) +h(1:3,3)*c(3)
c        write(6,'(a,3f10.5,f10.5)') 'ha=',ha(1:3)
c     &       ,sqrt(ha(1)**2+ha(2)**2+ha(3)**2)
c        write(6,'(a,3f10.5,f10.5)') 'hb=',hb(1:3)
c     &       ,sqrt(hb(1)**2+hb(2)**2+hb(3)**2)
c        write(6,'(a,3f10.5,f10.5)') 'hc=',hc(1:3)
c     &       ,sqrt(hc(1)**2+hc(2)**2+hc(3)**2)
        call cross(hb,hc,hbc)
        vtet(itet)= abs(dot(ha,hbc)/6)
      enddo

      end subroutine make_tetrahedra
c=======================================================================
      function lexists(i,nnmax,natm,lspr,j)
c
c Check whether atom-i exists in j-neighbor
c
      implicit none
      integer,intent(in):: i,j,nnmax,natm,lspr(0:nnmax,natm)

      integer:: l
      logical:: lexists
      
      lexists= .false.
      do l=1,lspr(0,j)
        if( i.eq.lspr(l,j) ) then
          lexists= .true.
          return
        endif
      enddo
c.....return false (there is no i in j-neighbor)
      return
      end function lexists
c=======================================================================
      function dot(v1,v2)
      implicit none
      real(8),intent(in):: v1(3),v2(3)
      real(8):: dot

      dot= v1(1)*v2(1) +v1(2)*v2(2) +v1(3)*v2(3)
      return
      end function dot
c=======================================================================
      subroutine cross(v1,v2,ans)
      implicit none
      real(8),intent(in):: v1(3),v2(3)
      real(8),intent(out):: ans(3)

      ans(1)= v1(2)*v2(3) -v1(3)*v2(2)
      ans(2)= v1(3)*v2(1) -v1(1)*v2(3)
      ans(3)= v1(1)*v2(2) -v1(2)*v2(1)
      return
      end subroutine cross
c=======================================================================
      subroutine invert(h,hi)
      implicit none
      real(8),intent(in):: h(3,3)
      real(8),intent(out):: hi(3,3)

      integer:: i,j,im,ip,jm,jp
      real(8):: sgm(3,3),hit(3,3),vol

c-----cofactor matrix, SGM
      do j=1,3
        jm=mod(j+1,3)+1
        jp=mod(j,  3)+1
        do i=1,3
          im=mod(i+1,3)+1
          ip=mod(i,  3)+1
          sgm(i,j)=h(ip,jp)*h(im,jm)-h(im,jp)*h(ip,jm)
        enddo
      enddo
c-----MD-box volume
      vol=h(1,1)*sgm(1,1)+h(2,1)*sgm(2,1)+h(3,1)*sgm(3,1)
      do j=1,3
        do i=1,3
          hit(i,j)= sgm(i,j)/vol
        enddo
      enddo
c-----transpose
      do j=1,3
        do i=1,3
          hi(i,j)= hit(j,i)
        enddo
      enddo
      
      end subroutine invert
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make intstl"
c     End:
