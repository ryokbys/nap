      program make_slices
c-----------------------------------------------------------------------
c  USAGE:
c    $ ./10mkslice <#-of-slice>
c  <#-of-slice> includes initial and final slices.
c-----------------------------------------------------------------------
      implicit none
      integer,parameter:: namax = 100000
      integer:: naini,nafin,i,j,ia,ib,nslc,natm
      real(8):: raini(3,namax),rafin(3,namax)
     &     ,tagini(namax),tagfin(namax),hini(3,3,0:1),hfin(3,3,0:1)
     &     ,hunit
      real(8):: rd(3,namax),h(3,3,0:1),tag(namax),va(3,namax)
      character:: cnum*3,argv*128
      integer:: nargc
      integer,external:: iargc
      real(8):: tempini = 100d0

c.....Get num of slices which includes initial and final slices
      nargc= iargc()
      if( nargc.ne.1 .and. nargc.ne.2 ) then
        stop ' [Usage] ./10mkslice <#-of-slice> (<tempmini>)'
      endif
      call getarg(1,argv)
      read(argv,*) nslc
      if( nslc.le.2 ) stop ' There is no need to run this program...'
c.....Read tempini if there is
      if( nargc.eq.2 ) then
        call getarg(2,argv)
        read(argv,*) tempini
      endif

c.....Get initial and final atom configurations
      call read_atom_conf(10,'nebini',namax,naini,hini,hunit
     &     ,raini,tagini)
      call read_atom_conf(10,'nebfin',namax,nafin,hfin,hunit
     &     ,rafin,tagfin)
c.....Check some
      if( naini.ne.nafin ) stop ' [Error] naini .ne. nafin !!!'
      natm= naini
      
c.....Make slices
      h(1:3,1:3,0:1)= hini(1:3,1:3,0:1)
      tag(1:natm)= tagini(1:natm)
      rd(1:3,1:natm)= (rafin(1:3,1:natm)-raini(1:3,1:natm)
     &     -anint(rafin(1:3,1:natm)-raini(1:3,1:natm)))/(nslc-1)
      do i=1,nslc
        call setv(natm,va,tag,tempini)
        write(cnum(1:3),'(i3.3)') i
        write(6,'(a)') ' writing neb'//cnum//'...'
        call write_atom_conf(20,'neb'//cnum,natm,h,hunit
     &       ,raini(1:3,1:natm)+rd(1:3,1:natm)*(i-1),va,tag)
c.....Write akr file
        open(21,file='akr'//cnum,status='replace')
        write(21,'(es15.7)') hunit
        write(21,'(3es11.3)') ((h(ia,ib,0)/hunit,ia=1,3),ib=1,3)
        write(21,'(i10,3i5)') natm, 3, 0, 0
        do j=1,natm
          write(21,'(i3,6es11.3)') int(tag(j))
     &         ,raini(1:3,j)+rd(1:3,j)*(i-1)
     &         ,0d0,0d0,0d0
        enddo
        close(21)
      enddo

      end program make_slices
c=======================================================================
      subroutine read_atom_conf(ionum,cfname,namax,natm,h,hunit,ra,tag)
      implicit none
      integer,intent(in):: namax,ionum
      character(len=*),intent(in):: cfname
      real(8),intent(out):: ra(3,namax),tag(namax),h(3,3,0:1),hunit
      integer,intent(out):: natm

      integer:: ia,ib,l,i

      open(ionum,file=trim(cfname),status='old')
      read(ionum,*) hunit
      read(ionum,*) (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
      h(1:3,1:3,0:1)= h(1:3,1:3,0:1) *hunit
      read(ionum,*) natm
      do i=1,natm
        read(ionum,*) tag(i),ra(1:3,i)
      enddo
      close(ionum)
      
      end subroutine read_atom_conf
c=======================================================================
      subroutine write_atom_conf(ionum,cfname,natm,h,hunit,ra,va,tag)
      implicit none
      integer,intent(in):: natm,ionum
      character(len=*),intent(in):: cfname
      real(8),intent(in):: h(3,3,0:1),ra(3,natm),va(3,natm),tag(natm)
     &     ,hunit

      integer:: ia,ib,l,i
      real(8):: ri(3)

      open(ionum,file=trim(cfname),status='replace')
      write(ionum,'(es23.14e3)') hunit
      write(ionum,'(3es23.14e3)') (((h(ia,ib,l)/hunit,ia=1,3)
     &     ,ib=1,3),l=0,1)
      write(ionum,'(i10)') natm
      do i=1,natm
        ri(1:3)= ra(1:3,i)
        if( ri(1).lt.0d0 ) ri(1)=ri(1) +1d0
        if( ri(2).lt.0d0 ) ri(2)=ri(2) +1d0
        if( ri(3).lt.0d0 ) ri(3)=ri(3) +1d0
        if( ri(1).ge.1d0 ) ri(1)=ri(1) -1d0
        if( ri(2).ge.1d0 ) ri(2)=ri(2) -1d0
        if( ri(3).ge.1d0 ) ri(3)=ri(3) -1d0
        write(ionum,'(7es23.14e3,11es12.4)') tag(i),ri(1:3)
     &       ,va(1:3,i)       ! va(1:3,i)
     &       ,0d0,0d0 ! eki(i),epi(i)
     &       ,0d0,0d0,0d0 ,0d0,0d0,0d0 ,0d0,0d0,0d0 ! strs(1:3,1:3,i)
      enddo
      close(ionum)
      
      end subroutine write_atom_conf
c=======================================================================
      subroutine setv(natm,va,tag,tempini)
      use variables, only: fkb,amu
      implicit none
      integer,intent(in):: natm
      real(8),intent(in):: tempini,tag(natm)
      real(8),intent(out):: va(3,natm)

      integer:: i,l,is
      real(8):: dseed,sumvx,sumvy,sumvz,rnd1,rnd2,tmp,facv,am
      real(8),parameter:: pi = 3.14159265358979d0
      real(8),parameter:: am_fe = 55.847d0 *amu

      facv=dsqrt(2d0*tempini*fkb/am_fe)
      am= am_fe

c-----velocities in Maxwell-Boltzmann distribution
      dseed=12345
      do i=1,natm
        is= int(tag(i))
        do l=1,3
          call myrnd(rnd1,dseed)
          call myrnd(rnd2,dseed)
          va(l,i)=facv*dsqrt(-dlog(rnd1))*dcos(2d0*pi*rnd2)
        enddo
      enddo
c-----set center of mass motion to zero
      sumvx=0d0
      sumvy=0d0
      sumvz=0d0
      do i=1,natm
        sumvx=sumvx+va(1,i)
        sumvy=sumvy+va(2,i)
        sumvz=sumvz+va(3,i)
      enddo
      do i=1,natm
        va(1,i)=va(1,i)-sumvx/dble(natm)
        va(2,i)=va(2,i)-sumvy/dble(natm)
        va(3,i)=va(3,i)-sumvz/dble(natm)
      enddo

      tmp=0d0
      do i=1,natm
        is= int(tag(i))
        tmp= tmp +0.5d0*am*(va(1,i)**2 +va(2,i)**2 +va(3,i)**2)
      enddo
      write(6,'(a,es12.4)') " ekin=",tmp
      write(6,'(a,es12.4)') " temp.=",tmp*2d0/3d0/fkb/natm

      end subroutine setv
c=======================================================================
      subroutine myrnd(rnd,dseed)
      real*8 rnd,dseed
      real*8 d2p31m,d2p31
      save d2p31m,d2p31
      data d2p31m/2147483647d0/
      data d2p31 /2147483648d0/
      
      dseed=dmod(16807d0*dseed,d2p31m)
      rnd=dseed/d2p31
      return
      end subroutine myrnd
c=======================================================================
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make 10mkslice"
c     End:
