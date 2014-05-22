      program slice2vnorm
c
c make normal vectors from NEB slices
c
c usage:
c   $ ./slc2vnrm
c   $ for num in `seq 1 ??`; do cnum=`print %03d $num`;\\
c      cp vnrm$cnum $cnum/qmclvnrm; done
c
c   - ??: number of slices in NEB
c
      use variables
      implicit none
      integer:: islc,ism,isp,i,ionum
      real(8),external:: sprod
      real(8),allocatable:: vn(:,:)
      real(8):: tmp
      character(len=3):: cslc,cnum 
      
      call read_input(10,'in.neb')
      allocate(ra(3,natm,nslc),va(3,natm,nslc),fa(3,natm,nslc)
     &     ,tag(natm,nslc),vn(3,natm))
c.....eV/A^2 ---> Hartree/Bohr^2
      scnst=scnst *ev2hrt /aa2bohr**2
c.....eV/A ---> Hartree/Bohr
      feps= feps *ev2hrt /aa2bohr
c.....A ---> Bohr
      deps= deps *aa2bohr

c.....Read atom configuration of each slice
      do islc=1,nslc
        ionum= 10+islc
        write(cslc(1:3),'(i3.3)') islc
        call read_atom_conf(ionum,'neb'//cslc,natm
     &       ,h,ra(1,1,islc),va(1,1,islc),tag(1,islc))
      enddo

      do islc=1,nslc
        ism=islc-1
        isp=islc+1
        if( islc.eq.1 ) then
          ism=1
          isp=islc+1
        elseif( islc.eq.nslc ) then
          ism=islc-1
          isp=islc
        endif
c.....tangential vector
        vn(1:3,1:natm)= ra(1:3,1:natm,isp)-ra(1:3,1:natm,ism)
c.....correct too big vectors
        do i=1,natm
          vn(1:3,i)= vn(1:3,i) -anint(vn(1:3,i))
        enddo
        do i=1,natm
          vn(1,i)= vn(1,i)*h(1,1,0)
          vn(2,i)= vn(2,i)*h(2,2,0)
          vn(3,i)= vn(3,i)*h(3,3,0)
        enddo
c.....normalize
        tmp= sprod(3*natm,vn,vn)
        vn(1:3,1:natm)= vn(1:3,1:natm)/tmp
c.....output
        write(cnum(1:3),'(i3.3)') islc
        open(10,file='vnrm'//cnum,status='replace')
        write(10,'(i10)') natm
        do i=1,natm
          write(10,'(3es22.14)') vn(1:3,i)
        enddo
        close(10)
        write(6,'(a)') ' vnrm'//cnum//' is created'
      enddo

      end program slice2vnorm
c=======================================================================
      function sprod(ndim,v1,v2)
      implicit none
      integer,intent(in):: ndim
      real(8),intent(in):: v1(ndim),v2(ndim)
      real(8):: sprod

      integer:: i

      sprod= 0d0
      do i=1,ndim
        sprod= sprod +v1(i)*v2(i)
      enddo
      return
      end function sprod
c=======================================================================
      subroutine read_atom_conf(ionum,cfname,natm,h,ra,va,tag)
      implicit none
      integer,intent(in):: ionum,natm
      character(len=*),intent(in):: cfname
      real(8),intent(out):: ra(3,natm),tag(natm),h(3,3,0:1),va(3,natm)

      integer:: ia,ib,l,i,itmp

      open(ionum,file=trim(cfname),status='old')
      read(ionum,*) itmp
      if( itmp.ne.natm) then
        write(6,'(a)') ' [Error] itmp.ne.natm !!!'
        write(6,'(a,2i6)') '  itmp,natm=',itmp,natm
        stop
      endif
      read(ionum,*) (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
      do i=1,natm
        read(ionum,*) tag(i),ra(1:3,i),va(1:3,i)
      enddo
      close(ionum)
      
      end subroutine read_atom_conf
c=======================================================================
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make slc2vnrm"
c     End:

