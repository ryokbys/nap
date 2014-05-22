      program convert_slices
c-----------------------------------------------------------------------
c  USAGE:
c    $ ./20conv <#-of-atoms> <#-old> <#-new>
c    <#-of-atoms>: number of atoms
c    <#-old>, <#-new>: old and new numbers of slices
c-----------------------------------------------------------------------
      implicit none
      integer:: i,j,ia,ib,islc,nslc,natm,inew,iold,nnew,nold
      real(8):: rd(3),h(3,3,0:1),pnew,pold,dnew,dold
      character:: cnum*3,argv*128,cslc*3
      integer:: nargc
      integer,external:: iargc
      real(8),allocatable:: ra0(:,:,:),va0(:,:,:),tag0(:,:)
     &     ,ra(:,:),va(:,:),tag(:)

c.....Get num of slices which includes initial and final slices
      nargc= iargc()
      if( nargc.ne.3 ) then
        stop ' [Usage] ./20conv 31104  2 9'
      endif
      call getarg(1,argv)
      read(argv,*) natm
      call getarg(2,argv)
      read(argv,*) nold
      allocate(ra0(3,natm,nold),va0(3,natm,nold),tag0(natm,nold)
     &     ,ra(3,natm),va(3,natm),tag(natm))
      call getarg(3,argv)
      read(argv,*) nnew
      if( nold.eq.1 .or. nnew.eq.1 ) then
        write(6,'(a)') ' [Error] nold & nnew must be larger than 1 !!!'
        stop
      endif

c.....Get slices of (from) configuration
      do islc=1,nold
        cslc='000'
        write(cslc,'(i3.3)') islc
        call read_atom_conf(10+islc,'neb'//cslc,natm,h,ra0(1,1,islc)
     &       ,va0(1,1,islc),tag0(1,islc))
        call system('mv neb'//cslc//' neb'//cslc//'.old')
      enddo

      va(1:3,1:natm)= 0d0
      tag(1:natm)= tag0(1:natm,1)

c.....Make slices by interpolating from old slices
      dold= 1d0/(nold-1)
      dnew= 1d0/(nnew-1)
      do inew=1,nnew
        pnew= (inew-1)*dnew
        do iold=2,nold
          pold= (iold-1)*dold
          if( pold.gt.pnew ) exit
        enddo
        if( iold.gt.nold ) iold= nold
        do i=1,natm
          rd(1:3)= ra0(1:3,i,iold)-ra0(1:3,i,iold-1)
          rd(1:3)= (rd(1:3)-anint(rd(1:3)))
     &         *(pnew-(pold-dold))/dold
          ra(1:3,i)=ra0(1:3,i,iold-1)+rd(1:3)
          call pbc(ra(1,i))
          call pbc(ra(2,i))
          call pbc(ra(3,i))
        enddo
        write(cnum(1:3),'(i3.3)') inew
        write(6,'(a)') ' writing neb'//cnum//'...'
        call write_atom_conf(20,'neb'//cnum,natm,h,ra,va,tag)
c.....Write akr file
        open(21,file='akr'//cnum,status='replace')
        write(21,'(i10,3i5)') natm, 3, 0, 0
        write(21,'(3es11.3)') ((h(ia,ib,0),ia=1,3),ib=1,3)
        do j=1,natm
          write(21,'(i3,6es11.3)') int(tag(j)),ra(1:3,j),va(1:3,j)
        enddo
        close(21)
      enddo


      end program convert_slices
c=======================================================================
      subroutine read_atom_conf(ionum,cfname,natm,h,ra,va,tag)
      implicit none
      integer,intent(in):: natm,ionum
      character(len=*),intent(in):: cfname
      real(8),intent(out):: ra(3,natm),va(3,natm),tag(natm),h(3,3,0:1)

      integer:: ia,ib,l,i,itmp

      open(ionum,file=trim(cfname),status='old')
      read(ionum,*) itmp
      if( itmp.ne.natm ) stop ' [Error] itmp.ne.natm !!!'
      read(ionum,*) (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
      do i=1,natm
        read(ionum,*) tag(i),ra(1:3,i),va(1:3,i)
      enddo
      close(ionum)
      
      end subroutine read_atom_conf
c=======================================================================
      subroutine write_atom_conf(ionum,cfname,natm,h,ra,va,tag)
      implicit none
      integer,intent(in):: natm,ionum
      character(len=*),intent(in):: cfname
      real(8),intent(in):: h(3,3,0:1),ra(3,natm),va(3,natm),tag(natm)

      integer:: ia,ib,l,i
      real(8):: ri(3)

      open(ionum,file=trim(cfname),status='replace')
      write(ionum,'(i10)') natm
      write(ionum,'(3es23.14e3)') (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
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
      subroutine pbc(x)
      implicit none
      real(8),intent(inout):: x

      if( x.lt.0d0 ) x=x+1d0
      if( x.ge.1d0 ) x=x-1d0
      end subroutine pbc
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make 20conv"
c     End:
