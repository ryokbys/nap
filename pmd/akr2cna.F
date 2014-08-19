      program akr2cna
c-----------------------------------------------------------------------
c Reduce atoms from akr file using Common Neighbor Analysis (CNA).
c   - assuming akr file to be redirect to this program as STDIN
c   - output akr format to STDOUT
c   - FCC crystal only?
c   - change cutoff radius, rc, to apply this to your system
c-----------------------------------------------------------------------
      implicit none
c-----max num of neighbors; 12 is enough
      integer,parameter:: nnmax = 100
cc.....Cutoff radius for LJ-Ar FCC, halfway between 1st and 2nd N.N.
cc      real(8),parameter:: rc    = 3.36d0
c      real(8),parameter:: rc    = 3.5d0
cc.....Cutoff radius for BCC Fe, halfway between 2nd and 3rd N.N.
c      real(8),parameter:: abcc = 2.835d0
c.....Cutoff radius for BCC W, halfway between 2nd and 3rd N.N.
cc.....W @1000K = 3.226d-10
c      real(8),parameter:: abcc = 3.226d0
c.....W @0K = 3.171d-10
      real(8),parameter:: abcc = 3.204d0
      real(8),parameter:: rc   = 0.5d0*(1d0+1.41421356d0)*abcc

      integer:: i,n,itmp,nauxdat,j,ix,iy,iz,m,ia,ib
      integer:: natm0,natm1,natm,nx,ny,nz,nlsa,ncpx,ncpy,ncpz
      real(8):: tmp,al(3),al1(3),xi(3),rij(3),r,hunit,h(3,3)
      integer,allocatable:: lsa(:),lspr(:,:),itype(:),itag(:),itag1(:)
      real(8),allocatable:: ra0(:,:),ra1(:,:),ra(:,:),auxdat0(:,:)
     &     ,auxdat(:,:)
      integer nargc
c      integer,external:: iargc
      character*128:: ciname,coname

c.....Read an argument, ipmd
c      nargc=iargc()
      nargc= command_argument_count()
      if( nargc.ne.2 ) then
        stop 'Usage: $ ./akr2cna akr0000 cna0000'
      endif
      call getarg(1,ciname)
      call getarg(2,coname)

c-----read akr file from STDIN
      open(10,file=trim(ciname),status="old")
      read(10,*) hunit
      read(10,*) ((h(ia,ib),ia=1,3),ib=1,3)
c      read(10,*) al(1),tmp,tmp
c      read(10,*) tmp,al(2),tmp
c      read(10,*) tmp,tmp,al(3)
c      al(1:3)= al(1:3) *hunit
      h(1:3,1:3)= h(1:3,1:3)*hunit
      al(1)= h(1,1)
      al(2)= h(2,2)
      al(3)= h(3,3)
      ncpx=1
      ncpy=1
      ncpz=1
      if( al(1).lt.3*rc ) then
        ncpx= int(3*rc/al(1)+1)
      endif
      if( al(2).lt.3*rc ) then
        ncpy= int(3*rc/al(2)+1)
      endif
      if( al(3).lt.3*rc ) then
        ncpz= int(3*rc/al(3)+1)
      endif
      al1(1)= al(1)*ncpx
      al1(2)= al(2)*ncpy
      al1(3)= al(3)*ncpz
      
c      write(6,'(a,7es12.4)') 'rc,al(1:3),al1(1:3)=',rc,al(1:3),al1(1:3)
c      write(6,'(a,3i6)') 'ncpx,ncpy,ncpz=',ncpx,ncpy,ncpz
      read(10,*) natm0,nauxdat
      natm1= natm0*ncpx*ncpy*ncpz
      allocate(ra0(3,natm0),ra(3,natm0),auxdat(nauxdat,natm0)
     &     ,auxdat0(nauxdat,natm0),itype(natm0),itag(natm0)
     &     ,ra1(3,natm1),itag1(natm1))
      do i=1,natm0
        read(10,*) itag(i),ra0(1:3,i),auxdat0(1:nauxdat,i)
      enddo
      close(10)
      m=0
      do ix=0,ncpx-1
        do iy=0,ncpy-1
          do iz=0,ncpz-1
            do i=1,natm0
              m=m+1
              ra1(1,m)= (ra0(1,i) +ix)/ncpx
              ra1(2,m)= (ra0(2,i) +iy)/ncpy
              ra1(3,m)= (ra0(3,i) +iz)/ncpz
              itag1(m)= itag(i)
            enddo
          enddo
        enddo
      enddo
      

c-----prepare for linked-list
      nx=max(1,int(al1(1)/rc))
      ny=max(1,int(al1(2)/rc))
      nz=max(1,int(al1(3)/rc))
      nlsa= natm1 +nx*ny*nz
c      write(6,'(a,4i10)') " nx,ny,nz,nlsa=",nx,ny,nz,nlsa
      allocate(lsa(nlsa))
c-----make linked-list
      call mkllist(natm1,nx,ny,nz,ra1,nlsa,lsa)

c-----make 1st nearest neighbor pair-list
      allocate(lspr(0:nnmax,natm1))
      call ll2lspr(nlsa,lsa,nx,ny,nz,natm1,nnmax,lspr,al1,ra1,rc,itag1)
      deallocate(lsa)

c-----common neighbor analysis
      call cna(natm0,natm1,nnmax,lspr,itag1,ra0,auxdat0,nauxdat
     &     ,natm,ra,auxdat,itype)

c-----write cna file 
      open(11,file=trim(coname),status="replace")
      write(11,'(es15.7)') hunit
      write(11,'(3es12.4)') ((h(ia,ib)/hunit,ia=1,3),ib=1,3)
c      write(11,'(3es12.4)') al(1)/hunit,0d0,0d0
c      write(11,'(3es12.4)') 0d0,al(2)/hunit,0d0
c      write(11,'(3es12.4)') 0d0,0d0,al(3)/hunit
      write(11,'(i10,3i4)') natm,nauxdat,0,0
      do i=1,natm
        write(11,'(i4,3es15.7,17es10.2)')
     &       itype(i),ra(1:3,i),auxdat(1:nauxdat,i)
      enddo
      close(11)

      end program akr2cna
c=======================================================================
      subroutine mkllist(natm,n1,n2,n3,ra,nlist,llist)
c-----------------------------------------------------------------------
c     make linked-list
c-----------------------------------------------------------------------
      implicit none 
      integer,intent(in):: natm,n1,n2,n3,nlist
      real(8),intent(inout):: ra(3,natm)
      integer,intent(out):: llist(nlist)
      integer:: i,j,n,nr1,nr2,nr3
      real(8):: a1,a2,a3
c.....To avoid nr?==n?, ra(?,i) has to be under 1.0d0
      real(8),parameter:: small=1d0-1d-8

      llist(1:nlist)=0

      a1= 1d0/n1
      a2= 1d0/n2
      a3= 1d0/n3
      do i=1,natm
        if( ra(1,i).ge.1d0 ) ra(1,i)=ra(1,i) -1d0
        if( ra(1,i).lt.0d0 ) ra(1,i)=ra(1,i) +1d0
        if( ra(2,i).ge.1d0 ) ra(2,i)=ra(2,i) -1d0
        if( ra(2,i).lt.0d0 ) ra(2,i)=ra(2,i) +1d0
        if( ra(3,i).ge.1d0 ) ra(3,i)=ra(3,i) -1d0
        if( ra(3,i).lt.0d0 ) ra(3,i)=ra(3,i) +1d0
c-------region index in order z,y,x
        nr3= int(ra(3,i)*small/a3)+1
        nr2= int(ra(2,i)*small/a2)
        nr1= int(ra(1,i)*small/a1)
        n=nr1*n2*n3 +nr2*n3 +nr3
        j=llist(natm+n)
        llist(i)=j
        llist(natm+n)=i
      enddo
      return
      end subroutine mkllist
c=======================================================================
      subroutine ll2lspr(nlls,lls,nx,ny,nz,natm,nnmax,lspr,al,ra,rc
     &     ,itag)
c-----------------------------------------------------------------------
c     make Verlet list from linked-list
c-----------------------------------------------------------------------
      implicit none 
      integer,intent(in):: nlls,lls(nlls),nx,ny,nz,natm,nnmax,itag(natm)
      real(8),intent(in):: al(3),ra(3,natm),rc
      integer,intent(out):: lspr(0:nnmax,natm)
      integer:: i,j,k,l,m,n,ir,jr,kx,ky,kz,irx,iry,irz,jrx,jry,jrz,ii
     &     ,is,js
      real(8):: dx,dy,dz,r

      real(8),allocatable:: da(:)

      allocate(da(nnmax))

      lspr(0:nnmax,1:natm)=0
      
c-----region (ir)
      do ir=1,nx*ny*nz
        irz=mod(ir-1,nz) +1
        iry=mod((ir-1)/nz,ny) +1
        irx=(ir-1)/(nz*ny) +1
c-------atoms (i) in rgn (ir)
        i=lls(natm+ir)
        do
          n=0
          if(i.eq.0) exit
          is= itag(i)
c---------region (jr)
          do kx=-1,1
            jrx= irx +kx
            if(jrx.gt.nx) jrx=jrx -nx
            if(jrx.le.0 ) jrx=jrx +nx
            do ky=-1,1
              jry= iry +ky
              if(jry.gt.ny) jry=jry -ny
              if(jry.le.0 ) jry=jry +ny
              do kz=-1,1
                jrz= irz +kz
                if(jrz.gt.nz) jrz=jrz -nz
                if(jrz.le.0 ) jrz=jrz +nz
                jr= (jrx-1)*nz*ny +(jry-1)*nz +jrz
                j=lls(natm+jr)
                do
                  if(j.eq.0) exit
                  js= itag(j)
                  dx=(ra(1,j)-ra(1,i)-anint(ra(1,j)-ra(1,i)))*al(1)
                  dy=(ra(2,j)-ra(2,i)-anint(ra(2,j)-ra(2,i)))*al(2)
                  dz=(ra(3,j)-ra(3,i)-anint(ra(3,j)-ra(3,i)))*al(3)
                  r= dsqrt(dx*dx +dy*dy +dz*dz)
                  if(r.lt.rc .and. j.ne.i .and. is.eq.js) then
c                  if(r.lt.rc .and. j.ne.i) then
c                    if(r.lt.1d-5) then
c                      write(6,'(a,3es12.4)') ' al(1:3)=',al(1:3)
c                      write(6,'(a,i8,3es12.4)') 'i,ra(1:3,i)='
c     &                     ,i,ra(1:3,i)
c                      write(6,'(a,i8,3es12.4)') 'j,ra(1:3,j)='
c     &                     ,j,ra(1:3,j)
c                      write(6,'(a,es12.4)') ' r=',r
c                      stop ' r.lt.1d-5 !!!'
c                    endif
                    n=n+1
                    if(n.gt.nnmax) then
                      write(6,'(a)') "mklspr: n.gt.nnmax"
                      write(6,'(a,2i10)') " n,nnmax=",n,nnmax
                      write(6,'(a,2i8,2es11.3)') "i,j,r,rc=",i,j,r,rc
                      stop "mklspr: n.gt.nnmax"
                    endif
                    da(n)= r
                    lspr(0,i)= lspr(0,i) +1
                    lspr(n,i)=j
                  endif
c-----------------next (j) in region (jr)
                  j=lls(j)
                enddo
              enddo
            enddo
          enddo
cc.....Sort
c          call heapsort_i(lspr(0,i),nnmax,da,lspr(1:nnmax,i))
c          if( i.eq.1 ) then
c            do l=1,lspr(0,i)
c              j=lspr(l,i)
c              write(6,'(a,2i10,es12.4)') " l,j,da=",l,j,da(l)
c            enddo
c          endif
c---------next (i) in region (ir)
          i=lls(i)
        enddo
      enddo

c      stop


c      write(6,'(a)') " checking lspr..."
c      n=0
c      do i=1,natm
cc        do l=1,lspr(0,i)
cc          j=lspr(l,i)
cc          write(6,'(2i10,3es12.4)') i,j,(ra(1:3,j)-ra(1:3,i))*al(1:3)
cc        enddo
c        if( lspr(0,i).gt.14 ) then
c          n=n+1
c          write(6,'(a,2i10)') "i,lspr=",i,lspr(0,i)
c        endif
c      enddo
c      write(6,'(a,i10)') " num of non 12..=",n
c      stop

      end subroutine ll2lspr
c=======================================================================
      subroutine cna(natm0,natm1,nnmax,lspr,itag,ra0,auxdat0,nauxdat
     &     ,natm,ra,auxdat,itype)
c-----------------------------------------------------------------------
c Reduce natm0 to natm by using Common Neighbor Analysis
c-----------------------------------------------------------------------
      implicit none
      integer,intent(in):: natm0,natm1,nnmax,lspr(0:nnmax,natm1),nauxdat
     &     ,itag(natm1)
      real(8),intent(in):: ra0(3,natm0),auxdat0(nauxdat,natm0)
      integer,intent(out):: itype(natm0),natm
      real(8),intent(out):: ra(3,natm0),auxdat(nauxdat,natm0)

      integer:: i,j,l,m,n,ii,iii,ni,jj,nj,il,jl,n1,n2,iil,nn1,im,iim
     &     ,ib1,ib2,iib1,iib2,n421,n422,n663,n443
      integer,allocatable:: icommon(:),ibond(:,:),nb(:),idc(:,:,:)
cc.....For FCC
c      integer,parameter:: lmax= 12
c      integer,parameter:: mmax= 12
c.....For BCC
c      integer,parameter:: lmax= 14
      integer,parameter:: lmax= 30
      integer,parameter:: mmax= lmax*(lmax-1)/2

      allocate(icommon(lmax),ibond(2,mmax),nb(mmax)
     &     ,idc(3,nnmax,natm1))

c-----init three indices
      idc(1:3,1:nnmax,1:natm1)= 0

c-----for each atom-i, store three indices (LMN)
      do i=1,natm1
        if( itag(i).ne.1 ) cycle

c-------for each 1st n.n.
        do ii=1,lspr(0,i)
          j=lspr(ii,i)
c---------j>i only
          if(j.le.i) cycle

c---------count num of common neighbors: L
          l= 0
          icommon(1:lmax)= 0
          do iii=1,lspr(0,i)
            ni=lspr(iii,i)
            if(ni.eq.j) cycle
            do jj=1,lspr(0,j)
              nj=lspr(jj,j)
              if(nj.eq.ni) then
                l=l+1
                if(l.gt.lmax) then
                  write(6,'(a,2i6)') ' i,lspr(0,i)=',i,lspr(0,i)
                  write(6,'(a,2i6)') ' j,lspr(0,j)=',j,lspr(0,j)
                  stop " l.gt.lmax!!!"
                endif
                icommon(l)= ni
                exit
              endif
            enddo
c---------end of counting L
          enddo
          
c---------count num of bonds between common neighbors: M
          m= 0
          ibond(1:2,1:mmax)= 0
c---------for each common neighbor-n1
          do il=1,l
            n1=icommon(il)
c-----------for common neighbor-n2 which must be larger than n1
            do jl=1,l
              n2=icommon(jl)
              if(n2.le.n1) cycle
c-------------scan 1st n.n. of n1
              do iil=1,lspr(0,n1)
                nn1=lspr(iil,n1)
                if(nn1.eq.n2) then
                  m=m+1
                  if(m.gt.mmax) stop " m.gt.mmax!!"
                  ibond(1:2,m)= (/ n1,n2 /)
                  exit
                endif
              enddo
            enddo
          enddo

c---------count max num of continuous bonds: N
          nb(1:mmax)= 1
c---------this does not distinguish star and chain connections
          do im=1,m-1
            ib1= ibond(1,im)
            ib2= ibond(2,im)
            do iim=im+1,m
              iib1= ibond(1,iim)
              iib2= ibond(2,iim)
c-------------if two bonds are connected, up nb
              if(iib1.eq.ib1 .or. iib2.eq.ib1
     &             .or. iib1.eq.ib2 .or. iib2.eq.ib2) then
                nb(im)=nb(im) +1
                nb(iim)=nb(iim) +1
              endif
            enddo
          enddo
c---------maximum nb
          n= 0
          do im=1,m
            n= max(nb(im),n)
          enddo

c---------store (LMN) to i
          idc(1:3,ii,i)= (/ l,m,n /)
c---------store (LMN) to j, too
          do jj=1,lspr(0,j)
            if(lspr(jj,j).eq.i) then
              idc(1:3,jj,j)= (/ l,m,n /)
              exit
            endif
          enddo
c-------end of 1st n.n. of i: j
        enddo
c-----end of atom-i
      enddo

c-----reduce atoms
      natm= 0
      do i=1,natm0
        if( itag(i).eq.2 ) then
          natm=natm+1
          itype(natm)=itag(i)
          goto 10
        endif
        n421= 0
        n422= 0
        n663= 0
        n443= 0
        do ii=1,lspr(0,i)
          l=idc(1,ii,i)
          m=idc(2,ii,i)
          n=idc(3,ii,i)
c          if( i.eq.41 ) then
c            write(6,'(a,4i10)') '  ii,l,m,n=',ii,l,m,n
c          endif
          if(l.eq.4 .and. m.eq.2 .and. n.eq.1 ) n421=n421 +1
          if(l.eq.4 .and. m.eq.2 .and. n.eq.2 ) n422=n422 +1
          if(l.eq.6 .and. m.eq.6 .and. n.eq.3 ) n663=n663 +1
          if(l.eq.4 .and. m.eq.4 .and. n.eq.3 ) n443=n443 +1
        enddo
c.....If perfect FCC structure, reduce the atom-i
        if(n421.eq.12 .and. n422.eq.0) cycle
c.....If perfect BCC structure, reduce the atom-i
        if( n663.eq.8 .and. n443.eq.6 ) cycle
        natm=natm +1
c.....If HCP structure, itype(i)= 2
        if(n421.eq.6 .and. n422.eq.6) then
          itype(natm)= 2
c.....Otherwise, itype(i)= 1
        else
          itype(natm)= 1
        endif

        itype(natm)= 1

c-------copy to ra,va
 10     ra(1:3,natm)= ra0(1:3,i)
        auxdat(1:nauxdat,natm)= auxdat0(1:nauxdat,i)
c        if( natm.eq.2444 ) then
c          do ii=1,lspr(0,i)
c            j=lspr(ii,i)
c            l=idc(1,ii,i)
c            m=idc(2,ii,i)
c            n=idc(3,ii,i)
c          enddo
c        endif
      enddo

      deallocate(icommon,ibond,nb,idc)
      end subroutine cna
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make akr2cna"
c     End:
