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
c-----cutoff radius
c      real(8),parameter:: rc    = 3.36d-10 /0.5291772d-10
      real(8),parameter:: rc    = 3.5d-10 /0.5291772d-10

      integer:: i,n,itmp,nauxdat,j
      integer:: natm0,natm,nx,ny,nz,nlsa
      real(8):: tmp,al(3),xi(3),rij(3),r
      integer,allocatable:: lsa(:),lspr(:,:),itype(:)
      real(8),allocatable:: ra0(:,:),ra(:,:),auxdat0(:,:)
     &     ,auxdat(:,:)

c-----read akr file from STDIN
      read(5,*) natm0,nauxdat
c      write(6,'(a,i10)') " natm =",natm0
      read(5,*) al(1),tmp,tmp
      read(5,*) tmp,al(2),tmp
      read(5,*) tmp,tmp,al(3)
c      write(6,'(a,3es12.4)') " al(1:3) =",al(1:3)
      allocate(ra0(3,natm0),ra(3,natm0),auxdat(nauxdat,natm0)
     &     ,auxdat0(nauxdat,natm0),itype(natm0))
      do i=1,natm0
        read(5,*) itmp,ra0(1:3,i),auxdat0(1:nauxdat,i)
      enddo

c-----prepare for linked-list
      nx=int(al(1)/(rc*1.05d0))
      ny=int(al(2)/(rc*1.05d0))
      nz=int(al(3)/(rc*1.05d0))
      nlsa= natm0 +nx*ny*nz
c      write(6,'(a,4i10)') " nx,ny,nz,nlsa=",nx,ny,nz,nlsa
      allocate(lsa(nlsa))
c-----make linked-list
      call mkllist(natm0,nx,ny,nz,ra0,nlsa,lsa)

c-----make 1st nearest neighbor pair-list
      allocate(lspr(0:nnmax,natm0))
      call ll2lspr(nlsa,lsa,nx,ny,nz,natm0,nnmax,lspr,al,ra0,rc)
      deallocate(lsa)
      
c-----common neighbor analysis
      call cna(natm0,nnmax,lspr,ra0,auxdat0,nauxdat
     &     ,natm,ra,auxdat,itype)

c-----write akr file to STDOUT
      write(6,'(i10,3i4)') natm,nauxdat,0,0
      write(6,'(3es12.4)') al(1),0d0,0d0
      write(6,'(3es12.4)') 0d0,al(2),0d0
      write(6,'(3es12.4)') 0d0,0d0,al(3)
      do i=1,natm
        write(6,'(i4,20es11.3)') itype(i),ra(1:3,i),auxdat(1:nauxdat,i)
      enddo

      end program akr2cna
c=======================================================================
      subroutine mkllist(natm,n1,n2,n3,ra,nlist,llist)
c-----------------------------------------------------------------------
c     make linked-list
c-----------------------------------------------------------------------
      implicit none 
      integer,intent(in):: natm,n1,n2,n3,nlist
      real(8),intent(in):: ra(3,natm)
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
      subroutine ll2lspr(nlls,lls,nx,ny,nz,natm,nnmax,lspr,al,ra,rc)
c-----------------------------------------------------------------------
c     make Verlet list from linked-list
c-----------------------------------------------------------------------
      implicit none 
      integer,intent(in):: nlls,lls(nlls),nx,ny,nz,natm,nnmax
      real(8),intent(in):: al(3),ra(3,natm),rc
      integer,intent(out):: lspr(0:nnmax,natm)
      integer:: i,j,k,l,m,n,ir,jr,kx,ky,kz,irx,iry,irz,jrx,jry,jrz,ii
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
                  dx=(ra(1,j)-ra(1,i)-anint(ra(1,j)-ra(1,i)))*al(1)
                  dy=(ra(2,j)-ra(2,i)-anint(ra(2,j)-ra(2,i)))*al(2)
                  dz=(ra(3,j)-ra(3,i)-anint(ra(3,j)-ra(3,i)))*al(3)
                  r= dsqrt(dx*dx +dy*dy +dz*dz)
c                  if(i.eq.1) write(6,'(a,4i4,2es11.3)')
c     &                 "ir,i,jr,j,r,rc=",ir,i,jr,j,r,rc*1.3d0
                  if(r.lt.rc .and. j.ne.i) then
c                    write(6,'(a,2i5,es11.3)') "i,j,r-rc*1.2=",i,j
c     &                   ,r-rc*1.2d0
                    n=n+1
                    if(n.gt.nnmax) then
                      write(6,'(a)') "mklspr: n.gt.nnmax"
                      write(6,'(a,2i10)') " n,nnmax=",n,nnmax
                      write(6,'(a,2i5,2es11.3)') "i,j,r,rc=",i,j,r,rc
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
c        if( lspr(0,i).gt.12 ) then
c          n=n+1
cc          write(6,'(a,2i10)') "i,lspr=",i,lspr(0,i)
c        endif
c      enddo
c      write(6,'(a,i10)') " num of non 12..=",n
c      stop

      end subroutine ll2lspr
c=======================================================================
      subroutine cna(natm0,nnmax,lspr,ra0,auxdat0,nauxdat
     &     ,natm,ra,auxdat,itype)
c-----------------------------------------------------------------------
c Reduce natm0 to natm by using Common Neighbor Analysis
c-----------------------------------------------------------------------
      implicit none
      integer,intent(in):: natm0,nnmax,lspr(0:nnmax,natm0),nauxdat
      real(8),intent(in):: ra0(3,natm0),auxdat0(nauxdat,natm0)
      integer,intent(out):: itype(natm0),natm
      real(8),intent(out):: ra(3,natm0),auxdat(nauxdat,natm0)

      integer:: i,j,l,m,n,ii,iii,ni,jj,nj,il,jl,n1,n2,iil,nn1,im,iim
     &     ,ib1,ib2,iib1,iib2,n421,n422
      integer,allocatable:: icommon(:),ibond(:,:),nb(:),idc(:,:,:)
      integer,parameter:: lmax= 12
      integer,parameter:: mmax= 12

      allocate(icommon(lmax),ibond(2,mmax),nb(mmax)
     &     ,idc(3,nnmax,natm0))

c-----init three indices
      idc(1:3,1:nnmax,1:natm0)= 0

c-----for each atom-i, store three indices (LMN)
      do i=1,natm0

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
                if(l.gt.lmax) stop " l.gt.lmax!!!"
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
        n421= 0
        n422= 0
        do ii=1,lspr(0,i)
          l=idc(1,ii,i)
          m=idc(2,ii,i)
          n=idc(3,ii,i)
          if(l.eq.4 .and. m.eq.2 .and. n.eq.1 ) n421=n421 +1
          if(l.eq.4 .and. m.eq.2 .and. n.eq.2 ) n422=n422 +1
        enddo
c-------if perfect FCC structure, reduce the atom-i
        if(n421.eq.12 .and. n422.eq.0) cycle
        natm=natm +1
c-------if HCP structure, itype(i)= 2
        if(n421.eq.6 .and. n422.eq.6) then
          itype(natm)= 2
c-------otherwise, itype(i)= 1
        else
          itype(natm)= 1
        endif
c-------copy to ra,va
        ra(1:3,natm)= ra0(1:3,i)
        auxdat(1:nauxdat,natm)= auxdat0(1:nauxdat,i)
        if( natm.eq.2444 ) then
          do ii=1,lspr(0,i)
            j=lspr(ii,i)
            l=idc(1,ii,i)
            m=idc(2,ii,i)
            n=idc(3,ii,i)
          enddo
        endif
      enddo

      deallocate(icommon,ibond,nb,idc)
      end subroutine cna
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make akr2cna"
c     End:
