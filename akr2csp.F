      program akr2csp
c-----------------------------------------------------------------------
c Reduce atoms from akr file using Centro-Symmetry Parameter (CSP).
c   - assuming akr file to be redirect to this program as STDIN
c   - output akr format to STDOUT
c   - change cutoff radius, rc, to apply this to your system
c-----------------------------------------------------------------------
      implicit none
c-----max num of neighbors; 12 is enough
      integer,parameter:: nnmax = 20
cc.....Cutoff radius for LJ-Ar FCC, halfway between 1st and 2nd N.N.
cc      real(8),parameter:: rc    = 3.36d-10 /0.5291772d-10
c      real(8),parameter:: rc    = 3.5d-10 /0.5291772d-10
c.....Cutoff radius for BCC Fe, halfway between 1st and 2nd N.N.
      real(8),parameter:: abcc = 2.8533d-10/0.5291772d-10
      real(8),parameter:: rc   = 0.5d0*(1d0+0.866025d0)*abcc

      integer:: i,n,itmp,nauxdat,j
      integer:: natm0,natm,nx,ny,nz,nlsa
      real(8):: tmp,h(3,3),xi(3),rij(3),r
      integer,allocatable:: lsa(:),lspr(:,:)
      real(8),allocatable:: ra0(:,:),ra(:,:),auxdat0(:,:)
     &     ,auxdat(:,:),cspi(:)

c-----read akr file from STDIN
      read(5,*) natm0,nauxdat
      read(5,*) h(1,1:3)
      read(5,*) h(2,1:3)
      read(5,*) h(3,1:3)
      allocate(ra0(3,natm0),ra(3,natm0),auxdat(nauxdat,natm0)
     &     ,auxdat0(nauxdat,natm0),cspi(natm0))
      do i=1,natm0
        read(5,*) itmp,ra0(1:3,i),auxdat0(1:nauxdat,i)
      enddo

c-----prepare for linked-list
      nx=int(h(1,1)/rc)
      ny=int(h(2,2)/rc)
      nz=int(h(3,3)/rc)
      nlsa= natm0 +nx*ny*nz
c      write(6,'(a,4i10)') " nx,ny,nz,nlsa=",nx,ny,nz,nlsa
      allocate(lsa(nlsa))
c-----make linked-list
      call mkllist(natm0,nx,ny,nz,ra0,nlsa,lsa)

c-----make 1st nearest neighbor pair-list
      allocate(lspr(0:nnmax,natm0))
      call ll2lspr(nlsa,lsa,nx,ny,nz,natm0,nnmax,lspr,h,ra0,rc)
      deallocate(lsa)

c-----common neighbor analysis
      call csp(natm0,nnmax,lspr,h,ra0,auxdat0,nauxdat
     &     ,natm,ra,auxdat,cspi)

c-----write akr file to STDOUT
      write(6,'(i10,3i4)') natm,nauxdat+1,0,0
      write(6,'(3es15.7)') h(1,1:3)
      write(6,'(3es15.7)') h(2,1:3)
      write(6,'(3es15.7)') h(3,1:3)
      do i=1,natm
        write(6,'(i4,20es11.3)') 1,ra(1:3,i)
     &       ,auxdat(1:nauxdat,i),cspi(i)
      enddo

      end program akr2csp
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
      subroutine ll2lspr(nlls,lls,nx,ny,nz,natm,nnmax,lspr,h,ra,rc)
c-----------------------------------------------------------------------
c     make Verlet list from linked-list
c-----------------------------------------------------------------------
      implicit none 
      integer,intent(in):: nlls,lls(nlls),nx,ny,nz,natm,nnmax
      real(8),intent(in):: h(3,3),ra(3,natm),rc
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
                  dx=(ra(1,j)-ra(1,i)-anint(ra(1,j)-ra(1,i)))*h(1,1)
                  dy=(ra(2,j)-ra(2,i)-anint(ra(2,j)-ra(2,i)))*h(2,2)
                  dz=(ra(3,j)-ra(3,i)-anint(ra(3,j)-ra(3,i)))*h(3,3)
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
      subroutine csp(natm0,nnmax,lspr,h,ra0,auxdat0,nauxdat
     &     ,natm,ra,auxdat,cspi)
c-----------------------------------------------------------------------
c  Reduce natm0 to natm based on the Centro-Symmetry Parameter.
c-----------------------------------------------------------------------
      implicit none
      integer,intent(in):: natm0,nnmax,nauxdat,lspr(0:nnmax,natm0)
      real(8),intent(in):: ra0(3,natm0),auxdat0(nauxdat,natm0),h(3,3)
      integer,intent(out):: natm
      real(8),intent(out):: ra(3,natm0),auxdat(nauxdat,natm0)
     &     ,cspi(natm0)

      integer:: i,j,k,l,m,kmin
      real(8):: xi(3),xj(3),xk(3),xij(3),xik(3),val,vmin,vjk

      natm= 0
      do i=1,natm0
        xi(1:3)= ra0(1:3,i)
        val= 0d0
        do l=1,lspr(0,i)
          j=lspr(l,i)
          
          xj(1:3)= ra0(1:3,j)-xi(1:3) -anint(ra0(1:3,j)-xi(1:3))
          xij(1:3)= h(1:3,1)*xj(1) +h(1:3,2)*xj(2) +h(1:3,3)*xj(3)

          vmin=1d10
          do m=1,lspr(0,i) 
            k=lspr(m,i)
            
            xk(1:3)= ra0(1:3,k)-xi(1:3) -anint(ra0(1:3,k)-xi(1:3))
            xik(1:3)= h(1:3,1)*xk(1) +h(1:3,2)*xk(2) +h(1:3,3)*xk(3)

            vjk= (xij(1)+xik(1))**2 +(xij(2)+xik(2))**2
     &           +(xij(3)+xik(3))**2
            if( vjk.lt.vmin ) then
              vmin= vjk
              kmin= k
            endif
          enddo

          val= val +vmin
        enddo
c.....Adopt half of val, because of double counting
        val= val/2
        if( val.gt.40.0d0 ) then
          natm=natm +1
          cspi(natm)= val
          ra(1:3,natm)= ra0(1:3,i)
          auxdat(1:nauxdat,natm)=auxdat0(1:nauxdat,i)
        endif
      enddo

      end subroutine csp
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make akr2csp"
c     End:
