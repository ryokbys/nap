      program akr2ddm
c-----------------------------------------------------------------------
c Create differential displacement map (DDM) in PS format
c   - akrref file needed
c   - assuming BCC and (111) as z-direction, which means this is
c     only applicable to screw dislocations
c   - change cutoff radius, rc, to apply this to your system
c-----------------------------------------------------------------------
c Usage:
c   $ ./akr2ddm akr000
c   - then, ddm000.ps will be created.
c-----------------------------------------------------------------------
      implicit none
c-----max num of neighbors; 12 is enough
      integer,parameter:: nnmax = 100
c.....Lattice constant for BCC W
      real(8),parameter:: abcc = 3.171d-10/0.5291772d-10
cc.....Lattice constant for BCC Fe
c      real(8),parameter:: abcc = 2.835d-10/0.5291772d-10
c.....Cutoff radius: halfway between 1st and 2nd N.N.
      real(8),parameter:: rc   = 0.5d0*(1d0+sqrt(3d0)/2)*abcc
c.....Cutoff radius: halfway between 2nd and 3rd N.N.
c      real(8),parameter:: rc   = 0.5d0*(1d0+sqrt(2d0))*abcc
c.....Scaling factor
      real(8),parameter:: sfac = 4
c.....radius of a circle
      real(8),parameter:: rcrc = 1.5d0
c.....width,height of the paper
      real(8),parameter:: alx  = 590d0
      real(8),parameter:: aly  = 840d0

      integer:: i,n,itmp,nauxdat,j
      integer:: natm,nx,ny,nz,nlsa,ndisp,nvec,isrc,idst
      real(8):: tmp,al(3),xi(3),rij(3),r,rlz,gray,x,y
     &     ,x1,y1,x2,y2,vmax,scale,xmax,ymax,sx,sy
      integer,allocatable:: lsa(:),lspr(:,:),itype(:),itag(:)
     &     ,lsza(:),lsdisp(:),ijdd(:,:)
      real(8),allocatable:: ra0(:,:),ra(:,:),auxdat0(:,:)
     &     ,auxdat(:,:),ga(:),vdd(:)

      integer,external:: iargc
      integer:: nargc
      character*128:: cin,cout

      nargc= iargc()
      if( nargc.ne.1 ) stop ' [Usage] ./ark2ddm akr000'
      call getarg(1,cin)


c-----read akr file from STDIN
      open(11,file=trim(cin),status='old')
      read(11,*) natm,nauxdat
      read(11,*) al(1),tmp,tmp
      read(11,*) tmp,al(2),tmp
      read(11,*) tmp,tmp,al(3)
      allocate(ra0(3,natm),ra(3,natm),auxdat(nauxdat,natm)
     &     ,auxdat0(nauxdat,natm),itype(natm),itag(natm)
     &     ,lsza(natm),lsdisp(natm),ga(natm)
     &     ,ijdd(2,8*natm),vdd(8*natm))
      do i=1,natm
        read(11,*) itag(i),ra(1:3,i),auxdat(1:nauxdat,i)
      enddo
      close(11)

c-----read akrref file
      open(10,file='akrref',status='old')
      read(10,*) itmp
      if( itmp.ne.natm ) stop ' [Error] itmp.ne.natm !!!'
      read(10,*) al(1),tmp,tmp
      read(10,*) tmp,al(2),tmp
      read(10,*) tmp,tmp,al(3)
      xmax= 0d0
      ymax= 0d0
      do i=1,natm
        read(10,*) itmp,ra0(1:3,i)
        xmax= max(xmax,ra0(1,i)*al(1)*sfac)
        ymax= max(ymax,ra0(2,i)*al(2)*sfac)
      enddo
      close(10)

      rlz= abcc*sqrt(3d0)/2
      call get_z_layer(natm,ra0,al,lsza,rlz,ga)

      write(6,'(a,2es12.4)') ' rlz,rlz/al(3)=',rlz,rlz/al(3)
      call get_display_atoms(natm,ra0,ndisp,lsdisp
     &     ,0d0,rlz/al(3))

c-----prepare for linked-list
      nx=max(int(al(1)/rc),1)
      ny=max(int(al(2)/rc),1)
      nz=max(int(al(3)/rc),1)
      nlsa= natm +nx*ny*nz
      write(6,'(a,4i10)') " nx,ny,nz,nlsa=",nx,ny,nz,nlsa
      allocate(lsa(nlsa))
c-----make linked-list
      call mkllist(natm,nx,ny,nz,ra0,nlsa,lsa)

c-----make 1st nearest neighbor pair-list
      allocate(lspr(0:nnmax,natm))
      call ll2lspr(nlsa,lsa,nx,ny,nz,natm,nnmax,lspr,al,ra0,rc)
      deallocate(lsa)

      call get_dd_vector(natm,ra,ra0,al,nvec,ijdd,vdd,nnmax,lspr
     &     ,abcc*sqrt(3d0)/2)

      sx= (alx-xmax)/2
      sy= (aly-ymax)/2
      write(6,'(a,2es12.4)') ' alx ,aly =',alx,aly
      write(6,'(a,2es12.4)') ' xmax,ymax=',xmax,ymax
      write(6,'(a,2es12.4)') ' sx  ,sy  =',sx,sy
c      stop

      cout=cin
      write(cout(1:3),'(a)') 'ddm'
      open(12,file=trim(cout)//'.ps',status='replace')
      call init_ps(12)
c.....Draw circles at atom positions
      do i=1,ndisp
        j= lsdisp(i)
        if( lsza(j).eq.1 ) then
          gray= 0.0d0 ! black
        elseif( lsza(j).eq.2 ) then
          gray= 0.5d0 ! gray
        else
          gray= 1.0d0 ! white
        endif
        x= ra0(1,j)*al(1)*sfac +sx
        y= ra0(2,j)*al(2)*sfac +sy
        call draw_circle(12,x,y,rcrc,gray)
      enddo
c.....Draw vectors of DDM
      vmax= 0d0
      do i=1,nvec
        vmax= max(vmax,vdd(i))
      enddo
      write(6,'(a,i6)') ' nvec=',nvec
      write(6,'(a,es12.4)') ' vmax=',vmax
      do i=1,nvec
        isrc= ijdd(1,i)
        idst= ijdd(2,i)
        x1= ra0(1,isrc)*al(1)*sfac +sx
        y1= ra0(2,isrc)*al(2)*sfac +sy
        x2= (ra0(1,idst)-anint(ra0(1,idst)-ra0(1,isrc)))*al(1)*sfac +sx
        y2= (ra0(2,idst)-anint(ra0(2,idst)-ra0(2,isrc)))*al(2)*sfac +sy
        scale= vdd(i)/vmax
        if( scale.lt.1d0/20 ) cycle
        call draw_arrow(12,x1,y1,x2,y2,scale,rcrc)
      enddo
      call final_ps(12)
      close(12)

      end program akr2ddm
c=======================================================================
      subroutine get_z_layer(natm,ra,al,lsza,rlz,ga)
c
c  Get layer number (1,2, or 3) of each atom looking at z-coordiate.
c
      implicit none
      integer,intent(in):: natm
      real(8),intent(in):: ra(3,natm),al(3),rlz
      integer,intent(out):: lsza(natm)
      real(8),intent(out):: ga(natm)
      
      integer:: i
      real(8):: z10
      
      do i=1,natm
        z10= ra(3,i)*al(3)
        ga(i)= (z10-int(z10/rlz)*rlz)/rlz
        if( ga(i).lt.1d0/6 .or. ga(i).gt.5d0/6 ) then
          lsza(i)= 1 ! to be black
        elseif( 1d0/6.le.ga(i) .and. ga(i).lt.3d0/6 ) then
          lsza(i)= 2 ! to be gray
        elseif( 3d0/6.le.ga(i) .and. ga(i).lt.5d0/6 ) then
          lsza(i)= 3 ! to be white
        endif
c        if( i.lt.9762 ) write(6,'(2i10,2es12.4)') i,lsza(i),z10,ga(i)
      enddo
      return
      end subroutine get_z_layer
c=======================================================================
      subroutine get_display_atoms(natm,ra,ndisp,lsdisp,z0,zw)
c
c  Extract atoms to be displayed within the z-width [Z0:Z0+ZW] in RA
c  and list atom-id in lsdisp.
c
      implicit none
      integer,intent(in):: natm
      real(8),intent(in):: ra(3,natm),z0,zw
      integer,intent(out):: ndisp,lsdisp(natm)
      real(8),parameter:: eps= 1d-8

      integer:: i

      ndisp= 0
      do i=1,natm
        if( ra(3,i).ge.z0 .and. ra(3,i).lt.z0+zw-eps ) then
          ndisp= ndisp +1
          lsdisp(ndisp)= i
        endif
      enddo
      write(6,'(a,i10)') ' num of atoms to be displayed=',ndisp
      return
      end subroutine get_display_atoms
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
      subroutine get_dd_vector(natm,ra,ra0,al,nvec,ijdd,vdd
     &     ,nnmax,lspr,bv)
      implicit none
      integer,intent(in):: natm,nnmax,lspr(0:nnmax,natm)
      real(8),intent(in):: ra(3,natm),ra0(3,natm),al(3),bv
      integer,intent(out):: nvec,ijdd(2,8*natm)
      real(8),intent(out):: vdd(8*natm)

      integer:: i,j,l,itmp
      real(8):: xi0(3),xi(3),ui(3),xj0(3),xj(3),uj(3),du(3)

      nvec= 0
      do i=1,natm-1
        xi0(1:3)= ra0(1:3,i)
        xi(1:3)= ra(1:3,i)
        ui(1:3)= (xi(1:3)-xi0(1:3)-anint(xi(1:3)-xi0(1:3)))*al(1:3)
c.....make z-displacement positive
        if( ui(3).lt.0d0 ) ui(3)=ui(3) +bv
c        write(6,'(a,2i6)') ' i,lspr(0,i)=',i,lspr(0,i)
        do l=1,lspr(0,i)
          j= lspr(l,i)
          if( j.le.i ) cycle
          xj0(1:3)= ra0(1:3,j)
          xj(1:3)= ra(1:3,j)
          if( nint((xi0(1)-xj0(1))*al(1)).eq.0 .and.
     &        nint((xi0(2)-xj0(2))*al(2)).eq.0 ) cycle
          uj(1:3)= (xj(1:3)-xj0(1:3)-anint(xj(1:3)-xj0(1:3)))*al(1:3)
c.....make z-displacement positive
          if( uj(3).lt.0d0 ) uj(3)=uj(3) +bv
          du(1:3)= abs(uj(1:3)-ui(1:3))
          nvec=nvec +1
          if( nvec.gt.8*natm ) stop ' [Error] nvec.gt.8*natm !!!'
          if( abs(ui(3)).ge.abs(uj(3)) ) then
            ijdd(1,nvec)= j
            ijdd(2,nvec)= i
          else
            ijdd(1,nvec)= i
            ijdd(2,nvec)= j
          endif
          if( du(3).gt.abs(du(3)-bv) ) then
            du(3)= abs(du(3) -bv)
c.....turn around the vector
            itmp= ijdd(1,nvec)
            ijdd(1,nvec)= ijdd(2,nvec)
            ijdd(2,nvec)= itmp
          endif
          vdd(nvec)= du(3)
        enddo
      enddo

      end subroutine get_dd_vector
c=======================================================================
      subroutine init_ps(ionum)
      implicit none
      integer,intent(in):: ionum

      write(ionum,'(a)') '%!PS'
      write(ionum,'(a)') '%'
c      write(ionum,'(a)') '/cm { 28.3 mul } def'
      write(ionum,'(a)') '0.5 setlinewidth'
      write(ionum,'(a)') '%'
      write(ionum,'(a)') 'gsave'
      return
      end subroutine init_ps
c=======================================================================
      subroutine final_ps(ionum)
      implicit none
      integer,intent(in):: ionum

      write(ionum,'(a)') 'showpage'
      write(ionum,'(a)') 'grestore'
      return
      end subroutine final_ps
c=======================================================================
      subroutine draw_circle(ionum,x,y,r,gray)
      implicit none
      integer,intent(in):: ionum
      real(8),intent(in):: x,y,r,gray
      write(ionum,'(3f8.2,a)') x,y,r,' 0  360  arc'
      write(ionum,'(a)') 'gsave'
      write(ionum,'(2x,f8.2,a)') gray,' setgray'
      write(ionum,'(a)') '  fill'
      write(ionum,'(a)') 'grestore'
      write(ionum,'(a)') 'stroke'
      return
      end subroutine draw_circle
c=======================================================================
      subroutine draw_triangle(ionum,x,y,al,cs,sn)
c
c  draw triangle of an edge-length 'al', orientation (cs,sn), 
c  center of mass position (x,y)
c
      implicit none
      integer,intent(in):: ionum
      real(8),intent(in):: x,y,al,cs,sn
      real(8),parameter:: pi = 3.14159265358979d0
      real(8),parameter:: ang= pi*5d0/6

      real(8):: ex,ey,ex1,ey1,ex2,ey2

      ex= cos(ang)
      ey= sin(ang)
      ex1= cs*ex -sn*ey
      ey1= sn*ex +cs*ey
      ex= cos(-pi/2)
      ey= sin(-pi/2)
      ex2= cs*ex -sn*ey
      ey2= sn*ex +cs*ey

      write(ionum,'(a)') ' newpath'
      write(ionum,'(2f8.2,a)')  x+cs*al/sqrt(3d0),y+sn*al/sqrt(3d0)
     &     ,' moveto'
      write(ionum,'(2f8.2,a)')  ex1*al, ey1*al,' rlineto'
      write(ionum,'(2f8.2,a)')  ex2*al, ey2*al,' rlineto'
      write(ionum,'(a)') ' closepath fill'

      return
      end subroutine draw_triangle
c=======================================================================
      subroutine draw_arrow(ionum,x1,y1,x2,y2,scale,rcrc)
      implicit none
      integer,intent(in):: ionum
      real(8),intent(in):: x1,y1,x2,y2,scale,rcrc

      real(8):: rx,ry,r,cs,sn,rl,dr,x0,y0,dx,dy,x3,y3,ex,ey
     &     ,ex1,ey1,ex2,ey2
      real(8),parameter:: pi = 3.14159265358979d0
      real(8),parameter:: angv = pi -pi/6
c.....arrow length, 3 points
      real(8),parameter:: alen = 2d0

      rx= x2-x1
      ry= y2-y1
      r= sqrt(rx*rx +ry*ry)
      cs= rx/r
      sn= ry/r

      rl= r-2d0*rcrc
      dr= rl*scale
      x0= x1+ (r-dr)/2 *cs
      y0= y1+ (r-dr)/2 *sn
      dx= dr*cs
      dy= dr*sn
      x3= x0+dx
      y3= y0+dy
      write(ionum,'(a)') ' newpath'
      write(ionum,'(2f8.2,a)') x0,y0,' moveto'
      write(ionum,'(2f8.2,a)') x3,y3,' lineto'
      write(ionum,'(a)') ' stroke'

c.....arrow decoration with triangle
      call draw_triangle(ionum,x3,y3,alen,cs,sn)

      return
      end subroutine draw_arrow
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make akr2ddm"
c     End:
