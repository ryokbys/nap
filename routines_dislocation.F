      subroutine init_shear_stress(shrst,shrfx,natm,tag,ra,h,sorg,acon
     &     ,myid_md,mpi_md_world)
      implicit none
      include "mpif.h"
      integer,intent(in):: natm,myid_md,mpi_md_world
      real(8),intent(in):: shrst,tag(natm),ra(3,natm)
     &     ,h(3,3,0:1),sorg(3),acon
      real(8),intent(out):: shrfx

      integer:: n,m,nt,mt,i,ierr,l
      real(8):: area

      n=0
      m=0
      do i=1,natm
        l= int(mod(tag(i)*10,10d0))
        if( l.eq.2 ) then
          if( ra(2,i)+sorg(2).gt.0.5d0 ) then
            n=n+1
          else
            m=m+1
          endif
        endif
      enddo
      nt=0
      mt=0
      call mpi_allreduce(n,nt,1,mpi_integer,mpi_sum
     &     ,mpi_md_world,ierr)
      call mpi_allreduce(m,mt,1,mpi_integer,mpi_sum
     &     ,mpi_md_world,ierr)

      area= h(1,1,0)*h(3,3,0)

c.....Change unit of shear stress (GPa ---> a.u.)
c.....Shear stress to shear force on each atom
      shrfx= shrst/2.94d+4 *area /nt

      if( myid_md.eq.0 ) then
        write(6,'(a,2i8)') " Num of atoms to be shifted=",nt,mt
        write(6,'(a,es12.4)') " Top & bottom area [Bohr^2]=",area
        write(6,'(a,es12.4)') " Shear stress [a.u.]=",shrst/2.94d+4
        write(6,'(a,es12.4)') " Shear force  [a.u.]=",shrfx
      endif

c.....Normalize shrfx by Lx, and multiply 0.5d0*dt**2/am(i)
      shrfx= shrfx/h(1,1,0) *acon

      end subroutine init_shear_stress
c=======================================================================
      subroutine add_shear_stress(natm,ra,aa,tag,sorg,shrfx)
      implicit none
      integer,intent(in):: natm
      real(8),intent(in):: ra(3,natm),tag(natm),sorg(3),shrfx
      real(8),intent(inout):: aa(3,natm)

      integer:: i,l
      
c.....Shear stress on top and bottom edge
      do i=1,natm
        l= int(mod(tag(i)*10,10d0))
        if( l.eq.2 ) then
          if( ra(2,i)+sorg(2).gt.0.5d0 ) then
            aa(1,i)= aa(1,i) +shrfx
          else
            aa(1,i)= aa(1,i) -shrfx
          endif
        endif
      enddo

      end subroutine add_shear_stress
c=======================================================================
      subroutine average_stress_around_core(rcore,namax,natm
     &     ,ra,h,myid_md,mpi_md_world,strs,istp)
c-----------------------------------------------------------------------
c  Compute average stress around dislocation core using strs(3,3,natm).
c  Results are written to the file #91, 'dat.strs'.
c-----------------------------------------------------------------------
      implicit none
      include 'mpif.h'
      include "./params_au.h"
      integer,intent(in):: namax,natm,myid_md,mpi_md_world,istp
      real(8),intent(in):: ra(3,natm),h(3,3,0:1),rcore(3),strs(3,3,natm)

c.....Cut off distance for average calculation
      real(8),parameter:: rcutave= 20.0d-10 /bohr
      real(8),parameter:: rcutave2= rcutave**2

      integer:: ia,ierr,n1,nr1,n2,nr2
      real(8):: al(3),strs1(2,2),strsr1(2,2),strs2(2,2),strsr2(2,2)
     &     ,xic(2),ric2

c.....Assuming orthogonal 2D system
      al(1)= h(1,1,0)
      al(2)= h(2,2,0)

c.....Find atoms within a certain distance from the disl core
      n1= 0
      strs1(1:2,1:2)= 0d0
      n2= 0
      strs2(1:2,1:2)= 0d0
      do ia=1,natm
        xic(1:2)= (ra(1:2,ia)-rcore(1:2))*al(1:2)
        ric2= xic(1)**2 +xic(2)**2
c.....Store stress data of atoms within the distance from the core
        if( ric2 .le. rcutave2 ) then
          n1=n1 +1
          strs1(1:2,1:2)=strs1(1:2,1:2) +strs(1:2,1:2,ia)
        endif
c.....Store stress data of atoms inside the band along x
        if( abs(xic(2)) .lt. rcutave ) then
          n2=n2+1
          strs2(1:2,1:2)=strs2(1:2,1:2) +strs(1:2,1:2,ia)
        endif
      enddo

c.....Gather stress data of the atoms to node-0
      call mpi_reduce(n1,nr1,1,mpi_integer,mpi_sum,0,mpi_md_world,ierr)
      call mpi_reduce(n2,nr2,1,mpi_integer,mpi_sum,0,mpi_md_world,ierr)
      call mpi_reduce(strs1,strsr1,4,mpi_double_precision,mpi_sum,0
     &     ,mpi_md_world,ierr)
      call mpi_reduce(strs2,strsr2,4,mpi_double_precision,mpi_sum,0
     &     ,mpi_md_world,ierr)

      if( myid_md .eq. 0 ) then
c.....Compute average of the stress
        strs1(1:2,1:2)= strsr1(1:2,1:2)/nr1
        strs2(1:2,1:2)= strsr2(1:2,1:2)/nr2

c.....Output the stress value only node-0
        write(91,'(i8,8es12.4,2i8)') istp,strs1(1:2,1:2),strs2(1:2,1:2)
     &       ,nr1,nr2
      endif
      
      end subroutine average_stress_around_core
c=======================================================================
      subroutine disl_pos_2Dtri(namax,natm,nb,nnmax,ra,tag,h,hi,epi
     &     ,myid_md,mpi_md_world,nodes_md,anxi,anyi,anzi
     &     ,sorg,rcore,rcpre,ndpbc)
c
c Find dislocation core by looking at number of nearest neighbours.
c Note this can be applicable only to 2D triangle lattice.
c
      implicit none
      include 'mpif.h'
      include './params_au.h'
      include './params_LJ_Ar.h'
c.....Arguments
      integer,intent(in):: namax,natm,nb,nnmax,myid_md,mpi_md_world
     &     ,nodes_md
      real(8),intent(in):: sorg(3),ra(3,natm),tag(natm),epi(natm)
     &     ,h(3,3,0:1),hi(3,3),anxi,anyi,anzi
      real(8),intent(inout):: rcore(3),rcpre
      integer,intent(inout):: ndpbc
c.....Constants
      integer,parameter:: ncmax = 100
      real(8),parameter:: rcnn  = 1.5d0*sgmlj
c.....Variables
      integer:: i,j,k,l,m,n,nc,ierr
      integer,allocatable,save:: lnc(:),nrecvs(:),idisps(:),lsnn(:,:)
      real(8),allocatable,save:: rci(:,:),rcis(:,:),add(:)
      logical,save:: l1st=.true.
      
      if(l1st) then
        allocate(lnc(0:nodes_md-1),nrecvs(0:nodes_md-1)
     &       ,idisps(0:nodes_md-1)
     &       ,rci(3,ncmax),rcis(3,ncmax),add(ncmax))
        rcpre= 0.5d0
c.....Nearest neighbor search
        allocate(lsnn(0:nnmax,namax))
        l1st=.false.
      endif

      lnc(0:nodes_md-1)= 0
      nrecvs(0:nodes_md-1)= 0
      idisps(0:nodes_md-1)= 0

c.....Count num of nearest neighbors
      call mk_lspr(namax,natm,nb,nnmax,tag,ra,rcnn,h,hi
     &     ,anxi,anyi,anzi,lsnn)
c.....Extract disl-core atoms according to num of nearest neighbors
      nc=0
      do i=1,natm
c.....Exclude top and bottom surface atoms
        l= int(mod(tag(i)*10,10d0))
        if( l.eq.2 ) cycle
        if( lsnn(0,i).lt.6 ) then
          nc=nc+1
          if( nc.gt.ncmax ) then
            write(6,'(a,i4)') '[Error] nc.gt.ncmax, myid_md=',myid_md
            stop 
          endif
          rci(1:3,nc)= ra(1:3,i)+sorg(1:3)
        endif
      enddo
      
c      nc=0
c      do i=1,natm
cc.....Exclude top and bottom surface atoms
c        l= int(mod(tag(i)*10,10d0))
c        if( l.eq.2 ) cycle
cc.....Extract disl-core atoms by looking at potential energy per atom
c        if( epi(i).gt.-8.5d-4 ) then
c          nc=nc+1
c          if( nc.gt.ncmax ) then
c            write(6,'(a,i4)') '[Error] nc.gt.ncmax, myid_md=',myid_md
c            stop 
c          endif
c          rci(1:3,nc)= ra(1:3,i)+sorg(1:3)
c        endif
c      enddo
      
c.....Gather disl-core information
      lnc(myid_md)= nc
      call mpi_allgather(nc,1,mpi_integer,lnc,1,mpi_integer
     &     ,mpi_md_world,ierr)
c.....Count num of core atoms
      n=0
      do i=0,nodes_md-1
        n= n +lnc(i)
      enddo
      if( n.gt.ncmax ) then
        write(6,'(a)') '[Error] n.gt.ncmax'
        call mpi_finalize(ierr)
        stop
      endif
      idisps(0:nodes_md)= 0
      do i=1,nodes_md-1
        idisps(i)= idisps(i-1) +lnc(i-1)
      enddo
      do i=0,nodes_md-1
        nrecvs(i)= lnc(i)
      enddo
      idisps(0:nodes_md-1)= idisps(0:nodes_md-1)*3
      nrecvs(0:nodes_md-1)= nrecvs(0:nodes_md-1)*3
      call mpi_gatherv(rci,3*nc,mpi_double_precision
     &     ,rcis,nrecvs,idisps,mpi_double_precision
     &     ,0,mpi_md_world,ierr)
      if( myid_md.eq.0 ) then
c.....Check pbc
        add(1:n)=0d0
        do i=1,n
          if( rcpre-rcis(1,i).gt.0.5d0 ) add(i)=1d0
          if( rcpre-rcis(1,i).lt.-0.5d0 ) add(i)=-1d0
        enddo
c.....Determine disl core as center of core atoms
        rcore(1:3)= 0d0
        do i=1,n
          rcore(1:3)= rcore(1:3) +rcis(1:3,i) +add(i)
        enddo
        rcore(1:3)= rcore(1:3)/n
        if( rcore(1).gt.1d0 ) rcore(1)=rcore(1)-1d0
        if( rcpre-rcore(1).gt.0.5d0 ) then
          ndpbc=ndpbc+1
        endif
        rcpre=rcore(1)
      endif

      end subroutine disl_pos_2Dtri
c=======================================================================
      subroutine disl_pos_by_pot(epith,natm,ra,h,epi,sorg
     &     ,nodes,myid_md,mpi_md_world,iocnt,ionum)
c
c Find dislocation core by looking at a potential energy value.
c Taking into account there are two disl cores (dissociated disl)
c in the system.
c
      implicit none
      include 'mpif.h'
      include './params_au.h'
      include './params_EAM_Fe-H.h'
      integer,intent(in):: natm,myid_md,mpi_md_world,nodes,iocnt,ionum
c.....Potential energy threshold
      real(8),intent(in):: epith
      real(8),intent(in):: ra(3,natm),h(3,3,0:1),epi(natm),sorg(3)

c.....Max num of core atoms
      integer,parameter:: mnc= 10000
c.....Max num of core groups
      integer,parameter:: mng= 2

      integer:: i,n,nc,ierr,ndx,ndy,ndz,ndxz,ndxy,mx,my,mz,m
     &     ,ig,ng,kx,ky,kz,m1x,m1y,m1z,m1,j
      real(8):: al(3),dx,dy,dz,avx
      integer,allocatable:: lscl(:),lshd(:)
      character:: cfname*9

c.....Cutoff radius
      real(8),parameter:: rc   = 10.0d-10/bohr
      real(8),parameter:: rcz  = rc_vphi
c.....Num of disl cores to be assumed
      integer,parameter:: ncrs = 2
      logical,save:: l1st=.true.
      integer,allocatable,save:: lnc(:),idisps(:)
      real(8),allocatable,save:: rci(:,:),rcis(:,:,:)
      integer,allocatable,save:: ngc(:),igc(:,:)

      if( l1st ) then
        allocate(lnc(0:nodes-1),idisps(0:nodes-1)
     &       ,rci(3,mnc),rcis(3,mnc,0:nodes-1)
     &       ,ngc(mng),igc(mnc,mng))
        l1st=.false.
      endif

c.....Search and list up disl core atoms
      nc=0
      do i=1,natm
        if( epi(i).gt.epith .and.
     &       ra(2,i)+sorg(2).lt.0.75d0 .and.
     &       ra(2,i)+sorg(2).gt.0.25d0 ) then
          nc=nc+1
          rci(1:3,nc)= ra(1:3,i)+sorg(1:3)
        endif
      enddo

c.....Gather lsdc from all nodes to node-0
      lnc(myid_md)= nc
      call mpi_allgather(nc,1,mpi_integer,lnc,1,mpi_integer
     &     ,mpi_md_world,ierr)

      rcis(1:3,1:nc,myid_md)= rci(1:3,1:nc)
      idisps(0)= 0
      do i=1,nodes-1
        idisps(i)=idisps(i-1) +lnc(i-1)*3
      enddo
      call mpi_gatherv(rci,3*nc,mpi_double_precision
     &     ,rcis,lnc*3,idisps
     &     ,mpi_double_precision,0,mpi_md_world,ierr)

      if( myid_md.eq.0 ) then
        nc=0
        do i=0,nodes-1
          nc=nc +lnc(i)
        enddo
c        write(6,'(a,i8)') " nc=",nc
        
c.....Assuming that there are two cores, seperate core atoms into
c.....  two groups using linked-list-cell technique
        al(1)= h(1,1,0)
        al(2)= h(2,2,0)
        al(3)= h(3,3,0)
        ndx= al(1)/rc
        ndy= al(2)/rc
        ndz= al(3)/rc
        dx = 1d0/ndx
        dy = 1d0/ndy
        dz = 1d0/ndz
        ndxy= ndx*ndy
        allocate(lscl(nc),lshd(ndxy))
c.....Register atoms to cells
        lshd(1:ndxy)= 0
        
        do i=1,nc
          mx= (rcis(1,i,0)+dx)/dx ! [1:ndx]
          my= (rcis(2,i,0)+dy)/dy
c          mz= (rcis(3,i,0)+dz)/dz
          m= (mx-1)*ndy +my
          lscl(i)=lshd(m)
          lshd(m)= i
        enddo
c.....Make groups
        ig=0
        ngc(1:mng)= 0
        igc(1:mnc,1:mng)=0
        do 10 mx=1,ndx
        do 10 my=1,ndy
          m= (mx-1)*ndy +my

          if( lshd(m).eq.0 ) goto 10
          i=lshd(m)
          ig=ig+1
          if(ig.gt.mng) goto 90
          ngc(ig)=ngc(ig) +1
          igc(ngc(ig),ig)= i
          do kx=-1,1
            m1x=mx +kx
            if( m1x.lt.1 ) m1x=m1x +ndx
            if( m1x.gt.ndx ) m1x=m1x -ndx
            do ky=-1,1
              m1y= my +ky
              if( m1y.lt.1 ) m1y=m1y +ndy
              if( m1y.gt.ndy ) m1y=m1y -ndy
              m1= (m1x-1)*ndy +m1y

              if( lshd(m1).eq.0 ) cycle

 1            continue
              j=lshd(m1)
              lshd(m1)= 0

 2            continue
              if( j.eq.i ) goto 3

              ngc(ig)=ngc(ig) +1
              igc(ngc(ig),ig)= j

 3            j=lscl(j)
              if( j.gt.0 ) goto 2

            enddo
          enddo
 10     continue

 90     ng=ig
        if( ng.gt.mng ) ng=mng

c        do ig=1,ng
c          do j=1,ngc(ig)
c            write(6,'(a,3i6)') " ig,j,igc(j,ig)=",ig,j,igc(j,ig)
c          enddo
c        enddo

c.....Get disl pos at each y-slice
        cfname='disl00000'
        write(cfname(5:9),'(i5.5)') iocnt
        open(ionum,file=cfname,status='replace')
        ndz= al(3)/rcz
        dz = 1d0/ndz
        do mz=1,ndz
          do ig=1,ng
            avx=0d0
            n=0
            do i=1,ngc(ig)
              j= igc(i,ig)
              if( rcis(3,j,0).ge.(mz-1)*dz .and.
     &             rcis(3,j,0).lt.mz*dz ) then
                avx=avx +rcis(1,j,0)
                n=n+1
              endif
            enddo
            avx=avx/n
            write(ionum,'(2es12.4,$)') (mz-1d0+0.5d0)*dz,avx
          enddo
          write(ionum,'(a)') ''
        enddo
        close(ionum)

        deallocate(lshd,lscl)
      endif
      
      end subroutine disl_pos_by_pot
c=======================================================================
      subroutine perf_disl_pos_by_pot(epith,natm,ra,h,epi,sorg
     &     ,nodes,myid_md,mpi_md_world,iocnt,ionum)
c
c  Find dislocation core of perfect dislocation by looking at a 
c  potential energy value of each atom.
c  Assuming disl line to be z-direction and there are surfaces 
c  at the top and bottom of y-direction.
c
      implicit none
      include 'mpif.h'
      include './params_au.h'
      include './params_EAM_Fe-H.h'
      integer,intent(in):: natm,myid_md,mpi_md_world,nodes,iocnt,ionum
c.....Potential energy threshold
      real(8),intent(in):: epith
      real(8),intent(in):: ra(3,natm),h(3,3,0:1),epi(natm),sorg(3)

c.....Max num of core atoms
      integer,parameter:: mnc= 1000

      integer:: i,n,nc,mz,ierr
      real(8):: avx,xref,xav
      character:: cfname*9

c.....Cutoff radius
      real(8),parameter:: rcz  = rc_vphi
      logical,save:: l1st=.true.
      integer,save:: ndz
      real(8),save:: dz
      integer,allocatable,save:: lnc(:),idisps(:)
      real(8),allocatable,save:: rci(:,:),rcis(:,:,:)

      if( l1st ) then
        allocate(lnc(0:nodes-1),idisps(0:nodes-1)
     &       ,rci(3,mnc),rcis(3,mnc,0:nodes-1))
        ndz= h(3,3,0)/rcz
        dz= 1d0/ndz
        l1st=.false.
      endif

c.....Search and list up disl core atoms
      nc=0
      do i=1,natm
        if( epi(i).gt.epith .and.
     &       ra(2,i)+sorg(2).lt.0.75d0 .and.
     &       ra(2,i)+sorg(2).gt.0.25d0 ) then
          nc=nc+1
          rci(1:3,nc)= ra(1:3,i)+sorg(1:3)
        endif
      enddo

c.....Gather lsdc from all nodes to node-0
      lnc(myid_md)= nc
      call mpi_allgather(nc,1,mpi_integer,lnc,1,mpi_integer
     &     ,mpi_md_world,ierr)

      rcis(1:3,1:nc,myid_md)= rci(1:3,1:nc)
      idisps(0)= 0
      do i=1,nodes-1
        idisps(i)=idisps(i-1) +lnc(i-1)*3
      enddo
      call mpi_gatherv(rci,3*nc,mpi_double_precision
     &     ,rcis,lnc*3,idisps
     &     ,mpi_double_precision,0,mpi_md_world,ierr)

      if( myid_md.eq.0 ) then
        nc= 0
        do i=0,nodes-1
          nc=nc +lnc(i)
        enddo

        cfname='disl00000'
        write(cfname(5:9),'(i5.5)') iocnt
        open(ionum,file=cfname,status='replace')
        xref= rcis(1,1,0)
        do mz=1,ndz
          n=0
          do i=1,nc
            if( rcis(3,i,0).ge.(mz-1)*dz .and.
     &           rcis(3,i,0).lt.mz*dz ) then
              xav=xav +rcis(1,i,0) -anint(rcis(1,i,0)-xref)
              n=n+1
            endif
          enddo
          xav=xav/n
          write(ionum,'(2es12.4)') (mz-1d0+0.5d0)*dz, xav
        enddo
        close(ionum)
      endif
      
      end subroutine perf_disl_pos_by_pot
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make pmd"
c     End:
