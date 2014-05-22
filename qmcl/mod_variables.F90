  module variables
    implicit none
    save

!.....max. num. of atoms
    integer,parameter:: namax = 200000
!.....max. num. of boundary atoms
    integer,parameter:: nbmax = 100000
!.....max. num. of neighbors
    integer,parameter:: nnmax = 20
!.....max. num. of species
    integer,parameter:: nismax= 9
!.....Max. num. of atoms in a cluster
    integer,parameter:: nacmax= 1000
!.....Max. num. of term-H atoms in a cluster
    integer,parameter:: nhcmax= 200

    integer,parameter:: ioerg   = 11
    character,parameter:: cferg*7 = 'out.erg'

    integer:: nstp,nouterg,noutpmd,istp,nerg &
         ,npmd,iocntpmd,iocnterg
    integer:: natm,nb,nis
    integer:: istp0 = 0
    real(8):: tcpu,tcpu1,tcpu2,tcom
    real(8):: dt,dmp,treq,trlx,temp,epot,ekin,epot0
    real(8):: tinit = -10d0
    real(8):: xc(3)
    real(8):: rcut  =  5.3d0
    real(8):: rclst = 10d0 ! 10 AA (initial value)
    real(8):: rcqm  =  5d0 !  5 AA (initial value)
    real(8):: rclx  = 10d0
    real(8):: rcly  = 10d0
    real(8):: rclz  =-10d0
    real(8):: rcqmx =  5d0
    real(8):: rcqmy =  5d0
    real(8):: rcqmz = -5d0
    real(8):: wvac(3)=(/ 0d0, 0d0, 0d0 /)
    character(len=6):: cqmshape=  'circle'
    character(len=6):: ciofmt  = 'ascii'
    character(len=4):: cfrct   = 'qmcl'
    character(len=7):: cfnorm  = 'vnrm000'
    character(len=4):: crtype  = 'MD' ! MD, CG, test
    logical:: lpmd    = .true.
    logical:: ltctl   = .false.
    logical:: ldmp    = .false.
    logical:: lcsdmp  = .false.
    logical:: lconv   = .false.
    logical:: lconst  = .false.
    logical:: lterm   = .false.
    real(8):: feps    = 1d-2
    real(8):: deps    = 1d-2

!.....simulation box
    real(8):: hunit,h(3,3,0:1),hi(3,3),vol,sgm(3,3),al(3),avol
!.....factors on each moving direction
    real(8):: fmv(3,0:9)=reshape([ &
         0d0, 0d0, 0d0, & ! 0: fixed
         1d0, 1d0, 1d0, & ! 1: free
         0d0, 1d0, 1d0, & ! 2: yz-plane
         1d0, 0d0, 1d0, & ! 3: xz-plane
         1d0, 1d0, 0d0, & ! 4: xy-plane
         1d0, 0d0, 0d0, & ! 5: x-direction
         0d0, 1d0, 0d0, & ! 6: y-direction
         0d0, 0d0, 1d0, & ! 7: z-direction
         1d0, 1d0, 1d0, & ! 8:
         1d0, 1d0, 1d0  & ! 9:
         ],[3,10])
!.....positions, velocities, and accelerations
    real(8):: ra(3,namax),va(3,namax),fa(3,namax) &
         ,strs(3,3,namax),stt(3,3,namax)
!.....real*8 identifier which includes species, index of FMV, total id
    real(8):: tag(namax)
!.....potential and kinetic energy per atoms
    real(8):: epi(namax),eki(namax),stp(3,3,namax),erg
!.....mass, prefactors
    real(8):: am(nismax),acon(nismax),fack(nismax)
!.....strain
    real(8):: stn(3,3,namax),ra0(3,namax),h0(3,3,0:1)
!.....Pair list
    integer:: lspr(0:nnmax,namax)
!.....for free energy calculation?
    real(8):: vn(3,namax),facc(namax)

!.....QM cluster
    integer:: nac,naqm
    real(8):: rac(3,nacmax),tagc(nacmax),fac(3,nacmax)
    real(8):: hc(3,3) ! box size of cluster calc.
    integer:: idclst(nacmax),idqm(nacmax)
!.....H termination
    integer:: nhc
    integer:: rhc(3,nhcmax)

  contains
!=======================================================================
    subroutine bacopy(rc)
      implicit none
      real(8),intent(in):: rc

      integer:: i,ixyz
      real(8):: rch

      nb=0

      do i=1,natm
        do ixyz=1,3
          rch= rc /h(ixyz,ixyz,0)
          if( ra0(ixyz,i).lt.rch ) then ! bottom
            nb=nb +1
            if( nb.gt.nbmax ) then
              write(6,'(a)') ' [bacopy] nb.gt.nbmax !!!'
              stop
            endif
            ra0(1:3,natm+nb)=ra0(1:3,i)
            ra0(ixyz,natm+nb)= ra0(ixyz,natm+nb) +1d0
            ra(1:3,natm+nb)=ra(1:3,i)
            if( abs(ra0(ixyz,i)-ra(ixyz,i)).lt.0.5d0 ) then
              ra(ixyz,natm+nb)= ra(ixyz,natm+nb) +1d0
            endif
          endif
          if( ra0(ixyz,i).gt.1d0-rch ) then ! top
            nb=nb +1
            if( nb.gt.nbmax ) then
              write(6,'(a)') ' [bacopy] nb.gt.nbmax !!!'
              stop
            endif
            ra0(1:3,natm+nb)=ra0(1:3,i)
            ra0(ixyz,natm+nb)= ra0(ixyz,natm+nb) -1d0
            ra(1:3,natm+nb)=ra(1:3,i)
            if( abs(ra0(ixyz,i)-ra(ixyz,i)).lt.0.5d0 ) then
              ra(ixyz,natm+nb)= ra(ixyz,natm+nb) -1d0
            endif
          endif
        enddo
      enddo
      
    end subroutine bacopy
!=======================================================================
    subroutine mklspr(rc)
!-----------------------------------------------------------------------
!  Make pair-list from linked-list
!-----------------------------------------------------------------------
      implicit none 
      real(8),intent(in):: rc

      integer:: i,ic1,ic2,ic3,ic,mx,my,mz,kux,kuy,kuz,m1x,m1y,m1z,m1,j,m
      real(8):: xij(3),rij2,xi(3)

      integer,save:: lcx,lcy,lcz,lcx2,lcy2,lcz2,lcyz2,ncl
      real(8),save:: a1,a2,a3,rc2
      integer,allocatable,save:: lshd(:),lscl(:)
      logical,save:: l1st = .true.

      if( l1st ) then
!.....Number of cells in each orientation, assuming orthorhombic
        lcx= int(h(1,1,0)/rc)
        lcy= int(h(2,2,0)/rc)
        lcz= int(h(3,3,0)/rc)
        if( lcx.eq.0 .or. lcy.eq.0 .or. lcz.eq.0 ) then
          write(6,'(a,3i5)') ' [mklsr] lc?.eq.0, lcx,lcy,lcz='&
               ,lcx,lcy,lcz
!          stop
          if(lcx.eq.0) lcx= 1
          if(lcy.eq.0) lcy= 1
          if(lcz.eq.0) lcz= 1
        endif
        a1= 1d0/lcx
        a2= 1d0/lcy
        a3= 1d0/lcz
        lcx2= lcx+2
        lcy2= lcy+2
        lcz2= lcz+2
        lcyz2= lcy2*lcz2
        ncl= lcx2*lcy2*lcz2
        rc2= rc**2
        write(6,'(a,3i8)') ' lcx,lcy,lcz       =',lcx,lcy,lcz
        write(6,'(a,4i8)') " lcx2,lcy2,lcz2,ncl=",lcx2,lcy2,lcz2,ncl
        allocate(lshd(ncl),lscl(natm+nbmax))
        l1st= .false.
      endif

      call bacopy(rc)

!-----Here we use ra0 for making list
      lshd(1:ncl)= 0
      do i=1,natm+nb
!-------region index in order z,y,x
        ic3= (ra0(3,i)+a3)/a3
        ic2= (ra0(2,i)+a2)/a2
        ic1= (ra0(1,i)+a1)/a1
        if( i.le.natm ) then
          ic3= min(max(ic3,1),lcz)
          ic2= min(max(ic2,1),lcy)
          ic1= min(max(ic1,1),lcz)
        else
          ic3= min(max(ic3,0),lcz+1)
          ic2= min(max(ic2,0),lcy+1)
          ic1= min(max(ic1,0),lcz+1)
        endif
        ic= ic1*lcyz2 +ic2*lcz2 +ic3+1
        lscl(i)= lshd(ic)
        lshd(ic)= i
      enddo
!      stop

!.....Only residences (1:natm) need neighbors
      lspr(0:nnmax,1:natm)=0
      
!.....Scan resident cells
      do mz=1,lcz
        do my=1,lcy
          do mx=1,lcx
            m= mx*lcyz2 +my*lcz2 +mz+1
            if( lshd(m).eq.0 ) cycle
            do kuz=-1,1
              do kuy=-1,1
                do kux=-1,1
                  m1x= mx +kux
                  m1y= my +kuy
                  m1z= mz +kuz
                  m1= m1x*lcyz2 +m1y*lcz2 +m1z+1
                  if( lshd(m1).eq.0 ) cycle
                  
                  i=lshd(m)
1                 continue
                  if( i.gt.natm ) goto 4

                  xi(1:3)= ra0(1:3,i)
                  j=lshd(m1)
2                 continue
                  if( j.eq.i ) goto 3
                  xij(1:3)= ra0(1:3,j) -xi(1:3)
                  xij(1:3)= h(1:3,1,0)*xij(1) &
                       +h(1:3,2,0)*xij(2) &
                       +h(1:3,3,0)*xij(3)
                  rij2= xij(1)**2 +xij(2)**2 +xij(3)**2

                  if( rij2.le.rc2 ) then
                    lspr(0,i)=lspr(0,i)+1
                    if( lspr(0,i).gt.nnmax ) then
                      write(6,'(a)') ' [mklspr] lspr(0,i).gt.nnmax !!!'
                      stop
                    endif
                    lspr(lspr(0,i),i)= j
                  endif

!.....Continue until j=0
3                 j=lscl(j)
                  if( j.gt.0 ) goto 2
!.....Continue until i=0
4                 i=lscl(i)
                  if( i.gt.0 ) goto 1
                enddo
              enddo
            enddo
          enddo !mx
        enddo !my
      enddo !mz

    end subroutine mklspr
!=======================================================================
    function is_tag(tag)
      implicit none
      real(8),intent(in):: tag
      integer:: is_tag
      
      is_tag= int(tag)
      return
    end function is_tag
!=======================================================================
    function ifmv_tag(tag)
      implicit none
      real(8),intent(in):: tag
      integer:: ifmv_tag
      
      ifmv_tag= int((tag -int(tag))*10)

      return
    end function ifmv_tag
!=======================================================================
    function itot_tag(tag)
      implicit none
      real(8),intent(in):: tag
      integer:: itot_tag
      real(8):: tmp

      tmp= tag -is_tag(tag) -ifmv_tag(tag)*1d-1
      itot_tag= nint(tmp*1d+14)
      return
    end function itot_tag
!=======================================================================
  end module variables
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make qmcl"
!     End:
