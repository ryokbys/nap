      program mkconf_Si_disl_scew_60
c-----------------------------------------------------------------------
c  Make atomic configuration of diamond structure with:
c  (in case of screw or 60-degree)
c    x: [-1, 2,-1]
c    y: [ 1, 1, 1]
c    z: [-1, 0, 1]
c  (in case of 30-degree)
c    x: [-1, 0, 1]
c    y: [ 1, 1, 1]
c    z: [ 1,-2, 1]/2
c  (in case of 41-degree)
c    x: [-5, 1, 2]/2
c    y: [ 1, 1, 1]
c    z: [-1, 3,-2]/2
c  and create dislocation at the center of the xy plane (z-direction).
c  This configuration allows us to make screw or 60 dislocation.
c-----------------------------------------------------------------------
c  OUTPUT
c  ------
c    * pm00000
c    * akr0000
c-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      include "../pmd/params_unit.h"
c-----max # of atoms
      integer,parameter::nmax=1000000
c.....2*pi
      real(8),parameter:: pi   = 3.14159265358979d0
      real(8),parameter:: twopi= 2d0 *pi
c.....Poisson ratio
      real(8),parameter:: anu = 0.25d0
c.....outfile
      character:: cfout*10='out.mkconf'
      real(8),allocatable:: ua(:,:)
      real(8):: ra(3,nmax)
      real(8):: h(3,3,0:1)
      real(8):: va(3,nmax)
      real(8):: tag(nmax)
      real(8):: eki(nmax),epi(nmax),strs(3,3,nmax)
      real(8):: xi(3),du(3)
      real(8):: radius,theta
      integer:: ifshffl,iangle,nuc(3)
      real(8),external:: get_theta

      open(10,file=cfout,status='replace')

      small=1d-7

 1    write(6,'(a)') ' Input nx ny nz:'
      read(5,*) nuc(1:3)
 2    write(6,'(a)') ' Input cylinder radius in Angstrom:'
      read(5,*) radius

c-----be: equilibrium bond length: 2.35 Ang
      be= 2.35d0
c.....al: length of an edge of conventional cubic cell
      al= 4d0*be/sqrt(3d0)


 20   write(6,'(a)') ' Choose (0) bulk, (1) screw, (2) 60 degree,'
     &     //' (3) 30 degree, or (4) 41 degree:'
      read(5,*) iangle
c.....Burgers vector
      bv= al *sqrt(2d0)/2
c.....screw or 60 dislocation
      if( iangle .eq. 1 ) then  ! screw
        bvx= 0d0
        bvz= bv
      else if( iangle .eq. 2 ) then ! 60 degree
        bvx= bv *sqrt(3d0)/2
        bvz= bv /2
      else if( iangle .eq. 3 ) then
        bvx= bv /2
        bvz= bv *sqrt(3d0)/2
c        bvx= 0d0
c        bvz= 0d0
      else if( iangle .eq. 4 ) then
        theta= 41d0/180d0 *pi
        bvx= bv *sin(theta)
        bvz= bv *cos(theta)
c        bvx= 0d0
c        bvz= 0d0
      else if( iangle .eq. 0 ) then
        bvx= 0d0
        bvz= 0d0
      else
        write(6,'(a)') ' [Error] iangle must be (1) screw or (2) 60 !!'
        goto 20
      endif

c.....Simulation box size
      h(1:3,1:3,0:1)= 0d0
      if( iangle.eq.3 ) then ! 30-degree
        h(1,1,0)= nuc(1) *al *sqrt(2d0)
        h(2,2,0)= nuc(2) *al *sqrt(3d0)
        h(3,3,0)= nuc(3) *al *sqrt(6d0)/2
      elseif( iangle.eq.4 ) then ! 41-degree
        h(1,1,0)= nuc(1) *al *3.2403703492d0
        h(2,2,0)= nuc(2) *al *sqrt(3d0)
        h(3,3,0)= nuc(3) *al *1.87082869339d0
      else
        h(1,1,0)= nuc(1) *al *sqrt(6d0)/2
        h(2,2,0)= nuc(2) *al *sqrt(3d0)
        h(3,3,0)= nuc(3) *al *sqrt(2d0)/2
      endif
      
      if( h(1,1,0).lt.2*radius .or. h(2,2,0).lt.2*radius )then
        write(6,'(a)') ' [Error] h(1,1,0) or h(2,2,0) is too small !!!'
        write(6,'(a,es12.4)') '  h(1,1,0) =',h(1,1,0)
        write(6,'(a,es12.4)') '  h(2,2,0) =',h(2,2,0)
        write(6,'(a,es12.4)') '  2*radius =',2*radius
        goto 1
      endif

 10   write(6,'(a)') ' Choose (1) shuffle or (2) glide:'
      read(5,*) ifshffl
      if( iangle.eq.3 ) then ! 30-degree
        if( ifshffl.eq.1 ) then ! shuffl
c.....Center position in normalized unit
          xc= (1d0/8 +nuc(1)/2)/nuc(1)
          yc= (5.5d0/12 +nuc(2)/2)/nuc(2)
        else if( ifshffl.eq.2 ) then ! glide
c.....Center position in normalized unit
          xc= (1d0/8 +nuc(1)/2)/nuc(1)
          yc= (3.5d0/12 +nuc(2)/2)/nuc(2)
        else
          write(6,'(a)') ' [Error] ifshffl must be (1) shuffle'
     &         //' or (2) glide !!!'
          goto 10
        endif
      elseif( iangle.eq.4 ) then ! 41-degree
        if( ifshffl.eq.1 ) then ! shuffl
c.....Center position in normalized unit
          xc= (0d0 +nuc(1)/2)/nuc(1)
          yc= (5.5d0/12 +nuc(2)/2)/nuc(2)
        else if( ifshffl.eq.2 ) then ! glide
c.....Center position in normalized unit
          xc= (0d0/8 +nuc(1)/2)/nuc(1)
          yc= (3.5d0/12 +nuc(2)/2)/nuc(2)
        else
          write(6,'(a)') ' [Error] ifshffl must be (1) shuffle'
     &         //' or (2) glide !!!'
          goto 10
        endif
      else
        if( ifshffl.eq.1 ) then ! shuffl
c.....Center position in normalized unit
          xc= (5d0/12 +nuc(1)/2)/nuc(1)
          yc= (11d0/24 +nuc(2)/2)/nuc(2)
        else if( ifshffl.eq.2 ) then ! glide
c.....Center position in normalized unit
          xc= (0.5d0 +nuc(1)/2)/nuc(1)
          yc= (15d0/24 +nuc(2)/2)/nuc(2)
        else
          write(6,'(a)') ' [Error] ifshffl must be (1) shuffle'
     &         //' or (2) glide !!!'
          goto 10
        endif
      endif

      write(6,'(a,3es12.4)') ' bond length=',be
      write(6,'(a,3es12.4)') ' cubic diamond edge length=',al
      write(10,'(a,3i5)') ' nx,ny,nz=',nuc(1:3)
      write(10,'(a,i5)') ' ifshffl=',ifshffl
      write(10,'(a,i5)') ' iangle =',iangle
      write(10,'(a,3es12.4)') ' bond length=',be
      write(10,'(a,3es12.4)') ' cubic diamond edge length=',al


      write(6,'(a)') ' h-matrix:'
      write(6,'(2x,3es12.4)') h(1,1:3,0)
      write(6,'(2x,3es12.4)') h(2,1:3,0)
      write(6,'(2x,3es12.4)') h(3,1:3,0)
      write(10,'(a)') ' h-matrix:'
      write(10,'(2x,3es12.4)') h(1,1:3,0)
      write(10,'(2x,3es12.4)') h(2,1:3,0)
      write(10,'(2x,3es12.4)') h(3,1:3,0)

      if( iangle.eq.3 ) then ! 30-degree
        nacell= 24
        allocate(ua(3,nacell))
        ua(1:3,1 )=(/2.50E-01, 0.000000E+00, 2.50000E-01 *2/)
        ua(1:3,2 )=(/2.50E-01, 2.500000E-01, 2.50000E-01 *2/)
        ua(1:3,3 )=(/5.00E-01, 0.000000E+00, 0.00000E+00 *2/)
        ua(1:3,4 )=(/5.00E-01, 2.500000E-01, 0.00000E+00 *2/)
        ua(1:3,5 )=(/0.00E+00, 0.000000E+00, 0.00000E+00 *2/)
        ua(1:3,6 )=(/0.00E+00, 3.333333E-01, 1.66667E-01 *2/)
        ua(1:3,7 )=(/0.00E+00, 2.500000E-01, 0.00000E+00 *2/)
        ua(1:3,8 )=(/0.00E+00, 5.833333E-01, 1.66667E-01 *2/)
        ua(1:3,9 )=(/7.50E-01, 0.000000E+00, 2.50000E-01 *2/)
        ua(1:3,10)=(/7.50E-01, 2.500000E-01, 2.50000E-01 *2/)
        ua(1:3,11)=(/2.50E-01, 3.333333E-01, 4.16667E-01 *2/)
        ua(1:3,12)=(/7.50E-01, 3.333333E-01, 4.16667E-01 *2/)
        ua(1:3,13)=(/2.50E-01, 5.833333E-01, 4.16667E-01 *2/)
        ua(1:3,14)=(/7.50E-01, 5.833333E-01, 4.16667E-01 *2/)
        ua(1:3,15)=(/5.00E-01, 3.333333E-01, 1.66667E-01 *2/)
        ua(1:3,16)=(/2.50E-01, 6.666667E-01, 8.33333E-02 *2/)
        ua(1:3,17)=(/5.00E-01, 6.666667E-01, 3.33333E-01 *2/)
        ua(1:3,18)=(/7.50E-01, 6.666667E-01, 8.33333E-02 *2/)
        ua(1:3,19)=(/5.00E-01, 5.833333E-01, 1.66667E-01 *2/)
        ua(1:3,20)=(/2.50E-01, 9.166667E-01, 8.33333E-02 *2/)
        ua(1:3,21)=(/5.00E-01, 9.166667E-01, 3.33333E-01 *2/)
        ua(1:3,22)=(/7.50E-01, 9.166667E-01, 8.33333E-02 *2/)
        ua(1:3,23)=(/0.00E+00, 6.666667E-01, 3.33333E-01 *2/)
        ua(1:3,24)=(/0.00E+00, 9.166667E-01, 3.33333E-01 *2/)
      elseif( iangle.eq.4 ) then ! 41-degree
        nacell= 84
        allocate(ua(3,nacell))
        ua(1:3,1 )=(/ 7.142857142E-02, 0.000000E+00, 6.42857143E-01/)
        ua(1:3,2 )=(/ 2.142857142E-01, 0.000000E+00, 9.28571429E-01/)
        ua(1:3,3 )=(/ 7.142857142E-02, 2.500000E-01, 6.42857143E-01/)
        ua(1:3,4 )=(/ 2.142857142E-01, 2.500000E-01, 9.28571429E-01/)
        ua(1:3,5 )=(/ 1.428571428E-01, 0.000000E+00, 2.85714286E-01/)
        ua(1:3,6 )=(/ 3.571428571E-01, 0.000000E+00, 2.14285714E-01/)
        ua(1:3,7 )=(/ 1.428571428E-01, 2.500000E-01, 2.85714286E-01/)
        ua(1:3,8 )=(/ 3.571428571E-01, 2.500000E-01, 2.14285714E-01/)
        ua(1:3,9 )=(/ 0.000000000E+00, 0.000000E+00, 0.00000000E+00/)
        ua(1:3,10)=(/ 1.190476190E-01, 3.333333E-01, 7.14285714E-02/)
        ua(1:3,11)=(/ 0.000000000E+00, 2.500000E-01, 0.00000000E+00/)
        ua(1:3,12)=(/ 1.190476190E-01, 5.833333E-01, 7.14285714E-02/)
        ua(1:3,13)=(/ 4.285714285E-01, 0.000000E+00, 8.57142857E-01/)
        ua(1:3,14)=(/ 5.000000000E-01, 0.000000E+00, 5.00000000E-01/)
        ua(1:3,15)=(/ 6.428571428E-01, 0.000000E+00, 7.85714286E-01/)
        ua(1:3,16)=(/ 4.285714285E-01, 2.500000E-01, 8.57142857E-01/)
        ua(1:3,17)=(/ 5.000000000E-01, 2.500000E-01, 5.00000000E-01/)
        ua(1:3,18)=(/ 6.428571428E-01, 2.500000E-01, 7.85714286E-01/)
        ua(1:3,19)=(/ 2.857142857E-01, 0.000000E+00, 5.71428571E-01/)
        ua(1:3,20)=(/ 1.904761904E-01, 3.333333E-01, 7.14285714E-01/)
        ua(1:3,21)=(/ 2.619047619E-01, 3.333333E-01, 3.57142857E-01/)
        ua(1:3,22)=(/ 4.047619047E-01, 3.333333E-01, 6.42857143E-01/)
        ua(1:3,23)=(/ 2.857142857E-01, 2.500000E-01, 5.71428571E-01/)
        ua(1:3,24)=(/ 1.904761904E-01, 5.833333E-01, 7.14285714E-01/)
        ua(1:3,25)=(/ 2.619047619E-01, 5.833333E-01, 3.57142857E-01/)
        ua(1:3,26)=(/ 4.047619047E-01, 5.833333E-01, 6.42857143E-01/)
        ua(1:3,27)=(/ 4.761904761E-02, 3.333333E-01, 4.28571429E-01/)
        ua(1:3,28)=(/ 2.380952380E-02, 6.666667E-01, 2.14285714E-01/)
        ua(1:3,29)=(/ 1.666666666E-01, 6.666667E-01, 5.00000000E-01/)
        ua(1:3,30)=(/ 4.761904761E-02, 5.833333E-01, 4.28571429E-01/)
        ua(1:3,31)=(/ 2.380952380E-02, 9.166667E-01, 2.14285714E-01/)
        ua(1:3,32)=(/ 1.666666666E-01, 9.166667E-01, 5.00000000E-01/)
        ua(1:3,33)=(/ 5.476190476E-01, 3.333333E-01, 9.28571429E-01/)
        ua(1:3,34)=(/ 5.476190476E-01, 5.833333E-01, 9.28571429E-01/)
        ua(1:3,35)=(/ 3.095238095E-01, 6.666667E-01, 7.85714286E-01/)
        ua(1:3,36)=(/ 3.095238095E-01, 9.166667E-01, 7.85714286E-01/)
        ua(1:3,37)=(/ 9.523809523E-02, 6.666667E-01, 8.57142857E-01/)
        ua(1:3,38)=(/ 9.523809523E-02, 9.166667E-01, 8.57142857E-01/)
        ua(1:3,39)=(/ 5.714285714E-01, 0.000000E+00, 1.42857143E-01/)
        ua(1:3,40)=(/ 7.857142857E-01, 0.000000E+00, 7.14285714E-02/)
        ua(1:3,41)=(/ 5.714285714E-01, 2.500000E-01, 1.42857143E-01/)
        ua(1:3,42)=(/ 7.857142857E-01, 2.500000E-01, 7.14285714E-02/)
        ua(1:3,43)=(/ 3.333333333E-01, 3.333333E-01, 0.00000000E+00/)
        ua(1:3,44)=(/ 3.333333333E-01, 5.833333E-01, 0.00000000E+00/)
        ua(1:3,45)=(/ 8.571428571E-01, 0.000000E+00, 7.14285714E-01/)
        ua(1:3,46)=(/ 9.285714285E-01, 0.000000E+00, 3.57142857E-01/)
        ua(1:3,47)=(/ 8.571428571E-01, 2.500000E-01, 7.14285714E-01/)
        ua(1:3,48)=(/ 9.285714285E-01, 2.500000E-01, 3.57142857E-01/)
        ua(1:3,49)=(/ 7.142857142E-01, 0.000000E+00, 4.28571429E-01/)
        ua(1:3,50)=(/ 6.190476190E-01, 3.333333E-01, 5.71428571E-01/)
        ua(1:3,51)=(/ 6.904761904E-01, 3.333333E-01, 2.14285714E-01/)
        ua(1:3,52)=(/ 8.333333333E-01, 3.333333E-01, 5.00000000E-01/)
        ua(1:3,53)=(/ 7.142857142E-01, 2.500000E-01, 4.28571429E-01/)
        ua(1:3,54)=(/ 6.190476190E-01, 5.833333E-01, 5.71428571E-01/)
        ua(1:3,55)=(/ 6.904761904E-01, 5.833333E-01, 2.14285714E-01/)
        ua(1:3,56)=(/ 8.333333333E-01, 5.833333E-01, 5.00000000E-01/)
        ua(1:3,57)=(/ 4.761904761E-01, 3.333333E-01, 2.85714286E-01/)
        ua(1:3,58)=(/ 3.809523809E-01, 6.666667E-01, 4.28571429E-01/)
        ua(1:3,59)=(/ 4.523809523E-01, 6.666667E-01, 7.14285714E-02/)
        ua(1:3,60)=(/ 5.952380952E-01, 6.666667E-01, 3.57142857E-01/)
        ua(1:3,61)=(/ 4.761904761E-01, 5.833333E-01, 2.85714286E-01/)
        ua(1:3,62)=(/ 3.809523809E-01, 9.166667E-01, 4.28571429E-01/)
        ua(1:3,63)=(/ 4.523809523E-01, 9.166667E-01, 7.14285714E-02/)
        ua(1:3,64)=(/ 5.952380952E-01, 9.166667E-01, 3.57142857E-01/)
        ua(1:3,65)=(/ 2.380952380E-01, 6.666667E-01, 1.42857143E-01/)
        ua(1:3,66)=(/ 2.380952380E-01, 9.166667E-01, 1.42857143E-01/)
        ua(1:3,67)=(/ 9.761904761E-01, 3.333333E-01, 7.85714286E-01/)
        ua(1:3,68)=(/ 9.761904761E-01, 5.833333E-01, 7.85714286E-01/)
        ua(1:3,69)=(/ 7.619047619E-01, 3.333333E-01, 8.57142857E-01/)
        ua(1:3,70)=(/ 7.380952380E-01, 6.666667E-01, 6.42857143E-01/)
        ua(1:3,71)=(/ 8.809523809E-01, 6.666667E-01, 9.28571429E-01/)
        ua(1:3,72)=(/ 7.619047619E-01, 5.833333E-01, 8.57142857E-01/)
        ua(1:3,73)=(/ 7.380952380E-01, 9.166667E-01, 6.42857143E-01/)
        ua(1:3,74)=(/ 8.809523809E-01, 9.166667E-01, 9.28571429E-01/)
        ua(1:3,75)=(/ 5.238095238E-01, 6.666667E-01, 7.14285714E-01/)
        ua(1:3,76)=(/ 5.238095238E-01, 9.166667E-01, 7.14285714E-01/)
        ua(1:3,77)=(/ 9.047619047E-01, 3.333333E-01, 1.42857143E-01/)
        ua(1:3,78)=(/ 8.095238095E-01, 6.666667E-01, 2.85714286E-01/)
        ua(1:3,79)=(/ 9.047619047E-01, 5.833333E-01, 1.42857143E-01/)
        ua(1:3,80)=(/ 8.095238095E-01, 9.166667E-01, 2.85714286E-01/)
        ua(1:3,81)=(/ 6.666666666E-01, 6.666667E-01, 0.00000000E+00/)
        ua(1:3,82)=(/ 6.666666666E-01, 9.166667E-01, 0.00000000E+00/)
        ua(1:3,83)=(/ 9.523809523E-01, 6.666667E-01, 5.71428571E-01/)
        ua(1:3,84)=(/ 9.523809523E-01, 9.166667E-01, 5.71428571E-01/)
      else
        nacell= 12
        allocate(ua(3,nacell))
c-----Unit cell of diamond structure with z-axis towrads (111): type2
        ua(1:3, 1)= (/ 0.00d0, 0.00d0, 0.00d0 /)
        ua(1:3, 2)= (/ 0.00d0, 1d0/4,  0.00d0 /)
        ua(1:3, 3)= (/ 1d0/6,  1d0/3,  1d0/2  /)
        ua(1:3, 4)= (/ 1d0/6,  7d0/12, 1d0/2  /)
        ua(1:3, 5)= (/ 1d0/3,  2d0/3,  0.00d0 /)
        ua(1:3, 6)= (/ 1d0/3, 11d0/12, 0.00d0 /)
        ua(1:3, 7)= (/ 1d0/2,  0.00d0, 1d0/2  /)
        ua(1:3, 8)= (/ 1d0/2,  1d0/4,  1d0/2  /)
        ua(1:3, 9)= (/ 2d0/3,  1d0/3,  0.00d0 /)
        ua(1:3,10)= (/ 2d0/3,  7d0/12, 0.00d0 /)
        ua(1:3,11)= (/ 5d0/6,  2d0/3,  1d0/2  /)
        ua(1:3,12)= (/ 5d0/6, 11d0/12, 1d0/2  /)
      endif


c.....Normalized Burgers vector in x- and z-direction
      bvxn= bvx /h(1,1,0)
      bvzn= bvz /h(3,3,0)
      tan04pi= tan(0.4d0*pi)
      write(6,'(a,es12.4)') ' bv   =',bv
      write(6,'(a,es12.4)') ' bvx  =',bvx
      write(6,'(a,es12.4)') ' bvz  =',bvz
      write(6,'(a,es12.4)') ' bvxn =',bvxn
      write(6,'(a,es12.4)') ' bvzn =',bvzn
      write(10,'(a,es12.4)') ' bv   =',bv
      write(10,'(a,es12.4)') ' bvx  =',bvx
      write(10,'(a,es12.4)') ' bvz  =',bvz
      write(10,'(a,es12.4)') ' bvxn =',bvxn
      write(10,'(a,es12.4)') ' bvzn =',bvzn

      inc=0 
      xmin= 1d10
      xmax= 0d0
      do ix=0,nuc(1)-1
        do iy=0,nuc(2)-1
          do iz=0,nuc(3)-1
            do m=1,nacell
              xi(1)= (ua(1,m)+dble(ix))/nuc(1) +small
              xi(2)= (ua(2,m)+dble(iy))/nuc(2) +small
              xi(3)= (ua(3,m)+dble(iz))/nuc(3) +small
              x= (xi(1)-xc)*h(1,1,0)
              y= (xi(2)-yc)*h(2,2,0)
              r= sqrt(x**2+y**2)
c.....Cut atoms outside the radius
              if( r.gt.radius ) cycle
c.....Introduce dislocation by subtracting one atomic layer from
c.....lower half of the system.
              if( xi(2).lt.yc .and.
c     &             xi(1).gt.xc-bvxn/2 .and. xi(1).lt.xc+bvxn/2 ) cycle
     &             xi(1).gt.xc .and. xi(1).le.xc+bvxn ) cycle
c.....Displace only bottom half atoms to remove the space where atoms
c.....are eliminated. (See, my note on 2014-01-30)
              if( xi(2).lt.yc ) then
                du(1:3)= 0d0
                dx= xi(1)-xc
                if( dx.gt. 0.5d0 ) then dx=dx -1d0
                if( dx.le.-0.5d0 ) then dx=dx +1d0
                if( dx.gt.0d0 ) du(1)= +bvxn*dx -bvxn/2
                if( dx.le.0d0 ) du(1)= +bvxn*dx +bvxn/2
                dy= (xi(2)-yc) *h(2,2,0)
                du(1)= du(1) *atan(dy/2/bv*(-tan04pi))/(pi/2)
                xi(1)= xi(1) +du(1)
              endif
c.....Displace atoms according to dislocation theory,
c.....only atoms far from the core
              du(1:3)= 0d0
              theta= get_theta(x,y,r)
              if( theta .ge. twopi*0.75d0 ) theta=theta-twopi
              du(3)= bvzn/twopi *theta
     &             *atan(r/2/bv*(-tan04pi))/(pi/2)
              xi(3)= xi(3) +du(3)
c.....Store xmin and xmax
              xmin= min(xmin,xi(1))
              xmax= max(xmax,xi(1))
c.....Add this atom
              inc=inc+1
              if(inc.gt.nmax)then
                write(*,*)'Error inc>nmax',inc,nmax
                stop
              endif
              call pbc(xi(1))
              call pbc(xi(2))
              call pbc(xi(3))
              ra(1:3,inc)= xi(1:3)
              is= 1
              ifmv= 1
              if( r.gt.radius-al ) ifmv= 0
              tag(inc)= 1d0*is +0.1d0*ifmv +1d-14*inc
            enddo
          enddo
        enddo
      enddo
      write(6,'(a,i10)') " natm=",inc
      write(10,'(a,i10)') " natm=",inc

      if( iangle .gt. 1 ) then ! if the disl has edge component
c.....Compress bvx/2 along x
        do i=1,inc
          ra(1,i)=ra(1,i) -bvxn/4*(ra(1,i)-xc)/(xmax-xc)
        enddo
      endif

      va(1:3,1:inc)= 0d0

      call write_pmd0_ascii(15,'pmd00000','replace',inc,tag
     &     ,ra,va,h,al,eki,epi,strs)
      
c-----output 'kvs000' for KVS visualization
      open(15,file='akr0000',form='formatted',status='replace')
      write(15,'(es15.7)') al
      write(15,'(3es11.3)') ((h(ia,ib,l)/al,ia=1,3),ib=1,3)
      write(15,'(i10,3i5)')inc,3,0,0
      do i=1,inc
        write(15,'(i3,6es11.3)') int(tag(i)),ra(1:3,i),va(1:3,i)
      enddo
      close(15)

      close(10)
      
      end program mkconf_Si_disl_scew_60
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
      subroutine pbc(x)
      implicit none
      real(8),intent(inout):: x

      if( x.lt.0d0 ) x= x+1d0
      if( x.ge.1d0 ) x= x-1d0
      return
      end subroutine pbc
c=======================================================================
      function get_theta(x,y,r)
c
c Returns [0:2*pi]
c
      implicit none
      real(8),intent(in):: x,y,r
      real(8):: get_theta
      real(8):: asn,acs
      real(8),parameter:: pi = 3.14159265358979d0

      get_theta= 0d0
      asn= asin(y/r)
      if( x.ge.0d0 .and. y.ge.0d0 ) then
        get_theta= asn
      else if( x.lt.0d0 .and. y.ge.0d0 ) then
        get_theta= pi -asn
      else if( x.lt.0d0 .and. y.lt.0d0 ) then
        get_theta= pi -asn
      else if( x.ge.0d0 .and. y.lt.0d0 ) then
        get_theta= 2d0*pi +asn
      endif

      return
      end function get_theta
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make mkconf_Si_disl"
c     End:
