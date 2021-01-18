module element
  implicit none
  save

  integer,parameter:: nelem = 118
  
  type atom
    character(len=3):: symbol
    real(8):: mass = 1.0d0
    real(8):: val = 0d0  ! valence in e
  end type atom
  type(atom):: elmts(nelem)

contains
!=======================================================================
  subroutine init_element()

    elmts(  1)%symbol='H'  ; elmts(  1)%mass= 1.008d0  ; elmts(  1)%val= 1d0
    elmts(  2)%symbol='He' ; elmts(  2)%mass= 4.0026d0 ; elmts(  2)%val= 0d0
    elmts(  3)%symbol='Li' ; elmts(  3)%mass= 6.94d0   ; elmts(  3)%val= 1d0
    elmts(  4)%symbol='Be' ; elmts(  4)%mass= 9.0122d0 ; elmts(  4)%val= 2d0
    elmts(  5)%symbol='B'  ; elmts(  5)%mass= 10.81d0  ; elmts(  5)%val= 3d0
    elmts(  6)%symbol='C'  ; elmts(  6)%mass= 12.011d0 ; elmts(  6)%val= 4d0
    elmts(  7)%symbol='N'  ; elmts(  7)%mass= 14.007d0 ; elmts(  7)%val= 5d0
    elmts(  8)%symbol='O'  ; elmts(  8)%mass= 15.999d0 ; elmts(  8)%val= 6d0
    elmts(  9)%symbol='F'  ; elmts(  9)%mass= 18.998d0 ; elmts(  9)%val= 7d0
    elmts( 10)%symbol='Ne' ; elmts( 10)%mass= 20.180d0 ; elmts( 10)%val= 0d0
    elmts( 11)%symbol='Na' ; elmts( 11)%mass= 22.990d0 ; elmts( 11)%val= 1d0
    elmts( 12)%symbol='Mg' ; elmts( 12)%mass= 24.305d0 ; elmts( 12)%val= 2d0
    elmts( 13)%symbol='Al' ; elmts( 13)%mass= 26.982d0 ; elmts( 13)%val= 3d0
    elmts( 14)%symbol='Si' ; elmts( 14)%mass= 28.085d0 ; elmts( 14)%val= 4d0
    elmts( 15)%symbol='P'  ; elmts( 15)%mass= 30.974d0 ; elmts( 15)%val= 5d0
    elmts( 16)%symbol='S'  ; elmts( 16)%mass= 32.06d0  ; elmts( 16)%val= 6d0
    elmts( 17)%symbol='Cl' ; elmts( 17)%mass= 35.45d0  ; elmts( 17)%val= 7d0
    elmts( 18)%symbol='Ar' ; elmts( 18)%mass= 39.948d0 ; elmts( 18)%val= 0d0
    elmts( 19)%symbol='K'  ; elmts( 19)%mass= 39.098d0 ; elmts( 19)%val= 1d0
    elmts( 20)%symbol='Ca' ; elmts( 20)%mass= 40.078d0 ; elmts( 20)%val= 2d0
    elmts( 21)%symbol='Sc' ; elmts( 21)%mass= 44.956d0 ; elmts( 21)%val= 3d0
    elmts( 22)%symbol='Ti' ; elmts( 22)%mass= 47.867d0 ; elmts( 22)%val= 4d0
    elmts( 23)%symbol='V'  ; elmts( 23)%mass= 50.942d0 ; elmts( 23)%val= 5d0
    elmts( 24)%symbol='Cr' ; elmts( 24)%mass= 51.996d0 ; elmts( 24)%val= 6d0
    elmts( 25)%symbol='Mn' ; elmts( 25)%mass= 54.938d0 ; elmts( 25)%val= 7d0
    elmts( 26)%symbol='Fe' ; elmts( 26)%mass= 55.845d0 ; elmts( 26)%val= 8d0
    elmts( 27)%symbol='Co' ; elmts( 27)%mass= 58.933d0 ; elmts( 27)%val= 8d0
    elmts( 28)%symbol='Ni' ; elmts( 28)%mass= 58.693d0 ; elmts( 28)%val= 8d0
    elmts( 29)%symbol='Cu' ; elmts( 29)%mass= 63.546d0 ; elmts( 29)%val= 9d0
    elmts( 30)%symbol='Zn' ; elmts( 30)%mass= 65.38d0  ; elmts( 30)%val=10d0
    elmts( 31)%symbol='Ga' ; elmts( 31)%mass= 69.723d0 ; elmts( 31)%val= 3d0
    elmts( 32)%symbol='Ge' ; elmts( 32)%mass= 72.631d0 ; elmts( 32)%val= 4d0
    elmts( 33)%symbol='As' ; elmts( 33)%mass= 74.922d0 ; elmts( 33)%val= 5d0
    elmts( 34)%symbol='Se' ; elmts( 34)%mass= 78.972d0 ; elmts( 34)%val= 6d0
    elmts( 35)%symbol='Br' ; elmts( 35)%mass= 79.904d0 ; elmts( 35)%val= 7d0
    elmts( 36)%symbol='Kr' ; elmts( 36)%mass= 84.798d0 ; elmts( 36)%val= 0d0
    elmts( 37)%symbol='Rb' ; elmts( 37)%mass= 85.468d0 ; elmts( 37)%val= 1d0
    elmts( 38)%symbol='Sr' ; elmts( 38)%mass= 87.62d0  ; elmts( 38)%val= 2d0
    elmts( 39)%symbol='Y'  ; elmts( 39)%mass= 88.906d0 ; elmts( 39)%val= 3d0
    elmts( 40)%symbol='Zr' ; elmts( 40)%mass= 91.224d0 ; elmts( 40)%val= 4d0
    elmts( 41)%symbol='Nb' ; elmts( 41)%mass= 92.906d0 ; elmts( 41)%val= 5d0
    elmts( 42)%symbol='Mo' ; elmts( 42)%mass= 95.95d0  ; elmts( 42)%val= 6d0
    elmts( 43)%symbol='Tc' ; elmts( 43)%mass= 98.907d0 ; elmts( 43)%val= 7d0
    elmts( 44)%symbol='Ru' ; elmts( 44)%mass= 101.07d0 ; elmts( 44)%val= 8d0
    elmts( 45)%symbol='Rh' ; elmts( 45)%mass= 102.90d0 ; elmts( 45)%val= 8d0
    elmts( 46)%symbol='Pd' ; elmts( 46)%mass= 106.42d0 ; elmts( 46)%val= 8d0
    elmts( 47)%symbol='Ag' ; elmts( 47)%mass= 107.86d0 ; elmts( 47)%val= 9d0
    elmts( 48)%symbol='Cd' ; elmts( 48)%mass= 112.41d0 ; elmts( 48)%val=10d0
    elmts( 49)%symbol='In' ; elmts( 49)%mass= 114.81d0 ; elmts( 49)%val= 3d0
    elmts( 50)%symbol='Sn' ; elmts( 50)%mass= 118.71d0 ; elmts( 50)%val= 4d0
    elmts( 51)%symbol='Sb' ; elmts( 51)%mass= 121.76d0 ; elmts( 51)%val= 5d0
    elmts( 52)%symbol='Te' ; elmts( 52)%mass= 127.60d0 ; elmts( 52)%val= 6d0
    elmts( 53)%symbol='I'  ; elmts( 53)%mass= 126.90d0 ; elmts( 53)%val= 7d0
    elmts( 54)%symbol='Xe' ; elmts( 54)%mass= 131.29d0 ; elmts( 54)%val= 0d0
    elmts( 55)%symbol='Cs' ; elmts( 55)%mass= 132.90d0 ; elmts( 55)%val= 1d0
    elmts( 56)%symbol='Ba' ; elmts( 56)%mass= 137.32d0 ; elmts( 56)%val= 2d0
    elmts( 57)%symbol='La' ; elmts( 57)%mass= 138.90d0 ; elmts( 57)%val= 3d0
    elmts( 58)%symbol='Ce' ; elmts( 58)%mass= 140.11d0 ; elmts( 58)%val= 3d0
    elmts( 59)%symbol='Pr' ; elmts( 59)%mass= 140.90d0 ; elmts( 59)%val= 3d0
    elmts( 60)%symbol='Nd' ; elmts( 60)%mass= 144.24d0 ; elmts( 60)%val= 3d0
    elmts( 61)%symbol='Pm' ; elmts( 61)%mass= 144.91d0 ; elmts( 61)%val= 3d0
    elmts( 62)%symbol='Sm' ; elmts( 62)%mass= 150.36d0 ; elmts( 62)%val= 3d0
    elmts( 63)%symbol='Eu' ; elmts( 63)%mass= 151.96d0 ; elmts( 63)%val= 3d0
    elmts( 64)%symbol='Gd' ; elmts( 64)%mass= 157.25d0 ; elmts( 64)%val= 3d0
    elmts( 65)%symbol='Tb' ; elmts( 65)%mass= 158.92d0 ; elmts( 65)%val= 3d0
    elmts( 66)%symbol='Dy' ; elmts( 66)%mass= 162.50d0 ; elmts( 66)%val= 3d0
    elmts( 67)%symbol='Ho' ; elmts( 67)%mass= 164.93d0 ; elmts( 67)%val= 3d0
    elmts( 68)%symbol='Er' ; elmts( 68)%mass= 167.25d0 ; elmts( 68)%val= 3d0
    elmts( 69)%symbol='Tm' ; elmts( 69)%mass= 168.93d0 ; elmts( 69)%val= 3d0
    elmts( 70)%symbol='Yb' ; elmts( 70)%mass= 173.05d0 ; elmts( 70)%val= 3d0
    elmts( 71)%symbol='Lu' ; elmts( 71)%mass= 174.96d0 ; elmts( 71)%val= 3d0
    elmts( 72)%symbol='Hf' ; elmts( 72)%mass= 178.49d0 ; elmts( 72)%val= 4d0
    elmts( 73)%symbol='Ta' ; elmts( 73)%mass= 180.94d0 ; elmts( 73)%val= 5d0
    elmts( 74)%symbol='W'  ; elmts( 74)%mass= 183.84d0 ; elmts( 74)%val= 6d0
    elmts( 75)%symbol='Re' ; elmts( 75)%mass= 186.20d0 ; elmts( 75)%val= 7d0
    elmts( 76)%symbol='Os' ; elmts( 76)%mass= 190.23d0 ; elmts( 76)%val= 8d0
    elmts( 77)%symbol='Ir' ; elmts( 77)%mass= 192.21d0 ; elmts( 77)%val= 8d0
    elmts( 78)%symbol='Pt' ; elmts( 78)%mass= 195.08d0 ; elmts( 78)%val= 8d0
    elmts( 79)%symbol='Au' ; elmts( 79)%mass= 196.96d0 ; elmts( 79)%val= 9d0
    elmts( 80)%symbol='Hg' ; elmts( 80)%mass= 200.59d0 ; elmts( 80)%val=10d0
    elmts( 81)%symbol='Tl' ; elmts( 81)%mass= 204.38d0 ; elmts( 81)%val= 3d0
    elmts( 82)%symbol='Pb' ; elmts( 82)%mass= 207.20d0 ; elmts( 82)%val= 4d0
    elmts( 83)%symbol='Bi' ; elmts( 83)%mass= 208.98d0 ; elmts( 83)%val= 5d0
    elmts( 84)%symbol='Po' ; elmts( 84)%mass= 208.98d0 ; elmts( 84)%val= 6d0
    elmts( 85)%symbol='At' ; elmts( 85)%mass= 209.98d0 ; elmts( 85)%val= 7d0
    elmts( 86)%symbol='Rn' ; elmts( 86)%mass= 222.01d0 ; elmts( 86)%val= 8d0
    elmts( 87)%symbol='Fr' ; elmts( 87)%mass= 223.20d0 ; elmts( 87)%val= 1d0 
    elmts( 88)%symbol='Ra' ; elmts( 88)%mass= 226.02d0 ; elmts( 88)%val= 1d0 
    elmts( 89)%symbol='Ac' ; elmts( 89)%mass= 227.02d0 ; elmts( 89)%val= 1d0 
    elmts( 90)%symbol='Th' ; elmts( 90)%mass= 232.03d0 ; elmts( 90)%val= 1d0 
    elmts( 91)%symbol='Pa' ; elmts( 91)%mass= 231.03d0 ; elmts( 91)%val= 1d0 
    elmts( 92)%symbol='U'  ; elmts( 92)%mass= 238.02d0 ; elmts( 92)%val= 1d0 
    elmts( 93)%symbol='Np' ; elmts( 93)%mass= 237.04d0 ; elmts( 93)%val= 1d0 
    elmts( 94)%symbol='Pu' ; elmts( 94)%mass= 244.06d0 ; elmts( 94)%val= 1d0 
    elmts( 95)%symbol='Am' ; elmts( 95)%mass= 243.06d0 ; elmts( 95)%val= 1d0 
    elmts( 96)%symbol='Cm' ; elmts( 96)%mass= 247.00d0 ; elmts( 96)%val= 1d0 
    elmts( 97)%symbol='Bk' ; elmts( 97)%mass= 247.07d0 ; elmts( 97)%val= 1d0 
    elmts( 98)%symbol='Cf' ; elmts( 98)%mass= 251.08d0 ; elmts( 98)%val= 1d0 
    elmts( 99)%symbol='Es' ; elmts( 99)%mass= 254.00d0 ; elmts( 99)%val= 1d0 
    elmts(100)%symbol='Fm' ; elmts(100)%mass= 257.09d0 ; elmts(100)%val= 1d0 
    elmts(101)%symbol='Md' ; elmts(101)%mass= 258.10d0 ; elmts(101)%val= 1d0 
    elmts(102)%symbol='No' ; elmts(102)%mass= 259.10d0 ; elmts(102)%val= 1d0 
    elmts(103)%symbol='Lr' ; elmts(103)%mass= 262.00d0 ; elmts(103)%val= 1d0 
    elmts(104)%symbol='Rf' ; elmts(104)%mass= 261.00d0 ; elmts(104)%val= 1d0 
    elmts(105)%symbol='Db' ; elmts(105)%mass= 262.00d0 ; elmts(105)%val= 1d0 
    elmts(106)%symbol='Sg' ; elmts(106)%mass= 266.00d0 ; elmts(106)%val= 1d0 
    elmts(107)%symbol='Bh' ; elmts(107)%mass= 264.00d0 ; elmts(107)%val= 1d0 
    elmts(108)%symbol='Hs' ; elmts(108)%mass= 269.00d0 ; elmts(108)%val= 1d0 
    elmts(109)%symbol='Mt' ; elmts(109)%mass= 268.00d0 ; elmts(109)%val= 1d0 
    elmts(110)%symbol='Ds' ; elmts(110)%mass= 269.00d0 ; elmts(110)%val= 1d0 
    elmts(111)%symbol='Rg' ; elmts(111)%mass= 272.00d0 ; elmts(111)%val= 1d0 
    elmts(112)%symbol='Cn' ; elmts(112)%mass= 277.00d0 ; elmts(112)%val= 1d0 
    elmts(113)%symbol='Uut'; elmts(113)%mass= 300.00d0 ; elmts(113)%val= 1d0 !! no data
    elmts(114)%symbol='Fl' ; elmts(114)%mass= 289.00d0 ; elmts(114)%val= 1d0 
    elmts(115)%symbol='Uup'; elmts(115)%mass= 300.00d0 ; elmts(115)%val= 1d0 !! no data
    elmts(116)%symbol='Lv' ; elmts(116)%mass= 298.00d0 ; elmts(116)%val= 1d0 
    elmts(117)%symbol='Uus'; elmts(117)%mass= 300.00d0 ; elmts(117)%val= 1d0 !! no data
    elmts(118)%symbol='Uuo'; elmts(118)%mass= 300.00d0 ; elmts(118)%val= 1d0 !! no data

    
    return
  end subroutine init_element
!=======================================================================
  function get_element(symbol)
    character(len=*),intent(in):: symbol
    type(atom):: get_element

    integer:: i

    do i=1,nelem
      if( trim(symbol).eq.trim(elmts(i)%symbol) ) then
        get_element = elmts(i)
        return
      endif
    enddo
    return
  end function get_element
!=======================================================================
  subroutine get_cube_info(nspmax,specorder,nums,vals)
    integer,intent(in):: nspmax
    character(len=*),intent(in):: specorder(nspmax)
    integer,intent(out):: nums(nspmax)
    real(8),intent(out):: vals(nspmax)

    integer:: is,iel
    character(len=3):: s 

    nums(:) = 0
    vals(:) = 0d0
    do is=1,nspmax
      if( trim(specorder(is)).eq.'x' ) exit
      s = specorder(is)
      do iel=1,nelem
        if( trim(s).eq.trim(elmts(iel)%symbol) ) then
          nums(is) = iel
          vals(is) = elmts(iel)%val
          exit
        endif
      enddo
    enddo
    return
  end subroutine get_cube_info
end module element
