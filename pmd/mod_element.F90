module element
  implicit none
  save

  integer,parameter:: nelem = 118
  
  type atom
    character(len=3):: symbol
    real(8):: mass = 1.0d0
  end type atom
  type(atom):: elements(nelem)

contains
!=======================================================================
  subroutine init_element()

    elements(  1)%symbol='H'  ; elements(  1)%mass= 1.008d0 
    elements(  2)%symbol='He' ; elements(  2)%mass= 4.0026d0
    elements(  3)%symbol='Li' ; elements(  3)%mass= 6.94d0  
    elements(  4)%symbol='Be' ; elements(  4)%mass= 9.0122d0
    elements(  5)%symbol='B'  ; elements(  5)%mass= 10.81d0 
    elements(  6)%symbol='C'  ; elements(  6)%mass= 12.011d0
    elements(  7)%symbol='N'  ; elements(  7)%mass= 14.007d0
    elements(  8)%symbol='O'  ; elements(  8)%mass= 15.999d0
    elements(  9)%symbol='F'  ; elements(  9)%mass= 18.998d0
    elements( 10)%symbol='Ne' ; elements( 10)%mass= 20.180d0
    elements( 11)%symbol='Na' ; elements( 11)%mass= 22.990d0
    elements( 12)%symbol='Mg' ; elements( 12)%mass= 24.305d0
    elements( 13)%symbol='Al' ; elements( 13)%mass= 26.982d0
    elements( 14)%symbol='Si' ; elements( 14)%mass= 28.085d0
    elements( 15)%symbol='P'  ; elements( 15)%mass= 30.974d0
    elements( 16)%symbol='S'  ; elements( 16)%mass= 32.06d0 
    elements( 17)%symbol='Cl' ; elements( 17)%mass= 35.45d0 
    elements( 18)%symbol='Ar' ; elements( 18)%mass= 39.948d0
    elements( 19)%symbol='K'  ; elements( 19)%mass= 39.098d0
    elements( 20)%symbol='Ca' ; elements( 20)%mass= 40.078d0
    elements( 21)%symbol='Sc' ; elements( 21)%mass= 44.956d0
    elements( 22)%symbol='Ti' ; elements( 22)%mass= 47.867d0
    elements( 23)%symbol='V'  ; elements( 23)%mass= 50.942d0
    elements( 24)%symbol='Cr' ; elements( 24)%mass= 51.996d0
    elements( 25)%symbol='Mn' ; elements( 25)%mass= 54.938d0
    elements( 26)%symbol='Fe' ; elements( 26)%mass= 55.845d0
    elements( 27)%symbol='Co' ; elements( 27)%mass= 58.933d0
    elements( 28)%symbol='Ni' ; elements( 28)%mass= 58.693d0
    elements( 29)%symbol='Cu' ; elements( 29)%mass= 63.546d0
    elements( 30)%symbol='Zn' ; elements( 30)%mass= 65.38d0 
    elements( 31)%symbol='Ga' ; elements( 31)%mass= 69.723d0
    elements( 32)%symbol='Ge' ; elements( 32)%mass= 72.631d0
    elements( 33)%symbol='As' ; elements( 33)%mass= 74.922d0
    elements( 34)%symbol='Se' ; elements( 34)%mass= 78.972d0
    elements( 35)%symbol='Br' ; elements( 35)%mass= 79.904d0
    elements( 36)%symbol='Kr' ; elements( 36)%mass= 84.798d0
    elements( 37)%symbol='Rb' ; elements( 37)%mass= 85.468d0
    elements( 38)%symbol='Sr' ; elements( 38)%mass= 87.62d0
    elements( 39)%symbol='Y'  ; elements( 39)%mass= 88.906d0
    elements( 40)%symbol='Zr' ; elements( 40)%mass= 91.224d0
    elements( 41)%symbol='Nb' ; elements( 41)%mass= 92.906d0
    elements( 42)%symbol='Mo' ; elements( 42)%mass= 95.95d0
    elements( 43)%symbol='Tc' ; elements( 43)%mass= 98.907d0
    elements( 44)%symbol='Ru' ; elements( 44)%mass= 101.07d0
    elements( 45)%symbol='Rh' ; elements( 45)%mass= 102.90d0
    elements( 46)%symbol='Pd' ; elements( 46)%mass= 106.42d0
    elements( 47)%symbol='Ag' ; elements( 47)%mass= 107.86d0
    elements( 48)%symbol='Cd' ; elements( 48)%mass= 112.41d0
    elements( 49)%symbol='In' ; elements( 49)%mass= 114.81d0
    elements( 50)%symbol='Sn' ; elements( 50)%mass= 118.71d0
    elements( 51)%symbol='Sb' ; elements( 51)%mass= 121.76d0
    elements( 52)%symbol='Te' ; elements( 52)%mass= 127.60d0
    elements( 53)%symbol='I'  ; elements( 53)%mass= 126.90d0
    elements( 54)%symbol='Xe' ; elements( 54)%mass= 131.29d0
    elements( 55)%symbol='Cs' ; elements( 55)%mass= 132.90d0
    elements( 56)%symbol='Ba' ; elements( 56)%mass= 137.32d0
    elements( 57)%symbol='La' ; elements( 57)%mass= 138.90d0
    elements( 58)%symbol='Ce' ; elements( 58)%mass= 140.11d0
    elements( 59)%symbol='Pr' ; elements( 59)%mass= 140.90d0
    elements( 60)%symbol='Nd' ; elements( 60)%mass= 144.24d0
    elements( 61)%symbol='Pm' ; elements( 61)%mass= 144.91d0
    elements( 62)%symbol='Sm' ; elements( 62)%mass= 150.36d0
    elements( 63)%symbol='Eu' ; elements( 63)%mass= 151.96d0
    elements( 64)%symbol='Gd' ; elements( 64)%mass= 157.25d0
    elements( 65)%symbol='Tb' ; elements( 65)%mass= 158.92d0
    elements( 66)%symbol='Dy' ; elements( 66)%mass= 162.50d0
    elements( 67)%symbol='Ho' ; elements( 67)%mass= 164.93d0
    elements( 68)%symbol='Er' ; elements( 68)%mass= 167.25d0
    elements( 69)%symbol='Tm' ; elements( 69)%mass= 168.93d0
    elements( 70)%symbol='Yb' ; elements( 70)%mass= 173.05d0
    elements( 71)%symbol='Lu' ; elements( 71)%mass= 174.96d0
    elements( 72)%symbol='Hf' ; elements( 72)%mass= 178.49d0
    elements( 73)%symbol='Ta' ; elements( 73)%mass= 180.94d0
    elements( 74)%symbol='W'  ; elements( 74)%mass= 183.84d0
    elements( 75)%symbol='Re' ; elements( 75)%mass= 186.20d0
    elements( 76)%symbol='Os' ; elements( 76)%mass= 190.23d0
    elements( 77)%symbol='Ir' ; elements( 77)%mass= 192.21d0
    elements( 78)%symbol='Pt' ; elements( 78)%mass= 195.08d0
    elements( 79)%symbol='Au' ; elements( 79)%mass= 196.96d0
    elements( 80)%symbol='Hg' ; elements( 80)%mass= 200.59d0
    elements( 81)%symbol='Tl' ; elements( 81)%mass= 204.38d0
    elements( 82)%symbol='Pb' ; elements( 82)%mass= 207.20d0
    elements( 83)%symbol='Bi' ; elements( 83)%mass= 208.98d0
    elements( 84)%symbol='Po' ; elements( 84)%mass= 208.98d0
    elements( 85)%symbol='At' ; elements( 85)%mass= 209.98d0
    elements( 86)%symbol='Rn' ; elements( 86)%mass= 222.01d0
    elements( 87)%symbol='Fr' ; elements( 87)%mass= 223.20d0
    elements( 88)%symbol='Ra' ; elements( 88)%mass= 226.02d0
    elements( 89)%symbol='Ac' ; elements( 89)%mass= 227.02d0
    elements( 90)%symbol='Th' ; elements( 90)%mass= 232.03d0
    elements( 91)%symbol='Pa' ; elements( 91)%mass= 231.03d0
    elements( 92)%symbol='U'  ; elements( 92)%mass= 238.02d0
    elements( 93)%symbol='Np' ; elements( 93)%mass= 237.04d0
    elements( 94)%symbol='Pu' ; elements( 94)%mass= 244.06d0
    elements( 95)%symbol='Am' ; elements( 95)%mass= 243.06d0
    elements( 96)%symbol='Cm' ; elements( 96)%mass= 247.00d0
    elements( 97)%symbol='Bk' ; elements( 97)%mass= 247.07d0
    elements( 98)%symbol='Cf' ; elements( 98)%mass= 251.08d0
    elements( 99)%symbol='Es' ; elements( 99)%mass= 254.00d0
    elements(100)%symbol='Fm' ; elements(100)%mass= 257.09d0
    elements(101)%symbol='Md' ; elements(101)%mass= 258.10d0
    elements(102)%symbol='No' ; elements(102)%mass= 259.10d0
    elements(103)%symbol='Lr' ; elements(103)%mass= 262.00d0
    elements(104)%symbol='Rf' ; elements(104)%mass= 261.00d0
    elements(105)%symbol='Db' ; elements(105)%mass= 262.00d0
    elements(106)%symbol='Sg' ; elements(106)%mass= 266.00d0
    elements(107)%symbol='Bh' ; elements(107)%mass= 264.00d0
    elements(108)%symbol='Hs' ; elements(108)%mass= 269.00d0
    elements(109)%symbol='Mt' ; elements(109)%mass= 268.00d0
    elements(110)%symbol='Ds' ; elements(110)%mass= 269.00d0
    elements(111)%symbol='Rg' ; elements(111)%mass= 272.00d0
    elements(112)%symbol='Cn' ; elements(112)%mass= 277.00d0
    elements(113)%symbol='Uut'; elements(113)%mass= 300.00d0  !! no data
    elements(114)%symbol='Fl' ; elements(114)%mass= 289.00d0  
    elements(115)%symbol='Uup'; elements(115)%mass= 300.00d0  !! no data
    elements(116)%symbol='Lv' ; elements(116)%mass= 298.00d0
    elements(117)%symbol='Uus'; elements(117)%mass= 300.00d0  !! no data
    elements(118)%symbol='Uuo'; elements(118)%mass= 300.00d0  !! no data

    
    return
  end subroutine init_element
!=======================================================================
  function get_element(symbol)
    character(len=*),intent(in):: symbol
    type(atom):: get_element

    integer:: i

    do i=1,nelem
      if( trim(symbol).eq.trim(elements(i)%symbol) ) then
        get_element = elements(i)
        return
      endif
    enddo
    return
  end function get_element
end module element
