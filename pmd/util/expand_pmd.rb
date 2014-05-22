#!/usr/bin/ruby
#
# Expand pmd if the system is smaller than the cutoff radius
#
# Usage:
#   $ ./expand_pmd.rb <pmd-file> <cutoff-radius>
#
#.....add the directory where this file exists to the search path
$: << File.dirname(__FILE__)
require 'MD_classes.rb'

################################################## subroutines #########
def read_pmd(filename='pmd00000')
  file= open(filename)
  # 1st: hunit
  hunit= file.gets.to_f
  # 2nd-4th: cell vectors
  a1=[]
  a2=[]
  a3=[]
  i=0
  (file.gets.split).each do |a|
    a1[i]= a.to_f*hunit
    i += 1
  end
  i=0
  (file.gets.split).each do |a|
    a2[i]= a.to_f*hunit
    i += 1
  end
  i=0
  (file.gets.split).each do |a|
    a3[i]= a.to_f*hunit
    i += 1
  end
  sys= MD_system.new(a1,a2,a3,hunit)
  # 5th-7th: skip 3 other lines
  tmp=file.gets
  tmp=file.gets
  tmp=file.gets
  # 8th: num of atoms
  natm= file.gets.to_i
  # 9th-: atoms positions from here
  (0..natm-1).each do |i|
    data= file.gets.split
    sid= data[0].to_i
    atom= MD_atom.new(data[1].to_f, data[2].to_f, data[3].to_f, sid)
    sys.add_atom(atom)
  end
  file.close
  return sys
end

def expand_system(sys0,rcut)
  a1= sys0.a1
  a2= sys0.a2
  a3= sys0.a3
  a12= cross_product(a1,a2)
  a23= cross_product(a2,a3)
  a31= cross_product(a3,a1)
  la1= dot_product(a1,a23) /Math.sqrt(dot_product(a23,a23))
  la2= dot_product(a2,a31) /Math.sqrt(dot_product(a31,a31))
  la3= dot_product(a3,a12) /Math.sqrt(dot_product(a12,a12))
  m1=1
  m2=1
  m3=1
  if la1 < 2.0*rcut then
    m1= (2.0*rcut/la1 +1.0).to_i
  end
  if la2 < 2.0*rcut then
    m2= (2.0*rcut/la2 +1.0).to_i
  end
  if la3 < 2.0*rcut then
    m3= (2.0*rcut/la3 +1.0).to_i
  end
  at1=[]
  at2=[]
  at3=[]
  (0..2).each do |i|
    at1[i]= a1[i]*m1
    at2[i]= a2[i]*m2
    at3[i]= a3[i]*m3
  end
  #.....create new system
  sys1= MD_system.new(at1,at2,at3,sys0.hunit)
  m1.times do |i1|
    m2.times do |i2|
      m3.times do |i3|
        sys0.natm.times do |i|
          atom0= sys0.atoms[i]
          x= (i1 +atom0.x)/m1
          y= (i2 +atom0.y)/m2
          z= (i3 +atom0.z)/m3
          atom1= MD_atom.new(x,y,z,atom0.species)
          sys1.add_atom(atom1)
        end
      end
    end
  end
  return sys1
end

def cross_product(v1,v2)
  v3= []
  v3[0]= v1[1]*v2[2] -v1[2]*v2[1]
  v3[1]= v1[2]*v2[0] -v1[0]*v2[2]
  v3[2]= v1[0]*v2[1] -v1[1]*v2[0]
  return v3
end

def dot_product(v1,v2)
  return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]
end

def write_pmd(sys)
  hunit= sys.hunit
  a1= sys.a1
  a2= sys.a2
  a3= sys.a3
  3.times do |i|
    a1[i]=a1[i]/hunit
    a2[i]=a2[i]/hunit
    a3[i]=a3[i]/hunit
  end
  printf("%15.7f\n" % hunit)
  printf(" %12.7f %12.7f %12.7f\n" % [a1[0],a1[1],a1[2]])
  printf(" %12.7f %12.7f %12.7f\n" % [a2[0],a2[1],a2[2]])
  printf(" %12.7f %12.7f %12.7f\n" % [a3[0],a3[1],a3[2]])
  printf(" %12.7f %12.7f %12.7f\n" % [0.0,0.0,0.0])
  printf(" %12.7f %12.7f %12.7f\n" % [0.0,0.0,0.0])
  printf(" %12.7f %12.7f %12.7f\n" % [0.0,0.0,0.0])
  printf("%8d\n" % sys.natm)
  i=0
  (sys.atoms).each do |atom|
    i += 1
    printf(" %22.14e %12.7f %12.7f %12.7f %6.3f %6.3f %6.3f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n"\
           % [species_to_tag(atom.species,i), atom.x, atom.y, atom.z, 0.0, 0.0, 0.0, \
           0.0, 0.0, \
           0.0, 0.0, 0.0, \
           0.0, 0.0, 0.0, \
           0.0, 0.0, 0.0] )
  end
end

def species_to_tag(species,num)
  return species.to_f() +1.0e-1*1 +1.0e-14*num
end

################################################## main routine ########
usage=' Usage: ./expand_pmd.rb <pmd-file> <cutoff-radius>'

if ARGV.length != 2 then
  puts usage
  exit 0
end

pmdfile=ARGV[0]
rcut= ARGV[1].to_f

sys0=read_pmd(pmdfile)
sys1= expand_system(sys0,rcut)
write_pmd(sys1)
