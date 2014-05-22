#!/usr/bin/ruby
#
# Make directories with POSCAR files which have atomic configurations
# with atoms of small displacements
#

require './MD_classes.rb'

#.....PARAMETERS
DISPLACEMENT= 0.1  # Angstrom

def read_POSCAR
  filename= "./POSCAR"
  file= open(filename)
  #.....1st line: comment
  $c1=file.gets
  #.....2nd line: multiply factor
  $afac= file.gets.to_f
  #.....3rd-5th lines: lattice vectors
  a1=[]
  i=0
  (file.gets.split).each do |a|
    a1[i]= a.to_f
    i += 1
  end
  a2=[]
  i=0
  (file.gets.split).each do |a|
    a2[i]= a.to_f
    i += 1
  end
  a3=[]
  i=0
  (file.gets.split).each do |a|
    a3[i]= a.to_f
    i += 1
  end
  $system= MD_system.new(a1,a2,a3)
  #.....6th line: num of atoms
  nums=[]
  i=0
  (file.gets.split).each do |n|
    nums[i]= n.to_i
    i += 1
  end
  #.....7th line: comment
  $c7= file.gets
  #.....8th--: atom positions
  pos=[]
  sid=0
  nums.each do |n|
    sid += 1
    (0..n-1).each do |j|
      pos= file.gets.split
      atom= MD_atom.new(pos[0].to_f, pos[1].to_f, pos[2].to_f, sid)
      $system.add_atom(atom)
    end
  end
  file.close
#  p afac
#  p a1
#  p a2
#  p a3
end

def read_idqm
  filename= "./out.idqm"
  file= open(filename)
  #.....naqm
  file.gets
  #.....idqm
  $idqm=[]
  iaqm=0
  while text= file.gets
    $idqm[iaqm]= text.to_i
    iaqm += 1
  end
  file.close
#  p $idqm
end

read_POSCAR
read_idqm

#.....displacement in relative scale
ux= DISPLACEMENT/$system.a1[0]

$idqm.each do |id|
  sid= sprintf("%04d",id)
  Dir.mkdir(sid) unless File.exist?(sid)
  filename=sid+"/POSCAR"
  file= File.open(filename,"w")
  #.....comment
  file.puts $c1
  #.....multiply factor
  file.puts $afac
  #.....box
  file.printf(" %10.5f %10.5f %10.5f\n",
              $system.a1[0],
              $system.a1[1],
              $system.a1[2])
  file.printf(" %10.5f %10.5f %10.5f\n",
              $system.a2[0],
              $system.a2[1],
              $system.a2[2])
  file.printf(" %10.5f %10.5f %10.5f\n",
              $system.a3[0],
              $system.a3[1],
              $system.a3[2])
  #.....num of atoms of species
  $system.species.each do |s|
    file.printf(" %5d",$system.num_atoms_of_species(s))
  end
  file.printf("\n")
  #.....comment
  file.puts $c7
  #.....atom positions
  i=0
  $system.atoms.each do |a|
    i += 1
    u= 0.0
    u= ux if i == id
    file.printf(" %12.7f %12.7f %12.7f\n", a.x+u, a.y, a.z)
  end
  file.close
  
end
