#!/usr/bin/ruby
#
# Gather forces of QM atoms identified in out.idqm
#

require './MD_classes.rb'

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

filename='dat.frc-id'
outfile=open(filename,"w")
$idqm.each do |id|
  cid= sprintf("%04d",id)
  if !File.exist?(cid)
    puts 'Dir '+cid+' does not exist !!!'
  end
  filename=cid+'/frc.vasp'
  frcfile= open(filename)
  natm= frcfile.gets.to_i
  afrc= Array.new
  iatm=0
  natm.times do
    atom= $system.atoms[iatm]
    iatm += 1
    text= frcfile.gets.split
    if iatm == id
      outfile.printf(" %5d %s %12.7f\n",id,text[0],atom.z)
      break
    end
  end
  frcfile.close
end
outfile.close

