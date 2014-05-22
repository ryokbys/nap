#!/usr/bin/ruby

require './MD_classes.rb'

a1= [10.0, 0, 0]
a2= [0, 10.0, 0]
a3= [0, 0, 20.0]

system= MD_system.new(a1,a2,a3)

system.p_box

10.times do |i|
  if i < 5
    j=1
  else
    j=2
  end
  a= MD_atom.new(i*1,i*2,i*3,j)
  system.add_atom(a)
end

system.natm.times do |i|
  p system.atoms[i]
end

p system.species
