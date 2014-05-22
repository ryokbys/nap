#!/usr/bin/ruby
#
# Class definitions for molecular dynamics simulation
#

class MD_atom
  attr_accessor :x, :y, :z, :species

  def initialize(x,y,z,species)
    @x= x
    @y= y
    @z= z
    @species= species
  end
end

class MD_system
  attr_reader :atoms, :natm, :a1, :a2, :a3, :species
  
  def initialize(a1,a2,a3)
    @a1=[]
    @a2=[]
    @a3=[]
    @a1[0]= a1[0]
    @a1[1]= a1[1]
    @a1[2]= a1[2]
    @a2[0]= a2[0]
    @a2[1]= a2[1]
    @a2[2]= a2[2]
    @a3[0]= a3[0]
    @a3[1]= a3[1]
    @a3[2]= a3[2]
    @atoms= []
    @natm= 0
    @species= []
  end

  def add_atom(atom)
    @atoms[@atoms.size]= atom
    @natm += 1
    @species[@species.size]= atom.species unless species_exist?(atom.species)
  end

  def species_exist?(sid)
    @species.each do |s|
      return true if s == sid
    end
    return false
  end

  def num_atoms_of_species(sid)
    n=0
    @atoms.each do |a|
      n += 1 if a.species == sid
    end
    return n
  end

  def p_box
    p @a1
    p @a2
    p @a3
  end

end

