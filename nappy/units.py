#!/usr/bin/env python
"""
Units
"""

Ry_to_eV = 13.605698066
eV_to_Ry = 1.0 /Ry_to_eV
eV_to_J = 1.60217657e-19


Bohr_to_Ang = 0.529177249
Ang_to_Bohr = 1.0 /Bohr_to_Ang
Ang_to_m = 1.0e-10
Ang_to_cm = Ang_to_m *100

amu_to_kg = 1.660539040e-27
kg_to_amu = 1.0 /amu_to_kg
amu_to_g = amu_to_kg *1000

kB = 8.6173303e-5  # eV/K
K_to_eV = kB

amu_prss = eV_to_J /(Ang_to_m**3)
amu_to_GPa = 1.0e-9 *amu_prss
