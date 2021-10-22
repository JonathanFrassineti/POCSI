# This is the main file, where we import the basis from
# 'Crystalline_structure.py' and the functions that calculate the EFG from point charge
# from 'EFG_functions.py'

import numpy as np
from EFG_functions import diagonalize_EFG, point_charge_EFG

atomic_EFG = - 9.71736166e21 # V/(m^2), atomic unit of EFG
angtom = 1.0e-10 # Angstrom unit of measure
elementary_charge = 1.6021766e-19 # Coulomb = ampere*second
epsilon0 = 8.8541878e-12 # dielectric constant of the vacuum, in ampere^2*kilogram^−1*meter^−3*second^4
h = 6.6260693e-34 # J*s
hbar = h/(2*np.pi) # J*s

lattice_par = 4.722 # Lattice parameter, in Angstrom
quadrupole_moment = -0.052e-28 # Quadrupolar moment of the atom, in m^2

EFG_original = point_charge_EFG(lattice_par,np.array([0.25, 0.00, 0.50]),10)
EFG_diag, eta, chi = diagonalize_EFG(EFG_original,quadrupole_moment)

print("The original point charge EFG is: ")
print(EFG_original)
print("V/(m^2)")

print("The diagonalized point charge EFG is: ")
print(EFG_diag)
print("V/(m^2)")

print("The EFG parameters are: ")
print("eta = ",eta)
print("chi = ",chi," MHz")