# -*- coding: utf-8 -*-
"""
@author: Jonathan Frassineti
"""

""" 
This is an example file, where we the functions that calculate the EFG from point charge
from pocsi.py'.
"""

import numpy as np
from pocsi import point_charge_EFG, diagonalize_EFG


"""
Now, in function 'EFG_original()' you have to put the atomic position
where you want to calculate the EFG,number of steps of the algorithm
(default, s = 15), the lattice parameters, the atomic base and the quadrupole moment.
"""
base = np.array([[0.25, 0.75, 0.00, 1, 1], # V, along a axis
                 [0.75, 0.75, 0.00, 1, 1], # //
                 [0.50, 0.50, 0.50, 1, 1], # V, along b axis
                 [0.50, 0.00, 0.50, 1, 1], # //
                 [0.00, 0.25, 0.75, 1, 1], # V, along c axis
                 [0.00, 0.25, 0.25, 1, 1], # //
                 [0.00, 0.75, 0.50, 1, -3],   # Si
                 [0.50, 0.25, 0.00, 1, -3]])  # Si

lattice_a = 4.722 # Lattice parameter 'a' of V3Si, in Angstrom
lattice_b = 4.722 # Lattice parameter 'b' of V3Si, in Angstrom
lattice_c = 4.722 # Lattice parameter 'c' of V3Si, in Angstrom

quadrupole_moment = -0.052e-28 # Quadrupolar moment of the atom (in this case, 51V), in m^2

site = np.array([0.00, 0.25, -0.25])
 
EFG_original = point_charge_EFG(site, base, lattice_a, lattice_b, lattice_c, quadrupole_moment, s = 10)
EFG_diag, eta, chi = diagonalize_EFG(EFG_original, quadrupole_moment)

print("The original point charge EFG is: ")
print(EFG_original)
print("V/(m^2)\n\n")

print("The diagonalized point charge EFG is: ")
print(EFG_diag)
print("V/(m^2)\n\n")

print("The EFG parameters at site {} are: ".format(site))
print("eta = ",eta)
print("chi = ",chi," MHz")

############################################################################################################
############################ PERTURBATION ##################################################################
############################################################################################################

"""
Now, the possibility to add an impurity, or vacancy, or an additional particle
(such as a positive muon) inside the crystal, to obtain the EFG with this perturbation.
First, define the position of the impurity, the atomic base of N nearest neighbours 
(for this version, only 4 considered) which are displaced due to the presence of the impurity, 
the displacement itself and the list of N nearest neighbours involved from the original atomic 
base (used for EFG calculations).
"""

base = np.array([[0.25, 0.75, 0.00, 1, 1], # V, along a axis
                 [0.75, 0.75, 0.00, 1, 1], # //
                 [0.50, 0.50, 0.50, 1, 1], # V, along b axis
                 [0.50, 0.00, 0.50, 1, 1], # //
                 [0.00, 0.25, 0.75, 1, 1], # V, along c axis
                 [0.00, 0.25, 0.25, 1, 1], # //
                 [0.00, 0.75, 0.50, 1, -3],   # Si
                 [0.50, 0.25, 0.00, 1, -3]])  # Si

lattice_a = 4.722 # Lattice parameter 'a' of V3Si, in Angstrom
lattice_b = 4.722 # Lattice parameter 'b' of V3Si, in Angstrom
lattice_c = 4.722 # Lattice parameter 'c' of V3Si, in Angstrom

quadrupole_moment = -0.052e-28 # Quadrupolar moment of the atom (in this case, 51V), in m^2

site = np.array([0.00, 0.25, -0.25])

xxx = np.array([0,0,0,1,1]) # Same as for V or Si atom, in this case a H atom at (0,0,0) position

shift_a = 1.05 # The displacement of N nearest neighbours 
               # due to the presence of 'xxx' impurity along 'a'
shift_b = 1.05 # Same, along 'b'
shift_c = 1.05 # Same, along 'c'

shifts = np.array([shift_a,shift_b,shift_c])

base_shifted = np.array([[ 0.25*shift_a, -0.25*shift_b,  0.00*shift_c, 1, 1], # V, along a, shifted in a and b
                         [-0.25*shift_a, -0.25*shift_b,	 0.00*shift_c, 1, 1], # V, along a, shifted in a and b
                         [ 0.00*shift_a,  0.25*shift_b,  0.25*shift_c, 1, 1], # V, along c, shifted in b and c
                         [ 0.00*shift_a,  0.25*shift_b, -0.25*shift_c, 1, 1]]) # V, along c, shifted in b and c
 
EFG_original_withPerturbation = point_charge_EFG(site, base, lattice_a, lattice_b, lattice_c, quadrupole_moment, s = 15, impurity=True, xxx = xxx, shifts = shifts, base_shifted = base_shifted)
EFG_diag_withPerturbation, eta_withPerturbation, chi_withPerturbation = diagonalize_EFG(EFG_original_withPerturbation, quadrupole_moment)

print("The original point charge EFG with the perturbation is: ")
print(EFG_original_withPerturbation)
print("V/(m^2)\n\n")

print("The diagonalized point charge EFG eith the perturbation is: ")
print(EFG_diag_withPerturbation)
print("V/(m^2)\n\n")

print("The EFG parameters at site {} are: ".format(site))
print("eta = ",eta_withPerturbation)
print("chi = ",chi_withPerturbation," MHz")