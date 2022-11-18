# -*- coding: utf-8 -*-
"""
@author: Jonathan Frassineti
"""

""" 
This is the main file, where we import the basis from
'Crystalline_structure.py' and the functions that calculate the EFG from point charge
from 'EFG_functions.py'.
"""

import numpy as np
from EFG_functions import point_charge_EFG, diagonalize_EFG


"""
Now, in function 'EFG_original()' you have to put the atomic position
where you want to calculate the EFG and the number of steps of the algorithm
(default, s = 15).
"""

site = np.array([0.00, 0.25, -0.25])
 
EFG_original = point_charge_EFG(site, s = 15)
EFG_diag, eta, chi = diagonalize_EFG(EFG_original)

print("The original point charge EFG is: ")
print(EFG_original)
print("V/(m^2)\n\n")

print("The diagonalized point charge EFG is: ")
print(EFG_diag)
print("V/(m^2)\n\n")

print("The EFG parameters at site {} are: ".format(site))
print("eta = ",eta)
print("chi = ",chi," MHz")