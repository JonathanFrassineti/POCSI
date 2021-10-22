# In this file, the user can add the basis that forms the crystalline structure of the material
# in which you can calculate the Electric Field Gradient (EFG) by the Point Charge approximation
# The following is an example with V3Si, a well known superconductor with cubic structure.

import numpy as np

atomic_EFG = - 9.71736166e21 # V/(m^2), atomic unit of EFG
angtom = 1.0e-10 # Angstrom unit of measure
elementary_charge = 1.6021766e-19 # Coulomb = ampere*second
epsilon0 = 8.8541878e-12 # dielectric constant of the vacuum, in ampere^2*kilogram^−1*meter^−3*second^4
h = 6.6260693e-34 # J*s
hbar = h/(2*np.pi) # J*s

# Now the basis is added. Each row consists of an atom with:
# [x position, y position, z position, occupation, electric charge]
# The user must add every atom that forms the basis.

base = np.array([[0.25, 0.00, 0.5, 1, 1], # V, along a axis
                 [0.75, 0.00, 0.5, 1, 1], # //
                 [0.50, 0.25, 0.00, 1, 1], # V, along b axis
                 [0.50, 0.75, 0.00, 1, 1], # //
                 [0.00, 0.50, 0.75, 1, 1], # V, along c axis
                 [0.00, 0.50, 0.25, 1, 1], # //
                 [0.00, 0.00, 0.00, 1, -3],   # Si
                 [0.50, 0.50, 0.50, 1, -3]])  # Si


