# -*- coding: utf-8 -*-
"""
@author: Jonathan Frassineti
"""

"""
In this file, the user can add the basis that forms the crystalline structure of the material
in which you can calculate the Electric Field Gradient (EFG) by the Point Charge approximation
The following is an example with V3Si, a well known superconductor with cubic structure.
"""

import numpy as np

""" 
Now, the atomic basis is added. Each row consists of an atom with:
[x position, y position, z position, occupation, electric charge]
The user must add every atom that forms the basis.
"""

base = np.array([[0.25, 0.75, 0.00, 1, 1], # V, along a axis
                 [0.75, 0.75, 0.00, 1, 1], # //
                 [0.50, 0.50, 0.50, 1, 1], # V, along b axis
                 [0.50, 0.00, 0.50, 1, 1], # //
                 [0.00, 0.25, 0.75, 1, 1], # V, along c axis
                 [0.00, 0.25, 0.25, 1, 1], # //
                 [0.00, 0.75, 0.50, 1, -3],   # Si
                 [0.50, 0.25, 0.00, 1, -3]])  # Si

lattice_par1 = 4.722 # Lattice parameter 'a' of V3Si, in Angstrom
lattice_par2 = 4.722 # Lattice parameter 'b' of V3Si, in Angstrom
lattice_par3 = 4.722 # Lattice parameter 'c' of V3Si, in Angstrom

quadrupole_moment = -0.052e-28 # Quadrupolar moment of the atom (in this case, 51V), in m^2