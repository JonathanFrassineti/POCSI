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

"""
Now, the possibility to add an impurity, or vacancy, or an additional particle
(such as a positive muon) inside the crystal, to obtain the EFG with this perturbation.
First, define the position of the impurity, the atomic base of N nearest neighbours 
(for this version, only 4 considered) which are displaced due to the presence of the impurity, the displacement itself
and the list of N nearest neighbours involved from the original atomic base
(used for EFG calculations).
"""
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

jump_atoms = np.array([0,1,4,5]) # List of atoms considered as N nearest neighbours,
                                 # from the original atomic base