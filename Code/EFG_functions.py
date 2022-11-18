# -*- coding: utf-8 -*-
"""
@author: Jonathan Frassineti
"""

"""
In this file, the funcions that calculate and diagonalize the EFG are defined.
"""

import numpy as np
from Crystalline_structure import base, lattice_par1, \
                                  lattice_par2, lattice_par3, \
                                  quadrupole_moment

atomic_EFG = - 9.71736166e21 # V/(m^2), atomic unit of EFG
angtom = 1.0e-10 # Angstrom unit of measure
elementary_charge = 1.6021766e-19 # Coulomb = ampere*second
epsilon0 = 8.8541878e-12 # dielectric constant of the vacuum, in ampere^2*kilogram^−1*meter^−3*second^4
h = 6.6260693e-34 # J*s
hbar = h/(2*np.pi) # J*s

"""
'lattice' refers to the lattice parameter, 'pos' to the position in which the user
wants to calculate the EFG from point charge, and 's' the number of real space vectors
from the origin that build the crystalline structure, starting from the basis; greater the 's',
more accurate the calculation, but also much slower.
"""

def point_charge_EFG(pos, s=15):

    a1 = np.array([lattice_par1, 0, 0])
    a2 = np.array([0, lattice_par2, 0])
    a3 = np.array([0, 0, lattice_par3])

    v = np.zeros((3, 3), dtype='float64')
        
    for u1 in range(-s, s + 1):

        for u2 in range(-s, s + 1):

            for u3 in range(-s, s + 1):

                for k in range(len(base)):
                                                    
                    D = np.zeros(3, dtype='float')
                    D[0] = u1 + base[k, 0]
                    D[1] = u2 + base[k, 1]
                    D[2] = u3 + base[k, 2]

                    if ((D[0] == pos[0]) and (D[1] == pos[1]) and (D[2] == pos[2])): continue

                    r_k = np.zeros(3, dtype='float')

                    r_k[0] = (D[0] - pos[0]) * a1[0] + (D[1] - pos[1]) * a2[0] + (D[2] - pos[2]) * a3[0]
                    r_k[1] = (D[0] - pos[0]) * a1[1] + (D[1] - pos[1]) * a2[1] + (D[2] - pos[2]) * a3[1]
                    r_k[2] = (D[0] - pos[0]) * a1[2] + (D[1] - pos[1]) * a2[2] + (D[2] - pos[2]) * a3[2]

                    rr = r_k[0] * r_k[0] + r_k[1] * r_k[1] + r_k[2] * r_k[2]

                    v[0, 0] += base[k, 3] * base[k, 4] * (3 * r_k[0] * r_k[0] - rr) / pow(rr, 2.5)
                    v[0, 1] += base[k, 3] * base[k, 4] * (3 * r_k[0] * r_k[1]) / pow(rr, 2.5)
                    v[0, 2] += base[k, 3] * base[k, 4] * (3 * r_k[0] * r_k[2]) / pow(rr, 2.5)
                    v[1, 1] += base[k, 3] * base[k, 4] * (3 * r_k[1] * r_k[1] - rr) / pow(rr, 2.5)
                    v[1, 2] += base[k, 3] * base[k, 4] * (3 * r_k[1] * r_k[2]) / pow(rr, 2.5)
                    v[2, 2] += base[k, 3] * base[k, 4] * (3 * r_k[2] * r_k[2] - rr) / pow(rr, 2.5)
           
    v *= (elementary_charge / (4 * np.pi * epsilon0)) * (1 / ((angtom) ** 3))

    EFG = np.array([[v[0, 0], v[0, 1], v[0, 2]],
                    [v[0, 1], v[1, 1], v[1, 2]],
                    [v[0, 2], v[1, 2], v[2, 2]]])
    
    max_EFG = EFG.max()
    
    threshold = max_EFG*0.01
    
    for ii in range(len(EFG)):
        for jj in range(len(EFG[0])):
            if np.abs(EFG[ii][jj]) < threshold:
                EFG[ii][jj] = 0.000
            else: continue

    return EFG

"""
This function diagonalize the EFG in the Principal Axis System (PAS) of the tensor, 
and finds the quadrupolar parameters.
"""

def diagonalize_EFG(tensor):
    eigenvalues, P = np.linalg.eig(tensor)
    P_inv = np.linalg.inv(P)
    X = np.dot(P_inv, tensor)
    D = np.dot(X, P)
    
    max_D = D.max()
    
    threshold = max_D*0.01
    
    for ii in range(len(D)):
        for jj in range(len(D[0])):
            if np.abs(D[ii][jj]) < threshold:
                D[ii][jj] = 0.000
            else: continue
    
    # Now, the diagonal components of the EFG tnesor
    Vzz = max(eigenvalues, key=abs)
    Vyy = min(eigenvalues, key=abs)
    Vxx = -(Vzz+Vyy)
    eta = np.abs((Vyy - Vxx) / Vzz) # Asymmetry parameter of the EFG tensor
    eta = round(eta, 3)
    chi = np.abs(Vzz * elementary_charge * quadrupole_moment / h )*1e-6 # Quadrupolar constant, in MHz
    return D, eta, chi