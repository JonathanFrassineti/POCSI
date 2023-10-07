# -*- coding: utf-8 -*-
"""
@author: Jonathan Frassineti
"""

"""
In this file, the funcions that calculate and diagonalize the EFG are defined.
"""

import numpy as np

atomic_EFG = - 9.71736166e21 # V/(m^2), atomic unit of EFG
angtom = 1.0e-10 # Angstrom unit of measure
elementary_charge = 1.6021766e-19 # Coulomb = ampere*second
epsilon0 = 8.8541878e-12 # dielectric constant of the vacuum, in ampere^2*kilogram^−1*meter^−3*second^4
h = 6.6260693e-34 # J*s
hbar = h/(2*np.pi) # J*s

def point_charge_EFG(pos, base, lattice_a, lattice_b, lattice_c, quadrupole_moment, s=15, impurity = False, xxx = [], shifts = [], base_shifted = []):
    """This function calculates the point charge EFG for a given nucleus in a specific position in the lattice structure.
       
    Parameters
        pos : position of the nucleus in which the user wants to calculate the EFG from point charge.

        s: the number of real space vectors from the origin that build the crystalline structure, starting from the basis; 
        greater the 's', more accurate the calculation, but also much slower.

        base: atomic base of the material under study.
        
        lattice_a/b/c : lattice parameters a, b and c of the crystal, in Angstrom.
        
        quadrupole_moment: the quadrupole moment of the nucleus in which calculating the EFG.
        
        impurity = False: if True, add the possibility to add an impurity, or vacancy, or an additional particle
        (such as a positive muon) inside the crystal, to obtain the EFG with this perturbation.
        This part is currently under improvement (use with caution!).
        
        xxx = []: if impurity = True, define the position of the impurity (as in the 'base' parameter').

        shifts = []: if impurity = True, defines the displacement of each nearest neighbour atom of the nucleus under study.

        base_shifted = []: if impurity = True, defines the atomic base of N nearest neighbours 
        (for this version, only 4 considered) as for the 'base' parameter

        Returns:
        The point charge EFG of the nucleus under study. 
    
    """
    a1 = np.array([lattice_a, 0, 0])
    a2 = np.array([0, lattice_b, 0])
    a3 = np.array([0, 0, lattice_c])

    if impurity == True: 
        for i in range(len(pos)):
            pos[i] *= shifts[i]
            
        v = np.zeros((3, 3), dtype='float64')

        for u1 in range(-s, s + 1):
    
            for u2 in range(-s, s + 1):
    
                for u3 in range(-s, s + 1):
    
                    for k in range(len(base)):
                        
                        if ((u1 ==  0) and (u2 == -1) and (u3 ==  0) and (k == 0)):
                            continue
                        if ((u1 == -1) and (u2 == -1) and (u3 ==  0) and (k == 1)):
                            continue
                        if ((u1 ==  0) and (u2 ==  0) and (u3 == -1) and (k == 4)):
                            continue
                        if ((u1 ==  0) and (u2 ==  0) and (u3 ==  0) and (k == 5)):
                            continue
                            
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
        
            
        for k in range(len(base_shifted)):
                
            D = np.zeros(3, dtype = 'float')
            D[0] = base_shifted[k,0];
            D[1] = base_shifted[k,1];
            D[2] = base_shifted[k,2];
                                        
            if ((D[0] == pos[0]) and (D[1] == pos[1]) and (D[2] == pos[2])): continue
            
            r_k = np.zeros(3, dtype = 'float')
                                                               
            r_k[0] = ( D[0] - pos[0] )*a1[0] + (D[1] - pos[1])*a2[0] + (D[2] - pos[2])*a3[0];
            r_k[1] = ( D[0] - pos[0] )*a1[1] + (D[1] - pos[1])*a2[1] + (D[2] - pos[2])*a3[1];
            r_k[2] = ( D[0] - pos[0] )*a1[2] + (D[1] - pos[1])*a2[2] + (D[2] - pos[2])*a3[2];
                            
            rr = r_k[0]*r_k[0]+r_k[1]*r_k[1]+r_k[2]*r_k[2]
                                               
            v[0,0] += base_shifted[k,3] * base_shifted[k,4] *(3*r_k[0]*r_k[0] - rr) / pow(rr,2.5);
            v[0,1] += base_shifted[k,3] * base_shifted[k,4] *(3*r_k[0]*r_k[1]) / pow(rr,2.5);
            v[0,2] += base_shifted[k,3] * base_shifted[k,4] *(3*r_k[0]*r_k[2]) / pow(rr,2.5);
            v[1,1] += base_shifted[k,3] * base_shifted[k,4] *(3*r_k[1]*r_k[1] - rr) / pow(rr,2.5);
            v[1,2] += base_shifted[k,3] * base_shifted[k,4] *(3*r_k[1]*r_k[2]) / pow(rr,2.5);
            v[2,2] += base_shifted[k,3] * base_shifted[k,4] *(3*r_k[2]*r_k[2] - rr) / pow(rr,2.5);
                 
        D = np.zeros(3, dtype = 'float')
        D[0] = xxx[0];
        D[1] = xxx[1];
        D[2] = xxx[2];
    
        r_k = np.zeros(3, dtype = 'float')
    
        r_k[0] = ( D[0] - pos[0] )*a1[0] + (D[1] - pos[1])*a2[0] + (D[2] - pos[2])*a3[0];
        r_k[1] = ( D[0] - pos[0] )*a1[1] + (D[1] - pos[1])*a2[1] + (D[2] - pos[2])*a3[1];
        r_k[2] = ( D[0] - pos[0] )*a1[2] + (D[1] - pos[1])*a2[2] + (D[2] - pos[2])*a3[2];
    
        rr = r_k[0]*r_k[0]+r_k[1]*r_k[1]+r_k[2]*r_k[2]
    
        v[0,0] += xxx[3] * xxx[4] *(3*r_k[0]*r_k[0] - rr) / pow(rr,2.5);
        v[0,1] += xxx[3] * xxx[4] *(3*r_k[0]*r_k[1]) / pow(rr,2.5);
        v[0,2] += xxx[3] * xxx[4] *(3*r_k[0]*r_k[2]) / pow(rr,2.5);
        v[1,1] += xxx[3] * xxx[4] *(3*r_k[1]*r_k[1] - rr) / pow(rr,2.5);
        v[1,2] += xxx[3] * xxx[4] *(3*r_k[1]*r_k[2]) / pow(rr,2.5);
        v[2,2] += xxx[3] * xxx[4] *(3*r_k[2]*r_k[2] - rr) / pow(rr,2.5);
                                 
    
    else:
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

def diagonalize_EFG(tensor, quadrupole_moment):
    """This function diagonalize the EFG in the Principal Axis System (PAS).
       
    Parameters
        tensor: the EFG tensor calculated from point charge.
        quadrupole_moment: the quadrupole_moment of the nucleus under study.
    
    Returns:
        D: the diagonal EFG from point charge, in the PAS.
        eta: the asymmetry parameter of the EFG.
        chi: the quadrupolar coupling constant of the nucleus, in MHz.
    
    """
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
    
    # Now, the diagonal components of the EFG tensor
    Vzz = max(eigenvalues, key=abs)
    Vyy = min(eigenvalues, key=abs)
    Vxx = -(Vzz+Vyy)
    eta = np.abs((Vyy - Vxx) / Vzz) # Asymmetry parameter of the EFG tensor
    eta = round(eta, 3)
    chi = np.abs(Vzz * elementary_charge * quadrupole_moment / h )*1e-6 # Quadrupolar constant, in MHz
    return D, eta, chi