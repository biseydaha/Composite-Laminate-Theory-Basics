#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  8 22:38:26 2024

@author: seydaakgun
"""

import numpy as np

# Given ply properties
E1 = 20010000  # lb/in^2
E2 = 1301000   # lb/in^2
G12 = 1001000  # lb/in^2
v12 = 0.3
v21 = 0.02
ply_thickness = 0.005      # in (ply thickness)

# Q matrix components
Q11 = E1 / (1 - v12 * v21)
Q22 = E2 / (1 - v12 * v21)
Q12 = v12 * E2 / (1 - v12 * v21)
Q66 = G12

Q = np.array([
    [Q11, Q12, 0],
    [Q12, Q22, 0],
    [0, 0, Q66]
])

    
# Function to calculate transformed reduced stiffnesses
def transformed_reduced_stiffness(theta):
    c = np.cos(np.radians(theta))
    s = np.sin(np.radians(theta))
    
    Q_bar_30 = np.zeros((3, 3))
    Q_bar_30[0, 0] = Q11 * c**4 + 2 * (Q12 + 2 * Q66) * s**2 * c**2 + Q22 * s**4
    Q_bar_30[0, 1] = (Q11 + Q22 - 4 * Q66) * s**2 * c**2 + Q12 * (c**4 + s**4)
    Q_bar_30[1, 1] = Q11 * s**4 + 2 * (Q12 + 2 * Q66) * s**2 * c**2 + Q22 * c**4
    Q_bar_30[0, 2] = (Q11 - Q12 - 2 * Q66) * c**3 * s - (Q22 - Q12 - 2 * Q66) * s**3 * c
    Q_bar_30[1, 2] = (Q11 - Q12 - 2 * Q66) * s**3 * c - (Q22 - Q12 - 2 * Q66) * c**3 * s
    Q_bar_30[2, 2] = (Q11 + Q22 - 2 * Q12 - 2 * Q66) * s**2 * c**2 + Q66 * (s**4 + c**4)
    
    Q_bar[1, 0] = Q_bar[0, 1]
    Q_bar[2, 0] = Q_bar[0, 2]
    Q_bar[2, 1] = Q_bar[1, 2]
    
    return Q_bar

# Ply orientations
orientations = [0, 45, 45, 0]
    
num_plies = len(orientations)
total_thickness = ply_thickness * num_plies


# Calculate z values (distances from midplane)
h = np.linspace(-total_thickness / 2, total_thickness / 2, num_plies + 1)

                      
# Calculate A matrix
A = np.zeros((3, 3))

for theta in orientations:
    Q_bar = transformed_reduced_stiffness(theta)
    A += Q_bar * ply_thickness


A11 =  A[0,0]
A12 =  A[0,1]
A22 =  A[1,1]
A16 =  A[0,2]
A26 =  A[1,2]
A66 =  A[2,2]
"""
Nx_dividedby_Ex0 = (A11 + 
                    (A12 * (A26*A16 - A12 * A66 / (A22 * A66 - A26 ** 2)))
                     + A16 * ((-A16/A66) +  ((A12 * A26 * A66 - (A26**2 * A16)) /
                                           (A22 * A66**2 - (A26 ** 2 * A66)))))


Ex = (A11/total_thickness + 
                    (A12/total_thickness  * (A26*A66 - A12 * A26 / (A22 * A66 - A26 ** 2)))
                     + (A16/total_thickness  * ((-A16/A66) +  ((A12 * A26 * A66 - (A26**2 * A16)) /
                                        (A22 * A66**2 - (A26 ** 2 * A66))))))

# Normalize by total thickness
Ex = Nx_dividedby_Ex0 / total_thickness                          

"""

Nx_dividedby_Ex0 = (A11 + 
                    A12 * (A26 * A16 - A12 * A66) / (A22 * A66 - A26 ** 2) +
                    A16 * (-A16 / A66 + (A12 * A26 * A66 - A26 ** 2 * A16) / 
                    (A22 * A66 ** 2 - A26 ** 2 * A66)))

# Normalize by total thickness
Ex = Nx_dividedby_Ex0 / total_thickness

print("Ex of the laminate is: {:.2f} lb/in^2".format(Ex))

