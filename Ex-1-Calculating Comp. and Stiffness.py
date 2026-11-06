#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 15:13:41 2024

@author: seydaakgun
"""

import numpy as np

# Given material properties

E1 = 181.0e9  # Pa
E2 = 10.3e9   # Pa
E3 = 10.3e9   # Pa
v12 = 0.28    # Dimensionless
v23 = 0.60    # Dimensionless
v13 = 0.27    # Dimensionless
G12 = 7.17e9  # Pa
G23 = 3.0e9   # Pa
G31 = 7.0e9   # Pa
#


# Calculating v21, v31, v32 from given Poisson's ratios and elastic moduli
v21 = v12 * E2 / E1
v31 = v13 * E3 / E1
v32 = v23 * E3 / E2

# Calculation of the determinant Δ
Delta = (1 - v12 * v21 - v23 * v32 - v13 * v31 - 2 * v21 * v32 * v13)/(E1*E2*E3)

# Calculating the compliance matrix [C] which is based on the determinant Δ
C_matrix = np.array([
    [(1-v23*v32)/(E2*E3*Delta), (v21+v23*v31)/(E2*E3*Delta), (v31+v21*v32)/(E2*E3*Delta), 0, 0, 0],
    [(v21+v23*v31)/(E2*E3*Delta), (1-v13*v31)/(E1*E3*Delta), (v32+v12*v31)/(E1*E3*Delta), 0, 0, 0],
    [(v31+v21*v32)/(E2*E3*Delta), (v32+v12*v31)/(E1*E3*Delta), (1-v12*v21)/(E1*E2*Delta), 0, 0, 0],
    [0, 0, 0, G23, 0, 0],
    [0, 0, 0, 0, G31, 0],
    [0, 0, 0, 0, 0, G12]
])

# Now we will solve for the stiffness matrix [S] which is the inverse of the compliance matrix [C]
S_matrix_fromtheinverse = np.linalg.inv(C_matrix)

S_matrix_simplified = np.array([
    [1/E1, -v12/E1, -v13/E1, 0, 0, 0],
    [-v21/E2, 1/E2, -v23/E2, 0, 0, 0],
    [-v31/E3, -v32/E3, 1/E3, 0, 0, 0],
    [0, 0, 0, 1/G23, 0, 0],
    [0, 0, 0, 0, 1/G31, 0],
    [0, 0, 0, 0, 0, 1/G12]
])



# Determine the width needed for each column
column_width = max(len(f"{element:.2e}") for row in S_matrix_simplified for element in row) + 2

# Loop through each row in the matrix
for row in S_matrix_simplified:
    # Loop through each element in the row
    for element in row:
        # Print the element with two decimal places and right-aligned
        print(f"{element:>{column_width}.3e}", end=" ")
    # Print a new line for the next row
    print()


# Use numpy's set_printoptions for setting the precision
np.set_printoptions(precision=3, suppress=True, formatter={'float': '{:0.3e}'.format})
print(C_matrix)
    
    
