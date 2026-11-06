#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 10:15:44 2024

@author: seydaakgun
"""

import numpy as np

# Given material properties for a unidirectional graphite/epoxy lamina
E1 = 181.0e9  # Longitudinal elastic modulus in Pa
E2 = 10.3e9   # Transverse elastic modulus in Pa
v12 = 0.28    # Major Poisson's ratio, dimensionless
G12 = 7.17e9  # Shear modulus in Pa

# Using the reciprocity relationship to calculate the minor Poisson's ratio v21
v21 = v12 * E2 / E1

# Given applied stresses (in Pa)
stress_matrix = np.array([1e3, 1e3, 0])  # σ1, σ2, τ12

# Calculating the compliance matrix [C]
S_matrix = np.array([
    [1/E1, -v12/E1, 0],
    [-v12/E1, 1/E2, 0],
    [0, 0, 1/G12]
])

# Reduced stiffness matrix [Q] is the inverse of the compliance matrix [S]
Q_reduced_stiffness = np.linalg.inv(S_matrix)

# Calculate the strains in the 1-2 coordinate system using the compliance matrix [S]
# ε = [C] * {σ}
strain_matrix = np.dot(S_matrix, stress_matrix)


# Printing the compliance matrix with three decimal places
print("Compliance matrix [S]:")
print(np.array2string(S_matrix, formatter={'float_kind':'{0:.3e}'.format}))

# Printing the reduced stiffness matrix [Q] with two decimal places
print("\nReduced stiffness matrix [Q]:")
print(np.array2string(Q_reduced_stiffness, formatter={'float_kind':'{0:.3e}'.format}))

# Printing the minor Poisson's ratio
print("\nMinor Poisson's ratio v21:")
print("{:.2e}".format(v21))

# Printing the strains with three decimal places
print("\nStrains in the 1-2 coordinate system:")
print(np.array2string(strain_matrix, formatter={'float_kind':'{0:.3e}'.format}))
