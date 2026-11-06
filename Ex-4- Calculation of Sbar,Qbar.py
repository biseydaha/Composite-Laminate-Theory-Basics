#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 15:24:29 2024

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

# Define the rotation angle
theta = np.deg2rad(60)  # Convert degrees to radians

# Given applied stresses (in Pa)
sigma_global = np.array([2e6, -3e6, 4e6])  # σ1, σ2, τ12

# Define the transformation matrix [T] , the inverse of the transformation matrix [T_inv] and T transpose
T = np.array([
    [np.cos(theta)**2, np.sin(theta)**2, 2 * np.sin(theta) * np.cos(theta)],
    [np.sin(theta)**2, np.cos(theta)**2, -2 * np.sin(theta) * np.cos(theta)],
    [-np.sin(theta) * np.cos(theta), np.sin(theta) * np.cos(theta), np.cos(theta)**2 - np.sin(theta)**2]
])

T_inv = np.linalg.inv(T)
T_transpose = T.T

# Define the compliance matrix [S] for the unidirectional lamina

S = np.array([
    [1/E1, -v21/E2, 0],
    [-v12/E1, 1/E2, 0],
    [0, 0, 1/G12]
])

# Calculate the transformed compliance matrix [S_bar] 
S_bar = T_transpose @ S @ T

# Calculate the reduced stiffness matrix [Q] which is the inverse of the compliance matrix [S]
Q = np.linalg.inv(S)

# Calculate the transformed reduced stiffness matrix [Q'] = [T_inv]T[Q][T_inv]
Q_bar = np.linalg.inv(S_bar)

# The matrix [R] as given in the image
R = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 2]])
R_inv = np.linalg.inv(R)

# Calculate the global strains using the transformed reduced stiffness matrix [Q']
epsilon_global = np.dot(S_bar, sigma_global)

# Calculate the local strains using the transformation matrix [T] and global strains
epsilon_local = R @ T @ R_inv @ epsilon_global

# Calculate the local stresses using the compliance matrix [S] and local strains
sigma_local = np.dot(T, sigma_global)

np.set_printoptions(precision=3, suppress=True, formatter={'float': '{:0.3e}'.format})

# Print out the results
print("\nCompliance matrix [S]: ")
print(S)
print("\nTransformed compliance matrix [S bar]:")
print(S_bar)
print("\nTransformed reduced stiffness matrix [Q bar]:")
print(Q_bar)
print("\nGlobal strains:")
print(epsilon_global)
print("\nLocal strains:")
print(epsilon_local)
print("\nLocal stresses:")
print(sigma_local)

