#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  9 08:01:45 2024

@author: seydaakgun
"""

import numpy as np

# Define material properties
E1 = 181e9  # Pa
E2 = 10.3e9  # Pa
G12 = 7.17e9  # Pa
nu12 = 0.28

# Calculate nu21
nu21 = (nu12 * E2) / E1

# Calculate the Q matrix for a unidirectional lamina
Q11 = E1 / (1 - nu12 * nu21)
Q12 = (nu12 * E2) / (1 - nu12 * nu21)
Q22 = E2 / (1 - nu12 * nu21)
Q66 = G12

Q = np.array([
    [Q11, Q12, 0],
    [Q12, Q22, 0],
    [0, 0, Q66]
])

# Function to transform Q matrix for a given angle theta
def transform_Q(Q, theta):
    m = np.cos(np.radians(theta))
    n = np.sin(np.radians(theta))
    
    T = np.array([
        [m**2, n**2, 2*m*n],
        [n**2, m**2, -2*m*n],
        [-m*n, m*n, m**2 - n**2]
    ])
    
    T_inv = np.linalg.inv(T)
    Q_bar = T_inv @ Q @ T_inv.T
    return Q_bar

# Define laminate layup
angles = [0, 30, -45]
thickness = 0.005  # m per ply
total_thickness = thickness * len(angles)

# Calculate z values (distances from midplane)
z = np.linspace(-total_thickness / 2, total_thickness / 2, len(angles) + 1)

# Initialize A, B, D matrices
A = np.zeros((3, 3))
B = np.zeros((3, 3))
D = np.zeros((3, 3))

# Calculate A, B, D matrices
for i, angle in enumerate(angles):
    Q_bar = transform_Q(Q, angle)
    A += Q_bar * (z[i+1] - z[i])
    B += 0.5 * Q_bar * (z[i+1]**2 - z[i]**2)
    D += (1/3) * Q_bar * (z[i+1]**3 - z[i]**3)

# Define applied loads
Nx = 1000  # N/mm
Ny = 1000  # N/mm
Nxy = 0

N = np.array([Nx, Ny, Nxy])

# Solve for midplane strains and curvatures
ABD = np.block([
    [A, B],
    [B, D]
])

strain_curvature = np.linalg.inv(ABD) @ np.concatenate((N, np.zeros(3)))

midplane_strains = strain_curvature[:3]
curvatures = strain_curvature[3:]

print("Midplane Strains:", midplane_strains)
print("Curvatures:", curvatures)

# Function to calculate global stresses in a ply
def global_stresses(Q_bar, midplane_strains, curvatures, z):
    epsilon = midplane_strains + z * curvatures
    sigma = Q_bar @ epsilon
    return sigma

# Calculate global and local stresses on the top surface of the 30° ply
z_top_30 = z[1]
Q_bar_30 = transform_Q(Q, 30)

global_stress_30 = global_stresses(Q_bar_30, midplane_strains, curvatures, z_top_30)

# Transform global stresses to local stresses
T_30 = np.array([
    [np.cos(np.radians(30))**2, np.sin(np.radians(30))**2, 2*np.sin(np.radians(30))*np.cos(np.radians(30))],
    [np.sin(np.radians(30))**2, np.cos(np.radians(30))**2, -2*np.sin(np.radians(30))*np.cos(np.radians(30))],
    [-1*np.sin(np.radians(30))*np.cos(np.radians(30)), np.sin(np.radians(30))*np.cos(np.radians(30)), (np.cos(np.radians(30))**2- (np.sin(np.radians(30))**2))]
])

local_stress_30 = T_30 @ global_stress_30

print("Global Stress on 30° Ply:", global_stress_30)
print("Local Stress on 30° Ply:", local_stress_30)
