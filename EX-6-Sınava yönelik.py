#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  9 09:15:57 2024

@author: seydaakgun
"""

import numpy as np
import pandas as pd

# Material properties for unidirectional graphite/epoxy
E1 = 181e9  # Longitudinal modulus (Pa)
E2 = 10.3e9  # Transverse modulus (Pa)
G12 = 7.17e9  # Shear modulus (Pa)
nu12 = 0.28  # Poisson's ratio

# Calculate Poisson's ratio nu21
nu21 = (nu12 * E2) / E1

# Q matrix for a unidirectional lamina
Q11 = E1 / (1 - nu12 * nu21)
Q12 = (nu12 * E2) / (1 - nu12 * nu21)
Q22 = E2 / (1 - nu12 * nu21)
Q66 = G12

Q = np.array([
    [Q11, Q12, 0],
    [Q12, Q22, 0],
    [0, 0, Q66]
])

def transform_Q(Q, theta):
    """Transform the Q matrix for a given angle theta."""
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

# Display Q matrix
Q_matrix = pd.DataFrame(Q, columns=[' ', ' ', ' '], index=['', ' ', ' '])
print("Qbar 0 Matrix:")
print(Q_matrix)
print() 

Q_matrix_30 = pd.DataFrame(transform_Q(Q, 30), columns=[' ', ' ', ' '], index=['', ' ', ' '])
print("Qbar 30 Matrix:")
print(Q_matrix_30)
print() 

Q_matrix_45 = pd.DataFrame(transform_Q(Q, -45), columns=[' ', ' ', ' '], index=['', ' ', ' '])
print("Qbar -45 Matrix:")
print(Q_matrix_45)
print() 
# Laminate layup definition
angles = [0, 30, -45]
thickness = 0.005  # Thickness of each ply (m)
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


A_matrix = pd.DataFrame(A, columns=[' ', ' ', ' '], index=['', ' ', ' '])
print("A:")
print(A_matrix)
print() 


B_matrix = pd.DataFrame(B, columns=[' ', ' ', ' '], index=['', ' ', ' '])
print("B:")
print(B_matrix)
print() 


D_matrix = pd.DataFrame(D, columns=[' ', ' ', ' '], index=['', ' ', ' '])
print("D:")
print(D_matrix)
print() 

# Applied loads
Nx = 1000  # Load in x-direction (N/m)
Ny = 1000  # Load in y-direction (N/m)
Nxy = 0  # Shear load (N/m)

N = np.array([Nx, Ny, Nxy])

# Solve for midplane strains and curvatures
ABD = np.block([
    [A, B],
    [B, D]
])

strain_curvature = np.linalg.inv(ABD) @ np.concatenate((N, np.zeros(3)))

midplane_strains = strain_curvature[:3]
curvatures = strain_curvature[3:]


midplane_strainss = pd.DataFrame(midplane_strains, columns=[' '], index=['', ' ', ' '])
print("Midplane Strains:")
print(midplane_strainss)
print() 


curvaturess = pd.DataFrame(curvatures, columns=[' '], index=['', ' ', ' '])
print("Curvatures:")
print(curvaturess)
print() 


def global_stresses(Q_bar, midplane_strains, curvatures, z):
    """Calculate global stresses in a ply."""
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
    [-np.sin(np.radians(30))*np.cos(np.radians(30)), np.sin(np.radians(30))*np.cos(np.radians(30)), (np.cos(np.radians(30))**2 - np.sin(np.radians(30))**2)]
])

local_stress_30 = T_30 @ global_stress_30

Global_stress_30_ply = pd.DataFrame(global_stress_30, columns=[' '], index=['', ' ', ' '])
print("Global Stress on 30° Ply:")
print(Global_stress_30_ply)
print() 


Local_stress_30_ply = pd.DataFrame(local_stress_30, columns=[' '], index=['', ' ', ' '])
print("Local Stress on 30° Ply:")
print(Local_stress_30_ply)
print() 

# Calculate percentage of load Nx taken by each ply

strains_0_middle_global = midplane_strains + (-0.0050) * curvatures
stresses_0_middle_global = transform_Q(Q, 0) @ strains_0_middle_global
portion_of_0_ply = stresses_0_middle_global[0]*0.005
Percentage_of_0_ply = portion_of_0_ply *100/1000

print("Global Strains on 0° Ply's midpoint':", strains_0_middle_global)
print("Global Stress on  0° Ply's midpoint':", stresses_0_middle_global)
print("Portion of 0 degree ply :", portion_of_0_ply )
print("Percentage of 0 degree ply % ':", Percentage_of_0_ply)

print() 

strains_30_middle_global = midplane_strains + (0) * curvatures
stresses_30_middle_global = transform_Q(Q, 30) @ strains_30_middle_global
portion_of_30_ply = stresses_30_middle_global[0]*0.005
Percentage_of_30_ply = portion_of_30_ply *100/1000


print("Global Strains on 30° Ply's midpoint':", strains_30_middle_global)
print("Global Stress on  30° Ply's midpoint':", stresses_30_middle_global)
print("Portion of 30 degree ply :", portion_of_30_ply )
print("Percentage of 30 degree ply % ':", Percentage_of_30_ply)

print() 

strains_45_middle_global = midplane_strains + (0.005) * curvatures
stresses_45_middle_global = transform_Q(Q, -45)@ strains_45_middle_global
portion_of_45_ply = stresses_45_middle_global[0]*0.005
Percentage_of_45_ply = portion_of_45_ply *100/1000


print("Global Strains on 45° Ply's midpoint':", strains_45_middle_global)
print("Global Stress on  45° Ply's midpoint':", stresses_45_middle_global)
print("Portion of -45 degree ply :", portion_of_45_ply )
print("Percentage of -45 degree ply % ':", Percentage_of_45_ply)


print() 

# Display percentage of load Nx taken by each ply
percentage_load_data = pd.DataFrame({
    'Ply Angle (°)': [0, 30, -45],
    'Portion of Ply (N/m)': [portion_of_0_ply, portion_of_30_ply, portion_of_45_ply],
    'Percentage of Ply (%)': [Percentage_of_0_ply, Percentage_of_30_ply, Percentage_of_45_ply]
})
print(percentage_load_data)
print()

# Material strength properties for maximum stress theory
Xt = 1500e6  # Tensile strength in fiber direction (Pa)
Xc = -1500e6  # Compressive strength in fiber direction (Pa)
Yt = 40e6  # Tensile strength transverse to fiber direction (Pa)
Yc = -246e6  # Compressive strength transverse to fiber direction (Pa)
S = 68e6  # Shear strength (Pa)

def failure_mode(stress_local):
    """Determine potential failure mode using maximum stress theory."""
    sigma1, sigma2, tau12 = stress_local
    if sigma1 > Xt:
        return "Fiber tension failure"
    elif sigma1 < Xc:
        return "Fiber compression failure"
    elif sigma2 > Yt:
        return "Matrix tension failure"
    elif sigma2 < Yc:
        return "Matrix compression failure"
    elif abs(tau12) > S:
        return "Shear failure"
    else:
        return "No failure"

failure_30 = failure_mode(local_stress_30)

print("Potential failure mode in 30° ply:", failure_30)