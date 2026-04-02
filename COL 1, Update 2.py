#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 10:24:35 2026

@author: maggie
"""

"""
Spyder Editor

This is a temporary script file.

@Author: Maggie Eckhart
"""


import math 
import numpy as np
import matplotlib.pyplot as plt

#COL 1

N = 5 #Number of beads
l = 100e-6 #Length of fiber's in between beads
d = 0.3e-6 #diameter of fiber 
E = 32e6 #Young's Modulus
A = (math.pi * d**2) / 4 # cross-sectional area of fiber
k = (E * A) / l #spring const
#MASS = 1e-16 #pressumed mass 
DT = 1e-6 #time step
N_steps = 1000000
Damping_C = 30 # Damping coefficient, (Pa * s) = (N/m^2)*s
K_bending = 1e-3 # Pa

#Drag 
tissue_viscosity = 2.0 #Pa * s
C_drag = 6 * np.pi * tissue_viscosity * (d/2)   #(Pa*s*m = ((N/m^2)*s)*m = (N/m)*s

#postioning the fibers
pos = np.zeros((N,2)) #creating array of size N with 2 entries each

pos[:,0] = np.arange(N) * l #updating the X coords (indepdent)

pos[:, 1] = np.sin((pos[:, 0] * 5e4)) *2e-6 #updating y coords in terms of how X changes (dependent)

initial_pos = pos.copy()

#velocites 
vel = np.zeros((N,2)) #velocity array

#Force Calc

    #Internal force, this will account for stretching
def compute_force(pos, vel, k, l_0, Damping_C):
    forces = np.zeros_like(pos) #creating an array, like pos
    for i in range(N - 1):
        diff = pos[i+1] - pos[i] #difference in position between 2 beads
        dist = np.linalg.norm(diff) #distance between beads
        unit_vec = diff / dist #the vector representaion of the bead to bead connection (i.e. the spring)

        # Spring Force
        spring_f = k * (dist - l_0) * unit_vec #spring const * (distance - the updated length between beads) * the unit vector
        #^^ this calculates the internal force created by the force on the spring

        # Damping Force
        rel_vel = (vel[i+1] - vel[i]) #the relative veloicty between two beads
        # Project velocity onto the fiber direction
        damping_mag = (Damping_C * A / l_0) * np.dot(rel_vel, unit_vec)
        damping_f = damping_mag * unit_vec #damping force = magnttiude and direction of force

        total_f = spring_f + damping_f #adding up the interal forces
        
        # Apply equal and opposite forces to the connected beads
        forces[i] += total_f
        forces[i+1] -= total_f
       

    #Linearization Logic, without this, the fiber is only streched.
    #The Y-coords need to be corrected for the fiber to linearize
    for i in range(1, N - 1):
        #This will make the middle bead move towards the midpoint of it's 2 neighboring beads
        midpoint = (pos[i-1] + pos[i+1]) / 2.0
        Bead_vec = midpoint - pos[i]
        #Adding on the bending force
        Bending_f = K_bending * Bead_vec
        
        forces[i] += Bending_f
        forces[i-1] -= 0.5 * Bending_f
        forces[i+1] -= 0.5 * Bending_f
        
    # Apply Boundary Conditions AFTER internal force summation
    forces[0] = 0  # Fixed anchor, so it's not floating through space (THIS WILL CHANGE WHEN SIMULATION SCALES)
    #^will end up making fiber look to it's 'neighbors' for a boundary condition
    return forces   

#creating a stressor
stress = 1e3 #1kPa, typical of healthy breast tissue
f_ext = stress * A #external stress = stress * cross-sectional area



for step in range(N_steps): #running the function
    forces = compute_force(pos, vel, k, l, Damping_C)
    forces[-1,0] += f_ext # (column, row) selecting where the force is being applied
    # EXAMPLE (0,0) = x coord, first entry
    forces[-1,1] += f_ext #also apply to last bead, y-coord
    
    #Trying the overdamping logic from the literature- 
    #This is used bc in something thick like tissue, there isn't a constant momentum
    #Therefore, vecloity is = to force applied over the drag
    #vel = force (N) / drag ((N/m)*s) = 1/(s/m) = m/s
    vel = forces / C_drag
    pos += vel * DT
    
    if step % 100000 == 0: 
        # Create a string for all (x, y) pairs
        print("(X,Y) positions")
        bead_data = " | ".join([f"({p[0]:.6e}, {p[1]:.6e})" for p in pos])
        print(f"Step {step:<5} | {bead_data}\n") 
        
        
scale = 1e6
plt.figure(figsize=(10, 5))
plt.plot(initial_pos[:, 0] * scale, initial_pos[:, 1] * scale, color = 'darkorange', label='Initial (Crimped State)', alpha=0.5)
plt.plot(pos[:, 0] * scale, pos[:, 1] * scale, color = 'blue', label='Final (Linearized)')
plt.ylim(-5, 10) 
plt.yticks(np.arange(-4, 11, 2))
plt.xlabel('Position X (μm)')
plt.ylabel('Position Y (μm)')
plt.title(f'Fiber Simulation: {N} Beads over {N_steps} Steps = 1 second')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.7)
plt.show()
    
    