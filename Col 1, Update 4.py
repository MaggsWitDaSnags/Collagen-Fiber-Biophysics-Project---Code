#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 15:41:10 2026

@author: maggie
"""


import math 
import numpy as np
import matplotlib.pyplot as plt

#COL 1


N = 5 #Number of beads
L = 100e-6 #Length of fiber (m)
l = 100e-6/4 #length of each 'spring', fiber split into 4 segments (m)
d = 0.3e-6 #diameter of fiber (m)
E = 32e6 #Young's Modulus (Pa)
A = (math.pi * d**2) / 4 # cross-sectional area of fiber (m^2)
k = (E * A) / l #spring const (N/m)
#MASS = 1e-16 #pressumed mass (kg)
DT = 2e-6 #time step (seconds)
N_steps = 1000000
Damping_C = 30 # Damping coefficient, (Pa * s) = (N/m^2)*s
I = (np.pi * d**4) / 64  # Area moment of inertia (m^4)
# Physical bending stiffness:
# k_b = (E * I) / l
K_bending  = (E * I) / l #(Pa*m^4)/m = Pa*m^3 = (N*m)

#Drag 
tissue_viscosity = 1.7 # (Pa*s) (CAN CHANGE DEPENDING ON ENVIROMENT)
C_drag = 6 * np.pi * tissue_viscosity * (d/2)   #(Pa*s*m = ((N/m^2)*s)*m = ((N/m)*s)

#postioning the fibers
pos = np.zeros((N,2)) #creating array of size N with 2 entries each

pos[:, 0] = np.linspace(0, L, N)#updating the X coords

pos[:, 1] = 5e-6 * np.sin(np.linspace(0, 2 * np.pi, N)) #Updating Y coords

#Crimped state for graph 
initial_pos = pos.copy()

#velocites 
vel = np.zeros((N,2)) #velocity array

#Force Calc

    #Internal force, this will account for stretching
def compute_force(pos, vel, k, l_0, Damping_C):
    forces = np.zeros_like(pos) #creating an array, like pos
    for i in range(N - 1):
        diff = pos[i+1] - pos[i] #difference in position between 2 beads (m)
        dist = np.linalg.norm(diff) #distance between beads (m)
        unit_vec = diff / dist #the vector representaion of the bead to bead connection (i.e. the spring) (unitless)

        # Spring Force
        spring_f = k * (dist - l_0 - (l_0 * 0.15)) * unit_vec 
        #spring const * (distance - the updated length between beads) * the unit vector = (N)
        #^^ this calculates the internal force created by the force on the spring 
        #'(l_0 * 0.15)' effectively stretches the bead 

        # Damping Force, internally
        #rel_vel = (vel[i+1] - vel[i]) #the relative veloicty between two beads (m/s)
        # Project velocity onto the fiber direction
        #damping_mag = (Damping_C * A / l_0) * np.dot(rel_vel, unit_vec) #in (N)
        #damping_f = damping_mag * unit_vec #damping force = magnttiude and direction of force (N)

        total_f = spring_f #+ damping_f #adding up the interal forces (N)
        
        # Apply equal and opposite forces to the connected beads
        forces[i] += total_f
        forces[i+1] -= total_f
       

    #Linearization Logic, without this, the fiber is only streched.
    #Currently the best way I can think of doing this
    #The Y-coords need to be corrected for the fiber to linearize
    for i in range(1, N - 1):
        #This will make the middle bead move towards the midpoint of it's 2 neighboring beads
        midpoint = (pos[i-1] + pos[i+1]) / 2.0 #(m)
        Bead_vec = midpoint - pos[i] #we want to arrange the center bead towards the midpoint (m)
        #Adding on the bending force
        Bending_f = K_bending * Bead_vec #(N)
        
        forces[i] += Bending_f #(N)
        forces[i-1] -= 0.5 * Bending_f
        forces[i+1] -= 0.5 * Bending_f
        
    # Apply Boundary Conditions AFTER internal force summation
    forces[0] = 0  # Fixed anchor, so it's not floating through space (THIS WILL CHANGE WHEN SIMULATION SCALES)
    #^will end up making fiber look to it's 'neighbors' for a boundary condition
    return forces   

#creating a stressor
stress = 3e3 #Pa (CAN CHANGE DEPENDING ON ENVIROMENT)
f_ext = stress * A #external stress = stress * cross-sectional area (N)



for step in range(N_steps): #running the function
    forces = compute_force(pos, vel, k, l, Damping_C)
    forces[-1,0] += f_ext # last bead in the x column
    forces[-1,1] += f_ext # last bead in the y column
    
    #Trying the overdamping logic from the literature- 
    #This is used bc in something thick like tissue, there isn't a constant momentum
    #Therefore, vecloity is = to force applied over the drag
    #vel = force (N) / drag ((N/m)*s) = 1/(s/m) = m/s
    vel = forces / C_drag #(m/s)
    pos += vel * DT # (m)
    
    if step % 100000 == 0: 
        # Create a string for all (x, y) pairs
        print("(X,Y) positions")
        bead_data = " | ".join([f"({p[0]:.6e}, {p[1]:.6e})" for p in pos])
        print(f"Step {step:<5} | {bead_data}\n") 
        
        
scale = 1e6
plt.figure(figsize=(10, 5))
plt.plot(initial_pos[:, 0] * scale, initial_pos[:, 1] * scale, color = 'darkorange', label='Initial (Crimped State)', alpha=0.5)
plt.plot(pos[:, 0] * scale, pos[:, 1] * scale, color = 'blue', label='Final (Linearized)')
plt.scatter(pos[:, 0] * scale, pos[:, 1] * scale, color='black', zorder=5, label='Beads')
plt.scatter(initial_pos[:, 0] * scale, initial_pos[:, 1] * scale, color='black', zorder=5)
plt.ylim(-5, 55) 
plt.xlabel('Position X (μm)')
plt.ylabel('Position Y (μm)')
plt.title(f'Fiber Simulation: {N} Beads over {N_steps} Steps = 1 second')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.7)
plt.show()

    