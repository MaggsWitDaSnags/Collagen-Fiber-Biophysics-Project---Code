# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.

@Author: Maggie Eckhart

3/31/26
"""


import math 
import numpy as np

#COL 1

N = 5 #Number of beads
l = 100e-6 #Length of fiber's in between beads
d = 0.3e-6 #diameter of fiber 
E = 32e6 #Young's Modulus
A = (math.pi * d**2) / 4 # cross-sectional area of fiber
k = -(E * A) / l #spring const
MASS = 1e-16 #pressumed mass 
DT = 0.00000001 #time step
N_steps = 200

#postioning the fibers
pos = np.zeros((N,2)) #creating array of size N with 2 entries each
print(pos)
pos[:,0] = np.arange(N) * l #updating the X coords (indepdent)
pos[:, 1] = np.sin(pos[:, 0]) #updating y coords in terms of how X changes (dependent)

#velocites 
vel = np.zeros((N,2)) #velocity array

#Force Calc
def compute_force(pos, vel, K, l_0):
    forces = np.zeros_like(pos) #creating a new array similar to pos
    forces[N-1] = 0 #anchor one bead
    
    #Internal force
    for i in range(N-1): #does this for all beads
        diff = pos[i+1] - pos[i] #calulates the difference in position of two beads
        dist = np.linalg.norm(diff) #calulates the distance between two beads
        
        #spring force (hooke's law)
        unit_vec = diff / dist #displacement
        force_mag = k * (dist - l) #spring const * (distance moved - length of fiber)
        spring_f = force_mag * unit_vec #force = -kx
        
        forces[i] += spring_f #adding this force to the bead
        forces[i+1] -= spring_f #subracting this force from the next bead (force propigates)

    return forces



#creating a stressor
stress = 1e6
f_ext = stress * A #external stress = stress * cross-sectional area

for step in range(N_steps): #running the function
    forces = compute_force(pos, vel, k, l)
    forces[0,0] += f_ext # (column, row) selecting where the force is being applied
    # EXAMPLE (0,0) = x coord, first entry
    accel = forces / MASS #acceleration
   
    
    vel += accel * DT #integral of acceleration += change in velocity 
    pos += vel * DT #integral of vecloity += change in position
    
    if step % 5 == 0: #for evey 50 steps
        # Create a string for all (x, y) pairs
        print("(X,Y) positions")
        bead_data = " | ".join([f"({p[0]:.2e}, {p[1]:.2e})" for p in pos])
        print(f"Step {step:<5} | {bead_data}\n") 
    
    
