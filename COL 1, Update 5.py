#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 13:37:10 2026

@author: maggie
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 15:41:10 2026

@author: maggie
"""


import math
import numpy as np
import matplotlib.pyplot as plt

# COL 1
N = 5 #Number of beads
L = 100e-6 #Length of fiber (m)
l = 100e-6/4 #length of each spring , fiber split into 4 segments (m)
d = 0.3e-6 #diameter of fiber (m)
E = 32e6 #Young s Modulus (Pa)
A = (math.pi * d**2) / 4 # cross-sectional area of fiber (m^2)
k = (E * A) / l #spring const (N/m)
DT = 2e-6 #time step (seconds)
N_steps = 1000000
Damping_C = 30 # Damping coefficient, (Pa * s) = (N/m^2)*s
I = (np.pi * d**4) / 64 # Area moment of inertia (m^4)
K_bending = (E * I) / l #(Pa*m^4)/m = Pa*m^3 = (N*m)

#Drag
tissue_viscosity = 1.7 # (Pa*s) (CAN CHANGE DEPENDING ON ENVIROMENT)
C_drag = 6 * np.pi * tissue_viscosity * (d/2) #(Pa*s*m = ((N/m^2)*s)*m = ((N/m)*s)

#Diffusion Constant 
D_Parrallel = ((1.380649e-23 * 309.65) / (2 * math.pi * tissue_viscosity * L)) * (math.log(L/d))
D_Perp = ((1.380649e-23 * 309.65) / (4 * math.pi * tissue_viscosity * L)) * (math.log(L/d))
D_avg = (D_Parrallel + (2*D_Perp))/3 # 2 forces perp, 1 parrallel


scale = 1e6
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

def fibers(count):
    for i in range(count):
        #postioning the fibers
        pos = np.zeros((N,3)) #creating array of size N with 2 entries each 
        pos[:, 0] = np.linspace(0, L, N) + (i*40e-6) #updating the X coords
        pos[:, 1] = 5e-6 * np.sin(np.linspace(0, 2 * np.pi, N)) + (i*40e-6) #Updating Y coords
        pos[:, 2] = 5e-6 * np.cos(np.linspace(0, 2 * np.pi, N)) + (i*40e-6) #Updating Z coords
        
        #Crimped state for graph
        initial_pos = pos.copy() 
        #velocites
        vel = np.zeros((N,3)) #velocity array

        def compute_force(pos, vel, k, l_0, Damping_C):
            forces = np.zeros_like(pos) #creating an array, like pos
            for j in range(N - 1):
                diff = pos[j+1] - pos[j] #difference in position between 2 beads (m)
                dist = np.linalg.norm(diff) #distance between beads (m)
                unit_vec = diff / (dist + 1e-20) #the vector representaion of the bead to bead connection
                
                # Spring Force
                spring_f = k * (dist - l_0 - (l_0 * 0.15)) * unit_vec 
                forces[j] += spring_f
                forces[j+1] -= spring_f

            #Linearization Logic
            for j in range(1, N - 1):
                midpoint = (pos[j-1] + pos[j+1]) / 2.0 #(m)
                Bead_vec = midpoint - pos[j] 
                Bending_f = K_bending * Bead_vec #(N)
                forces[j] += Bending_f #(N)
                forces[j-1] -= 0.5 * Bending_f
                forces[j+1] -= 0.5 * Bending_f

            forces[0] = 0 # Fixed anchor
            return forces 

        #creating a stressor
        stress = 0.2e3 #Pa (CAN CHANGE DEPENDING ON ENVIROMENT)
        f_ext = stress * A #external stress = stress * cross-sectional area (N)

        for step in range(N_steps):
            forces = compute_force(pos, vel, k, l, Damping_C)
            forces[-1,0] += f_ext # last bead in the x column
            forces[-1,1] += f_ext # last bead in the y column
            forces[-1,2] += f_ext # last bead in the z column

            #Thermal noise
            noise_mag = np.sqrt(2 * D_avg * DT)
            thermal_displacement = noise_mag * np.random.normal(0, 1, size=pos.shape)
            thermal_displacement[0] = 0 

            vel = forces / C_drag 
            pos += (vel * DT) + thermal_displacement # (m)

            if step % 100000 == 0:
                bead_data = " | ".join([f"({p[0]:.6e}, {p[1]:.6e}, {p[2]:.6e})" for p in pos])
                print(f"Fiber {i+1} Step {step:<5} | {bead_data}\n")
            
        #inital position
        ax.plot(initial_pos[:, 0] * scale, initial_pos[:, 1] * scale, initial_pos[:, 2] * scale, color='darkorange', linestyle='--', alpha=0.4, label=f'Fiber {i+1} Initial')
        ax.scatter(initial_pos[:, 0] * scale, initial_pos[:, 1] * scale, initial_pos[:, 2] * scale, color='orange', s=20, alpha=0.4)

        #final position
        ax.plot(pos[:, 0] * scale, pos[:, 1] * scale, pos[:, 2] * scale, color='blue', label=f'Fiber {i+1} Final')
        ax.scatter(pos[:, 0] * scale, pos[:, 1] * scale, pos[:, 2] * scale, color='black', s=30, zorder=5)


    #plotting
    ax.set_xlabel('X (μm)')
    ax.set_ylabel('Y (μm)')
    ax.set_zlabel('Z (μm)')
    ax.set_title(f'3D Fiber Simulation: {count} Fibers')
    ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1)) # Moved legend so it doesn't block the view
    plt.tight_layout()
    plt.show()


fibers(2)
