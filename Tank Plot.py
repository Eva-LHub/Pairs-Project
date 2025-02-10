#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 16:35:16 2025

@author: evaloughridge
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

def cylinder_cross_section_area(R, h):
    """
    Compute the area of the liquid portion of a horizontal cylindrical tank.
    
    """
    if h == 0:
        return 0
    else:
        segment_area = R**2 * np.arccos((R - h) / R) - (R - h) * np.sqrt(2 * R * h - h**2)
        return segment_area

def volume_at_height(R, h, L):
    """
    Compute the volume of the tank at a given liquid height.
    
    """
    area = cylinder_cross_section_area(R, h)
    return area * L

def plot_tank_fill(R, L):
    """
    Plot the volume of the tank filled for varying liquid height.
    
    """
    heights = np.linspace(0, 2*R, 1000)  # Heights from 0 to 2R (full tank)
    volumes = [volume_at_height(R, h, L) for h in heights]
    
    plt.plot(volumes, heights, label="Volume vs Height", color="m")
    plt.title(f"Volume of a Horizontal Cylindrical Tank\n(Radius = {R:.2f}m, Length = {L:.2f}m)")
    plt.xlabel("Volume Filled (m³)")  
    plt.ylabel("Height of Liquid (m)") 
    

# Parameters for the tank 
R = 0.5    
V = 1
I= 0.05 #intervals
L= V/(np.pi*(R**2))
print("L is:", L, "m") 
plot_tank_fill(R, L)


def find_height_for_percent_full(R, L, target_percent):
    """
    Find the liquid height corresponding to a specific percentage of the tank's volume.
    
    """
    V_total = np.pi * R**2 * L  
    target_volume = target_percent * V_total  
    
    def equation(h):
        return volume_at_height(R, h, L) - target_volume

    h_initial_guess = R / 2 
    h_solution = fsolve(equation, h_initial_guess)
    
    return h_solution[0]

h_list = [] 

for t in np.arange(I, 1 + I,I):
    height_percent = find_height_for_percent_full(R, L, t)
    h_list.append(height_percent)
    print(f"The height for {t*100:.0f}% full tank is: {height_percent:.3f} meters")
    plt.hlines(height_percent,0,V + I,linestyle = ':', colors= 'teal', dashes=(0, (1,1)))
    plt.vlines(t*V,0,2*R + I,linestyle = ':', colors= 'teal', dashes=(0, (1,1)))
    
print('-------------------------------------------------------------')
h2 = [h/R for h in h_list]
for t, val in zip(np.arange(I, 1+ I, I), h2):
    print(f"Height is {val:.3f}R at {t*1000:.0f}L")
    
plt.hlines(height_percent,0,V + I,linestyle = ':', colors= 'teal', dashes=(0, (1,1)), label= f'Height of Liquid at {I} m³ Intervals')
plt.legend(loc='upper left')
plt.show()



