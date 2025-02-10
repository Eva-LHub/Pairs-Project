import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

def cylinder_cross_section_area(a,b, h):
    if h == 0:
        return 0     
    else:
        segment_area = (a*b)*(np.arccos(1 - h/b) - (1-h/b)*np.sqrt(((2*h)/b)- (h**2)/b**2))
        return segment_area

def volume_at_height(a,b, h, L):
    area = cylinder_cross_section_area(a,b, h)
    return area * L

def plot_tank_fill(a,b, L):
    heights = np.linspace(0, 2*b, 1000)  # Heights from 0 to 2b (full tank)
    volumes = [volume_at_height(a,b, h, L) for h in heights]
    
    plt.plot(volumes, heights, label="Volume vs Height", color="m")
    plt.title(f"Volume of a Horizontal Cylindrical Tank\n(Semi Minor Axis = {b:.2f}m, Semi Major Axis = {a:.2f}m, Length = {L:.2f}m)")
    plt.xlabel("Volume Filled (m³)")  
    plt.ylabel("Height of Liquid (m)") 
    

# Parameters for the tank 
a= 0.3
b= 0.8   
V = 1
I= 0.05 #intervals
L= V/(np.pi*a*b)
print("L is:", L, "m") 

plot_tank_fill(a,b, L)


def find_height_for_percent_full(a,b, L, target_percent):

    V_total = np.pi *a *b * L  
    target_volume = target_percent * V_total  

    def equation(h):
        return volume_at_height(a,b, h, L) - target_volume

    h_initial_guess = b / 2 
    h_solution = fsolve(equation, h_initial_guess)
    
    return h_solution[0]

h_list = []

for t in np.arange(I, 1 + I,I):
    height_percent = find_height_for_percent_full(a,b, L, t)
    h_list.append(height_percent)
    print(f"The height for {t*100:.0f}% full tank is: {height_percent:.3f} meters")
    plt.hlines(height_percent,0,V + I,linestyle = ':', colors= 'teal', dashes=(0, (1,1)))
    plt.vlines(t*V,0,2*b + I,linestyle = ':', colors= 'teal', dashes=(0, (1,1)))
    
print('-------------------------------------------------------------')


h2 = [h/(a*b) for h in h_list]
for t, val in zip(np.arange(I, 1+ I, I), h2):
    print(f"Height is {val:.3f}ab at {t*1000:.0f}L")

    
plt.hlines(height_percent,0,V + I,linestyle = ':', colors= 'teal', dashes=(0, (1,1)), label= f'Height of Liquid at {I} m³ Intervals')
plt.legend(loc='upper left')
plt.show()


