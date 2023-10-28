import matplotlib.pyplot as plt
import numpy as np

import sys
sys.path.append("./main_method")

import srp_phat

# Build the cube of microphones 15.5 cm wide.
centre = np.array([[5, 5, 1]]).T
r_m = np.array([
    [1,1,1],    # NE top
    [1,1,-1],   # NE bottom
    [-1,1,1],   # NW top
    [-1,1,-1],  # NW bottom
    [1,-1,1],   # SE top
    [1,-1,-1],  # SE bottom
    [-1,-1,1],  # SW top
    [-1,-1,-1]  # SW bottom
]).T * 91.2*10**(-3) / 2

f_s = 24*10**3

"""
The Beamforming Method
"""

# Get the number of samples in the frame.
F = 1024

# Make the object for the SRP-PHAT.
estimator = srp_phat.srp_phat(r_m, f_s, F)
# Compute the grid and all TDOAs.
grid = estimator.make_spherical_grid(4)
# grid = make_circular_grid()
τ = estimator.compute_tdoa_grid(grid)
# τ_k = (1024 + τ.round().astype(int)) % 1024
τ_k = τ.round().astype(int)

table = np.zeros((28,2562))

key = 0
for i in range(8):
    for j in range(i+1,8):
        table[key+(j-i-1),:] = τ_k[:,i,j]
        print(key+(j-i-1))
    key = key + (8-i-1)

string = ""
for i in range(28):
    string = string + "__attribute__((section(\"tables_data_{}\")))\nconst volatile int8_t tau_table_{}[] = {{\n".format(i//20,i)
    for j in range(2562):
        string = string + str(int(table[i,j]))
        if j != 2561:
            string = string + ","
        if j%64 == 63:
            string = string + "\n"
    string = string + "};\n"

f = open("./main_method/table/tau.txt", "w")
f.write(string)
f.close()