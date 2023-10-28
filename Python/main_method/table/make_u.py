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

θ = np.arctan2(grid[1,:], grid[0,:])
ψ = np.arctan2(np.linalg.norm(grid[0:2,:], axis=0), grid[2,:])

string = "const volatile float32_t grid[] = {\n"

string = string + "{\n"
for i in range(2562):
    string = string + str(θ[i]) + ","
    if i%64 == 63:
        string = string + "\n"
string = string + "}\n"

string = string + "{\n"
for i in range(2562):
    string = string + str(ψ[i]) + ","
    if i%64 == 63:
        string = string + "\n"
string = string + "}\n"

string = string + "};\n"

print(string)

print(θ[264]*180/np.pi)
print(ψ[264]*180/np.pi)
