"""
This script makes the plots for the testing against noise.
"""

# Import modules.
import matplotlib.pyplot as plt
import numpy as np

# Set the font.
ssfont = {"fontname":"Times New Roman"}

# Build the rectangular pyramid of microphones 25 cm wide and 12.5 tall.
centre = np.array([[5, 5, 1]]).T
# r_m = np.array([
#  ,  [0,0,1],
#  ,  [1,1,0],
#  ,  [1,-1,0],
#  ,  [-1,-1,0],
#  ,  [-1,1,0]
# ]).T * 12.5*10**(-2)
r_m = np.array([
    [1,1,1],
    [1,1,-1],
    [1,-1,1],
    [1,-1,-1],
    [-1,1,1],
    [-1,1,-1],
    [-1,-1,1],
    [-1,-1,-1]
]).T * 15.5*10**(-2) / 2

p_m = r_m + centre

# The position of the source:
p_true = np.array([[5.5, 3, 1]]).T

fig = plt.figure()
ax1 = fig.add_subplot(projection = "3d")

ax1.scatter(p_m[0,:], p_m[1,:], p_m[2,:], c = "r", marker = ".", s = 4)
ax1.scatter(p_true[0], p_true[1], p_true[2], c = "b", marker = ".", s = 4)
ax1.plot([0,10],[0,0],[0,0], "k")
ax1.plot([10,10],[0,10],[0,0], "k")
ax1.plot([10,0],[10,10],[0,0], "k")
ax1.plot([0,0],[10,0],[0,0], "k")

ax1.plot([0,10],[0,0],[3,3], "k")
ax1.plot([10,10],[0,10],[3,3], "k")
ax1.plot([10,0],[10,10],[3,3], "k")
ax1.plot([0,0],[10,0],[3,3], "k")

ax1.plot([0,0],[0,0],[0,3], "k")
ax1.plot([0,0],[10,10],[0,3], "k")
ax1.plot([10,10],[10,10],[0,3], "k")
ax1.plot([10,10],[0,0],[0,3], "k")

ax1.set_xlim([-2, 12])
ax1.set_ylim([-2, 12])
ax1.set_zlim([-1, 4])

ax1.set_xlabel("x-coordinate (m)", **ssfont)
ax1.set_ylabel("y-coordinate (m)", **ssfont)
ax1.set_zlabel("z-coordinate (m)", **ssfont)

# ax1.set_title("Positions of the Rectangular Pyramid Array and the Source for Testing Against Noise", **ssfont)
ax1.legend(["microphones", "source"])

plt.savefig("./pyramid_robot/room_3d.png")
# plt.savefig("./main_method/room_3d.png")

plt.show()

fig, ax1 = plt.subplots()

ax1.set_aspect("equal")

ax1.plot([0,0],[0,10], "k")
ax1.plot([0,10],[10,10], "k")
ax1.plot([10,10],[10,0], "k")
ax1.plot([10,0],[0,0], "k")

ax1.plot(p_m[0,:], p_m[1,:], "r.", marker = ".", markersize = 4)
ax1.plot(p_true[0], p_true[1], "b.", marker = ".", markersize = 4)

ax1.set_xlim([-2, 12])
ax1.set_ylim([-2, 12])

ax1.set_xlabel("x-coordinate (m)", **ssfont)
ax1.set_ylabel("y-coordinate (m)", **ssfont)

ax1.grid(visible = True)

# ax1.set_title("Positions of the Rectangular Pyramid Array and the Source for Testing Against Noise", **ssfont)
# ax1.legend(["microphones", "source"])

plt.savefig("./pyramid_robot/room_2d.png")
# plt.savefig("./main_method/room_2d.png")

plt.show()