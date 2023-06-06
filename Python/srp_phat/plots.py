"""
This script makes the plots for the testing against noise.
"""

# Import modules.
import matplotlib.pyplot as plt
import numpy as np

# Set the font.
ssfont = {"fontname":"Times New Roman"}

# Import the data.
data = np.loadtxt("./srp_phat/noise_logs/log_angles.csv", delimiter=',')

# Set the true position of the source.
p_true = np.array([5.5, 3, 1])
centre = np.array([5, 5, 1])
r_true = p_true - centre
θ_true = np.arctan2(r_true[1], r_true[0])
φ_true = np.arctan2(np.linalg.norm(r_true[0:2]), r_true[2])

# Calculate the mean and the variance for each noise-level.
n_noises = 10
n_runs = 2000
θ_mean = np.zeros(n_noises)
φ_mean = np.zeros(n_noises)
θ_var = np.zeros(n_noises)
φ_var = np.zeros(n_noises)
for i_noise in range(n_noises):
    θ = data[i_noise*n_runs:(1+i_noise)*n_runs,2]
    φ = data[i_noise*n_runs:(1+i_noise)*n_runs,3]
    θ_mean[i_noise] = np.mean(θ)
    φ_mean[i_noise] = np.mean(φ)
    θ_var[i_noise] = np.var(θ)
    φ_var[i_noise] = np.var(φ)

noises = np.arange(25, -25, -5)

# Plot the results for the azimuth.
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
ax1.plot(noises, np.rad2deg(θ_mean), "kx-")
ax1.plot([-20,25], np.rad2deg([θ_true,θ_true]), "r-")
ax1.set_ylabel("Mean azimuth (degree)", **ssfont)
ax2.plot(noises, np.rad2deg(θ_var), "kx-")
ax2.set_ylabel("Variance of azimuth (degree^2)", **ssfont)
ax2.set_xlabel("SNR (dB)", **ssfont)
ax1.set_title("The Error in Azimuth", **ssfont)
plt.savefig("./srp_phat/noise_azimuth.png")
plt.show()

# Plot the results for the azimuth.
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
ax1.plot(noises, np.rad2deg(φ_mean), "kx-")
ax1.plot([-20,25], np.rad2deg([φ_true,φ_true]), "r-")
ax1.set_ylabel("Mean elevation (degree)", **ssfont)
ax2.plot(noises, np.rad2deg(φ_var), "kx-")
ax2.set_ylabel("Variance of elevation (degree^2)", **ssfont)
ax2.set_xlabel("SNR (dB)", **ssfont)
ax1.set_title("The Error in Elevation", **ssfont)
plt.savefig("./srp_phat/noise_elevation.png")
plt.show()

# # Build the rectangular pyramid of microphones 25 cm wide and 12.5 tall.
# centre = np.array([[5, 5, 1]]).T
# r_m = np.array([
#     [0,0,1],
#     [1,1,0],
#     [1,-1,0],
#     [-1,-1,0],
#     [-1,1,0]
# ]).T * 12.5*10**(-2)
# # r_m = np.array([
# #     [1,1,1],
# #     [1,1,-1],
# #     [1,-1,1],
# #     [1,-1,-1],
# #     [-1,1,1],
# #     [-1,1,-1],
# #     [-1,-1,1],
# #     [-1,-1,-1]
# # ]).T * 15.5*10**(-2) / 2

# p_m = r_m + centre

# # The position of the source:
# p_true = np.array([[5.5, 3, 1]]).T

# fig = plt.figure()
# ax1 = fig.add_subplot(projection = "3d")

# ax1.scatter(p_m[0,:], p_m[1,:], p_m[2,:], c = "r", marker = ".", s = 4)
# ax1.scatter(p_true[0], p_true[1], p_true[2], c = "b", marker = ".", s = 4)
# ax1.plot([0,10],[0,0],[0,0], "k")
# ax1.plot([10,10],[0,10],[0,0], "k")
# ax1.plot([10,0],[10,10],[0,0], "k")
# ax1.plot([0,0],[10,0],[0,0], "k")

# ax1.plot([0,10],[0,0],[3,3], "k")
# ax1.plot([10,10],[0,10],[3,3], "k")
# ax1.plot([10,0],[10,10],[3,3], "k")
# ax1.plot([0,0],[10,0],[3,3], "k")

# ax1.plot([0,0],[0,0],[0,3], "k")
# ax1.plot([0,0],[10,10],[0,3], "k")
# ax1.plot([10,10],[10,10],[0,3], "k")
# ax1.plot([10,10],[0,0],[0,3], "k")

# ax1.set_xlim([-2, 12])
# ax1.set_ylim([-2, 12])
# ax1.set_zlim([-1, 4])

# ax1.set_xlabel("x-coordinate (m)", **ssfont)
# ax1.set_ylabel("y-coordinate (m)", **ssfont)
# ax1.set_zlabel("z-coordinate (m)", **ssfont)

# # ax1.set_title("Positions of the Rectangular Pyramid Array and the Source for Testing Against Noise", **ssfont)
# ax1.legend(["microphones", "source"])

# plt.savefig("./pyramid_robot/room_3d.png")
# # plt.savefig("./srp_phat/room_3d.png")

# plt.show()

# fig, ax1 = plt.subplots()

# ax1.set_aspect("equal")

# ax1.plot([0,0],[0,10], "k")
# ax1.plot([0,10],[10,10], "k")
# ax1.plot([10,10],[10,0], "k")
# ax1.plot([10,0],[0,0], "k")

# ax1.plot(p_m[0,:], p_m[1,:], "r.", marker = ".", markersize = 4)
# ax1.plot(p_true[0], p_true[1], "b.", marker = ".", markersize = 4)

# ax1.set_xlim([-2, 12])
# ax1.set_ylim([-2, 12])

# ax1.set_xlabel("x-coordinate (m)", **ssfont)
# ax1.set_ylabel("y-coordinate (m)", **ssfont)

# ax1.grid(visible = True)

# # ax1.set_title("Positions of the Rectangular Pyramid Array and the Source for Testing Against Noise", **ssfont)
# # ax1.legend(["microphones", "source"])

# plt.savefig("./pyramid_robot/room_2d.png")
# # plt.savefig("./srp_phat/room_2d.png")

# plt.show()