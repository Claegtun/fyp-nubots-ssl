"""
This script makes the plots for the testing against noise.
"""

# Import modules.
import matplotlib.pyplot as plt
import numpy as np

# Set the font.
ssfont = {"fontname":"Times New Roman"}

# Import the data.
p_data = np.loadtxt("./pyramid_robot/noise_logs/log_r.csv", delimiter=',')
n_failures_data = np.loadtxt("./pyramid_robot/noise_logs/log_n_failures.csv", delimiter=',')

# Ignore the first hundred. Something weird happened there.
p_data = p_data[100:,:]

# Find where the failures happened, i.e. where the position is just the origin.
i_failure = ((p_data[:,0] != 0) | (p_data[:,1] != 0) | (p_data[:,2] != 0))\
    .reshape((100000,1))
i_failure = np.repeat(i_failure, 3, axis = 1)

# Set the true position of the source.
p_true = np.array([5.5, 3, 1])
centre = np.array([5, 5, 1])
r_true = p_true - centre
d_true = np.linalg.norm(r_true)
θ_true = np.arctan2(r_true[1], r_true[0])
φ_true = np.arctan2(np.linalg.norm(r_true[0:2]), r_true[2])

# Calculate the distance, the azimuth, and the elevation.
r_data = (p_data - centre) * i_failure
# d_data = np.linalg.norm(r_data, axis = 1)
# θ_data = np.arctan2(r_data[:,1], r_data[:,0])
# φ_data = np.arctan2(np.linalg.norm(r_data[:,0:2], axis = 1), r_data[:,2])

# Calculate the mean and the variance for each noise-level.
n_noises = 10
n_runs = 10000
d_mean = np.zeros(n_noises)
θ_mean = np.zeros(n_noises)
φ_mean = np.zeros(n_noises)
d_var = np.zeros(n_noises)
θ_var = np.zeros(n_noises)
φ_var = np.zeros(n_noises)
for i_noise in range(n_noises):
    # Choose the next hundred.
    r_noise = r_data[i_noise*n_runs:(1+i_noise)*n_runs]
    i_clean = i_failure[i_noise*n_runs:(1+i_noise)*n_runs,:]
    # Clean out any samples from failures.
    r_clean = r_noise[i_clean[:,0]]
    d = np.linalg.norm(r_clean, axis = 1)
    θ = np.arctan2(r_clean[:,1], r_clean[:,0])
    φ = np.arctan2(np.linalg.norm(r_clean[:,0:2], axis = 1), r_clean[:,2])
    d_mean[i_noise] = np.mean(d)
    θ_mean[i_noise] = np.mean(θ)
    φ_mean[i_noise] = np.mean(φ)
    d_var[i_noise] = np.var(d)
    θ_var[i_noise] = np.var(θ)
    φ_var[i_noise] = np.var(φ)

# Build the rectangular pyramid of microphones 25 cm wide and 12.5 tall.
centre = np.array([[5, 5, 1]]).T
r_m = np.array([
    [0,0,1],
    [1,1,0],
    [1,-1,0],
    [-1,-1,0],
    [-1,1,0]
]).T * 12.5*10**(-2)
# r_m = np.array([
#     [1,1,1],
#     [1,1,-1],
#     [1,-1,1],
#     [1,-1,-1],
#     [-1,1,1],
#     [-1,1,-1],
#     [-1,-1,1],
#     [-1,-1,-1]
# ]).T * 15.5*10**(-2) / 2

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
# plt.savefig("./srp_phat/room_3d.png")

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
# plt.savefig("./srp_phat/room_2d.png")

plt.show()

# noises = np.arange(25, -25, -5)

# # Plot the results for the distance.
# fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
# ax1.plot(noises, d_mean, "kx-")
# ax1.plot([-20,25], [d_true,d_true], "r-")
# ax1.set_ylabel("Mean distance (m)", **ssfont)
# ax2.plot(noises, d_var, "kx-")
# ax2.set_ylabel("Variance of distance (m^2)", **ssfont)
# ax3.plot(noises, n_failures_data, "kx-")
# ax3.set_ylabel("Number of failures", **ssfont)
# ax3.set_xlabel("SNR (dB)", **ssfont)
# ax1.set_title("The Error in Distance", **ssfont)
# plt.savefig("./pyramid_robot/noise_distance.png")
# plt.show()

# # Plot the results for the azimuth.
# fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
# ax1.plot(noises, np.rad2deg(θ_mean), "kx-")
# ax1.plot([-20,25], np.rad2deg([θ_true,θ_true]), "r-")
# ax1.set_ylabel("Mean azimuth (degree)", **ssfont)
# ax2.plot(noises, np.rad2deg(θ_var), "kx-")
# ax2.set_ylabel("Variance of azimuth (degree^2)", **ssfont)
# ax3.plot(noises, n_failures_data, "kx-")
# ax3.set_ylabel("Number of failures", **ssfont)
# ax3.set_xlabel("SNR (dB)", **ssfont)
# ax1.set_title("The Error in Azimuth", **ssfont)
# plt.savefig("./pyramid_robot/noise_azimuth.png")
# plt.show()

# # Plot the results for the azimuth.
# fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
# ax1.plot(noises, np.rad2deg(φ_mean), "kx-")
# ax1.plot([-20,25], np.rad2deg([φ_true,φ_true]), "r-")
# ax1.set_ylabel("Mean elevation (degree)", **ssfont)
# ax2.plot(noises, np.rad2deg(φ_var), "kx-")
# ax2.set_ylabel("Variance of elevation (degree^2)", **ssfont)
# ax3.plot(noises, n_failures_data, "kx-")
# ax3.set_ylabel("Number of failures", **ssfont)
# ax3.set_xlabel("SNR (dB)", **ssfont)
# ax1.set_title("The Error in Elevation", **ssfont)
# plt.savefig("./pyramid_robot/noise_elevation.png")
# plt.show()