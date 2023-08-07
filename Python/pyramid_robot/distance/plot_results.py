"""
This script makes the plots for the testing against noise.
"""

# Import modules.
import matplotlib.pyplot as plt
import numpy as np

# # Set the font.
# ssfont = {"fontname":"Times New Roman"}

# Import the data.
r_data = np.loadtxt("./pyramid_robot/distance/logs/r.csv", delimiter=',')
n_failures_data = np.loadtxt("./pyramid_robot/distance/logs/n_failures.csv", delimiter=',')

# # Ignore the first hundred. Something weird happened there.
# r_data = r_data[100:,:]

# Find where the failures happened, i.e. where the position is just the origin.
i_failure = ((r_data[:,0] != 0) | (r_data[:,1] != 0) | (r_data[:,2] != 0))\
    .reshape((1000,1))
i_failure = np.repeat(i_failure, 3, axis = 1)

# Set the true position of the source.
p_true = np.array([5.5, 3, 1])
centre = np.array([5, 5, 1])
r_true = p_true - centre
d_true = np.linalg.norm(r_true)
θ_true = np.arctan2(r_true[1], r_true[0])
φ_true = np.arctan2(np.linalg.norm(r_true[0:2]), r_true[2])

# Calculate the mean and the variance for each noise-level.
n_distances = 10
n_runs = 100
d_mean = np.zeros(n_distances)
θ_mean = np.zeros(n_distances)
φ_mean = np.zeros(n_distances)
d_var = np.zeros(n_distances)
θ_var = np.zeros(n_distances)
φ_var = np.zeros(n_distances)
for i_distance in range(n_distances):
    # Choose the next hundred.
    r_noise = r_data[i_distance*n_runs:(1+i_distance)*n_runs]
    i_clean = i_failure[i_distance*n_runs:(1+i_distance)*n_runs,:]
    # Clean out any samples from failures.
    r_clean = r_noise[i_clean[:,0]]
    d = np.linalg.norm(r_clean, axis = 1)
    θ = np.arctan2(r_clean[:,1], r_clean[:,0])
    φ = np.arctan2(np.linalg.norm(r_clean[:,0:2], axis = 1), r_clean[:,2])
    d_mean[i_distance] = np.mean(d)
    θ_mean[i_distance] = np.mean(θ)
    φ_mean[i_distance] = np.mean(φ)
    d_var[i_distance] = np.var(d)
    θ_var[i_distance] = np.var(θ)
    φ_var[i_distance] = np.var(φ)

distances = np.linspace(1, 5, 10)

# Plot the results for the distance.
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
ax1.plot(distances, d_mean, "kx-")
ax1.plot([1,5], [distances[0],distances[-1]], "r-")
ax1.set_ylabel("Mean distance (m)")#, **ssfont)
ax2.plot(distances, d_var, "kx-")
ax2.set_ylabel("Variance of \ndistance (m^2)")#, **ssfont)
ax3.plot(distances, n_failures_data, "kx-")
ax3.set_ylabel("Number of failures")#, **ssfont)
ax3.set_xlabel("Distance (m)")#, **ssfont)
ax1.set_title("The Error in Distance")#, **ssfont)
plt.savefig("./pyramid_robot/distance/distance.png", dpi = 512)
plt.show()

# Plot the results for the azimuth.
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
ax1.plot(distances, np.rad2deg(θ_mean), "kx-")
ax1.plot([1,5], np.rad2deg([θ_true,θ_true]), "r-")
ax1.set_ylabel("Mean azimuth \n(degree)")#, **ssfont)
ax2.plot(distances, np.rad2deg(θ_var), "kx-")
ax2.set_ylabel("Variance of \nazimuth (degree^2)")#, **ssfont)
ax3.plot(distances, n_failures_data, "kx-")
ax3.set_ylabel("Number of failures")#, **ssfont)
ax3.set_xlabel("Distance (m)")#, **ssfont)
ax1.set_title("The Error in Azimuth")#, **ssfont)
plt.savefig("./pyramid_robot/distance/azimuth.png", dpi = 512)
plt.show()

# Plot the results for the azimuth.
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
ax1.plot(distances, np.rad2deg(φ_mean), "kx-")
ax1.plot([1,5], np.rad2deg([φ_true,φ_true]), "r-")
ax1.set_ylabel("Mean elevation \n(degree)")#, **ssfont)
ax2.plot(distances, np.rad2deg(φ_var), "kx-")
ax2.set_ylabel("Variance of \nelevation (degree^2)")#, **ssfont)
ax3.plot(distances, n_failures_data, "kx-")
ax3.set_ylabel("Number of failures")#, **ssfont)
ax3.set_xlabel("Distance (m)")#, **ssfont)
ax1.set_title("The Error in Elevation")#, **ssfont)
plt.savefig("./pyramid_robot/distance/elevation.png", dpi = 512)
plt.show()