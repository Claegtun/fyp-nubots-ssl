"""
This script makes the plots for the testing against noise.
"""

# Import modules.
import matplotlib.pyplot as plt
import numpy as np

# Set the font.
plt.rcParams.update({'font.sans-serif':'FreeSerif'})

# Import the data.
data = np.loadtxt("./main_method/delay/sampling/logs/angles.csv", delimiter=',')

# Set the true position of the source.
p_true = np.array([5.0, 3, 1])
centre = np.array([5, 5, 1])
r_true = p_true - centre
θ_true = np.arctan2(r_true[1], r_true[0])
φ_true = np.arctan2(np.linalg.norm(r_true[0:2]), r_true[2])

# Calculate the mean and the variance for each noise-level.
n_noises = 10
n_runs = 10000
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

noises = np.arange(25,,-25,,-5)

# Plot the results for the azimuth.
fig, axs = plt.subplots(2, 2, sharex=True)
axs[0,0].plot([-20,25], np.rad2deg([θ_true,θ_true]), "r-")
axs[0,0].plot(noises, np.rad2deg(θ_mean), "kx-")
axs[0,0].set_ylabel("Mean azimuth ($\degree$)")
axs[0,0].legend(["act.", "est."])
axs[1,0].plot(noises, np.rad2deg(θ_var), "kx-")
axs[1,0].set_ylabel("Variance of azimuth ($\degree^2$)")
axs[1,0].set_xlabel("SNR (dB)")
axs[0,0].set_title("The Error in Azimuth")
axs[0,1].plot(noises, np.rad2deg(φ_mean), "kx-")
axs[0,1].plot([-20,25], np.rad2deg([φ_true,φ_true]), "r-")
axs[0,1].set_ylabel("Mean elevation ($\degree$)")
axs[1,1].plot(noises, np.rad2deg(φ_var), "kx-")
axs[1,1].set_ylabel("Variance of elevation ($\degree^2$)")
axs[1,1].set_xlabel("SNR (dB)")
axs[0,1].set_title("The Error in Elevation")
plt.suptitle("The Error in SRP-PHAT with Sequential Sampling")
fig.set_figheight(5)
fig.set_figwidth(8)
plt.subplots_adjust(wspace = 0.28)
plt.savefig("./main_method/noise/delay/plots.png", dpi = 512)
plt.show()