"""
This script makes the plots for the testing against frequency.
"""

# Import modules.
import matplotlib.pyplot as plt
import numpy as np

# Set the font.
plt.rcParams.update({'font.sans-serif':'FreeSerif'})

# Import the data.
data = np.loadtxt("./main_method/frequency/clap/logs/angles.csv", delimiter=',')

# Set the true position of the source.
p_true = np.array([5.5, 3, 2])
centre = np.array([5, 5, 1])
r_true = p_true - centre
θ_true = np.arctan2(r_true[1], r_true[0])
φ_true = np.arctan2(np.linalg.norm(r_true[0:2]), r_true[2])

# Calculate the mean and the variance for each frequency.
n_frequencies = 6
n_runs = 10000
θ_mean = np.zeros(n_frequencies)
φ_mean = np.zeros(n_frequencies)
θ_var = np.zeros(n_frequencies)
φ_var = np.zeros(n_frequencies)
for i_frequency in range(n_frequencies):
    θ = data[i_frequency*n_runs:(1+i_frequency)*n_runs,2]
    φ = data[i_frequency*n_runs:(1+i_frequency)*n_runs,3]
    θ_mean[i_frequency] = np.mean(θ)
    φ_mean[i_frequency] = np.mean(φ)
    θ_var[i_frequency] = np.var(θ)
    φ_var[i_frequency] = np.var(φ)

frequencies = np.array([16, 24, 48, 72, 96, 120])

# Plot the results for the azimuth.
fig, axs = plt.subplots(2, 2, sharex=True)
axs[0,0].plot([0,120], np.rad2deg([θ_true,θ_true]), "r-")
axs[0,0].plot(frequencies, np.rad2deg(θ_mean), "kx-")
axs[0,0].set_ylabel("Mean azimuth ($\degree$)")
axs[0,0].legend(["act.", "est."])
axs[1,0].plot(frequencies, np.rad2deg(θ_var), "kx-")
axs[1,0].set_ylabel("Variance of azimuth ($\degree^2$)")
axs[1,0].set_xlabel("Sampling frequency (kHz)")
axs[0,0].set_title("The Error in Azimuth")
axs[0,1].plot(frequencies, np.rad2deg(φ_mean), "kx-")
axs[0,1].plot([0,120], np.rad2deg([φ_true,φ_true]), "r-")
axs[0,1].set_ylabel("Mean elevation ($\degree$)")
axs[1,1].plot(frequencies, np.rad2deg(φ_var), "kx-")
axs[1,1].set_ylabel("Variance of elevation ($\degree^2$)")
axs[1,0].set_xlabel("Sampling frequency (kHz)")
axs[0,1].set_title("The Error in Elevation")
plt.xticks(frequencies)
plt.suptitle("The Simulated Error against Sampling Frequency")
fig.set_figheight(5)
fig.set_figwidth(8)
plt.subplots_adjust(wspace = 0.28)
plt.savefig("./main_method/frequency/clap/plots.png", dpi = 1024)
plt.show()