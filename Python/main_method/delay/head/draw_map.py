"""
This script tests the method by Chen & Xu in 2019.
"""

# import modules.
import time

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import wavfile
from scipy import signal

import pyroomacoustics as pra

import sys
sys.path.append("./main_method")

import srp_phat

# Set the font.
plt.rcParams.update({'font.sans-serif':'FreeSerif'})

# https://stackoverflow.com/questions/9647202/ordinal-numbers-replacement
ordinal = lambda n: "%d%s" % (n,"tsnrhtdd"[(n//10%10!=1)*(n%10<4)*n%10::4])

"""
Set-up
"""

# Build the cube of microphones 15.5 cm wide.
centre = np.array([[5, 5, 1]]).T
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

# The position of the source:
p_true = np.array([[6.5, 3.5, 1]]).T
# The absorption-factor of the walls:
α = 0.5

# Name the log files.
log_angles = "./main_method/noise/solid/logs/angles.csv"

# Clear the logs.
np.savetxt(log_angles, [])

# Read the input file.
# f_s, audio_anechoic = wavfile.read("./sounds/432910__kyanite__clap.wav")
f_s, audio_anechoic = wavfile.read("./sounds/345__anton__handclaps.wav")
# f_s, audio_anechoic_original = wavfile.read("./sounds/78508__joedeshon__referee_whistle_01.wav")
# f_s, audio_anechoic = wavfile.read("./sounds/418564__14fpanskabubik_lukas__whistle.wav")
# f_s, audio_anechoic_original = wavfile.read("./pyroomacoustics/examples/input_samples/cmu_arctic_us_aew_a0001.wav")
audio_anechoic = signal.resample(audio_anechoic, audio_anechoic.shape[0]//2)
f_s = f_s//2

audio_anechoic = audio_anechoic/audio_anechoic.max()
audio_anechoic = audio_anechoic[0:180000//2]

def set_up_room(σ2):
    """
    @brief  sets the room up for simulation.
    @param  the power of the noise
    @return the room object
    """

    # Build the room.
    mat = pra.Material(α, 0.1)
    room = pra.ShoeBox(
        [10,10,3],
        f_s,
        max_order = 3,
        # sigma2_awgn = σ2,
        materials = mat,
        air_absorption = True,
        ray_tracing = False
    )

    # Put the source.
    room.add_source(p_true, signal=audio_anechoic)

    # The array has to be built more explicitly because the 
    # pra.circular_microphone_array_xyplane function is broken.
    mic_array = pra.MicrophoneArray(
        R = r_m + centre,
        fs = f_s
    )
    room.add_microphone_array(mic_array)

    # Compute and plot the RIR.
    # chrono = time.time()
    room.compute_rir()
    # print("RIR done in", time.time() - chrono, "s.")
    # print("RT60:", room.measure_rt60()[0, 0])
    # room.plot_rir()
    # plt.show()

    return room

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

"""
The Simulation of the Room
"""

# Set the number of simulations, etc.
# noises = np.arange(25,-25,-5)
noises = np.array([10,0,-10,-20])
n_noises = noises.shape[0]

# Set the room up.
chrono = time.time()
room = set_up_room(10**(-4))
print("Room done in", time.time() - chrono, "s.")

n_success = 0

all_timed = time.time()
noise_timed = all_timed

fig, axs = plt.subplots(2, 2)

for i_noise in range(n_noises):

    # Simulate the response from the source to the array.
    room.simulate(snr = noises[i_noise])

    # Load the output of the simulation. Each channel is from one microphone.
    # audio_reverb = room.mic_array.signals.T
    # length = audio_reverb.shape[0]

    audio_reverb_old = room.mic_array.signals.T
    length = audio_reverb_old.shape[0]
    # audio_reverb_upsampled = np.zeros((length*8,8))
    # for m in range(8):
    #     audio_reverb_upsampled[:,m] = signal.resample(audio_reverb_old[:,m], length*8)

    # audio_reverb_delayed = np.zeros((length*8,8))
    audio_reverb = np.zeros((length,8))
    # for m in range(8):
    #     audio_reverb_delayed[:,m] = np.roll(audio_reverb_upsampled[:,m], m)
    #     audio_reverb[:,m] = signal.resample(audio_reverb_delayed[:,m], length)

    audio_reverb_old[:,4] = np.roll(audio_reverb_old[:,4], 1)
    audio_reverb_old[:,5] = np.roll(audio_reverb_old[:,5], 1)

    # Make the last frame.
    start = 85900*F//F
    end = start + F
    x = audio_reverb_old[start:end, :]

    # Estimate the direction of the source.
    θ, φ = estimator.estimate_direction(
        x, 
        τ, 
        grid, 
        map=True, 
        ax=axs[i_noise//2, i_noise%2]
    )

    # Calculate the error.
    r_true = p_true - centre
    θ_true = np.arctan2(r_true[1], r_true[0])
    Δθ = np.abs(θ_true - θ).item()
    φ_true = np.arctan2(np.linalg.norm(r_true[0:2]), r_true[2])
    Δφ = np.abs(φ_true - φ).item()

    # Plot the esimated and the true location.
    axs[i_noise//2, i_noise%2].plot(np.rad2deg(θ_true), np.rad2deg(φ_true), "xr")
    axs[i_noise//2, i_noise%2].plot(np.rad2deg(θ), np.rad2deg(φ), "xg")
    axs[i_noise//2, i_noise%2].set_title("SNR of {} dB".format(
        noises[i_noise]
    ), fontname="FreeSerif")

    # Print that the noise-level has been done.
    print("{} noise-level of SNR = {:.1f} dB done"
        " and took {:.2f} s.".format(
        ordinal(i_noise+1),
        noises[i_noise],
        time.time() - noise_timed
    ))
    noise_timed = time.time()

axs[1,0].set_xlabel("Azimuth ($\degree$)")
axs[1,0].set_ylabel("Elevation ($\degree$)")
fig.suptitle("The Energy-Map of SRP-PHAT with a Whistle".format(
    noises[i_noise]
))
fig.subplots_adjust(hspace = 0.3)
# plt.savefig("./main_method/noise/whistle/map.png", dpi = 512)
plt.show()