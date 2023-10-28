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
]).T * 91.2*10**(-3) / 2

# The position of the source:
p_true = np.array([[5.5, 3, 2]]).T
# The absorption-factor of the walls:
α = 0.5

# Name the log files.
log_angles = "./main_method/frequency/clap/logs/angles.csv"

# Clear the logs.
np.savetxt(log_angles, [])

# Read the input file.
# f_s, audio_anechoic = wavfile.read("./sounds/432910__kyanite__clap.wav")
f_s, audio_anechoic = wavfile.read("./sounds/345__anton__handclaps.wav")
# f_s, audio_anechoic_original = wavfile.read("./sounds/78508__joedeshon__referee_whistle_01.wav")
# f_s, audio_anechoic_original = wavfile.read("./sounds/418564__14fpanskabubik_lukas__whistle.wav")
# f_s, audio_anechoic_original = wavfile.read("./pyroomacoustics/examples/input_samples/cmu_arctic_us_aew_a0001.wav")
# audio_anechoic = signal.resample(audio_anechoic, audio_anechoic.shape[0]//2)
audio_anechoic = signal.decimate(audio_anechoic, 2)
f_s = f_s//2

audio_anechoic = audio_anechoic/audio_anechoic.max()
audio_anechoic = audio_anechoic[0:180000//2]

audio_anechoic_16 = signal.decimate(audio_anechoic, 3)
audio_anechoic_24 = signal.decimate(audio_anechoic, 2)
audio_anechoic_48 = audio_anechoic
audio_anechoic_72 = signal.resample(audio_anechoic, audio_anechoic.shape[0]*3//2)
audio_anechoic_96 = signal.resample(audio_anechoic, audio_anechoic.shape[0]*2)
audio_anechoic_120 = signal.resample(audio_anechoic, audio_anechoic.shape[0]*5//2)

frequencies = np.array([ f_s//2, f_s, f_s*3//2, f_s*2])

def set_up_room(audio, f_s):
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
    room.add_source(p_true, signal=audio)

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
The Simulation of the Room
"""

# Set the number of simulations, etc.
n_freqencies = frequencies.shape[0]

n_success = 0

all_timed = time.time()
frequency_timed = all_timed

fig, axs = plt.subplots(2, 2)

for i_frequency in range(n_freqencies):

    # Get the number of samples in the frame.
    F = 1024

    # Make the object for the SRP-PHAT.
    estimator = srp_phat.srp_phat(r_m, frequencies[i_frequency], F)

    # Compute the grid and all TDOAs.
    grid = estimator.make_spherical_grid(4)
    # grid = make_circular_grid()
    τ = estimator.compute_tdoa_grid(grid)

    # Set the room up.
    chrono = time.time()
    if i_frequency == 0:
        room = set_up_room(audio_anechoic_24, frequencies[i_frequency])
    elif i_frequency == 1:
        room = set_up_room(audio_anechoic_48, frequencies[i_frequency])
    elif i_frequency == 2:
        room = set_up_room(audio_anechoic_72, frequencies[i_frequency])
    else:
        room = set_up_room(audio_anechoic_96, frequencies[i_frequency])
    print("Room done in", time.time() - chrono, "s.")

    # Simulate the response from the source to the array.
    room.simulate(snr = 0)

    # Load the output of the simulation. Each channel is from one microphone.
    audio_reverb = room.mic_array.signals.T
    length = audio_reverb.shape[0]

    # Make the frame from all microphones.
    start = np.clip(np.argmax(audio_reverb[:,0]) - F, 0, None) + F//2
    end = start + F
    x = audio_reverb[start:end, :]

    # # For the first few frames, only 
    # start = 0
    # end = 0
    # for i_start in range(0, 85900*F//F, F):
    #     start = i_start
    #     end = start + F

    #     # Make the frame from all microphones.
    #     x = audio_reverb[start:end, :]

    #     # Calculate the noise in this frame.
    #     estimator.calculate_noise(x)

    # # Make the last frame.
    # start = start + F
    # end = end + F
    # x = audio_reverb[start:end, :]

    # Estimate the direction of the source.
    θ, φ = estimator.estimate_direction(
        x, 
        τ, 
        grid, 
        map=True, 
        ax=axs[i_frequency//2, i_frequency%2]
    )

    # Calculate the error.
    r_true = p_true - centre
    θ_true = np.arctan2(r_true[1], r_true[0])
    Δθ = np.abs(θ_true - θ).item()
    φ_true = np.arctan2(np.linalg.norm(r_true[0:2]), r_true[2])
    Δφ = np.abs(φ_true - φ).item()

    # Plot the esimated and the true location.
    axs[i_frequency//2, i_frequency%2].plot(np.rad2deg(θ_true), np.rad2deg(φ_true), "xr")
    axs[i_frequency//2, i_frequency%2].plot(np.rad2deg(θ), np.rad2deg(φ), "xg")
    axs[i_frequency//2, i_frequency%2].set_title("Frequency of {} kHz".format(
        frequencies[i_frequency]/1000
    ), fontname="FreeSerif")

    # Print that the frequency has been done.
    print("{} frequency of f_s = {:.0f} kHz done"
        " and took {:.2f} s.".format(
        ordinal(i_frequency+1),
        frequencies[i_frequency]/1000,
        time.time() - frequency_timed
    ))
    frequency_timed = time.time()

axs[1,0].set_xlabel("Azimuth ($\degree$)")
axs[1,0].set_ylabel("Elevation ($\degree$)")
fig.suptitle("The Energy-Map against Sampling Frequency with a Clap".format(
    frequencies[i_frequency]
))
fig.subplots_adjust(hspace = 0.3)
plt.savefig("./main_method/frequency/clap/map.png", dpi = 1024)
plt.show()