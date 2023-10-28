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
log_angles = "./main_method/noise/whistle/logs/angles.csv"

# Clear the logs.
np.savetxt(log_angles, [])

# Read the input file.
# f_s, audio_anechoic = wavfile.read("./sounds/432910__kyanite__clap.wav")
# f_s, audio_anechoic = wavfile.read("./sounds/345__anton__handclaps.wav")
# f_s, audio_anechoic = wavfile.read("./sounds/78508__joedeshon__referee_whistle_01.wav")
f_s, audio_anechoic = wavfile.read("./sounds/418564__14fpanskabubik_lukas__whistle.wav")
# f_s, audio_anechoic_original = wavfile.read("./pyroomacoustics/examples/input_samples/cmu_arctic_us_aew_a0001.wav")
# audio_anechoic = signal.resample(audio_anechoic, audio_anechoic.shape[0]//2)
f_s = f_s//2

audio_anechoic = audio_anechoic/audio_anechoic.max()
audio_anechoic = audio_anechoic[:40000//2,0]

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
n_noises = 10
n_runs = 10**4

noises = np.arange(25,-25,-5)

# Declare arrays to log data.
angles_log = np.zeros((100,4))
Δθ_log = np.zeros((n_noises, n_runs))
Δφ_log = np.zeros((n_noises, n_runs))

# Set the room up.
chrono = time.time()
room = set_up_room(10**(-4))
print("Room done in", time.time() - chrono, "s.")

n_success = 0

all_timed = time.time()
noise_timed = all_timed

for i_noise in range(n_noises):

    log_timed = time.time()

    for i_run in range(n_runs):

        # Simulate the response from the source to the array.
        room.simulate(snr = noises[i_noise])

        # Load the output of the simulation. Each channel is from one 
        # microphone.
        audio_reverb = room.mic_array.signals.T
        length = audio_reverb.shape[0]

        # Make the frame from all microphones.
        start = np.clip(np.argmax(audio_reverb[:,0]) - F, 0, None) + F//2
        end = start + F
        x = audio_reverb[start:end, :]

        # # For the first few frames, only 
        # start = 0
        # end = 0
        # for i_start in range(0, 34000//F*F, F):
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
        θ, φ = estimator.estimate_direction(x, τ, grid)

        # Calculate the error.
        r_true = p_true - centre
        θ_true = np.arctan2(r_true[1], r_true[0])
        Δθ = np.abs(θ_true - θ).item()
        φ_true = np.arctan2(np.linalg.norm(r_true[0:2]), r_true[2])
        Δφ = np.abs(φ_true - φ).item()

        # Store the logs.
        angles_log[i_run % 100, :] = np.array([i_noise, i_run, θ, φ])
        Δθ_log[i_noise, i_run] = Δθ
        Δφ_log[i_noise, i_run] = Δφ

        # Log progress for every hundredth run.
        progress = i_noise*100.0/n_noises + (i_run+1)*100.0/n_runs/n_noises
        if (i_run + 1) % 100 == 0:
            # Print an update.
            print("{:.0f} dB, {:.0f}, {:.2f} %:"
                " {:.2f} deg, {:.2f} deg;"
                " time taken {:.2f} s, {:.2f} s, {:.2f} s;"
                " time left {:.2f} h"
                .format(
                noises[i_noise],
                i_run,
                progress,
                np.mean(np.rad2deg(Δθ_log[i_noise, (i_run-99):(i_run+1)])),
                # np.mean(Δθ_log[0, i_run]),
                np.mean(np.rad2deg(Δφ_log[i_noise, (i_run-99):(i_run+1)])),
                # np.mean(Δφ_log[0, i_run]),
                time.time() - log_timed,
                time.time() - noise_timed,
                time.time() - all_timed,
                (time.time() - all_timed)*(1/progress*100.0 - 1)/60/60
            ))
            
            # Log the results so far.
            with open(log_angles, "a") as f:
                np.savetxt(
                    f, 
                    angles_log, 
                    delimiter = ","
                )

            # Time the next log.
            log_timed = time.time()

    # Print that the noise-level has been done.
    print("{} noise-level of SNR = {:.1f} dB done"
        " and took {:.2f} s.".format(
        ordinal(i_noise+1),
        noises[i_noise],
        time.time() - noise_timed
    ))
    noise_timed = time.time()
