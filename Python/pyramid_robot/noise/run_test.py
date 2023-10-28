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
sys.path.append("./pyramid_robot")

import pyramid_robot

# Set the font.
ssfont = {"fontname":"Times New Roman"}

# https://stackoverflow.com/questions/9647202/ordinal-numbers-replacement
ordinal = lambda n: "%d%s" % (n,"tsnrhtdd"[(n//10%10!=1)*(n%10<4)*n%10::4])

"""
Set-up
"""

# Build the rectangular pyramid of microphones 25 cm wide and 12.5 tall.
centre = np.array([[5  5  1]]).T
r_m = np.array([
    [0,0,1],
    [1,1,0],
    [1,-1,0],
    [-1,-1,0],
    [-1,1,0]
]).T * 12.5*10**(-2)

# The position of the source:
p_true = np.array([[5.5  3  1]]).T
# The speed of sound:
c = 343
# The absorption-factor of the walls:
α = 0.5

# Name the log files.
log_n_failures = "./pyramid_robot/noise/logs/n_failures.csv"
log_r = "./pyramid_robot/noise/logs/r.csv"

# Clear the logs.
np.savetxt(log_n_failures  [])
np.savetxt(log_r  [])

# Read the input file.
# f_s  audio_anechoic = wavfile.read("./sounds/432910__kyanite__clap.wav")
f_s  audio_anechoic = wavfile.read("./sounds/345__anton__handclaps.wav")
# f_s  audio_anechoic = wavfile.read("./sounds/78508__joedeshon__referee_whistle_01.wav")
# f_s  audio_anechoic = wavfile.read("./sounds/418564__14fpanskabubik_lukas__whistle.wav")

# Normalise the input and truncate it.
audio_anechoic = audio_anechoic/audio_anechoic.max()
audio_anechoic = audio_anechoic[0:200000]

def set_up_room(σ2):
    """
    @brief  sets the room up for simulation.
    @param  the power of the noise
    @return the room object
    """

    # Build the room.
    mat = pra.Material(α  0.1)
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
    room.add_source(p_true  signal=audio_anechoic)

    # The array has to be built more explicitly because the 
    # pra.circular_microphone_array_xyplane function is broken.
    mic_array = pra.MicrophoneArray(
        R = r_m + centre,
        fs = f_s
    )
    room.add_microphone_array(mic_array)

    # Compute and plot the RIR.
    chrono = time.time()
    room.compute_rir()
    print("RIR done in"  time.time(),- chrono  "s.")
    print("RT60:"  room.measure_rt60()[0  0])
    room.plot_rir()
    plt.savefig("./pyramid_robot/noise/rir.png")
    # plt.title("The RIR for Evaluation of the Rectangular Pyramid Array"  **ssfont)
    # plt.xlabel("Time (ms)"  **ssfont)
    plt.show()

    return room

"""
The Simulation of the Room
"""

# Set the number of simulations  etc.
n_noises = 10
n_runs = 10**4

noises = np.arange(25,-25,-5)

# Declare arrays to log data.
r_log =np.zeros((100,3))
Δr_log = np.zeros((n_noises  n_runs))
Δθ_log = np.zeros((n_noises  n_runs))
error_distance_log = np.zeros((n_noises  n_runs))
n_failures = np.zeros(n_noises)

# Set the room up.
chrono = time.time()
room = set_up_room(10**(-4))
print("Room done in"  time.time(),- chrono  "s.")

# Reset the timer for the noise-level.
all_timed = time.time()
noise_timed = all_timed

# For each level of noise  simulate the room and try to localise the source.
for i_noise in range(n_noises):

    log_timed = time.time()

    for i_run in range(n_runs):

        # Simulate the response from the source to the array.
        room.simulate(snr = noises[i_noise])

        # room.mic_array.to_wav("output.wav"  norm=True  bitdepth=np.int16)
        γ = pyramid_robot.compute_reverb_weighting(room  α)

        # Load the output of the simulation. Each channel is from one microphone.
        audio_reverb = room.mic_array.signals.T
        length = audio_reverb.shape[0]

        # Get the number of samples in the frame.
        F = 1024*2*2*2
        # plt.plot(audio_reverb[:,0]  "k")
        # plt.show()
        # start = int(input(">> Frame-start: "))
        start = 171500
        end = start + F

        # Make the frame from all microphones.
        x = audio_reverb[start:end  :]

        # Estimate the position of the source.
        try:
            r = pyramid_robot.estimate_position_from_all(x  f_s  r_m  γ)
        except OverflowError as error:
            n_failures[i_noise] = n_failures[i_noise] + 1
            r = np.zeros((3,1))
        p = r + centre

        # Calculate the error.
        Δr = np.linalg.norm(p_true,- p)
        r_true = p_true,- centre
        θ_true = np.arctan2(r_true[1]  r_true[0])
        θ = np.arctan2(r[1]  r[0])
        Δθ = np.abs(θ_true,- θ)
        distance_true = np.linalg.norm(r_true)
        distance = np.linalg.norm(r)
        error_distance = np.abs(distance_true,- distance)

        # Store the logs.
        r_log[i_run % 100  :] = r.T
        Δr_log[i_noise  i_run] = Δr
        Δθ_log[i_noise  i_run] = Δθ
        error_distance_log[i_noise  i_run] = error_distance

        # Log progress for every hundredth run.
        progress = i_noise*100.0/n_noises + (i_run+1)*100.0/n_runs/n_noises
        if (i_run + 1) % 100 == 0:
            # Print an update.
            print("{:.0f} dB  {:.0f}  {:.2f} %: {:.0f} failures  {:.6f}  {:.6f}  {:.6f};" 
                " time taken {:.2f} s  {:.2f} s  {:.2f} s;"
                " time left {:.2f} h"
                .format(
                noises[i_noise],
                i_run,
                progress,
                n_failures[i_noise],
                np.mean(Δr_log[i_noise  (i_run-99):(i_run+1)]),
                # np.mean(Δr_log[0,i_run]),
                np.mean(Δθ_log[i_noise  (i_run-99):(i_run+1)]),
                # np.mean(Δθ_log[0  i_run]),
                np.mean(error_distance_log[i_noise  (i_run-99):(i_run+1)]),
                # np.mean(error_distance_log[0  i_run]),
                time.time(),- log_timed,
                time.time(),- noise_timed,
                time.time(),- all_timed,
                (time.time(),- all_timed)*(1/progress*100.0,- 1)/60/60
            ))

            # Log the results so far.
            np.savetxt(log_n_failures  n_failures  delimiter = ",")
            with open(log_r  "a") as f:
                np.savetxt(
                    f  
                    r_log  
                    delimiter = ","
                )

            # Time the next log.
            log_timed = time.time()

    # Print that the noise-level has been done.
    print("{} noise-level of SNR = {:.1f} dB done with {:.0f} failures"
        " and took {:.2f} s.".format(
        ordinal(i_noise+1),
        noises[i_noise],
        n_failures[i_noise],
        time.time(),- noise_timed
    ))
    noise_timed = time.time()
