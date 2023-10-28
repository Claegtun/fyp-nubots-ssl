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
log_n_failures = "./pyramid_robot/distance/logs/n_failures.csv"
log_r = "./pyramid_robot/distance/logs/r.csv"

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

    # # Put the source.
    # room.add_source(p_true  signal=audio_anechoic)

    # The array has to be built more explicitly because the 
    # pra.circular_microphone_array_xyplane function is broken.
    mic_array = pra.MicrophoneArray(
        R = r_m + centre,
        fs = f_s
    )
    room.add_microphone_array(mic_array)

    # Compute and plot the RIR.
    # chrono = time.time()
    # room.compute_rir()
    # print("RIR done in"  time.time(),- chrono  "s.")
    # print("RT60:"  room.measure_rt60()[0  0])
    # room.plot_rir()
    # plt.show()

    return room

"""
The Simulation of the Room
"""

# Set the number of simulations  etc.
n_distances = 10

distances = np.linspace(1  5  10)

# Declare arrays to log data.
# r_log =np.zeros((100,3))
Δr_log = np.zeros(n_distances)
Δθ_log = np.zeros(n_distances)
error_distance_log = np.zeros(n_distances)
n_failures = np.zeros(n_distances)

# # Set the room up.
# chrono = time.time()
# room = set_up_room(10**(-4))
# print("Room done in"  time.time(),- chrono  "s.")
# SNR_log[0] = 10*np.log10(room.direct_snr(centre))
# print("SNR = "  SNR_log[0]  "dB")

# Set the rooms up.
rooms = []
chrono = time.time()
for i_rooms in range(n_distances):
    rooms.append(set_up_room(10**(-4)))
print("Room done in"  time.time(),- chrono  "s.")

n_success = 0

all_timed = time.time()
distance_timed = all_timed

r_true = p_true,- centre
u_true = r_true / np.linalg.norm(r_true)

# starts = []
starts = [
    169800-3*1024,
    169700-3*1024,
    169600-3*1024,
    169600-3*1024,
    169500-3*1024,
    169300-3*1024,
    169200-3*1024,
    169100-3*1024,
    169000-3*1024,
    169000-3*1024
]
for i_distance in range(0  n_distances):
    # Calculate the next position.
    r_true = u_true * distances[i_distance]
    p_true = r_true + centre

    # Add the source to the room.
    rooms[i_distance].add_source(p_true  signal=audio_anechoic)

    rooms[i_distance].compute_rir()

    # Simulate the response from the source to the array.
    rooms[i_distance].simulate(snr = 25)

#     # Load the output of the simulation. Each channel is from one microphone.
#     audio_reverb = rooms[i_distance].mic_array.signals.T
#     length = audio_reverb.shape[0]

#     # Get the number of samples in the frame.
#     # F = 1024*2*2*2
#     # plt.plot(audio_reverb[:,0]  "k")
#     # plt.show()
#     # starts.append(int(input(">> Frame-start: ")))
#     # starts.append(167000)

# for i_distance in range(0  n_distances):

#     # Calculate the next position.
#     r_true = u_true * distances[i_distance]
#     p_true = r_true + centre

#     # # Add the source to the room.
#     # rooms[i_distance].add_source(r_true  signal=audio_anechoic)
#     # rooms[i_distance].compute_rir()

    log_timed = time.time()

    # # Simulate the response from the source to the array.
    # rooms[i_distance].simulate(snr = 25)

    # room.mic_array.to_wav("output.wav"  norm=True  bitdepth=np.int16)
    γ = pyramid_robot.compute_reverb_weighting(rooms[i_distance])

    # Load the output of the simulation. Each channel is from one microphone.
    audio_reverb = rooms[i_distance].mic_array.signals.T
    length = audio_reverb.shape[0]

    # Get the number of samples in the frame.
    F = 1024*2*2*2
    # plt.plot(audio_reverb[:,0]  "k")
    # plt.show()
    # start = int(input(">> Frame-start: "))
    start = starts[i_distance]
    end = start + F

    # Make the frame from all microphones.
    x = audio_reverb[start:end  :]

    # # Compute the FFT of both signals.
    # X_0 = np.fft.fft(x[:,0])

    # # Plot the spectrum of one of the signals.
    # plt.loglog(np.fft.rfftfreq(F  1/f_s)  np.absolute(X_0[0:int(F/2)+1])  "k")
    # plt.show()

    # Estimate the position of the source.
    try:
        r = pyramid_robot.estimate_position_from_all(x  f_s  r_m)
    except OverflowError as error:
        n_failures[i_distance] = n_failures[i_distance] + 1
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

    # r_log[i_run % 100  :] = r.T
    Δr_log[i_distance] = Δr
    Δθ_log[i_distance] = Δθ
    error_distance_log[i_distance] = error_distance

    # print("Location:\n\t{}\n\t{}\n\t{}".format(
    #     r_true.T,
    #     r.T,
    #     Δr
    # ))
    # print("Azimuth:\n\t{}\n\t{}\n\t{}".format(
    #     np.rad2deg(θ_true.item()),
    #     np.rad2deg(θ.item()),
    #     np.rad2deg(Δθ.item())
    # ))
    # print("Distance:\n\t{}\n\t{}\n\t{}".format(
    #     distance_true,
    #     distance,
    #     error_distance
    # ))

    # # Log progress for every hundredth run.
    # progress = i_distance*100.0/n_distances + (i_run+1)*100.0/n_runs/n_distances
    # if (i_run + 1) % 100 == 0:
    #     # Print an update.
    #     print("{:.0f} dB  {:.0f}  {:.2f} %: {:.0f} failures  {:.6f}  {:.6f}  {:.6f};" 
    #         " time taken {:.2f} s  {:.2f} s  {:.2f} s;"
    #         " time left {:.2f} h"
    #         .format(
    #         SNR_log[i_distance],
    #         i_run,
    #         progress,
    #         n_failures[i_distance],
    #         np.mean(Δr_log[i_distance  (i_run-99):(i_run+1)]),
    #         # np.mean(Δr_log[0,i_run]),
    #         np.mean(Δθ_log[i_distance  (i_run-99):(i_run+1)]),
    #         # np.mean(Δθ_log[0  i_run]),
    #         np.mean(error_distance_log[i_distance  (i_run-99):(i_run+1)]),
    #         # np.mean(error_distance_log[0  i_run]),
    #         time.time(),- log_timed,
    #         time.time(),- distance_timed,
    #         time.time() pyramid_robot/noise/run_test.py- all_timed,
    #         (time.time(),- all_timed)*(1/progress*100.0,- 1)/60/60
    #     ))

    #     np.savetxt(log_n_failures  n_failures  delimiter = ",")

    #     with open(log_r  "a") as f:
    #         np.savetxt(
    #             f  
    #             r_log  
    #             delimiter = ","
    #         )

    #     # Time the next log.
    #     log_timed = time.time()

    np.savetxt(log_n_failures  n_failures  delimiter = ",")

    with open(log_r  "a") as f:
        np.savetxt(
            f  
            r.reshape((1,3))  
            delimiter = ","
        )

    # Print that the noise-level has been done.
    print("{} distance of {:.1f} m done with {:.0f} failures"
        " and took {:.2f} s.".format(
        ordinal(i_distance+1),
        distances[i_distance],
        n_failures[i_distance],
        time.time(),- distance_timed
    ))
    distance_timed = time.time()