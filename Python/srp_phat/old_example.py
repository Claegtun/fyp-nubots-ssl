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

# https://stackoverflow.com/questions/9647202/ordinal-numbers-replacement
ordinal = lambda n: "%d%s" % (n,"tsnrhtdd"[(n//10%10!=1)*(n%10<4)*n%10::4])

# Set the font.
ssfont = {"fontname":"Times New Roman"}

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
r_true = np.array([[5.5, 3, 1]]).T
# The speed of sound:
c = 343
# The absorption-factor of the walls:
α = 0.5

# Name the log files.
log_n_failures = "log_n_failures.csv"
log_r = "log_r.csv"
log_error = "log_error"

# # Clear the logs.
# np.savetxt(log_n_failures, [])
# with open(log_r) as f:
#     np.savetxt(f, [])

# Read the input file.
# f_s, audio_anechoic = wavfile.read("./sounds/432910__kyanite__clap.wav")
f_s, audio_anechoic = wavfile.read("./sounds/345__anton__handclaps.wav")
# f_s, audio_anechoic = wavfile.read("./sounds/78508__joedeshon__referee_whistle_01.wav")
# f_s, audio_anechoic = wavfile.read("./sounds/418564__14fpanskabubik_lukas__whistle.wav")

audio_anechoic = audio_anechoic/audio_anechoic.max()
audio_anechoic = audio_anechoic[0:200000]

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
    room.add_source(r_true, signal=audio_anechoic)

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

M = 8

def make_circular_grid():
    G = 144
    θ = np.deg2rad(2.5)
    v = np.zeros((3,G))
    v[:,0] = np.array([1,0,0])
    R = np.array([
        [np.cos(θ), -np.sin(θ), 0],
        [np.sin(θ), np.cos(θ), 0],
        [0, 0, 0]
    ])
    for g in range(1,G):
        v[:, g] = np.matmul(R, v[:,g-1])

    # fig = plt.figure()
    # ax = fig.add_subplot(projection='3d')
    # ax.scatter(v[0,:], v[1,:], v[2,:])
    # plt.show()

    return v


def make_spherical_grid(N_ev):
    """
    @brief  Makes the spherical grid as described in the papers by the French 
            Canadians, namely by tesselating an icosahedron.
    @param  N_ev: the number of evolutions,
    @return ndarray (3,?): an array of the points,
    """

    # Make an icosahedron with which to being.
    φ_r = (1 + np.sqrt(5))/2
    v = np.array([
        [0,     1,      φ_r ],
        [1,     φ_r,    0   ],
        [φ_r,   0,      1   ],
        [0,     -1,     φ_r ],
        [-1,    φ_r,    0   ],
        [φ_r,   0,      -1  ],
        [0,     1,      -φ_r],
        [1,     -φ_r,   0   ],
        [-φ_r,  0,      1   ],
        [0,     -1,     -φ_r],
        [-1,    -φ_r,   0   ],
        [-φ_r,  0,      -1  ]
    ]).T

    # Normalise the points.
    v = v / np.linalg.norm(v, axis = 0).reshape((1,12))

    # Calculate the shortest distance between two points.
    d_span = np.linalg.norm(v[:,0] - v[:,3])

    # Tesselate the grid for however many evolutions.
    for k in range(N_ev):
        # For each point and for each of its nearest neighbours and make 
        # another point between them.
        n_v = v.shape[1]
        for i in range(n_v):
            for j in range(i+1,n_v):
                # See if the two points are nearest neighbours by seeing 
                # whether they are withing the shortest distance.
                d = np.linalg.norm(v[:,i] - v[:,j])
                if d >= 1.5*d_span/(2**k):
                    continue
                # Make the new point and add it to the array.
                u = (v[:,i] + v[:,j]).reshape((3,1))
                u = u / np.linalg.norm(u)
                v = np.hstack((v,u))

    # print("{} points".format(v.shape[1]))

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(v[0,:], v[1,:], v[2,:], c="k")
    ax.set_xlabel("x-coordinate (m)", **ssfont)
    ax.set_ylabel("y-coordinate (m)", **ssfont)
    ax.set_zlabel("z-coordinate (m)", **ssfont)
    plt.show()

    return v

def compute_tdoa(u):
    """
    @brief  Computes the TDOA for every pair given the direction.
    @param  u (3,1): the direction,
    @return ndarray (M,M): the TDOAs in number of samples,
    """
    τ = np.zeros((M,M))

    # For each pair of microphones, compute the TDOA given by the formula in 
    # Valin et al. 2007.
    for i in range(M):
        for j in range(M):
            if i == j:
                continue
            τ[i,j] = np.dot(r_m[:,j] - r_m[:,i], u) * f_s / c
    return τ

def compute_tdoa_grid(grid):
    """
    @brief  Computes the TDOA for every pair and for every direction in the 
            grid.
    @param  grid (3,G): the sperical grid,
    @return ndarray (G,M,M): the TDOAs in number of samples,
    """

    G = grid.shape[1]
    τ = np.zeros((G,M,M))

    # For every direction in the grid, compute the TDOAs.
    for g in range(G):
        τ[g] = compute_tdoa(grid[:,g])

    return τ

# Compute all TDOAs.
grid = make_spherical_grid(4)
# grid = make_circular_grid()
G = grid.shape[1]
τ = compute_tdoa_grid(grid)

def compute_gcc_phat(x_0, x_1):
    """
    @brief  Computes the GCC-PHAT from two frames from two microphones.
    @param  x_0 (F, 1): the first channel
    @param  x_1 (F, 1): the second channel
    @return ndarray (F, 1): the GCC-PHAT
    """

    # Work out the frame-length.
    F = x_0.shape[0]

    # Compute the FFT of both signals.
    X_0 = np.fft.fft(x_0)
    X_1 = np.fft.fft(x_1)

    # Compute the mean across the whole spectrum and make a weighting based 
    # thereon.
    X_m = (np.absolute(X_0) + np.absolute(X_1))/2
    X_n = np.mean(X_m)
    w = (X_m >= X_n) + (X_m < X_n)*(X_m/X_n)**0.3

    # Compute the spectral weighting, i.e. PHAT.
    φ = 1 / (np.absolute(X_0) * np.absolute(X_1))
    # φ = w**2 / (γ * np.absolute(X_0) * np.absolute(X_1))
    # φ = 1

    # Compute the GCC.
    χ = np.multiply(φ, np.multiply(X_0, np.conj(X_1)))
    R_01 = np.real(np.fft.ifft(χ))

    return R_01

def compute_beam_energy(x, g, R):
    """
    @brief  Computes the beamformer's energy along the g-th direction in the 
            grid.
    @param  x (F, M): the frame with all channels,
    @param  g: the index of the direction along which to compute,
    @return float: the energy,
    """
    
    E = 0
    for i in range(M):
        for j in range(M):
            if i == j:
                continue
            E = E + R[i,j,int(np.rint(τ[g,i,j]))]
    return E

def estimate_direction(x):
    """
    @brief  Estimates the position of the source from one set of equations given 
            a microphone as a reference.
    @param  x (F, M): the frame with all channels,
    @return tuple: the estimated azimuth and elevation in radians,
    """

    # Compute the GCC-PHATs for all pairs.
    F = x.shape[0]
    R = np.zeros((M,M,F))
    for i in range(M):
        for j in range(M):
            if i == j:
                continue
            R[i,j,:] = compute_gcc_phat(x[:,i], x[:,j])

    # For each direction in the grid, compute the energy and see whether it is 
    # the largest energy.
    g_max = 0
    E_max = 0
    for g in range(G):
        E = compute_beam_energy(x, g, R)
        if E >= E_max:
            g_max = g
            E_max = E

    u = grid[:, g_max]
    θ = np.arctan2(u[1], u[0])
    φ = np.arctan2(np.linalg.norm(u[0:2]), u[2])

    return (θ, φ)

# Set the room up.
chrono = time.time()
room = set_up_room(10**(-4))
print("Room done in", time.time() - chrono, "s.")

# Simulate the response from the source to the array.
room.simulate(snr = 25)

# room.mic_array.to_wav("output.wav", norm=True, bitdepth=np.int16)

# Load the output of the simulation. Each channel is from one microphone.
audio_reverb = room.mic_array.signals.T
length = audio_reverb.shape[0]

# Get the number of samples in the frame.
F = 1024
plt.plot(audio_reverb[:,0], "k")
plt.show()
start = int(input(">> Frame-start: "))
# start = 172000
end = start + F

# Make the frame from all microphones.
x = audio_reverb[start:end, :]

θ, φ = estimate_direction(x)

print("θ = {:.1f} deg".format(np.rad2deg(θ)))
print("φ = {:.1f} deg".format(np.rad2deg(φ)))


# """
# The Simulation of the Room
# """

# # Set the number of simulations, etc.
# n_noises = 10
# n_runs = 10**4

# noises = np.arange(25,-25,-5)

# # Declare arrays to log data.
# r_log =np.zeros((100,3))
# Δr_log = np.zeros((n_noises, n_runs))
# Δθ_log = np.zeros((n_noises, n_runs))
# error_distance_log = np.zeros((n_noises, n_runs))
# n_failures = np.zeros(n_noises)
# SNR_log = np.zeros((n_noises))

# # Set the room up.
# chrono = time.time()
# room = set_up_room(10**(-4))
# print("Room done in", time.time() - chrono, "s.")
# # SNR_log[0] = 10*np.log10(room.direct_snr(centre))
# print("SNR = ", SNR_log[0], "dB")

# n_success = 0

# all_timed = time.time()
# noise_timed = all_timed

# for i_noise in range(n_noises):

#     SNR_log[i_noise] = noises[i_noise]

#     log_timed = time.time()

#     for i_run in range(n_runs):

#         # Simulate the response from the source to the array.
#         room.simulate(snr = noises[i_noise])

#         # room.mic_array.to_wav("output.wav", norm=True, bitdepth=np.int16)
#         γ = compute_reverb_weighting(room)

#         # Load the output of the simulation. Each channel is from one microphone.
#         audio_reverb = room.mic_array.signals.T
#         length = audio_reverb.shape[0]

#         # Get the number of samples in the frame.
#         F = 1024*2*2*2
#         # plt.plot(audio_reverb[:,0], "k")
#         # plt.show()
#         # start = int(input(">> Frame-start: "))
#         start = 171500
#         end = start + F

#         # Make the frame from all microphones.
#         x = audio_reverb[start:end, :]

#         # # Compute the FFT of both signals.
#         # X_0 = np.fft.fft(x[:,0])

#         # # Plot the spectrum of one of the signals.
#         # plt.loglog(np.fft.rfftfreq(F, 1/f_s), np.absolute(X_0[0:int(F/2)+1]), "k")
#         # plt.show()

#         # Estimate the position of the source.
#         try:
#             r = estimate_position_from_all(x, f_s) + centre
#         except OverflowError as error:
#             n_failures[i_noise] = n_failures[i_noise] + 1
#             r = np.zeros((3,1))

#         # Calculate the error.
#         Δr = np.linalg.norm(r_true - r)
#         θ_true = np.arctan2(r_true[1]-centre[0], r_true[0]-centre[0])
#         θ = np.arctan2(r[1]-centre[1], r[0]-centre[0])
#         Δθ = np.abs(θ_true - θ)
#         distance_true = np.linalg.norm(r_true - centre)
#         distance = np.linalg.norm(r - centre)
#         error_distance = np.abs(distance_true - distance)

#         r_log[i_run % 100, :] = r.T
#         Δr_log[i_noise, i_run] = Δr
#         Δθ_log[i_noise, i_run] = Δθ
#         error_distance_log[i_noise, i_run] = error_distance

#         # print("Location:\n\t{}\n\t{}\n\t{}".format(
#         #     r_true.T,
#         #     r.T,
#         #     Δr
#         # ))
#         # print("Azimuth:\n\t{}\n\t{}\n\t{}".format(
#         #     np.rad2deg(θ_true.item()),
#         #     np.rad2deg(θ.item()),
#         #     np.rad2deg(Δθ.item())
#         # ))
#         # print("Distance:\n\t{}\n\t{}\n\t{}".format(
#         #     distance_true,
#         #     distance,
#         #     error_distance
#         # ))

#         # Log progress for every hundredth run.
#         progress = i_noise*100.0/n_noises + (i_run+1)*100.0/n_runs/n_noises
#         if (i_run + 1) % 100 == 0:
#             # Print an update.
#             print("{:.0f} dB, {:.0f}, {:.2f} %: {:.0f} failures, {:.6f}, {:.6f}, {:.6f};" 
#                 " time taken {:.2f} s, {:.2f} s, {:.2f} s;"
#                 " time left {:.2f} h"
#                 .format(
#                 SNR_log[i_noise],
#                 i_run,
#                 progress,
#                 n_failures[i_noise],
#                 np.mean(Δr_log[i_noise, (i_run-99):(i_run+1)]),
#                 # np.mean(Δr_log[0,i_run]),
#                 np.mean(Δθ_log[i_noise, (i_run-99):(i_run+1)]),
#                 # np.mean(Δθ_log[0, i_run]),
#                 np.mean(error_distance_log[i_noise, (i_run-99):(i_run+1)]),
#                 # np.mean(error_distance_log[0, i_run]),
#                 time.time() - log_timed,
#                 time.time() - noise_timed,
#                 time.time() - all_timed,
#                 (time.time() - all_timed)*(1/progress*100.0 - 1)/60/60
#             ))

#             np.savetxt(log_n_failures, n_failures, delimiter = ",")

#             with open(log_r, "a") as f:
#                 np.savetxt(
#                     f, 
#                     r_log, 
#                     delimiter = ","
#                 )

#             # Time the next log.
#             log_timed = time.time()

#     # Print that the noise-level has been done.
#     print("{} noise-level of SNR = {:.1f} dB done with {:.0f} failures"
#         " and took {:.2f} s.".format(
#         ordinal(i_noise+1),
#         SNR_log[i_noise],
#         n_failures[i_noise],
#         time.time() - noise_timed
#     ))
#     noise_timed = time.time()

# # Save the errors in a file.
# np.savez(log_error, r=Δr_log, azimuth=Δθ_log, distance=error_distance_log)
