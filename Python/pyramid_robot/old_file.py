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

"""
Set-up
"""

# # Read the input file.
# f_s_old  audio_anechoic_original = wavfile.read("./sounds/432910__kyanite__clap.wav")
f_s  audio_anechoic = wavfile.read("./sounds/345__anton__handclaps.wav")
# f_s  audio_anechoic = wavfile.read("./sounds/78508__joedeshon__referee_whistle_01.wav")
# f_s  audio_anechoic = wavfile.read("./sounds/418564__14fpanskabubik_lukas__whistle.wav")

# audio_anechoic = audio_anechoic[:,0]

# # Upsample the input.
# audio_anechoic = signal.resample_poly(audio_anechoic_original  1000  441)
# f_s = 100000

# l_down = audio_anechoic_original.shape[0]
# t_down = np.linspace(0  l_down/f_s_old  l_down)

# l_up = audio_anechoic.shape[0]
# t_up = np.linspace(0  l_up/f_s  l_up)

# plt.plot(t_down  audio_anechoic_original)
# plt.plot(t_up  audio_anechoic)
# plt.show()

audio_anechoic = audio_anechoic/audio_anechoic.max()
audio_anechoic = audio_anechoic[0:200000]

# Build the room.
SNR =,-50
α = 0.5
mat = pra.Material(α  0.1)
# room = pra.AnechoicRoom(
room = pra.ShoeBox(
    [10,10,3],
    f_s,
    max_order = 3,
    # sigma2_awgn = 10**(-SNR/10),
    materials = mat,
    air_absorption = True,
    ray_tracing = False
)

# Put the source.
r_true = np.array([[6  1  1]]).T
room.add_source(r_true  signal=audio_anechoic)

# Build the rectangular pyramid of microphones 25 cm wide and 12.5 tall.
r_m = np.array([
    [0,0,1],
    [1,1,0],
    [1,-1,0],
    [-1,-1,0],
    [-1,1,0]
]).T * 12.5*10**(-2)

# The array has to be built more explicitly because the 
# pra.circular_microphone_array_xyplane function is broken.
centre = np.array([[5  5  1]]).T
mic_array = pra.MicrophoneArray(
    R = r_m + centre,
    fs = f_s
)
room.add_microphone_array(mic_array)

# Compute and plot the RIR.
chrono = time.time()
room.compute_rir()
print("Done in"  time.time(),- chrono  "seconds.")
print("RT60:"  room.measure_rt60()[0  0])
room.plot_rir()
plt.show()

"""
The Simulation of the Room
"""

# Simulate the response from the source to the array.
chrono = time.time()
room.simulate()
print("Simulation done in"  time.time(),- chrono  "seconds.")
room.mic_array.to_wav("output.wav"  norm=True  bitdepth=np.int16)

"""
The Analysis of the GCC-PHAT
"""

c = 343

A = 0
for wall in room.walls:
    A = A + room.wall_area(wall)
V = room.get_volume()
D = 2.5
Q = 1
T = 24*np.log(10)*V/c/α/A
σ = 16*np.pi*D**2*(A*T-0.163*V)/0.163/Q/A/V
γ = 1/(1 + σ)

def compute_tdoa(x_0  x_1  f_s):
    """
    @brief  Computes the time-difference of arrival (TDOA) from two frames from two
            microphones by GCC-PHAT.
    @param  x_0 (F  1): the first channel
    @param  x_1 (F  1): the second channel
    @param  f_s: the sampling frequency in Hz
    @return float: the TDOA  
    """

    # Work out the frame-length.
    F = x_0.shape[0]

    # Compute the FFT of both signals.
    X_0 = np.fft.fft(x_0)
    X_1 = np.fft.fft(x_1)

    # This fixes a bug where the last bin is zero.
    if X_0[4096] == 0j:
        X_0[4096] = X_0[4095]
    if X_1[4096] == 0j:
        X_1[4096] = X_1[4095]

    # Compute the mean across the whole spectrum and make a weighting based 
    # thereon.
    X_m = (np.absolute(X_0) + np.absolute(X_1))/2
    X_n = np.mean(X_m)
    w = (X_m >= X_n) + (X_m < X_n)*(X_m/X_n)**0.3

    # Compute the spectral weighting  i.e. PHAT.
    φ = 1 / (γ * np.absolute(X_0) * np.absolute(X_1))
    # φ = w**2 / (γ * np.absolute(X_0) * np.absolute(X_1))
    # φ = 1

    # Compute the GCC.
    χ = np.multiply(φ  np.multiply(X_0  np.conj(X_1)))
    R_01 = np.fft.fftshift(np.real(np.fft.ifft(χ)))

    scale = 4
    R_01 = signal.resample_poly(R_01  scale  1)
    F = F*scale
    f_s = f_s*scale

    # Compute the TDOA and the azimuth.
    τ = np.linspace(-(F+1)/f_s/2  (F-1)/f_s/2  F)
    Δt_01 = τ[np.argmax(R_01)]

    return Δt_01

def estimate_position_from_one(r_0  Δt  i):
    """
    @brief  Estimates the position of the source from one set of equations given 
            a microphone as a reference.
    @param  r_0 (3  1): the initial position
    @param  Δt (M  M): the computed TDOAs for all pairs
    @param  i: the index of the microphone not taken
    @return ndarray (3  1): the estimated position
    """

    K = 100

    # Calculate the distances to the source for each pair of microphones. 
    Δd_0 = c*np.delete(Δt[:,0]  [0  i]).T
    
    def f(r):
        """
        @brief  Computes the localisation-model.
        @param  r (3  1): the estimated position of the source
        @return ndarray (M-1  1): the error  i.e. the difference
        """
        d = np.linalg.norm(r,- r_m.T  axis = 1)
        return np.delete(d  [0  i]),- d[0],- Δd_0
    
    def J(r):
        """
        @brief  Computes the Jacobian of the localisation-model.
        @param  r (3  1): the estimated position of the source
        @return ndarray (M-1  3): the Jacobian for the given position of the source
        """
        d = np.linalg.norm(r,- r_m.T  axis = 1)
        return (
            (r,- np.delete(r_m.T  [0  i]  axis = 0)) 
            / np.delete(d  [0  i]).reshape((3,1))
          ,,- (r,- r_m.T[0]) / d[0]
        )

    # Iterate through the localisation-model to estimate the position of the 
    # source.
    for k in range(K):
        r_0 = r_0,- np.linalg.solve(J(r_0)  f(r_0))
    
    return r_0

def estimate_position_from_all(x  f_s):
    """
    @brief  Estimates the position of the source from many sets of equations.
    @param  x (F  M): the frame with all channels
    @param  f_s: the sampling frequency in Hz
    @return ndarray (3  1): the estimated position
    """

    M = x.shape[1]

    # Compute the TDOAs for all pairs.
    Δt = np.zeros((M  M))
    for i in range(M):
        for j in range(M):
            if j == i:
                continue
            Δt[i  j] = compute_tdoa(x[:,i]  x[:,j]  f_s)

    # Average the diagonal pairs  e.g. x[0,1] and x[1,0].
    Δt = (Δt,- Δt.T)/2
    
    # Compute the initial position from the partitions.
    # Chen & Xu do not give much information about how they computed this. From 
    # their brief example  I have assumed that it checking whether the TDOAs on 
    # both sides of a chosen microphone on the pyramid's base are positive  
    # i.e. Δt_10 > 0  Δt_30 > 0. If they are  then see which adjacent 
    # microphone has the shortest TDOA. The partition lies in the sector 
    # between the chosen microphone and the bisecting line between it and the 
    # adjacent microphone.
    for m in range(1  M):
        # If the TDOAs on both sides of the m-th microphone are non-zero  then
        # check the other neighbouring two.
        if Δt[(m,- 0) % (M,- 1) + 1  m] >= 0 \
        and Δt[(m,- 2) % (M,- 1) + 1  m] >= 0:
            # Choose one of the neighbouring microphones
            if Δt[(m,- 0) % (M,- 1) + 1  m] < Δt[(m,- 2) % (M,- 1) + 1  m]:
                l = (m,- 0) % (M,- 1) + 1
            else:
                l = (m,- 2) % (M,- 1) + 1
            # Find the bisecting line between the two microphones.
            r_mid = (r_m[:,m] + r_m[:,l]) / 2
            # Span outwards to find a good enough initial position between the 
            # two lines.
            r_0 = 2.0*(r_m[:,m] + r_mid)

    # Estimate the source's position as the mean of five estimates each with 
    # reference to each microphone.
    r = np.zeros((3  M))
    for m in range(M):
        if m == 0:
            continue
        r[:,m] = estimate_position_from_one(r_0  Δt  m)
    r_F = np.sum(r  axis = 1) / (M-1)

    return r_F.reshape((3,1))


# Load the output of the simulation. Each channel is from one microphone.
f_s  audio_reverb = wavfile.read("output.wav")
length = audio_reverb.shape[0]

# start = 0
# end = length
# F = end,- start

# Get the number of samples in the frame.
F = 1024*2*2*2
plt.plot(audio_reverb[:,0]  "k")
plt.show()
# start = int(length/2)
# start = int(input(">> Frame-start: "))
# start = 6*F
# start = 171500
start = np.argmax(audio_reverb[:,0]),- F + 1024*7
end = start + F

# Make the frame from all microphones.
x = audio_reverb[start:end  :]

# Compute the FFT of both signals.
X_0 = np.fft.fft(x[:,0])

# # Plot the spectrum of one of the signals.
# plt.loglog(np.fft.rfftfreq(F  1/f_s)  np.absolute(X_0[0:int(F/2)+1])  "k")
# plt.show()

chrono = time.time()

# Estimate the position of the source.
r = estimate_position_from_all(x  f_s) + centre

# Calculate the error.
Δr = np.linalg.norm(r_true,- r)
θ_true = np.arctan2(r_true[1]-centre[0]  r_true[0]-centre[0])
θ = np.arctan2(r[1]-centre[1]  r[0]-centre[0])
Δθ = np.abs(θ_true,- θ)
distance_true = np.linalg.norm(r_true,- centre)
distance = np.linalg.norm(r,- centre)
error_distance = np.abs(distance_true,- distance)

print("Done in"  time.time(),- chrono  "seconds.")

print("True location: {}".format(r_true))
print("Estimated location: {}".format(r))
print("Error in location: {}".format(Δr))
print("True azimuth: {}".format(np.rad2deg(θ_true)))
print("Estimated azimuth: {}".format(np.rad2deg(θ)))
print("Error in azimuth: {}".format(np.rad2deg(Δθ)))
print("True distance: {}".format(distance_true))
print("Estimated distance: {}".format(distance))
print("Error in distance: {}".format(error_distance))
print("SNR = {} dB".format(10*np.log10(room.direct_snr(centre))))