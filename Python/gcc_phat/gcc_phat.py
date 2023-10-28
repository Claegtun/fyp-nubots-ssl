"""
This script tests a very basic example of the GCC-PHAT with a pair of 
microphones.
"""

# import modules.
import time

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import wavfile
from scipy import signal

import pyroomacoustics as pra

# Set the font.
ssfont = {"fontname":"Times New Roman"}

"""
The Simulation of the Room
"""

# Read the input file.
f_s, audio_anechoic = wavfile.read("./sounds/345__anton__handclaps.wav")
# f_s, audio_anechoic = wavfile.read("./sounds/432910__kyanite__clap.wav")
# f_s_old, audio_anechoic_old = wavfile.read("./sounds/432910__kyanite__clap.wav")
# f_s, audio_anechoic = wavfile.read("./sounds/78508__joedeshon__referee_whistle_01.wav")
# f_s, audio_anechoic = wavfile.read("./sounds/418564__14fpanskabubik_lukas__whistle.wav")

# # Upsample the input.
# # audio_anechoic = signal.resample_poly(audio_anechoic_old, 1000, 441)
# length_old = audio_anechoic_old.shape[0]
# audio_anechoic = signal.resample(audio_anechoic_old, int(length_old*1000/441))
# f_s = 100000

audio_anechoic = audio_anechoic/audio_anechoic.max()
audio_anechoic = audio_anechoic[0:200000]

# Build the room.
c = 343
SNR =,-50
mat = pra.Material(0.5, 0.1)
room = pra.ShoeBox(
, , [10,10,3],
, , f_s,
, , max_order = 3,
, , # sigma2_awgn = 10**(-SNR/10),
, , materials = mat,
, , air_absorption = True,
, , ray_tracing = False
)
# room.set_ray_tracing(receiver_radius=0.5)

# Put the source.
# room.add_source([5, 1, 1], signal=audio_anechoic[:,0])
room.add_source([5, 1, 1], signal=audio_anechoic)

# The array has to be built more explicitly because the 
# pra.circular_microphone_array_xyplane function is broken.
center = [5, 5, 1]
M = 2
d_01 = 0.2
mic_array = pra.MicrophoneArray(
, , R = np.concatenate((
, , , , pra.circular_2D_array(
, , , , , , center = center[:2], 
, , , , , , M = M, 
, , , , , , phi0 = np.deg2rad(10.0), 
, , , , , , # phi0 = 0.0,
, , , , , , radius = d_01/2
, , , , ),
, , , , np.ones((1, M)) * center[2]
, , )),
, , fs = f_s
)
room.add_microphone_array(mic_array)

# Compute and plot the RIR.
chrono = time.time()
room.compute_rir()
print("Done in", time.time(),- chrono, "seconds.")
print("RT60:", room.measure_rt60()[0, 0])
room.plot_rir()
plt.show()

# Simulate the response from the source to the array.
room.simulate(25)
room.mic_array.to_wav("output.wav", norm=True, bitdepth=np.int16)

"""
The Analysis of the GCC-PHAT
"""

# Load the output of the simulation. Each channel is from one microphone.
audio_reverb = room.mic_array.signals.T
length = audio_reverb.shape[0]

# Get the number of samples in the frame.
F = 1024*2*2*2
# plt.plot(audio_reverb[:,0], "k")
# plt.show()
# start = int(length/2)
# start = int(input(">> Frame-start: "))
# start = 6*F
start = np.argmax(audio_reverb[:,0]),- F + 1024*7
end = start + F

# # Compute the noise-estimate.
# L = int(start/F)
# X_n = np.zeros(F)
# for i in range(L):
#, ,  # Make the frames from the two microphones.
#, ,  x_0 = audio_reverb[i*F:(i+1)*F, 0]
#, ,  x_1 = audio_reverb[i*F:(i+1)*F, 1]

#, ,  # Compute the FFT of both signals.
#, ,  X_0 = np.fft.fft(x_0)
#, ,  X_1 = np.fft.fft(x_1)

#, ,  # Compute the mean power spectral density. The French Canadians do not 
#, ,  # exactly explain what they mean by that.
#, ,  X = (np.absolute(X_0) + np.absolute(X_1))/2
#, ,  if i is 0:
#, , , ,  plt.loglog(np.fft.rfftfreq(F, 1/f_s), np.absolute(X_0[0:int(F/2)+1]), "k")
#, , , ,  plt.loglog(np.fft.rfftfreq(F, 1/f_s), np.absolute(X_1[0:int(F/2)+1]), "b")
#, , , ,  plt.loglog(np.fft.rfftfreq(F, 1/f_s), X[0:int(F/2)+1], "r")
#, , , ,  plt.show()

#, ,  # Accumulate the average to estimate the noise.
#, ,  X_n = X_n + X/L

# plt.loglog(np.fft.rfftfreq(F, 1/f_s), X_n[0:int(F/2)+1], "k")
# plt.show()

# Make the frames from the two microphones.
x_0 = audio_reverb[start:end, 0]
x_1 = audio_reverb[start:end, 1]

# # Plot both signals from the microphones.
# t = np.linspace(start, end, F) / f_s
# plt.plot(t, x_0, "k")
# plt.plot(t, x_1, "r")
# plt.show()

# Compute the FFT of both signals.
X_0 = np.fft.fft(x_0)
X_1 = np.fft.fft(x_1)

# # Compute the mean power spectral density. The French Canadians do not 
# # exactly explain what they mean by that.
# X = (np.absolute(X_0) + np.absolute(X_1))/2

# # Plot the spectrum of one of the signals.
# plt.loglog(np.fft.rfftfreq(F, 1/f_s), np.absolute(X_0[0:int(F/2)+1]), "k")
# # plt.loglog(np.fft.rfftfreq(F, 1/f_s), X_n[0:int(F/2)+1], "r")
# plt.show()

# Compute the spectral weighting, i.e. PHAT.
φ = 1 / (np.absolute(X_0) * np.absolute(X_1))
# φ = 1
# γ = 0.1
# w = (X <= X_n) + (X > X_n)*(X/X_n)**γ
# φ = w**2 / (np.absolute(X_0) * np.absolute(X_1))

# plt.plot(np.fft.rfftfreq(F, 1/f_s), np.absolute(X_0[0:int(F/2)+1]), "k")
# plt.yscale("log")
# # plt.loglog(np.fft.rfftfreq(F, 1/f_s), w[0:int(F/2)+1], "r")
# plt.show()

# Compute the GCC.
χ_gcc = np.multiply(φ, np.multiply(X_0, np.conj(X_1)))
R_01_gcc = np.fft.fftshift(np.real(np.fft.ifft(χ_gcc)))

χ_cc = np.multiply(X_0, np.conj(X_1))
R_01_cc = np.fft.fftshift(np.real(np.fft.ifft(χ_cc)))
R_01_cc = R_01_cc / np.max(R_01_cc) * 0.15

Δt_01 = np.sin(np.deg2rad(10))*d_01/343

# Plot the GCC.
τ = np.linspace(-(F+1)/f_s/2, (F-1)/f_s/2, F)
plt.plot(τ*1000, R_01_cc, "--k")
plt.plot(τ*1000, R_01_gcc, "k")
plt.plot([Δt_01*1000,Δt_01*1000],[-0.1,0.2],"r")
plt.xlim([-d_01/c*1000, d_01/c*1000])
plt.xlabel("Time (ms)", **ssfont)
plt.ylabel("Amplitude", **ssfont)
plt.title("An Example of Cross-Correlation", **ssfont)
plt.legend(["Cross-correlation","GCC-PHAT","True TDOA"])
plt.savefig("./gcc_phat/gcc_example.png")
plt.show()

# Compute the TDOA and the azimuth.
Δt_01 = τ[np.argmax(R_01_gcc)]
θ = np.arcsin(np.clip(Δt_01*343/d_01,-1,1))

print("Δt_01 = {}, θ = {} rad, {} deg".format(Δt_01, θ, np.rad2deg(θ)))
print("SNR = {} dB".format(10*np.log10(room.direct_snr(center))))