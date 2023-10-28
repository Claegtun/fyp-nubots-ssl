import matplotlib.pyplot as plt
import numpy as np
from scipy.io import wavfile
from scipy import signal

# Read the input file.
# f_s, audio_anechoic = wavfile.read("./sounds/432910__kyanite__clap.wav")
# f_s, audio_anechoic = wavfile.read("./sounds/345__anton__handclaps.wav")
# f_s, audio_anechoic = wavfile.read("./sounds/78508__joedeshon__referee_whistle_01.wav")
f_s, audio_anechoic = wavfile.read("./sounds/418564__14fpanskabubik_lukas__whistle.wav")
# f_s, audio_anechoic = wavfile.read("./pyroomacoustics/examples/input_samples/cmu_arctic_us_aew_a0001.wav")
# audio_anechoic = signal.resample(audio_anechoic, audio_anechoic.shape[0]*2)
# f_s = f_s*2

audio_anechoic = audio_anechoic/audio_anechoic.max()
# audio_anechoic = audio_anechoic[160000*2:180000*2],  
# audio_anechoic = audio_anechoic[160000:180000]    

print("f_s = {} Hz".format(f_s))
# print("{} channels".format(audio_anechoic.shape[1]))

plt.plot(audio_anechoic)
plt.show()

# x = audio_anechoic[:,0]
x = audio_anechoic[:]
X = np.fft.fft(x)
F = X.shape[0]
f = np.fft.fftfreq(F, 1/f_s)
plt.plot(f, np.abs(X), "k")
# plt.yscale("log")
plt.show()