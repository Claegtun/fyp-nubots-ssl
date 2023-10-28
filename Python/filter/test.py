import matplotlib.pyplot as plt
import numpy as np
from scipy.io import wavfile
from scipy import signal

plt.rcParams.update({'font.sans-serif':'FreeSerif'})

f_s = 192*10**3
n_taps = 64
f_c = 10/192*2

h = signal.firwin(n_taps, f_c)
H = np.fft.fft(h)

string = ""
for i in range(64):
    string = string + "{:.6e}".format(h[i]) + ","
    if i%8 == 7:
        string = string + "\n"
print(string)

fig, ax = plt.subplots()
ax.plot(h, "k.-")
ax.set_xlabel("$t$")
ax.set_ylabel("$h(t)$")
ax.set_title("The Filter's Impulse-Response")
plt.savefig("./filter/impulse",dpi=1024)
plt.show()

fig, ax = plt.subplots()
ax.set_ylim([-80,5])
ax.semilogx(np.fft.fftfreq(n_taps)*f_s/10**3, 20*np.log10(np.abs(H)), "kx-")
ax.plot([8.4,8.4], [-80,5], "k--")
ax.plot([12.4,12.4], [-80,5], "k", linestyle="dotted")
ax.fill_betweenx([-80,5],0,12,color="k",alpha=0.1)
ax.legend(["Response","approx. -3 dB Cut-off", "approx. -20 db Cut-off", "Band of Interest"])
ax.set_xlabel("Frequency (kHz)")
ax.set_ylabel("Gain $|H|$ (dB)")
ax.set_title("The Digital Filter's Frequency-Response")
plt.savefig("./filter/frequency",dpi=1024)
plt.show()