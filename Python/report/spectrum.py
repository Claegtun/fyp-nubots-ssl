import matplotlib.pyplot as plt
import matplotlib

import numpy as np

# Set the font.
plt.rcParams.update({'font.sans-serif':'FreeSerif'})

f = np.linspace(0, 800, 1000)

f_adc = 15
G_adc = 1/np.sqrt(1+f**2/f_adc**2)

f_s = 192

f_nq = f_s/2

fig, ax = plt.subplots()

ax.set_xlim([0.1,800])
ax.set_ylim([-40,5])

ax.semilogx(f, 20*np.log10(G_adc), "k")

ax.plot([f_adc,f_adc],[-40,5],"k", linestyle="dotted")
ax.text(f_adc,-20,"Analogue -3dB cut-off",ha="right",va="bottom",rotation="vertical")

ax.fill_betweenx([-40,5],0.1,12,color="k",alpha=0.1)
ax.text(1,-30,"Band of Interest",ha="center",va="top")

ax.plot([f_s,f_s],[-40,5],"k", linestyle="dotted")
ax.text(f_s,-20,"Sampling frequency",ha="right",va="bottom",rotation="vertical")
ax.plot([f_nq,f_nq],[-40,5],"k", linestyle="dotted")
ax.text(f_nq,-10,"Nyquist frequency",ha="right",va="bottom",rotation="vertical")

f_fold_0 = np.linspace(10,f_nq,100)
G_fold_0 = 1/np.sqrt(1+(f_s-f_fold_0)**2/f_adc**2)
ax.semilogx(f_fold_0, 20*np.log10(G_fold_0), "k--")
f_fold_1 = np.linspace(1,10,100)
G_fold_1 = 1/np.sqrt(1+(f_s-f_fold_1)**2/f_adc**2)
ax.semilogx(f_fold_1, 20*np.log10(G_fold_1), "k--")
ax.text(1,-23,"Aliasing",ha="left",va="top")

ax.set_xlabel("Frequency (kHz)")
ax.set_ylabel("Gain (dB)")
ax.set_title("The Anti-Alias Filter and Sampling")

plt.grid(True)
plt.savefig("./report/spectrum.png", dpi = 1024)
plt.show()

