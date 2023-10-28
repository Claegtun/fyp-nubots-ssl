import matplotlib.pyplot as plt
import matplotlib

# Set the font.
plt.rcParams.update({'font.sans-serif':'FreeSerif'})

fig, ax_0 = plt.subplots()
ax_1 = ax_0.twinx()
ax_2 = ax_0.twinx()
ax_3 = ax_0.twinx()
ax_4 = ax_0.twinx()

ax_0.spines["top"].set_color("none")
ax_0.spines["bottom"].set_color("none")
ax_1.spines["top"].set_color("none")
ax_1.spines["bottom"].set_color("none")
ax_2.spines["top"].set_color("none")
ax_2.spines["bottom"].set_color("none")
ax_3.spines["top"].set_color("none")
ax_3.spines["bottom"].set_color("none")
ax_4.spines["top"].set_color("none")
ax_4.spines["bottom"].set_color("none")

ax_1.spines["right"].set_position(("axes",0.25))
ax_2.spines["right"].set_position(("axes",0.5))
ax_3.spines["right"].set_position(("axes",0.75))

ax_0.set_ylabel("Sound-level (dB_SPL)")
ax_1.set_ylabel("Microphone (dBV)")
ax_2.set_ylabel("Preamplifier (V)")
ax_3.set_ylabel("Difference-amplifier (V)")
ax_4.set_ylabel("ADC")

ax_0.set_xlim([0,1])
ax_0.set_xticks([])

ax_2.set_yscale("log")
ax_3.set_yscale("log")
ax_4.set_yscale("log",base=2)

ax_0.set_ylim([0,120])
ax_1.set_ylim([-132,-12])
ax_2.set_ylim([10**(-99/20),10**(21/20)])
ax_3.set_ylim([10**(-82/20),10**(38/20)])
ax_4.set_ylim([10**(-82/20)/5*2**16,10**(38/20)/5*2**16])

ax_0.plot([0,1],[94-65,94-65], "k--")
ax_0.fill_between([0,1],94-65,0,color="k",alpha=0.1)
ax_0.set_yticks(list(ax_0.get_yticks()) + [29])
ax_0.text(0.05,29,"Noise-floor",ha="left",va="bottom")

ax_0.annotate(
    "",
    xytext=(0.22,94), 
    xy=(0.22,29),
    arrowprops=dict(arrowstyle="->, head_length = 1, head_width = .2", lw=1)
)
ax_0.text(0.22,40,"SNR = 65 dB",ha="right",va="bottom",rotation="vertical")

ax_0.plot([0,1],[90,90], "k--")
ax_0.set_yticks(list(ax_0.get_yticks()) + [90])
ax_0.text(0.05,90,"Chosen Ceiling",ha="left",va="top")

ax_0.plot([0,1],[94,94], "k--")
ax_0.set_yticks(list(ax_0.get_yticks()) + [94])
ax_1.set_yticks(list(ax_1.get_yticks()) + [-38])
ax_0.text(0.30,94,"Sensitivity",ha="left",va="bottom")

ax_3.fill_between([0.75,1],2.5, 10**(38/20),color="k",alpha=0.1)
# ax_3.get_xaxis().get_major_formatter().labelOnlyBase = False
# ax_3.set_yticks(list(ax_3.get_yticks()) + [2.5], list(ax_3.get_yticks()) + [2.5])
# ax_3.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax_3.text(0.8,2.5,"Rail-to-rail limit",ha="left",va="top")

ax_4.set_yticks([2**4,2**8,2**12,2**15,2**16])

ax_0.plot([0,1],[8.3,8.3],"k--")
ax_0.set_yticks(list(ax_0.get_yticks()) + [8.3])
ax_0.text(0.8,8.3,"ADC-floor",ha="left",va="bottom")

ax_0.annotate(
    "",
    xytext=(0.95,96), 
    xy=(0.95,8.3),
    arrowprops=dict(arrowstyle="->, head_length = 1, head_width = .2", lw=1)
)
ax_0.text(0.95,40,"SINAD = 87.7 dB",ha="right",va="bottom",rotation="vertical")

ax_0.set_ylim([0,120])
ax_1.set_ylim([-132,-12])
ax_2.set_ylim([10**(-99/20),10**(21/20)])
ax_3.set_ylim([10**(-82/20),10**(38/20)])
ax_4.set_ylim([10**(-82/20)/5*2**16,10**(38/20)/5*2**16])

# plt.savefig("./report/gains.png", dpi = 1024)
plt.show()