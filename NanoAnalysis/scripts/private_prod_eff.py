import matplotlib.pyplot as plt
import numpy as np

# categories = ("Gridpack", "NanoAOD", "Preselection", "Total Selected")
# states = dict(
#     ZLZL = (100, 100, 75.225, 27.939),
#     ZLZT = (100, 100, 74.527, 27.063),
#     ZTZT = (100, 100, 74.744, 25.048),
#     ZZ   = (100, 100, 75.361, 27.637)
# )

categories = ("SR", "OS_HM", "OS_MM", "OS_LM", "SS_HM", "SS_MM", "SS_LM")
states = dict(
    ZLZL = (21.810, 4.098, 1.514, 0.291, 0.177, 0.043, 0.008),
    ZLZT = (21.110, 4.078, 1.377, 0.270, 0.196, 0.025, 0.007),
    ZTZT = (19.821, 4.132, 0.625, 0.177, 0.265, 0.022, 0.007),
    ZZ   = (21.945, 4.303, 0.938, 0.181, 0.242, 0.022, 0.006)
)

x = np.arange(len(categories))
width = 0.15
multiplier = 0

fig, ax = plt.subplots(layout='constrained')
for state, percent in states.items():
    offset = width * multiplier
    ax.bar(x + offset, percent, width, label=state, log=True)
    multiplier +=1

ax.set_ylabel('Percent')
ax.set_title('Private MC Production Efficiency')
ax.set_xticks(x+width, categories)
ax.legend()
ax.set_ylim(0, 25)

plt.grid(axis = 'y')

fig.savefig("PrivateMC_Efficiency_SelectedEvents_log.png")