import matplotlib.pyplot as plt
import numpy as np

# categories = ("Gridpack", "NanoAOD", "Preselection", "Total Selected")
# states = dict(
#     ZLZL = (100, 100, 75.225, 27.939),
#     ZLZT = (100, 100, 74.527, 27.063),
#     ZTZT = (100, 100, 74.744, 25.048),
#     ZZ   = (100, 100, 75.361, 27.637)
# )

# categories = ("SR", "OS_HM", "OS_MM", "OS_LM", "SS_HM", "SS_MM", "SS_LM")
# states = dict(
#     ZLZL = (21.810, 4.098, 1.514, 0.291, 0.177, 0.043, 0.008),
#     ZLZT = (21.110, 4.078, 1.377, 0.270, 0.196, 0.025, 0.007),
#     ZTZT = (19.821, 4.132, 0.625, 0.177, 0.265, 0.022, 0.007),
#     ZZ   = (21.945, 4.303, 0.938, 0.181, 0.242, 0.022, 0.006)
# )

categories = ("4e", "4mu", "2e2mu", "4l")
states = dict(
    ZLZL = (3.473, 8.184, 10.776, 22.433),
    ZLZT = (14.609, 75.696, 40.162, 87.748),
    ZTZT = (46.217, 66.203, 120.247, 232.666),
    ZZ   = (75.696, 123.303, 173.891, 372.891)
)
errs = dict(
    ZLZL = (0.48, 0.61, 0.65, 0.8),
    ZLZT = (1.07, 1.31, 1.39, 1.71),
    ZTZT = (2.14, 2.42, 2.92, 3.48),
    ZZ   = (2.65, 3.04, 3.27, 4.08)
)

#errs = [[0.48, 0.61, 0.65, 0.8], [1.07, 1.31, 1.39, 1.71], [2.14, 2.42, 2.92, 3.48], [2.65, 3.04, 3.27, 4.08]]

x = np.arange(len(categories))
width = 0.2
multiplier = 0

fig, ax = plt.subplots(layout='constrained')
for state, percent in states.items():
    offset = width * multiplier
    #ax.bar(x + offset, percent, width, label=state, log=True)
    ax.bar(x + offset, percent, width, yerr=errs[state], label=state, log=True)
    multiplier +=1

ax.set_ylabel('Counts')
ax.set_title('Private MC Production Counts')
ax.set_xticks(x+width, categories)
ax.legend()
ax.set_ylim(0, 400)

plt.grid(axis = 'y')

fig.savefig("plot_log.png")