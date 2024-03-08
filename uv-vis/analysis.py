import numpy as np
import matplotlib.pyplot as plt

part1 = np.genfromtxt('part1.csv', delimiter=',')
part2_i = np.genfromtxt('part2_i.csv', delimiter=',')
part2_ii = np.genfromtxt('part2_ii.csv', delimiter=',')

####################################################
#PART 2
####################################################


#in session analysis
# 1. Correct each absorption spectrum using the recorded baseline spectrum (this only needs to be done
# if the absorbance of the blank is measured separately from that of the sample) and

#fuck off data labels
part2_i = part2_i[2:,:]
part2_ii = part2_ii[2:,:]

#correct Fe3+/Cu2+ ratio UV-VIS plots w recorded baseline
i = 1
while i < part2_i.shape[1]:
    part2_i[:,i] -= part2_i[:,23]
    i += 1

#Correct 5:5 Fe3+:Cu2+ dilution UV-VIS w blank (labeled sample zero due to linear dilution)
i = 1
while i < part2_ii.shape[1]:
    part2_ii[:,i] -= part2_ii[:,11]
    i += 1

# (a) plot the spectra for all mixture ratios at the same total Fe3+ + sal – concentration on the same
# graph. Identify the isosbestic point (note that it only occurs for mixture ratios in which Fe3+ is
# in excess).

plt.rcParams["font.family"] = "Comic Sans MS"
plt.rcParams['axes.facecolor'] = "#FFE1EF"
plt.figure(facecolor="#FFE1EF")
i = 0
while i < (part2_i.shape[1] - 2):
    plt.plot(part2_i[i], part2_i[i+1], color = "#fda0cc")
    i += 2
    
plt.axhline(y=0, color='k', linestyle='-', linewidth=0.5)
# plt.ylim([-1, max_rate + 3])
# plt.xlim([0, fattest_v])
plt.xlabel("Wavelength (nm)")
plt.ylabel("Absorbance")
# label = " ({}, {}) @ {} degrees.".format( max_rate_v, max_rate, max_rate_angle)
# plt.text(max_rate_v, max_rate, label , fontsize = "small")
plt.title("UV-Vis spectrum @ 300-800nm for varying Fe3+:Cu2+ ratio solutions")
plt.grid(color = 'w', linestyle = '-', linewidth = 0.5)
plt.show()

# (b) plot the spectra for the 1:1 mixtures at different dilutions on another graph


####################################################
#OLD
####################################################
# #plot results to visually affirm
# velocities = np.array(bruh[1])
# rates = np.array(bruh[0])
# max_rate = (bruh[2])
# max_rate_v = (bruh[3])
# max_rate_angle = (bruh[4])

# plt.rcParams["font.family"] = "Comic Sans MS"
# plt.rcParams['axes.facecolor'] = "#FFE1EF"
# plt.figure(facecolor="#FFE1EF")
# plt.plot(velocities, rates, color = "#fda0cc")
# plt.plot(max_rate_v, max_rate, marker ='.',color = "w")
# plt.axhline(y=0, color='k', linestyle='-', linewidth=0.5)
# plt.ylim([-1, max_rate + 3])
# plt.xlim([0, fattest_v])
# plt.xlabel("Airspeed (ms^{-1})")
# plt.ylabel("Rate of Climb (ms^{-1})")
# label = " ({}, {}) @ {} degrees.".format( max_rate_v, max_rate, max_rate_angle)
# plt.text(max_rate_v, max_rate, label , fontsize = "small")
# plt.title("Maximum Rate of Climb")
# plt.grid(color = 'w', linestyle = '-', linewidth = 0.5)
# plt.show()