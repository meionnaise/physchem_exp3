import numpy as np
import matplotlib.pyplot as plt

part1 = np.genfromtxt('part1.csv', delimiter=',')
part2_i = np.genfromtxt('part2_i.csv', delimiter=',')
part2_ii = np.genfromtxt('part2_ii.csv', delimiter=',')

####################################################
#PART 2
####################################################
print(part1[0, :])
print(part1[1:3, :])

#first bit


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