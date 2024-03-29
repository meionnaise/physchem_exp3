import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from matplotlib import font_manager

def mM_2_M(conc):
    results = conc / 1000
    return results

#assume length is the same for each array
def common_point(wavelengths, array1, array2):
    i = 0
    k = 0
    dy = 10 #arbitrarily large for initialising variable
    while i < len(wavelengths):
        while k < len(array1):
            difference = abs(array1[k] - array2[k])
            if difference < dy:
                dy = difference
                wl = wavelengths[i]
                absorbance = np.average((array1[k], array2[k]))
            k += 1
            i += 1

    return wl, absorbance

def avg_value(array):
    total = 0
    for i in array:
        total += i

    result = total / len(array)
    return result

part2_i = np.genfromtxt('uv_vis/part2_i.csv', delimiter=',')
part2_ii = np.genfromtxt('uv_vis/part2_ii.csv', delimiter=',')

font_path = 'Comic Sans MS.ttf'  
font_manager.fontManager.addfont(font_path)
prop = font_manager.FontProperties(fname=font_path)

####################################################
#PART 2
####################################################


#in session analysis
# 1. Correct each absorption spectrum using the recorded baseline spectrum (this only needs to be done
# if the absorbance of the blank is measured separately from that of the sample) and

#part2 i labels
part2_i_labels = ["0:10","","1:9","","2:8","","3:7","","4:6","","5:5","","6:4","","7:3","","8:2","","9:1","","10:0","","blank", ""]

#part2 ii labels
part2_ii_labels = ["1.3 mM","","1.04 mM","","0.78 mM","","0.52 mM","","0.26 mM","","0.00 mM", ""]
#fuck off data labels
part2_i = part2_i[2:,:]
part2_ii = part2_ii[2:,:]

#correct Fe3+/Cu2+ ratio UV-VIS plots w recorded baseline
i = 1
while i < part2_i.shape[1]:
    part2_i[:,i] -= part2_i[:,23]
    i += 2
part2_i = part2_i[:, :22] #get rid of blank

#Correct 5:5 Fe3+:Cu2+ dilution UV-VIS w blank (labeled sample zero due to linear dilution)
i = 1
while i < part2_ii.shape[1]:
    part2_ii[:,i] -= part2_ii[:,11]
    i += 2
#part2_ii = part2_ii[:, :10] #get rid of blank

# (a) plot the spectra for all mixture ratios at the same total Fe3+ + sal – concentration on the same
# graph. Identify the isosbestic point (note that it only occurs for mixture ratios in which Fe3+ is
# in excess).
    
#finding isosbestic point - from 6:4 - 10:1
i = 11
tot_common_wl = (0)
tot_common_abs = (0)
indexfrom = 300
indexto = 500
while i < 19:
    point = common_point(part2_i[indexfrom:indexto,0], part2_i[indexfrom:indexto, i], part2_i[indexfrom:indexto, i + 2])
    point =  np.array([point])
    if i == 11:
        common_points = point
        i += 2
        continue
    common_points = np.append(common_points, point, axis = 0)
    i += 2
isosbestic_point = (np.average(common_points[:,0]), np.average(common_points[:,1]))

print(f"\nQ1. a)\nIsosbestic point: {isosbestic_point}\n")

#plotting 2.1.a) with isosbestic point
plt.rcParams["font.family"] = "Comic Sans MS"
plt.rcParams['axes.facecolor'] = "#FFE1EF"
plt.figure(facecolor="#FFE1EF")
i = 0
while i < (part2_i.shape[1] - 2):
    plt.plot(part2_i[:, i], part2_i[:, i+1], label = part2_i_labels[i])
    i += 2
plt.xlabel("Wavelength (nm)")
plt.ylabel("Absorbance")
label = "Isosbestic point\n({}, {})".format( round((isosbestic_point[0]), 3), round((isosbestic_point[1]), 3))
plt.text(isosbestic_point[0], isosbestic_point[1], label , fontsize = "small")
plt.plot(isosbestic_point[0], isosbestic_point[1], 'ro')
plt.title("UV-Vis spectrum @ 300-800nm for varying [Fe3+]:[sal-] ratio solutions")
plt.grid(color = 'w', linestyle = '-', linewidth = 0.5)
plt.legend()
plt.savefig('part2_q1a.png')


#(b) plot the spectra for the 1:1 mixtures at different dilutions on another graph
#plotting 2.1.b)
plt.rcParams["font.family"] = "Comic Sans MS"
plt.rcParams['axes.facecolor'] = "#FFE1EF"
plt.figure(facecolor="#FFE1EF")
i = 0
while i < part2_ii.shape[1] - 2:
    plt.plot(part2_ii[:, i], part2_ii[:, i+1], label = part2_ii_labels[i])
    i += 2
plt.xlabel("Wavelength (nm)")
plt.ylabel("Absorbance")
plt.title("UV-Vis spectrum @ 300-800nm for 1:1 mixtures at different dilutions")
plt.grid(color = 'w', linestyle = '-', linewidth = 0.5)
plt.legend(title = "[Fe3+] + [sal-]")
plt.savefig('part2_q1_b.png')


###### 2.2. job plot & get empirical formula
#get index for abs @ 530nm
i = 0
absorbances_530nm = ()
while i < len(part2_i[:,0]):
    if part2_i[i,0] == 530:
        index = i
        break
    i += 1

#get absorbances @ 530nm
i = 1
while i < part2_i.shape[1]:
    point_wanted = part2_i[index, i]
    point_wanted = np.array([point_wanted])
    if i == 1:
        job_abs = point_wanted
        i += 2
        continue
    job_abs = np.append(job_abs, point_wanted, axis = 0)
    i += 2

# get corresp [Fe 3+] concs
i = 0
x_data = []
sal_conc = []
data_labels = []
while i < len(job_abs):
    label = part2_i_labels[i * 2]
    x = mM_2_M(1.3) * i / 10
    sal = mM_2_M(1.3) - x
    x = np.array([x])
    sal = np.array([sal])
    label = np.array([label])
    x_data = np.append(x_data, x, axis = 0)
    sal_conc = np.append(sal_conc, sal, axis = 0)
    data_labels = np.append(data_labels, label, axis = 0)
    i += 1

#find empirical formula
max_abs_a = np.amax(job_abs)
i = 0
while i < len(job_abs):
    if job_abs[i] == max_abs_a:
        empirical = data_labels[i]
        emp_x = x_data[i]
        break
    i += 1

print(f"\nQ2\nEmpirical formula in ratio {empirical}\n")

#plot
plt.rcParams["font.family"] = "Comic Sans MS"
plt.rcParams['axes.facecolor'] = "#FFE1EF"
plt.figure(facecolor="#FFE1EF")
plt.plot(x_data, job_abs, marker = 'o',label = data_labels, color = "#fda0cc")
plt.text(emp_x, max_abs_a, empirical , fontsize = "small")
plt.plot(emp_x, max_abs_a, 'o', color = "#fda0cc")
plt.xlabel("[Fe3+] (mM)")
plt.ylabel("Absorbance")
plt.title("Job Plot @ 530nm for different ratios of Fe3+ and sal-")
plt.grid(color = 'w', linestyle = '-', linewidth = 0.5)
plt.legend(title = "[Fe3+]")
plt.savefig('part2_q2_job_plot.png')

# 2.3. plot 1:1 abs @ 530nm and find extinction coeff using beer's law

#get index for abs @ 530nm
i = 0
absorbances_530nm = ()
while i < len(part2_ii[:,0]):
    if part2_ii[i,0] == 530:
        index = i
        break
    i += 1

#get absorbances @ 530nm
i = 0
absorbance = []
while i < part2_ii.shape[1]:
    point_wanted = part2_ii[index, i+1]
    point_wanted = np.array([point_wanted])
    absorbance = np.append(absorbance, point_wanted, axis = 0)
    i += 2
np.linalg.norm(absorbance)
# get corresp [Fe 3+] concs
i = 0
x_data_ii = []
data_labels = []
while i < len(absorbance):
    label = part2_ii_labels[i * 2]
    if i == 5:
        x = 0
    else:
        x = 1.3 * (2 * (i+1)) / 10
    x = np.array([x])
    label = np.array([label])
    x_data_ii = np.append(x_data_ii, x, axis = 0)
    data_labels = np.append(data_labels, label, axis = 0)
    i += 1

x_data_ii /= 2 #conc of complex is half of what was added
x_data_ii /= 1000 #get in mM

#plot
plt.rcParams["font.family"] = "Comic Sans MS"
plt.rcParams['axes.facecolor'] = "#FFE1EF"
plt.figure(facecolor="#FFE1EF")
plt.scatter(x_data_ii, absorbance, marker = 'o', color = "#fda0cc" )
m, b, *_ = stats.linregress(x_data_ii, absorbance)
q3_ext = m 
print(f"\nQ3\nExtinction coefficient: {q3_ext} M^(-1)cm^(-1)\n")
plt.axline(xy1=(0, b), slope=m, label=f'$y = {m:.1f}x {b:+.1f}$', color = "#fda0cc", linestyle = "--")
plt.xlabel("[Fe3+(sal-)] (M)")
plt.legend()
plt.ylabel("Absorbance")
plt.title("Absorbance @ 530nm for different concentrations of 1:1 [Fe3]+ :[ sal-]")
plt.grid(color = 'w', linestyle = '-', linewidth = 0.5)
plt.savefig('part2_q3.png')


#q4 find stability constant, gibbs free energy for formation, and extinction coefficient using methods A and B
#a (varying ratio) THIS IS FUCKED
i = 0
fe_conc_data_4a = x_data[1:-1] #remove zero conc, and get original concs, not for complex
sal_conc_data_4a = sal_conc [1:-1]
xaxis_4a = job_abs[1:-1] #remove abs data for zero conc
yaxis_4a = []
for concentration in fe_conc_data_4a:
    y4a = np.array([fe_conc_data_4a[i] * (sal_conc_data_4a[i]) / xaxis_4a[i]])
    yaxis_4a = np.append(yaxis_4a, y4a, axis = 0)
    i += 1


#plot
plt.rcParams["font.family"] = "Comic Sans MS"
plt.rcParams['axes.facecolor'] = "#FFE1EF"
plt.figure(facecolor="#FFE1EF")
plt.scatter(xaxis_4a, yaxis_4a, marker = 'o', color = "#fda0cc" )
xerra = np.std(xaxis_4a)
yerra = np.std(yaxis_4a)
plt.errorbar(xaxis_4a, yaxis_4a, xerr = xerra , yerr = yerra, color = "#fda0cc", ls = "none") #this line fucked
m, b, *_ = stats.linregress(xaxis_4a, yaxis_4a)
plt.axline(xy1=(0, b), slope=m, label=f'$y = {m:.8f}x {b:+.8f}$', color = "#fda0cc", linestyle = "--")
plt.ylabel("[Fe3+]0 x [sal-]0 / A (M^2)")
plt.legend()
plt.xlabel("Absorbance")
plt.title("Method A")
plt.grid(color = 'w', linestyle = '-', linewidth = 0.5)
plt.savefig('part2_q4a.png')
#get stuff
R = 8.3145 #J mol^-1 K^-1
T = 298.15 # K (assume 25 deg)
ext_coeff_a = 1 / np.sqrt(- m) 
#use point (0, y-intercept) to get values along line of best fit
k_a = -1 / ((b * (ext_coeff_a)) - (mM_2_M(1.3)))
Ko_a = 1 * k_a #see eqn 20, since previous stuff is in mM, multiply by 1000 mM to get unitless
dGo_a = - R * T * np.log(Ko_a)
print(f"\nQ4\na)\nMethod A:\nExtinction coefficient: {ext_coeff_a} M^(-1)cm^(-1)\nStability constant: {Ko_a} M^(-1)\ndG^o = {dGo_a} Jmol^(-1)\n")


#b
i = 0
yaxis_4b = []
xaxis_4b = []
while i < len(x_data_ii):
    x4b = np.sqrt(absorbance[i])
    if x4b == 0:
        i += 1
        continue
    x4b = np.array([x4b])
    y4b = x_data_ii[i] / (np.sqrt(absorbance[i]))
    y4b = np.array([y4b])
    yaxis_4b = np.append(yaxis_4b, y4b, axis = 0)
    xaxis_4b = np.append(xaxis_4b, x4b, axis = 0)

    i += 1

#plot
plt.rcParams["font.family"] = "Comic Sans MS"
plt.rcParams['axes.facecolor'] = "#FFE1EF"
plt.figure(facecolor="#FFE1EF")
plt.scatter(xaxis_4b, yaxis_4b, marker = 'o', color = "#fda0cc" )
xerrb = np.std(xaxis_4b)
yerrb = np.std(yaxis_4b)
plt.errorbar(xaxis_4b, yaxis_4b, xerr = xerrb , yerr = yerrb, color = "#fda0cc", ls = "none")
m, b, *_ = stats.linregress(xaxis_4b, yaxis_4b)
plt.axline(xy1=(0, b), slope=m, label=f'$y = {m:.8f}x {b:+.8f}$', color = "#fda0cc", linestyle = "--")
plt.plot(xaxis_4b, yaxis_4b, marker = 'o', color = "#fda0cc" )
plt.ylabel("[Fe3+] / A  ^ (1 / 2) (M^2)")
plt.legend()
plt.xlabel("(Absorbance)^(1/2)")
plt.title("Method B")
plt.grid(color = 'w', linestyle = '-', linewidth = 0.5)
plt.savefig('part2_q4b.png')
#get stuff
ext_coeff_b = 1 / m #extinction coefficient
k_b = (1 / (b ** 2)) / ext_coeff_b
Ko_b = 1 * k_b #see eqn 20
dGo_b = - R * T * np.log(Ko_b)
print(f"\nb)\nMethod B:\nExtinction coefficient: {ext_coeff_b} M^(-1)cm^(-1)\nStability constant: {Ko_b} M^(-1)\ndG^o = {dGo_b} Jmol^(-1)")
    
######
#write data to files

f = open("part2_q1a.txt", "w")
f.write(str(isosbestic_point))
f.close()

f = open("part2_q2.txt", "w")
f.write(str(empirical))
f.close()

f = open("part2_q3.txt", "w")
f.write(str(q3_ext))
f.close()

f = open("part2_q4_a_ext.txt", "w")
f.write(str(ext_coeff_a))
f.close()

f = open("part2_q4_a_stab.txt", "w")
f.write(str(Ko_a))
f.close()

f = open("part2_q4_a_gibbs.txt", "w")
f.write(str(dGo_a))
f.close()

f = open("part2_q4_b_ext.txt", "w")
f.write(str(ext_coeff_b))
f.close()

f = open("part2_q4_b_stab.txt", "w")
f.write(str(Ko_b))
f.close()

f = open("part2_q4_b_gibbs.txt", "w")
f.write(str(dGo_b))
f.close()