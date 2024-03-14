import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from matplotlib import font_manager

def normalise_col(x):
    i = 0
    while i < x.shape[1]:
        max_val = np.amax(x[:,i+1])
        x[:,i+1] /= max_val
        i += 2
    return x

def stern_volmer_y(x, noq):

    #find wl for max intensity
    wl_at_peak, peak = wl_4_peak(noq)
    i0 = peak
    
    result = []
   
    for sample in x:
        k = 0
        for wavelength in sample[:,0]:
            if wavelength != wl_at_peak:
                k += 1
                continue
            i = sample[k, 1]
            y = np.array([i0 / i])
            print(i, i0, y)
            result = np.append(result, y, axis = 0)
            break
    return result, wl_at_peak

def rate_coeff(slope):
    kd = 1.7 * (10 ** 6)
    kq = slope * kd
    return kq


def wl_4_peak(sample):
    pk = np.amax(sample[:,1])
    i = 0
    for wl in sample[:,0]:
        if sample[i,1] == pk:
            wavelength = sample[i,0]
            break
        i += 1
    return wavelength, pk

def kET(k22, x_red_potential, z):
    r_red_potential = -0.84
    r_ox = -r_red_potential
    R = 8.3145 #J mol^-1 K^-1
    T = 298.15 # K (assume 25 deg)
    F = 96485.3399 #Cmol^(-1) #faraday constant
    Z12 = 10 ** 11 #factor related to collision frequency
    k11 = 10 ** 8 #(Rubpy kii)

    #get cell potential
    cell_potential = x_red_potential + r_ox 
    # cell_potential *= 1.602 * (10 ** (-19)) #convert units
    # use eqn 13 to get from cell potential - K12 (equilibrium constant for cross rxn)
    K12 = np.exp(cell_potential * ( (z * F) / (R * T)))
    #use K12 --> f12 (eqn 12)
    f12 = np.exp(((np.log(K12)) ** 2) / (4 * np.log((k11 * k22) / (Z12 ** 2))))
    # use f12 , k11, k22, k12 to get kET
    kET = np.sqrt(k11 * k22 * K12 * f12)

    return kET


font_path = 'Comic Sans MS.ttf'  # Your font path goes here
font_manager.fontManager.addfont(font_path)
prop = font_manager.FontProperties(fname=font_path)




uv_vis = np.genfromtxt('uv_vis/part1.csv', delimiter=',')
rubpy_emission = np.genfromtxt('fluorescence/emission_rubpy.csv', delimiter=',')
rubpy_emission_norm = np.genfromtxt('fluorescence/emission_rubpy.csv', delimiter=',')
rubpy_excitation = np.genfromtxt('fluorescence/excitation Rubpy.csv', delimiter=',')
sample1 = np.genfromtxt('fluorescence/blank.csv', delimiter=',')
sample2 = np.genfromtxt('fluorescence/sample2_emission_q.csv', delimiter=',')
sample3 = np.genfromtxt('fluorescence/sample3_emission_q.csv', delimiter=',')
sample4 = np.genfromtxt('fluorescence/sample4_emission_q.csv', delimiter=',')
sample5 = np.genfromtxt('fluorescence/sample5_emission_q.csv', delimiter=',')
sample6 = np.genfromtxt('fluorescence/sample6_emission_q.csv', delimiter=',')
sample7 = np.genfromtxt('fluorescence/sample7_emission_q.csv', delimiter=',')
sample8 = np.genfromtxt('fluorescence/sample8_emission_q.csv', delimiter=',')
sample9 = np.genfromtxt('fluorescence/sample9_emission_q.csv', delimiter=',')
sample10 = np.genfromtxt('fluorescence/sample10_emission_q.csv', delimiter=',')
sample11 = np.genfromtxt('fluorescence/sample11_emission_q.csv', delimiter=',')
sample12 = np.genfromtxt('fluorescence/sample12_emission_q.csv', delimiter=',')
sample13 = np.genfromtxt('fluorescence/sample13_emission_q.csv', delimiter=',')

#fuck off data labels
rubpy_emission = rubpy_emission[2:,:]
rubpy_emission_norm = rubpy_emission_norm[2:,:]
rubpy_excitation = rubpy_excitation[2:,:]
uv_vis = uv_vis[2:,:]

sample1 = sample1[1:,:]
sample2 = sample2[1:,:]
sample3 = sample3[1:,:]
sample4 = sample4[1:,:]
sample5 = sample5[1:,:]
sample6 = sample6[1:,:]
sample7 = sample7[1:,:]
sample8 = sample8[1:,:]
sample9 = sample9[1:,:]
sample10 = sample10[1:,:]
sample11 = sample11[1:,:]
sample12 = sample12[1:,:]
sample13 = sample13[1:,:]

samples = np.array([sample1, sample2, sample3, sample4, sample5, sample6, sample7, sample8, sample9, sample10, sample11, sample12, sample13])

#correct w blanks
i = 0
while i < uv_vis.shape[1]:
    uv_vis[:,i] -= uv_vis[:,7]
    i += 2


rubpy_noq = uv_vis[:,0:2]
sample6_noq = uv_vis[:,2:4]
sample10_noq = uv_vis[:,4:6]
blank = uv_vis[:,6:]
rubpy_noq[:,1] -= blank[:,1]
sample6_noq[:,1] -= blank[:,1]
sample10_noq[:,1] -= blank[:,1]

# #normalise
# sample1 = normalise_col(sample1)
# sample2 = normalise_col(sample2)
# sample3 = normalise_col(sample3)
# sample4 = normalise_col(sample4)
# sample5 = normalise_col(sample5)
# sample6 = normalise_col(sample6)
# sample7 = normalise_col(sample7)
# sample8 = normalise_col(sample8)
# sample9 = normalise_col(sample9)
# sample10 = normalise_col(sample10)
# sample11 = normalise_col(sample11)
# sample12 = normalise_col(sample12)
# sample13 = normalise_col(sample13)

sample1, sample2, sample3, sample4, sample5, sample6, sample7, sample8, sample9, sample10, sample11, sample12, sample13 = samples

#no quencher

#normalise
rubpy_noq = normalise_col(rubpy_noq)
print(rubpy_emission.shape)
print(sample11.shape)
rubpy_emission_norm = normalise_col(rubpy_emission_norm)
rubpy_excitation = normalise_col(rubpy_excitation)


noq_labels = ["Absorption", ]
plt.rcParams["font.family"] = "Comic Sans MS"
plt.rcParams['axes.facecolor'] = "#FFE1EF"
plt.figure(facecolor="#FFE1EF")
plt.plot(rubpy_noq[:, 0], rubpy_noq[:,1], label = "Absorption")
plt.plot(rubpy_emission_norm[:, 0], rubpy_emission_norm[:,1], label = "Emission")
plt.plot(rubpy_excitation[:, 0], rubpy_excitation[:,1], label = "Excitation")
plt.xlabel("Wavelength (nm)")
plt.ylabel("Emission intensity")
plt.title("Spectra w/out quencher")
plt.grid(color = 'w', linestyle = '-', linewidth = 0.5)
plt.legend()
plt.savefig('part1_noq.png')

#plotting quenched ones

#fe
fe_labels = [0, 0.2, 0.4, 0.8, 1.2, 1.6, 1.8]
fe = np.array([ rubpy_emission, sample2, sample3, sample4, sample5, sample6, sample7])
plt.rcParams["font.family"] = "Comic Sans MS"
plt.rcParams['axes.facecolor'] = "#FFE1EF"
plt.figure(facecolor="#FFE1EF")
i = 0
while i < (fe.shape[0]):
    plt.plot(fe[i][:,0], fe[i][:, 1], label = fe_labels[i])
    i += 1
plt.xlabel("Wavelength (nm)")
plt.ylabel("Emission intensity")
plt.title("Emission intensity of Fe3+ solutions")
plt.grid(color = 'w', linestyle = '-', linewidth = 0.5)
plt.legend(title = "[Fe3+] (M)")
plt.savefig('part1_fe.png')

cu_labels = [0, 0.02, 0.04, 0.08, 0.12, 0.16, 0.18]
cu = np.array([rubpy_emission, sample8, sample9, sample10, sample11, sample12, sample13])
plt.rcParams["font.family"] = "Comic Sans MS"
plt.rcParams['axes.facecolor'] = "#FFE1EF"
plt.figure(facecolor="#FFE1EF")
i = 0
while i < (cu.shape[0]):
    plt.plot(cu[i][:,0], cu[i][:, 1], label = cu_labels[i])
    i += 1
plt.xlabel("Wavelength (nm)")
plt.ylabel("Emission intensity")
plt.title("Emission intensity of Cu2+ solutions")
plt.grid(color = 'w', linestyle = '-', linewidth = 0.5)
plt.legend(title = "[Cu2+] (M)")
plt.savefig('part1_cu.png')

#part 1

#q1

#fe
q1_feyaxis = stern_volmer_y(fe, rubpy_emission)
plt.rcParams["font.family"] = "Comic Sans MS"
plt.rcParams['axes.facecolor'] = "#FFE1EF"
plt.figure(facecolor="#FFE1EF")
q1_feyaxis, fe_peak_wl = stern_volmer_y(fe, rubpy_emission)
q1_fexaxis = fe_labels
plt.plot(q1_fexaxis, q1_feyaxis, color = "#fda0cc", marker = 'o', )
xerrb = np.std(q1_fexaxis) #multiply by two or nah bc theyre on either side already
yerrb = np.std(q1_feyaxis)
plt.errorbar(q1_fexaxis, q1_feyaxis, xerr = xerrb , yerr = yerrb, color = "#fda0cc")
m, b, *_ = stats.linregress(q1_fexaxis, q1_feyaxis)
plt.axline(xy1=(0, b), slope=m, label=f'$y = {m:.8f}x {b:+.8f}$', color = "w", linestyle = "--")
plt.xlabel("[Fe3+] (M)")
plt.legend()
plt.ylabel("I0/I")
# plt.xlim(0, q1_fexaxis[5])
plt.title(f"Stern-Volmer plot for Fe3+ @ {fe_peak_wl} nm")
plt.grid(color = 'w', linestyle = '-', linewidth = 0.5)
plt.savefig('part1_q1_Fe.png')
fe_slope = m
print(f"\nQ1\n\nFe3+:\nSlope = {fe_slope}\nwavelength = {fe_peak_wl} nm\n")

#cu
plt.rcParams["font.family"] = "Comic Sans MS"
plt.rcParams['axes.facecolor'] = "#FFE1EF"
plt.figure(facecolor="#FFE1EF")
q1_cuyaxis, cu_peak_wl = stern_volmer_y(cu, rubpy_emission)
q1_cuxaxis = cu_labels
plt.plot(q1_cuxaxis, q1_cuyaxis, color = "#fda0cc", marker = 'o', )
xerrb = np.std(q1_cuxaxis) #multiply by two or nah bc theyre on either side already
yerrb = np.std(q1_cuyaxis)
plt.errorbar(q1_cuxaxis, q1_cuyaxis, xerr = xerrb , yerr = yerrb, color = "#fda0cc")
m, b, *_ = stats.linregress(q1_cuxaxis, q1_cuyaxis)
plt.axline(xy1=(0, b), slope=m, label=f'$y = {m:.8f}x {b:+.8f}$', color = "w", linestyle = "--")
plt.xlabel("[Cu2+] (M)")
plt.legend()
# plt.xlim(0, q1_cuxaxis[5])
plt.ylabel("I0/I")
plt.title(f"Stern-Volmer plot for Cu2+ @ {cu_peak_wl} nm")
plt.grid(color = 'w', linestyle = '-', linewidth = 0.5)
plt.savefig('part1_q1_Cu.png')
cu_slope = m
print(f"\nCu3+:\nSlope = {cu_slope}\nwavelength = {cu_peak_wl} nm\n")


#q2 find rate coeffs for each
kq_fe = rate_coeff(fe_slope)
kq_cu = rate_coeff(cu_slope)
print(f"\nQ2\nFe3+ kq = {kq_fe} s^(-2)\nCu2+ kq = {kq_cu} s^(-2)\n")

#q3 find rate coeff for ET for Rubpy2+* to Fe3+ and Cu2+
#reduction potentials
fe_red_potential = 0.77
cu_red_potential = 0.16

#rates
fe_k22 = 4
cu_k22 = 10 ** (-5)

#electrons exchanged
fe_z = 3
cu_z = 2


fe_kET = kET(fe_k22, fe_red_potential, fe_z)
cu_kET = kET(cu_k22, cu_red_potential, cu_z)

print(f"\nQ3\nRate coefficient for electron transfer:\nFe: {fe_kET} M^(-1)s^(-1)\nCu: {cu_kET} M^(-1)s^(-1)\n")

f = open("part1_q1_fe_slope.txt", "w")
f.write(str(fe_slope))
f.close()

f = open("part1_q1_fe_wl.txt", "w")
f.write(str(fe_peak_wl))
f.close()


f = open("part1_q1_cu_slope.txt", "w")
f.write(str(cu_slope))
f.close()

f = open("part1_q1_cu_wl.txt", "w")
f.write(str(cu_peak_wl))
f.close()

f = open("part1_q2_fe.txt", "w")
f.write(str(kq_fe))
f.close()

f = open("part1_q2_cu.txt", "w")
f.write(str(kq_cu))
f.close()

f = open("part1_q3_fe.txt", "w")
f.write(str(fe_kET))
f.close()

f = open("part1_q3_cu.txt", "w")
f.write(str(cu_kET))
f.close()
