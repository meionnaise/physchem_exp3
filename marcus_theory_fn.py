import numpy as np
def kET(k22, x_red_potential, z):
    r_red_potential = -0.84
    r_ox = -r_red_potential
    R = 8.3145 #J mol^-1 K^-1
    T = 298.15 # K (assume 25 deg)
    F = 96485.3399 #Cmol^(-1) #faraday constant
    Z12 = 10 ** 11 #factor related to collision frequency
    k11 = 10 ** 8 #(Rubpy kii)

    cell_potential = x_red_potential + r_ox #get cell potential
    K12 = np.exp(cell_potential * ( (z * F) / (R * T))) # use eqn 13 to get from cell potential - K12
    f12 = np.exp(((np.log(K12)) ** 2) / (4 * np.log((k11 * k22) / (Z12 ** 2)))) #use K12 --> f12 (eqn 12)
    kET = np.sqrt(k11 * k22 * K12 * f12) # use f12 , k11, k22, k12 to get kET

    return kET