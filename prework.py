
ratios = [[0,10] , [1,9], [2, 8], [3, 7], [4, 6], [5,5], [6, 4], [7, 3], [8, 2], [9, 1], [10, 0]]

fe_conc = 1.3 * (10 ** (-3)) #M
sal_conc = 1.3 * (10 ** (-3)) #M

#number of moles in 10mL
moles = 0.1 * 1.3 * (10 ** (-3))

calcs = []
for i in ratios:
    ratio = str(ratios[i])
    fe = (ratio[0] * moles) / 10
    sal = (ratio[1] * moles) / 10
    calcs.append([ratio, fe, sal])

print(f"ratios\tfe\tsal\n")
print(calcs)
