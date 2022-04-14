# pTHatBins = [5, 10, 20, 40, 60, 80, 100]
pTHatBins = [1., 2., 3., 4., 5., 7., 9., 11., 13., 15., 17., 20., 25., 30., 35., 40., 45., 50., 55., 60., 70., 80., 90., 100.]

print('#!/usr/bin/env bash')
print()

for i in range(40): 
    for j in range(len(pTHatBins)-1): 
        pTHatMin = pTHatBins[j]
        pTHatMax = pTHatBins[j+1]
        print("cd /global/homes/t/td115/JETSCAPE3.4/JETSCAPE/build && python3 AA_wrapper_dp8.py --task %d --pTHatMin %.6f --pTHatMax %.6f" %(i, pTHatMin, pTHatMax))
