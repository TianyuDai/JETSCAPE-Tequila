# pTHatBins = [5, 10, 20, 40, 60, 80, 100]
pTHatBins = [10., 20., 50., 80., 110., 160., 210., 260., 310., 400., 500., 600., 800., 1000., 1380.]

for i in range(10): 
    for j in range(len(pTHatBins)-1): 
        pTHatMin = pTHatBins[j]
        pTHatMax = pTHatBins[j+1]
        print("cd /global/homes/t/td115/JETSCAPE3.4/JETSCAPE/build && python3 AA2760_wrapper.py --task %d --pTHatMin %.6f --pTHatMax %.6f" %(i, pTHatMin, pTHatMax))
