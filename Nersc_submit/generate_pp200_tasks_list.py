# pTHatBins = [5, 10, 20, 40, 60, 80, 100]
# pTHatBins = [1., 2., 3., 4., 5., 7., 9., 11., 13., 15., 17., 20., 25., 30., 35., 40., 45., 50., 55., 60., 70., 80., 90., 100.]
# pTHatBins = [10., 20., 50., 80., 110., 160., 210., 260., 310., 400., 500., 600., 800., 1000., 1380.]

pTHatBins = [1., 3., 5., 8., 11., 15., 20., 25., 30., 35., 40., 45., 50., 60., 70., 80., 100.]

for i in range(1): 
    for j in range(len(pTHatBins)-1): 
        pTHatMin = pTHatBins[j]
        pTHatMax = pTHatBins[j+1]
        print("cd /global/homes/t/td115/running_coupling/JETSCAPE-Tequila/build && python3 pp_wrapper.py --task %d --pTHatMin %.6f --pTHatMax %.6f" %(i, pTHatMin, pTHatMax))
