import argparse

# pTHatBins = [5, 10, 20, 40, 60, 80, 100]
# pTHatBins = [1., 3., 5., 8., 11., 15., 20., 25., 30., 35., 40., 45., 50., 60., 70., 80., 100.]
# pTHatBins = [3., 5., 8., 11.]
# pTHatBins = [60., 70.]
# pTHatBins = [3., 4., 5., 7., 9., 11., 13., 15., 17., 20., 25., 30., 35., 40., 45., 50., 55., 60., 70., 80., 90., 100.]
# pTHatBins = [4., 6., 10., 15., 20., 30., 45., 60., 80., 110., 160., 210., 260., 310., 400., 500., 600., 800., 1000., 1380.]
# pTHatBins = [3., 5., 8., 11.]
# pTHatBins = [10., 20., 50., 80., 110., 160., 210., 260., 310., 400., 500., 600., 800., 1000., 1380.]
pTHatBins = [4., 6., 10., 15., 20., 30., 45., 60., 80., 110., 160., 210., 260., 310., 400., 500., 600., 800., 1000., 1380.]
# pTHatBins = [4., 6., 10., 15., 20., 30., 45., 60., 80.]

description = 'generate xml input files, input design point index'
parser = argparse.ArgumentParser(description=description)
parser.add_argument('--dp', type=int, default=0, help='design point index')
args = parser.parse_args()

dp = args.dp

print('#!/usr/bin/env bash\n')

for j in range(len(pTHatBins)-1): 
    for i in range(40): 
        pTHatMin = pTHatBins[j]
        pTHatMax = pTHatBins[j+1]
        print("cd /global/homes/t/td115/running_coupling/JETSCAPE-Tequila/build && python3 PbPb30-40_wrapper.py --dp %d --task %d --pTHatMin %.6f --pTHatMax %.6f" %(dp, i, pTHatMin, pTHatMax))
