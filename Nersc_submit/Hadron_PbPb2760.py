#!/usr/bin/env python

import os

print('#!/usr/bin/env bash\n')

# pTHat_list = [5., 10., 20., 40., 60., 80]
# pTHat_list = [1., 3., 5., 8., 11., 15., 20., 25., 30., 35., 40., 45., 50., 60., 70., 80., 100.]
# pTHat_list = [1., 3., 5., 8., 11., 15., 20., 25., 80., 100.]
# pTHat_list = [3., 5., 8., 11.]
# pTHat_list = [3., 4., 5., 7., 9., 11., 13., 15., 17., 20., 25., 30., 35., 40., 45., 50., 55., 60., 70., 80., 90., 100.]
pTHat_list = [4., 6., 10., 15., 20., 30., 45., 60., 80., 110., 160., 210., 260., 310., 400., 500., 600., 800., 1000., 1380.]
# pTHat_list = [4., 6., 10., 15., 20., 30., 45., 60., 80.]

# dp_list = range(2, 3, 1)
# dp_list = [16, 17, 21]

for dp in range(50, 60): 
    for pTHat in pTHat_list[:-1]: 
        for i in range(40): 
            # print('cd /global/homes/t/td115/running_coupling/JETSCAPE-Tequila/build && ./FinalStateHadrons $SCRATCH/output/Tequila/running_coupling/PbPb2760/centrality0-5/dps/dp%d/%.6f_i%d.dat $SCRATCH/output/Tequila/running_coupling/PbPb2760/centrality0-5/dps/dp%d/hadrons/%.6f_i%d_chargedHadrons.txt' %(dp, pTHat, i, dp, pTHat, i))
            print('cd /global/homes/t/td115/running_coupling/JETSCAPE-Tequila/build && ./FinalStateHadrons $SCRATCH/output/Tequila/running_coupling/PbPb2760/centrality30-40/dps/dp%d/%.6f_i%d.dat $SCRATCH/output/Tequila/running_coupling/PbPb2760/centrality30-40/dps/dp%d/hadrons/%.6f_i%d_chargedHadrons.txt' %(dp, pTHat, i, dp, pTHat, i))
            # print('cd /global/homes/t/td115/running_coupling/JETSCAPE-Tequila/build && ./FinalStateHadrons $SCRATCH/output/Tequila/running_coupling/PbPb2760/centrality30-40/val/%.6f_i%d.dat $SCRATCH/output/Tequila/running_coupling/PbPb2760/centrality30-40/val/hadrons/%.6f_i%d_chargedHadrons.txt' %(pTHat, i, pTHat, i))
            # print('cd /global/homes/t/td115/running_coupling/JETSCAPE-Tequila/build && ./FinalStateHadrons $SCRATCH/output/Tequila/running_coupling/PbPb2760/centrality0-5/val/%.6f_i%d.dat $SCRATCH/output/Tequila/running_coupling/PbPb2760/centrality0-5/val/hadrons/%.6f_i%d_chargedHadrons.txt' %(pTHat, i, pTHat, i))
