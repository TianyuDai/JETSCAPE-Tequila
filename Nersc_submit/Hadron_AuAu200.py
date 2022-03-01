#!/usr/bin/env python

import os

# pTHat_list = [5., 10., 20., 40., 60., 80]
pTHat_list = [1., 2., 3., 4., 5., 7., 9., 11., 13., 15., 17., 20., 25., 30., 35., 40., 45., 50., 55., 60., 70., 80., 90., 100.]

for pTHat in pTHatBin: 
    for i in range(10): 
        os.system('./FinalStateHadrons $SCRATCH/output/STAT/AuAu200/20220106/%.6f_i%d.dat $SCRATCH/output/STAT/AuAu200/20220106/hadrons/%.6f_i%d_chargedHadrons.txt' %(pTHat, i, pTHat, i))
