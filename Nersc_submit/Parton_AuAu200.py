#!/usr/bin/env python

import os

pTHatBin = [5., 10., 20., 40., 60., 80]

for pTHat in pTHatBin: 
    os.system('./FinalStatePartons $SCRATCH/output/STAT/AuAu200/20211105/%.6f_i0.dat $SCRATCH/output/STAT/AuAu200/20211105/partons/%.6f_i0_partons.txt' %(pTHat, pTHat))
