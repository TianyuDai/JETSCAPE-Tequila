#!/usr/bin/env python

import os
import xml.etree.ElementTree as ET
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--task', type=int, help='parse task number to output')
parser.add_argument('--pTHatMin', type=float, help='parse minimum of pTHatBin')
parser.add_argument('--pTHatMax', type=float, help='parse maximum of pTHatBin')
parser.add_argument('--dp', type=int, help='parse maximum of pTHatBin')
# parser.add_argument('--coef', type=float, help='parse coefficient of qhat')
# parser.add_argument('--cut', type=float, help='parse cut of momentum transfer')
args = parser.parse_args()

new_pT_hat_min = args.pTHatMin
new_pT_hat_max = args.pTHatMax
task = args.task

os.system('./runJetscape ../dps/dp%d/Tequila_AuAu200_%.6f_i%d.xml' %(dp, new_pT_hat_min, task))
