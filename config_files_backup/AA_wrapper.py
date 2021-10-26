#!/usr/bin/env python

import os
import xml.etree.ElementTree as ET
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--task', type=int, help='parse task number to output')
parser.add_argument('--pTHatMin', type=float, help='parse minimum of pTHatBin')
parser.add_argument('--pTHatMax', type=float, help='parse maximum of pTHatBin')
parser.add_argument('--cut', type=float, help='parse hard-soft cutoff of momentum exchange')
parser.add_argument('--coef', type=float, help='parse coefficient of qhat')
args = parser.parse_args()

new_pT_hat_min = args.pTHatMin
new_pT_hat_max = args.pTHatMax
task = args.task
new_cut = args.cut
new_coef = args.coef

os.system('./runJetscape config_files/jetscape_AA_cut%.1f_coef%.1f_%.6f_i%d.xml' %(new_cut, new_coef, new_pT_hat_min, task))
