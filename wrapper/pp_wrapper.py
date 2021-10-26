#!/usr/bin/env python

import os
import xml.etree.ElementTree as ET
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--task', type=int, help='parse task number to output')
parser.add_argument('--pTHatMin', type=float, help='parse minimum of pTHatBin')
parser.add_argument('--pTHatMax', type=float, help='parse maximum of pTHatBin')
args = parser.parse_args()

new_pT_hat_min = args.pTHatMin
new_pT_hat_max = args.pTHatMax
task = args.task

os.system('./runJetscape config_files/jetscape_pp_%.6f_i%d.xml' %(new_pT_hat_min, task))
