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

with open('/global/homes/t/td115/JETSCAPE3.0-Tequila/config/jetscape_AA.xml', 'rb') as xml_file: 
    tree = ET.parse(xml_file)
    root = tree.getroot()
        
    name = root.find('outputFilename') 
    file_name = '/global/cscratch1/sd/td115/output/AuAu200/data_cut2_coef2/%.6f_i%d.dat' %(new_pT_hat_min, task)
    name.text = file_name
    name.set('updated', 'yes')

    hard = root.find('Hard')
    pythia = hard.find('PythiaGun')
    pT_hat_min = pythia.find('pTHatMin')
    pT_hat_max = pythia.find('pTHatMax')

    pT_hat_min.text = str(new_pT_hat_min)
    pT_hat_min.set('updated', 'yes')
    pT_hat_max.text = str(new_pT_hat_max)
    pT_hat_max.set('updated', 'yes')

    tree.write('/global/homes/t/td115/JETSCAPE3.0-Tequila/config/jetscape_AA.xml')

os.system('/global/homes/t/td115/JETSCAPE3.0-Tequila/build/runJetscape /global/homes/t/td115/JETSCAPE3.0-Tequila/config/jetscape_AA.xml')
