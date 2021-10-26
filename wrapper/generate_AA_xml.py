#!/usr/bin/env python

import os
import xml.etree.ElementTree as ET
import argparse

pTHat_list = [5., 10., 20., 40., 60., 80., 100.]
for i, new_pT_hat_min in enumerate(pTHat_list[:-1]): 
    for task in range(100): 
        with open('jetscape_AA.xml', 'rb') as xml_file: 
            tree = ET.parse(xml_file)
            root = tree.getroot()
        
            name = root.find('outputFilename') 
            file_name = '/global/cscratch1/sd/td115/output/AuAu200/data_cut2_coef05/%.6f_i%d' %(new_pT_hat_min, task)
            name.text = file_name
            name.set('updated', 'yes')

            hard = root.find('Hard')
            pythia = hard.find('PythiaGun')
            pT_hat_min = pythia.find('pTHatMin')
            pT_hat_max = pythia.find('pTHatMax')

            pT_hat_min.text = str(new_pT_hat_min)
            pT_hat_min.set('updated', 'yes')
            new_pT_hat_max = pTHat_list[i+1]
            pT_hat_max.text = str(new_pT_hat_max)
            pT_hat_max.set('updated', 'yes')

            tree.write('config_files/jetscape_AA_cut2_coef05_%.6f_i%d.xml' %(new_pT_hat_min, task), xml_declaration=True, encoding='utf-8')

