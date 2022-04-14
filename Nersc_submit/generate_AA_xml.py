#!/usr/bin/env python

import os
import xml.etree.ElementTree as ET
import argparse

pTHat_list = [1., 2., 3., 4., 5., 7., 9., 11., 13., 15., 17., 20., 25., 30., 35., 40., 45., 50., 55., 60., 70., 80., 90., 100.]
# pTHat_list = [10., 20., 50., 80., 110., 160., 210., 260., 310., 400., 500., 600., 800., 1000., 1380.]
for i, new_pT_hat_min in enumerate(pTHat_list[:-1]): 
    for task in range(80): 
        with open('../config/Tequila_AuAu200.xml', 'rb') as xml_file: 
            tree = ET.parse(xml_file)
            root = tree.getroot()
        
            name = root.find('outputFilename') 
            file_name = '/global/cscratch1/sd/td115/output/Tequila/AuAu200/centrality40-50/dp8/%.6f_i%d' %(new_pT_hat_min, task)
            name.text = file_name
            name.set('updated', 'yes')

            random = root.find('Random')
            seed = random.find('seed')
            seed.text = str(task)
            seed.set('updated', 'yes')

            hard = root.find('Hard')
            pythia = hard.find('PythiaGun')
            pT_hat_min = pythia.find('pTHatMin')
            pT_hat_max = pythia.find('pTHatMax')

            pT_hat_min.text = str(new_pT_hat_min)
            pT_hat_min.set('updated', 'yes')
            new_pT_hat_max = pTHat_list[i+1]
            pT_hat_max.text = str(new_pT_hat_max)
            pT_hat_max.set('updated', 'yes')

            eloss = root.find('Eloss')
            tequila = eloss.find('Tequila')
            alphas_soft = tequila.find('alphas_soft')
            alphas_soft.text = str(0.0738)
            alphas_soft.set('updated', 'yes')

            alphas_hard_elas = tequila.find('alphas_hard_elas')
            alphas_hard_elas.text = str(0.0613)
            alphas_hard_elas.set('updated', 'yes')
            
            alphas_hard_inel = tequila.find('alphas_hard_inel')
            alphas_hard_inel.text = str(0.106)
            alphas_hard_inel.set('updated', 'yes')
            

            tree.write('../config_files/Tequila_AuAu200_dp8_%.6f_i%d.xml' %(new_pT_hat_min, task), xml_declaration=True, encoding='utf-8')

