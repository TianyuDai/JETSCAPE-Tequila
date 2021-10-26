#!/usr/bin/env python

import os
import xml.etree.ElementTree as ET
import argparse

pTHat_list = [5., 10., 20., 40., 60., 80., 100.]
coef_list = [1.]
cut_list = [1.]
point_label = [1]

for i_point, cut, coef in zip(point_label, cut_list, coef_list): 
    for i, new_pT_hat_min in enumerate(pTHat_list[:-1]): 
        for task in range(10): 
            with open('jetscape_AA.xml', 'rb') as xml_file: 
                tree = ET.parse(xml_file)
                root = tree.getroot()
        
                name = root.find('outputFilename') 
                file_name = '/home/td115/research/Result/JETSCAPE3.0/test_for_Nersc/data_cut%d_coef1/%.6f_i%d' %(cut, new_pT_hat_min, task)
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

                eloss = root.find('Eloss')
                tequila = eloss.find('Tequila')
                mu_w = tequila.find('muomega_over_T')
                mu_q = tequila.find('muqperp_over_T')
                qhat_coef = tequila.find('qhat_coef')

                mu_w.text = str(cut)
                mu_w.set('updated', 'yes')
                mu_q.text = str(cut)
                mu_q.set('updated', 'yes')
                qhat_coef.text = str(coef)
                qhat_coef.set('updated', 'yes')

                tree.write('config_files/jetscape_AA_cut%.1f_coef%.1f_%.6f_i%d.xml' %(cut, coef, new_pT_hat_min, task), xml_declaration=True, encoding='utf-8')

