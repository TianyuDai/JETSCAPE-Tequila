#!/usr/bin/env python

import os
import xml.etree.ElementTree as ET
import argparse


# pTHat_list = [1., 3., 5., 8., 11., 15., 20., 25., 30., 35., 40., 45., 50., 60., 70., 80., 100.]
# pTHat_list = [5., 8., 11.]
pTHat_list = [3., 4., 5., 7., 9., 11., 13., 15., 17., 20., 25., 30., 35., 40., 45., 50., 55., 60., 70., 80., 90., 100.]
# pTHat_list = [10., 20., 50., 80., 110., 160., 210., 260., 310., 400., 500., 600., 800., 1000., 1380.]
for i, new_pT_hat_min in enumerate(pTHat_list[:-1]):
    for task in range(20): 
        with open('../config/Tequila_AuAu200.xml', 'rb') as xml_file: 
            tree = ET.parse(xml_file)
            root = tree.getroot()
        
            name = root.find('outputFilename') 
            file_name = '/global/cscratch1/sd/td115/output/Tequila/running_coupling/AuAu200/centrality0-10/tau0_0.8/%.6f_i%d' %(new_pT_hat_min, task)
            name.text = file_name
            name.set('updated', 'yes')

            random = root.find('Random')
            seed = random.find('seed')
            seed.text = str(task)
            seed.set('updated', 'yes')

            IS = root.find('IS')
            is_path = IS.find('initial_profile_path')
            is_path.text = '/global/cscratch1/sd/td115/AuAu200_hydro/centrality0-10/part-%d' %task
            is_path.set('updated', 'yes')

            hydro = root.find('Hydro')
            hydro_file = hydro.find('hydro_from_file')
            hydro_path = hydro_file.find('hydro_files_folder')
            hydro_path.text = '/global/cscratch1/sd/td115/AuAu200_hydro/centrality0-10/part-%d' %task
            hydro_path.set('updated', 'yes')

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
            matter = eloss.find('Matter')
            Q0 = matter.find('Q0')
            Q0.text = str(1.6)
            Q0.set('updated', 'yes')

            tequila = eloss.find('Tequila')
            beta1 = tequila.find('beta_para')
            beta1.text = str(0.)
            beta1.set('updated', 'yes')

            beta2 = tequila.find('beta_perp')
            beta2.text = str(0.)
            beta2.set('updated', 'yes')

            T_star = tequila.find('T_star')
            T_star.text = str(0.3)
            T_star.set('updated', 'yes')
            
            Q0 = tequila.find('Q0')
            Q0.text = str(1.6)
            Q0.set('updated', 'yes')
            
            alphas_hard_inel = tequila.find('alphas_hard_inel')
            alphas_hard_inel.text = str(0.25)
            alphas_hard_inel.set('updated', 'yes')
            
            tree.write('../config_files/Tequila_AuAu200_centrality0-10_tau0_0.8_%.6f_i%d.xml' %(new_pT_hat_min, task), xml_declaration=True, encoding='utf-8')

