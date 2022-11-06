#!/usr/bin/env python

import os
dp_list = [x for x in range(43)]
dp_list.append(60)
dp_list.append(72)
dp_list.append(73)
dp_list.append(88)

for dp in dp_list: 
    os.system('cp /global/cscratch1/sd/td115/output/Tequila/running_coupling/AuAu200/centrality0-10/dp0/hadrons/cross_section.py /global/cscratch1/sd/td115/output/Tequila/running_coupling/AuAu200/centrality0-10/dp%d/hadrons' %dp)
    os.system('cp /global/cscratch1/sd/td115/output/Tequila/running_coupling/AuAu200/centrality0-10/dp0/hadrons/sigmaGenErr.txt /global/cscratch1/sd/td115/output/Tequila/running_coupling/AuAu200/centrality0-10/dp%d/hadrons' %dp)

    os.system('python /global/cscratch1/sd/td115/output/Tequila/running_coupling/AuAu200/centrality0-10/dp%d/hadrons/cross_section.py --dp %d' %(dp, dp))

