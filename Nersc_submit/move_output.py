#!/usr/bin/env python

import os
import xml.etree.ElementTree as ET
import argparse

dp_list = [x for x in range(43)]
dp_list.append(60)
dp_list.append(72)
dp_list.append(73)
dp_list.append(88)
dp_list.append(85)
dp_list.append(92)
dp_list.append(128)
dp_list.append(134)

for dp in dp_list: 
    os.system('cp /global/cscratch1/sd/td115/output/Tequila/running_coupling/AuAu200/centrality0-10/dp%d/hadrons/AA200_dp%d_pion_cs.txt ../output/' %(dp, dp))
