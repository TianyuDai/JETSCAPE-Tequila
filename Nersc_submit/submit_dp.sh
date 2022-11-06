#!/usr/bin/env bash

# rm /global/homes/t/td115/running_coupling/JETSCAPE-Tequila/build/*.tfin

for i in "$@";do
mkdir /global/cscratch1/sd/td115/output/Tequila/running_coupling/AuAu200/centrality0-10/dps/dp$i
mkdir /global/cscratch1/sd/td115/output/Tequila/running_coupling/AuAu200/centrality0-10/dps/dp$i/hadrons

mkdir /global/homes/t/td115/running_coupling/JETSCAPE-Tequila/dps/AuAu200/centrality0-10/dp$i

cd /global/homes/t/td115/running_coupling/JETSCAPE-Tequila/Nersc_submit

python generate_AA_xml.py --dp $i
done

python generate_tasks_list.py --dp1 $1 --dp2 $4 > ../build/task_list_dps/AuAu200/centrality0-10/Tequila_AuAu200_dp$1_tasks_list.txt

sbatch Tequila_AuAu200_dp $1 
