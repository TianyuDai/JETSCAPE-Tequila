#!/usr/bin/env bash

# rm /global/homes/t/td115/running_coupling/JETSCAPE-Tequila/build/*.tfin

for i in "$@";do
mkdir /global/cscratch1/sd/td115/output/Tequila/running_coupling/PbPb2760/centrality0-5/dps/dp$i
mkdir /global/cscratch1/sd/td115/output/Tequila/running_coupling/PbPb2760/centrality0-5/dps/dp$i/hadrons

mkdir /global/homes/t/td115/running_coupling/JETSCAPE-Tequila/dps/PbPb2760/centrality0-5/dp$i

cd /global/homes/t/td115/running_coupling/JETSCAPE-Tequila/Nersc_submit

python generate_PbPb_xml.py --dp $i
done

python generate_PbPb_tasks_list.py --dp $1 > ../build/task_list_dps/PbPb2760/centrality0-5/Tequila_PbPb2760_dp$1_tasks_list.txt

sbatch Tequila_PbPb2760_dp $1 
