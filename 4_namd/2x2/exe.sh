#!/bin/bash
#rm -rf *namd_res_* avg_deco_* run_namd_* slurm-*
# igeo values
for i in {0..1000..10}; do
  # power values
  for j in {0..0}; do
    echo $j
    cp run_namd.py "run_namd_$((i))_$((j)).py"
    sed -i "s/igeo =.*/igeo =$i/g;s/power =.*/power =$((j))/g" "run_namd_$((i))_$((j)).py"
    sed -i "s/python run.*/python run_namd_$((i))_$((j)).py/g" submit.slm
    sbatch submit.slm
    echo $i
  done
done
