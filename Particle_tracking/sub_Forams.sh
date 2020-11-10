#!/bin/bash
##SBATCH -t 5-00:00:00
#SBATCH -t 24:00:00
###SBATCH -N 6 --ntasks-per-node=9
##SBATCH --cpus-per-task=3
##SBATCH --mem-per-cpu=16000
#SBATCH -p broadwell
###SBATCH -p fat

for dd in 50 300 -1
do
  python atsf.py $dd & 
  sleep 5
done

wait
