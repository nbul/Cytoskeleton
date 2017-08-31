#!/bin/bash                               
#$ -l rmem=10G
#$ -l h_rt=24:00:00                      
#$ -M n.bulgakova@sheffield.ac.uk
#$ -m bea

cd /data/md1nbu/Density/Patronin_homo

module load apps/matlab            

matlab -nodesktop -r 'ActinClusterV1'
