#!/bin/bash                               
#$ -l rmem=10G
#$ -l mem=25G
#$ -l h_rt=24:00:00                      
#$ -M n.bulgakova@sheffield.ac.uk
#$ -m bea

cd /data/md1nbu/Density/Khc_homo

module load apps/matlab            

matlab -nodesktop -r 'ActinClusterV1'
