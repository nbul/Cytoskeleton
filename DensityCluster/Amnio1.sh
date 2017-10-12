#!/bin/bash                               
#$ -l rmem=10G
#$ -l mem=25G
#$ -l h_rt=48:00:00                      
#$ -M n.bulgakova@sheffield.ac.uk
#$ -m bea

cd /data/md1nbu/Density/Amnioserosa1

module load apps/matlab            

matlab -nodesktop -r 'ActinClusterV1'
