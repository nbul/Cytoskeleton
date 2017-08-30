#!/bin/bash                               
#$ -l rmem=48G                      
#$ -M n.bulgakova@sheffield.ac.uk
#$ -m bea

cd /data/md1nbu/Density/noEB1DN

module load apps/matlab            

matlab -nodesktop -r ‘actin_cluster_v1’
