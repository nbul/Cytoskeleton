#!/bin/bash                               
#$ -l rmem=48G                      
#$ -M n.bulgakova@sheffield.ac.uk
#$ -m bea

module load apps/matlab            

R CMD BATCH correlate_counts_qsub.R correlate_counts_qsub.R.o$JOB_ID
matlab -nodesktop -r 'validation_cluster'
