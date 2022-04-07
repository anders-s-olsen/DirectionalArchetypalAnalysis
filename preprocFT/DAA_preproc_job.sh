#!/bin/sh
#BSUB -J DAAjob
#BSUB -q hpc
# -- specify that we need 2GB of memory per core/slot -- 
#BSUB -R "rusage[mem=128GB]"
# -- Notify me by email when execution begins --
#BSUB -B
# -- Notify me by email when execution ends   --
#BSUB -N
#BSUB -o DAAjob_out_%J.txt
# -- Error File --
#BSUB -e DAAjob_err_%J.txt
# -- estimated wall clock time (execution time): hh:mm -- 
#BSUB -W 06:00 
# -- Number of cores requested -- 
#BSUB -n 1 
# -- end of LSF options -- 

# -- commands you want to execute -- 
module load matlab/R2020a
matlab -nodisplay -batch DAA_preproc
