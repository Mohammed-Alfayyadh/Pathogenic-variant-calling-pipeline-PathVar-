#!/bin/bash -l
#PBS -N PathVar-Analysis
#PBS -l select=1:ncpus=16:mem=512g
#PBS -l walltime=050:00:00
cd $PBS_O_WORKDIR

conda activate PathVar
python /home/n9261834/phd/data/data-hg38/GATK.py
python /home/n9261834/phd/data/data-hg38/VEP.py
python /home/n9261834/phd/data/data-hg38/ACMG.py
python /home/n9261834/phd/data/data-hg38/Concordant-SNVs.py
