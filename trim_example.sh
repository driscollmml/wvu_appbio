#! /bin/bash



#PBS -q comm_mmem_day

#PBS -lwalltime=24:00:00

#PBS -lnodes=1:ppn=1,pvmem=54gb

#PBS -N JOB_NAME

#PBS -m abe -M YOUR_EMAIL_HERE




module load genomics/bioconda
source /shared/software/miniconda3/etc/profile.d/conda.sh
conda activate tpd0001


cd $SCRATCH

trimmomatic PE -threads 6 SRR254178_1.fastq SRR254178_2.fastq -baseout ec_cropped_70.fastq CROP:70 MINLEN:70


conda deactivate


