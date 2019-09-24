#! /bin/bash


# - Spruce Knob utilizes the TORQUE/PBS resource manager to manage compute resources.
#
#    To submit a batch job, issue:       qsub job_script.sh
#    To show all queued jobs, issue:     showq or qstat
#    To show YOUR queued jobs, issue:     qstat -au <userId>
#
#    To check a job status:		checkjob <jobId>
#    To kill a queued job, issue:        canceljob <jobId>
#
#    Example PBS job scripts are located at:
#    http://wiki.hpc.wvu.edu/hpc_wiki/index.php/Sample_Job_Scripts
#
#    The following man pages provide helpful pbs information:
#	man pbs_resources_linux (How-to request cluster resources)
#	man qsub
#	man qstat
#
# - Spruce Knob utilizes Environment Modules to help manage different software
#   packages available on the cluster.  "module avail" shows the available
#   modules.
#
# - Spruce has two file systems available to users:
#   - $HOME (permanent storage that is backed up via snapshots, 10GB Limit)
#   - $SCRATCH (temporary storage that is NOT backed up, current allocation 130 TB)



# Set the queue that will get your job
#
#PBS -q debug
#
# for longer jobs, you can also use the comm_mmem_week queue
# to test a job, use the debug queue (this has a max walltime of 5 min)


# Set the walltime, which is the maximum time that your job can run in HH:MM:SS
# this can not exceed the maximum walltime determined by the queue that you set above
#
#PBS -lwalltime=00:05:00

# Set the number of nodes, and the number of processors per node (up to 12), that you want to use
#
#PBS -lnodes=1:ppn=12


## Give your job a name
#
#PBS -N newfastqc

## Provide your email address, to receive notification when your job starts and ends
#
#PBS -m abe -M njs0012




# load the bioconda environment
# this provides access to all of the software that you will need for class
# DO NOT CHANGE THESE THREE LINES
#
module load genomics/bioconda
source /shared/software/miniconda3/etc/profile.d/conda.sh
conda activate tpd0001


# cd into your working dir on the scratch drive
# your home dir is limited to 5GB so it is usually advisable to work on the 'scratch'
# you have essentially unlimited space on scratch, but files can be deleted to make room for other users
#
# the absolute path to your top-level scratch dir is /scratch/USERNAME (eg., /scratch/tpd0001)
# you can use the environmental variable $SCRATCH as a shortcut
#
cd $SCRATCH


# put your commmand(s) in here
# you will need to change the strings in CAPITAL LETTERS
# the 'j' argument tells abyss-pe how many threads to run (each thread runs on a different processor)
# the manual for abyss-pe can be found at https://github.com/bcgsc/abyss#assembling-a-paired-end-library
# the publication describing abyss-pe can be found at https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2694472/
#
fastqc SRR8245541.1.fastq




# deactivate the conda session when you have finished with your pipeline
conda deactivate
