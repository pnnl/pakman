#!/bin/bash -l

#SBATCH -q debug
#SBATCH -N 64
#SBATCH -t 0:10:00
#SBATCH -J human_t1k
#SBATCH -e p_io.error
#SBATCH -o p_io.out
#SBATCH -C haswell

export STRIPE_DATA_PATH=<Path to input FASTA dataset>
#lfs setstripe -c 50 -S 1m $STRIPE_DATA_PATH
#cp $DATA_PATH $STRIPE_DATA_PATH

module unload darshan
export MPICH_MPIIO_HINTS="*:romio_cb_write=enable:cb_nodes=256:romio_ds_write=disable"
#export MPICH_MPIIO_HINTS_DISPLAY=1

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

export BINARY=<Path to pakman binary>

export JOB_OUTPUT=$PWD
export JOB_PREFIX=cori_haswell_p1k_Human_100x
export MAX_KMER_COUNT=100000000
export COVERAGE=100
export BUCKET_CUTOFF=21
export MERGE_CUTOFF=100000
export AVG_READ_LEN=100

export DATE_WITH_TIME=`date "+%Y%m%d-%H%M%S"`

for procs in 1024
do

#echo $DATE_WITH_TIME >> $JOB_OUTPUT/$JOB_PREFIX.out

srun -n $procs -c 4 --cpu_bind=cores $BINARY -f $STRIPE_DATA_PATH -b $MAX_KMER_COUNT -r $AVG_READ_LEN -c $COVERAGE -t $BUCKET_CUTOFF -n $MERGE_CUTOFF >> $JOB_OUTPUT/$JOB_PREFIX.txt

echo "----------------------------------- job's done -----------------------------------------" >> $JOB_OUTPUT/$JOB_PREFIX.out

done

