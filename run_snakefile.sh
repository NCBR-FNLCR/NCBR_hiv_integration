#! /bin/bash

###################
#
# Launching shell script for NIAID CSI batch processing of WES data
#
###################
module load snakemake/5.7.4
module load python/3.7

mkdir -p snakejobs
mkdir -p BATCH_QC

raw_dir=$1

if ! test -d "rawdata"; then
    if test -d $raw_dir; then
        echo "Linking rawdata subdirectory to $raw_dir"
        ln -s $raw_dir rawdata
    else
        echo "Unable to locate raw data directory $raw_dir"
        echo "Exiting"
        exit
    fi
else
    echo "input directory rawdata already exists"
fi


##
## Run snakemake
##
echo "Run snakemake"

CLUSTER_OPTS="sbatch --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name={params.rname} -e snakejobs/slurm-%j_{params.rname}.out -o snakejobs/slurm-%j_{params.rname}.out --chdir=$batchdir"

if [ "$2" == "npr" ]
then
    snakemake -npr --snakefile Snakefile
fi

if [ "$2" == "process" ]
then
    snakemake --stats snakemake.stats --rerun-incomplete -j 150 --cluster "$CLUSTER_OPTS" --cluster-config cluster.json --keep-going --snakefile Snakefile 2>&1|tee -a batch_processing.log
fi
