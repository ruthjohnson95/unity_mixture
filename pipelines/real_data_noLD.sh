#!/usr/bin/env sh

GWAS_FILE=$1 

# path to folder 
MASTER_PATH=/u/home/r/ruthjohn/ruthjohn/bayesM_noLD/unity_mixture/
SCRIPT_DIR=${MASTER_PATH}/scripts 
SRC_DIR=${MASTER_PATH}/src 
DATA_DIR=${MASTER_PATH}/data 

# simulation params 
SIM_NAME=height 
BINS=10
N=100000
SEED=2018 # can replace with SGE_TASK_ID
ITS=10

DATE=`date '+%Y-%m-%d %H:%M:%S'`
echo $DATE" Starting simulation for unity-mixture: "${SIM_NAME}

# Hoffman paths 
source /u/local/Modules/default/init/modules.sh
module load python/2.7

# data will be output to DATA_DIR 
mkdir -p $DATA_DIR 

python ${SRC_DIR}/mixture_em_noLD.py --name $SIM_NAME --gwas_file $GWAS_FILE --bins $BINS --N $N --seed $SEED --outdir $DATA_DIR  --its $ITS 

