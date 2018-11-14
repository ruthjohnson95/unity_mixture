#!/usr/bin/env sh


# select steps of the pipeline 
STEPS=$1

# default is all steps 
if [ -z "$STEPS"]
	STEPS="1,2,3,4"

# simulation params 
SIM_NAME=test_identity 
P_VEC=".25 .25 .50"
MU_VEC="-.10 0 .10"
SIGMA_VEC=".001 .001 .001"
LD_FILE=${DATA_DIR}/identity.100.ld 
M=100 
N=100000
SEED=2018 # can replace with SGE_TASK_ID 

DATE=`date '+%Y-%m-%d %H:%M:%S'`
echo $DATE" Starting simulation for unity-mixture: "${SIM_NAME}

# global paths 

# path to folder 
MASTER_PATH=/Users/ruthiejohnson/Development/unity_mixture
SCRIPT_DIR=${MASTER_PATH}/scripts 
DATA_DIR=${MASTER_PATH}/data 

# Hoffman paths 
#source /u/local/Modules/default/init/modules.sh
#module load python/2.7

# data will be output to DATA_DIR 
mkdir -p $DATA_DIR 

# STEP 1: simulate gwas 
if [[ "$STEPS" =~ "1" ]]
then
	DATE=`date '+%Y-%m-%d %H:%M:%S'`
	echo $DATE" Simulting GWAS effect sizes"
	python ${SCRIPT_DIR}/simulate.py --name $SIM_NAME --p_vec $P_VEC --mu_vec $MU_VEC --sigma_vec $SIGMA_VEC --ld_file $LD_FILE --M $M --N $N --seed $SEED --outdir $DATA_DIR
fi

# STEP 2: transform betas 
if [[ "$STEPS" =~ "2" ]]
then
	DATE=`date '+%Y-%m-%d %H:%M:%S'`
	echo $DATE" Transforming GWAS effect sizes"
	GWAS_FILE=${DATA_DIR}/${SIM_NAME}.${SEED}.txt 
	python ${SCRIPT_DIR}/transform_betas.py --gwas_file --ld_file $LD_FILE
fi 


# STEP 3: take 1/2 power of LD 
if [[ "$STEPS" =~ "3" ]]
then
	python ${SCRIPT_DIR}/half_ld.py --ld_file $LD_FILE
fi

# STEP 4: run inference 

