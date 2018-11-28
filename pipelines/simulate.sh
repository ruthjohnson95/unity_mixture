#!/usr/bin/env sh


# select steps of the pipeline
STEPS=$1
SEED=$2

# default is all steps
if [ -z "$STEPS" ]
then
	STEPS="1,2,3,4,8"
fi

# default is all steps
if [ -z "$SEED" ]
then
        SEED=2018
fi

# path to folder
#MASTER_PATH=/u/home/r/ruthjohn/ruthjohn/bayesM_noLD/unity_mixture
MASTER_PATH=/Users/ruthiejohnson/Development/unity_mixture
SCRIPT_DIR=${MASTER_PATH}/scripts
SRC_DIR=${MASTER_PATH}/src
DATA_DIR=${MASTER_PATH}/data

# simulation params
# simulation params
SIM_NAME=test_identity
P_VEC="0.01,.01,0.05,.05,.76,.05,.05,0.01,0.01"
BINS=9
LDSC_H2=.50
MU_VEC="-.10,-0.05,0,.05,.10"
SIGMA_VEC=".000001,.000001,.0000001,.000001,.000001"
M=100
LD_FILE=${DATA_DIR}/ukbb.${M}.ld
COEF=0.9
#LD_FILE=${DATA_DIR}/simulated_${COEF}.${M}.txt
N=100000
SEED=$SEED # can replace with SGE_TASK_ID
ITS=200

DATE=`date '+%Y-%m-%d %H:%M:%S'`
echo $DATE" Starting simulation for unity-mixture: "${SIM_NAME}

# global paths

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
	#python ${SCRIPT_DIR}/simulate.py --name $SIM_NAME --p_vec $P_VEC --mu_vec $MU_VEC --sigma_vec $SIGMA_VEC --ld_file $LD_FILE --M $M --N $N --seed $SEED --outdir $DATA_DIR --sigma_g $LDSC_H2
	python ${SCRIPT_DIR}/simulate.py --name $SIM_NAME --p_vec $P_VEC --bins $BINS --ld_file $LD_FILE --M $M --N $N --seed $SEED --outdir $DATA_DIR --sigma_g $LDSC_H2
fi

# STEP 2: transform betas
GWAS_FILE=${DATA_DIR}/${SIM_NAME}.${SEED}.txt
if [[ "$STEPS" =~ "2" ]]
then
	DATE=`date '+%Y-%m-%d %H:%M:%S'`
	echo $DATE" Transforming GWAS effect sizes"
	python ${SCRIPT_DIR}/transform_betas.py --gwas_file $GWAS_FILE --ld_file $LD_FILE
fi


# STEP 3: take 1/2 power of LD
if [[ "$STEPS" =~ "3" ]]
then
	python ${SCRIPT_DIR}/half_ld.py --ld_file $LD_FILE
fi

# STEP 7: Calculate true log-likelihood
if [[ "$STEPS" =~ "8" ]]
then
	LD_HALF_FILE=${LD_FILE%.*}.half_ld
	python ${SCRIPT_DIR}/true_likelihood.py --name $SIM_NAME --gwas_file $GWAS_FILE --ld_half_file $LD_HALF_FILE --N $N --seed $SEED --outdir $DATA_DIR --ldsc_h2 $LDSC_H2
fi


# STEP 4: run inference (with LD)
if [[ "$STEPS" =~ "4" ]]
then
	LD_HALF_FILE=${LD_FILE%.*}.half_ld
	#python ${SRC_DIR}/mixture_gibbs.py --name $SIM_NAME --gwas_file $GWAS_FILE --mu_vec $MU_VEC --sigma_vec $SIGMA_VEC --ld_half_file $LD_HALF_FILE --N $N --seed $SEED --outdir $DATA_DIR --precompute 'y' --its $ITS --ldsc_h2 $LDSC_H2
	python ${SRC_DIR}/mixture_gibbs.py --name $SIM_NAME --gwas_file $GWAS_FILE --bins $BINS --ld_half_file $LD_HALF_FILE --N $N --seed $SEED --outdir $DATA_DIR --precompute 'y' --its $ITS --ldsc_h2 $LDSC_H2

fi


# STEP 5: run inference (with NO LD)
if [[ "$STEPS" =~ "5" ]]
then
	python ${SRC_DIR}/mixture_em_noLD.py --name $SIM_NAME --gwas_file $GWAS_FILE --bins $BINS --N $N --seed $SEED --outdir $DATA_DIR  --its $ITS --ldsc_h2 $SIGMA_G
fi


# Step 6: plot gwas effect size histogram
if [[ "$STEPS" =~ "6" ]]
then
	python ${SCRIPT_DIR}/plot_histogram.py --name $SIM_NAME --gwas_file $GWAS_FILE --outdir ${DATA_DIR}
fi


# Step 7: plot binned effect size histogram
if [[ "$STEPS" =~ "7" ]]
then
	RESULTS_FILE=${DATA_DIR}/${SIM_NAME}.${SEED}.results
	python ${SCRIPT_DIR}/plot_EM_histogram.py --name $SIM_NAME --results_file $RESULTS_FILE --outdir ${DATA_DIR}
fi
