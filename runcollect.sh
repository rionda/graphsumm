# This script contains common variables for the grid scripts which
# need to be sourced on before the dataset specific variables.
# It also contains the code to handle command line options

usage() {
  echo "Usage: $0 -[abdgnprsy]" > /dev/stderr
  echo "\t-a: approx computation of error" > /dev/stderr
  echo "\t-b: mini-batch algorithm for summ" > /dev/stderr
  echo "\t-d: dry run" > /dev/stderr
  echo "\t-g: do grass" > /dev/stderr
  echo "\t-n: random alg for summ" > /dev/stderr
  echo "\t-p: approx algorithm for summ" > /dev/stderr
  echo "\t-r: do random" > /dev/stderr
  echo "\t-s: do summ" > /dev/stderr
  echo "\t-y: greedy algorithm for summ" > /dev/stderr
}

#BASE_DIR="/home/matteo/myres/yahoo/graphsumm/exper"
BASE_DIR="/home/matteo/yahoo/graphsumm/exper"
DATA_DIR="${BASE_DIR}/data"
SCRIPT_DIR="${BASE_DIR}/gridscripts"
RES_DIR="${BASE_DIR}/res"
TOP_SCRIPT_GRID="${BASE_DIR}/top_grid_script.sh" # contains the "header" for the base scripts
#PYTHON="/home/matteo/bin/python3"
PYTHON="python3"
QUERY_ERROR_SCRIPT="${BASE_DIR}/summary_errors.py"
STATISTICS_SCRIPT="${BASE_DIR}/get_summary_statistics.py"
COLLECT_SCRIPT="${BASE_DIR}/collect_stats.py"

USE_APPROX_ERR="false" # if true, use approx computation of error. Activate with -a
DO_MINI="false"  # if true, activate mini-batch alg for summ. Activate with -b
DRY_RUN="false"  # if true, don't run qsub, but just cat the scripts. Activate with -d
DO_GRASS="false" # if true, create and run scripts for grass. Activate with -g
DO_SUMM_RAND="false" # if true, activate random alg for summ. Activate with -n
DO_APPROX="false" # if true, activate approx alg for summ. Activate with -p
DO_RAND="false"  # if true, create and run scripts for rand. Activate with -r
DO_SUMM="false"  # if true, create and run scripts for summ. Activate with -s
DO_GREEDY="false" # if true, activate greedy alg for summ. Activate with -y

while getopts ":abdghnprsy" opt; do
    case ${opt} in
      a)
	USE_APPROX_ERR="true"
	;;
      b)
	DO_MINI="true"
	;;
      d)
	DRY_RUN="true"
	;;
      g)
	DO_GRASS="true"
	;;
      h)
	usage
	exit 0
	;;
      n)
	DO_SUMM_RAND="true"
	;;
      p)
	DO_APPROX="true"
	;;
      r)
	DO_RAND="true"
	;;
      s)
	DO_SUMM="true"
	;;
      y) 
	DO_GREEDY="true"
	;;
      \?) 
	echo "Invalid option: -$OPTARG" >&2
	usage
	exit 1
	;;
    esac
done

shift $(( OPTIND - 1 ));

if [ $# -lt 2 ]; then
  echo "Wrong number of parameters ({grid|sequential|collect {stats_csv}} datasetname)" > /dev/stderr
  exit 1
fi

case $1 in
  grid)
    COMMAND="create_and_run_script_grid"
    GRASS="\${GRASS_TMP}"
    SUMM="\${SUMM_TMP}"
    RAND="\${RAND_TMP}"
    ;;
  collect)
    COMMAND="collect"
    if [ $# -lt 4 ]; then
      echo "Wrong number of parameters ({grid|sequential|collect {stats_csv suffix}} variablesfile)" > /dev/stderr
      exit 1
    fi
    COLLECT_ARGS=$2
    SUFFIX=$3
    ;;
  sequential)
    COMMAND="create_and_run_script_sequential"
    GRASS="${BASE_DIR}/grass"
    SUMM="${BASE_DIR}/summ"
    RAND="${BASE_DIR}/random"
    ;;
  *)
    echo "Wrong command. Must be 'grid', 'sequential', or 'collect'" > /dev/stderr
    exit 1
    ;;
esac

# Go to the last parameter. "for" iterates over arguments if you don't
# tell it what to do
for SOURCE_FILE; do :; done

if [ ! -r ${SOURCE_FILE} ]; then
  echo "Variable file ${SOURCE_FILE} not readable" > /dev/stderr
  exit 1
fi

. ${BASE_DIR}/${SOURCE_FILE}
#echo $DATASET
#echo $INPUT
#echo $ERRS
#echo $KS
#echo $CS
#echo $DIMS
#echo $QUEUE
#echo $REPS

if [ ! -r ${DATA_DIR}/${INPUT} ]; then
  echo Input file ${DATA_DIR}/${INPUT} not readable > /dev/stderr
  exit 1
fi

# 1st argument is useless, second is result file base name
collect()
{
  if [ ${DRY_RUN:-false} != "true" ]; then
    ${PYTHON} ${COLLECT_SCRIPT} "${COLLECT_ARGS}" ${RES_DIR}/$2.*.err > ${RES_DIR}/collect/$2_${SUFFIX}.csv
  else
    echo "${PYTHON} ${COLLECT_SCRIPT} \"${COLLECT_ARGS}\" ${RES_DIR}/$2.*.err > ${RES_DIR}/collect/$2_${SUFFIX}.csv"
  fi
}

# 1st argument is command line (command + arguments), 2nd is script base name
create_and_run_script_sequential()
{
  SCRIPT_NAME="$2.sh"
  cat > ${SCRIPT_DIR}/${SCRIPT_NAME} << EOF
#! /bin/sh
if [ -n "\${REP}" ]; then
  REP="empty"
fi
$1 ${DATA_DIR}/${INPUT} > ${RES_DIR}/$2.\${REP}.out 2> ${RES_DIR}/$2.\${REP}.err
EOF
  for REP in $(seq 1 $REPS); do
    if [ ${DRY_RUN:-false} != "true" ]; then 
      echo "$1 ${INPUT}"
      $1 ${DATA_DIR}/${INPUT} > ${RES_DIR}/$2.${REP}.out 2> ${RES_DIR}/$2.${REP}.err
    else
      echo "$1 ${DATA_DIR}/${INPUT} > ${RES_DIR}/$2.${REP}.out 2> ${RES_DIR}/$2.${REP}.err"
    fi
  done
}

# 1st argument is command line (command + arguments), 2nd is script base name
#${BASE_DIR}/$1 ${DATA_DIR}/${INPUT} > >(tee ${RES_DIR}/$2.\${TASK_ID}.out) 2> >(tee ${RES_DIR}/$2.\${TASK_ID}.err > /dev/stderr)
create_and_run_script_grid()
{
  SCRIPT_NAME="$2.sh"
  cp ${TOP_SCRIPT_GRID} ${SCRIPT_DIR}/${SCRIPT_NAME}
  cat >> ${SCRIPT_DIR}/${SCRIPT_NAME} << EOF
$1 ${DATA_DIR}/${INPUT} > >(tee ${RES_DIR}/$2.\${TASK_ID}.out) 2> >(tee ${RES_DIR}/$2.\${TASK_ID}.err > /dev/stderr)
EOF
  if [ ${DRY_RUN:-false} != "true" ]; then 
    qsub -t 1-${REPS} -l ${QUEUE} ${SCRIPT_DIR}/${SCRIPT_NAME}
  else
    cat ${SCRIPT_DIR}/${SCRIPT_NAME}
    echo "qsub -t 1-${REPS} -l ${QUEUE} ${SCRIPT_DIR}/${SCRIPT_NAME}"
  fi
}

if [ ${USE_APPROX_ERR:-false} = "true" ]; then
  ERR_APPROX_FLAG="-a"
  ERR_APPROX_FLAG_STR="eapp"
else
  ERR_APPROX_FLAG=""
  ERR_APPROX_FLAG_STR="neap"
fi

for K in `echo ${KS}`; do
  # create grid script for ${RAND}
  if [ ${DO_RAND:-false} = "true" ]; then
    BASE_NAME="${DATASET}_rand_${K}_${ERR_APPROX_FLAG_STR}"
    ${COMMAND} "${RAND} -k ${K} ${ERR_APPROX_FLAG}" ${BASE_NAME}
  fi
  for ERR in `echo ${ERRS}`; do
    # create grid script for ${SUMM}
    if [ ${DO_SUMM:-false} = "true" ]; then
      for DIM in `echo ${DIMS}`; do
	if [ ${DO_MINI:-false} = "true" ]; then
	  ALG_FLAG="-b"
	  ALG_FLAG_STR="mini"
	  BASE_NAME="${DATASET}_summ_${K}_${ERR}_${DIM}_${ALG_FLAG_STR}_${ERR_APPROX_FLAG_STR}"
	  ${COMMAND} "${SUMM} -t ${ERR} -k ${K} -d ${DIM} ${ERR_APPROX_FLAG} ${ALG_FLAG} ${ADDITIONAL_FLAGS}" ${BASE_NAME}
	fi
	if [ ${DO_APPROX:-false} = "true" ]; then
	  ALG_FLAG="-p"
	  ALG_FLAG_STR="appr"
	  BASE_NAME="${DATASET}_summ_${K}_${ERR}_${DIM}_${ALG_FLAG_STR}_${ERR_APPROX_FLAG_STR}"
	  ${COMMAND} "${SUMM} -t ${ERR} -k ${K} -d ${DIM} ${ERR_APPROX_FLAG} ${ALG_FLAG} ${ADDITIONAL_FLAGS}" ${BASE_NAME}
	fi
	if [ ${DO_SUMM_RAND:-false} = "true" ]; then
	  ALG_FLAG="-r"
	  ALG_FLAG_STR="rand"
	  BASE_NAME="${DATASET}_summ_${K}_${ERR}_${DIM}_${ALG_FLAG_STR}_${ERR_APPROX_FLAG_STR}"
	  ${COMMAND} "${SUMM} -t ${ERR} -k ${K} -d ${DIM} ${ERR_APPROX_FLAG} ${ALG_FLAG} ${ADDITIONAL_FLAGS}" ${BASE_NAME}
	fi
	if [ ${DO_GREEDY:-false} = "true" ]; then
	  ALG_FLAG=""
	  ALG_FLAG_STR="gree"
	  BASE_NAME="${DATASET}_summ_${K}_${ERR}_${DIM}_${ALG_FLAG_STR}_${ERR_APPROX_FLAG_STR}"
	  ${COMMAND} "${SUMM} -t ${ERR} -k ${K} -d ${DIM} ${ERR_APPROX_FLAG} ${ALG_FLAG} ${ADDITIONAL_FLAGS}" ${BASE_NAME}
	fi
      done
    fi
    # create grid scripts for ${GRASS}
    if [ ${DO_GRASS:-false} = "true" ]; then
      for C in `echo ${CS}`; do
	BASE_NAME="${DATASET}_grass_${K}_${ERR}_${C}_${ERR_APPROX_FLAG_STR}"
	${COMMAND} "${GRASS} -t ${ERR} -k ${K} -c 0.${C} ${ERR_APPROX_FLAG}" ${BASE_NAME}
      done
    fi
  done
done

