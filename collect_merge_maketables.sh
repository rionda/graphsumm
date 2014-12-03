#! /bin/bash

if [ $# -ne 2 ]; then
  echo "USAGE: $0 {2|3} dataset" > /dev/stderr
  exit 1
fi

L2_ALL_STRING="in_file,k,err_type,approx_alg,k,dims,size_avg,size_stddev,size_max,size_min,int_dens_avg,int_dens_stddev,int_dens_max,int_dens_min,ext_dens_avg,ext_dens_stddev,ext_dens_max,ext_dens_min,L2_err"
LC_ALL_STRING="in_file,k,err_type,approx_alg,k,dims,size_avg,size_stddev,size_max,size_min,int_dens_avg,int_dens_stddev,int_dens_max,int_dens_min,ext_dens_avg,ext_dens_stddev,ext_dens_max,ext_dens_min,cut-norm_err"
GRASS_ALL_STRING="in_file,k,err_type,k,size_avg,size_stddev,size_max,size_min,int_dens_avg,int_dens_stddev,int_dens_max,int_dens_min,ext_dens_avg,ext_dens_stddev,ext_dens_max,ext_dens_min,L2_err"
L2_SIZEDENS_STRING="in_file,k,approx_alg,dims,size_avg,size_stddev,size_max,size_min,int_dens_avg,int_dens_stddev,int_dens_max,int_dens_min,ext_dens_avg,ext_dens_stddev,ext_dens_max,ext_dens_min"
LC_SIZEDENS_STRING="in_file,k,size_avg,size_stddev,size_max,size_min,int_dens_avg,int_dens_stddev,int_dens_max,int_dens_min,ext_dens_avg,ext_dens_stddev,ext_dens_max,ext_dens_min"
GRASS_SIZEDENS_STRING=${LC_SIZEDENS_STRING}
L2_SIZE_STRING="in_file,k,approx_alg,dims,size_avg,size_stddev,size_max,size_min"
LC_SIZE_STRING="in_file,k,size_avg,size_stddev,size_max,size_min"
GRASS_SIZE_STRING=${LC_SIZE_STRING}
L2_DENS_STRING="in_file,k,approx_alg,dims,int_dens_avg,int_dens_stddev,int_dens_max,int_dens_min,ext_dens_avg,ext_dens_stddev,ext_dens_max,ext_dens_min"
LC_DENS_STRING="in_file,k,int_dens_avg,int_dens_stddev,int_dens_max,int_dens_min,ext_dens_avg,ext_dens_stddev,ext_dens_max,ext_dens_min"
GRASS_DENS_STRING=${LC_DENS_STRING}
L2_ERR_STRING="in_file,k,approx_alg,dims,L2_err"
LC_ERR_STRING="in_file,k,cut-norm_err"
GRASS_ERR_STRING="in_file,k,L2_err"
L2_RUNTIME_STRING="in_file,k,approx_alg,dims,total_time"
LC_RUNTIME_STRING="in_file,k,total_time"
GRASS_RUNTIME_STRING=${LC_RUNTIME_STRING}
L2_QUERIES_STRING="in_file,k,approx_alg,dims,triangles,adj_err_avg,adj_err_stdev,adj_err_max,adj_err_min,reldeg_err_avg,reldeg_err_stdev,reldeg_err_max,reldeg_err_min,deg_err_avg,deg_err_stdev,deg_err_max,deg_err_min,absdeg_err_avg,absdeg_err_max,absdeg_err_min"
LC_QUERIES_STRING="in_file,k,dims,triangles,adj_err_avg,adj_err_stdev,adj_err_max,adj_err_min,reldeg_err_avg,reldeg_err_stdev,reldeg_err_max,reldeg_err_min,deg_err_avg,deg_err_stdev,deg_err_max,deg_err_min,absdeg_err_avg,absdeg_err_max,absdeg_err_min"
GRASS_QUERIES_STRING=${LC_QUERIES_STRING}

declare -a suffixes
declare -a stats_strings

suffixes=("all" "sizedens" "size" "dens" "err" "time" "queries")

if [ $1 = 2 ]; then
  FLAGS="-a -s -y -p"
  SUFFIX_PREFIX="L2"
  stats_strings=("${L2_ALL_STRING}" "${L2_SIZEDENS_STRING}" "${L2_SIZE_STRING}" "${L2_DENS_STRING}" "${L2_ERR_STRING}" "${L2_RUNTIME_STRING}" "${L2_QUERIES_STRING}")
  VAR_FILE="$2_variables.sh"
  ERR=2
elif [ $1 = y2 ]; then
  FLAGS="-a -s -y"
  SUFFIX_PREFIX="L2"
  stats_strings=("${L2_ALL_STRING}" "${L2_SIZEDENS_STRING}" "${L2_SIZE_STRING}" "${L2_DENS_STRING}" "${L2_ERR_STRING}" "${L2_RUNTIME_STRING}" "${L2_QUERIES_STRING}")
  VAR_FILE="$2_variables.sh"
  ERR=2
elif [ $1 = 3 ]; then
  FLAGS="-a -s -y"
  SUFFIX_PREFIX="L3"
  stats_strings=("${LC_ALL_STRING}" "${LC_SIZEDENS_STRING}" "${LC_SIZE_STRING}" "${LC_DENS_STRING}" "${LC_ERR_STRING}" "${LC_RUNTIME_STRING}" "${LC_QUERIES_STRING}")
  VAR_FILE="$2_variables_3.sh"
  ERR=3
elif [ $1 = g2 ]; then
  FLAGS="-a -s -p"
  SUFFIX_PREFIX="L2"
  stats_strings=("${L2_ALL_STRING}" "${L2_SIZEDENS_STRING}" "${L2_SIZE_STRING}" "${L2_DENS_STRING}" "${L2_ERR_STRING}" "${L2_RUNTIME_STRING}" "${L2_QUERIES_STRING}")
  VAR_FILE="$2_variables.sh"
  ERR=2
elif [ $1 = gg ]; then
  FLAGS="-a -g"
  SUFFIX_PREFIX="GRASS"
  stats_strings=("${GRASS_ALL_STRING}" "${GRASS_SIZEDENS_STRING}" "${GRASS_SIZE_STRING}" "${GRASS_DENS_STRING}" "${GRASS_ERR_STRING}" "${GRASS_RUNTIME_STRING}" "${GRASS_QUERIES_STRING}")
  VAR_FILE="$2_variables.sh"
  ERR=2
else
  echo "ERROR: Unrecognized error type '$1'. Must be in {2,3,y2,g2,gg}." > /dev/stderr
  exit 1
fi

if [ ! -r ${VAR_FILE} ]; then
  echo "ERROR: Variable file ${VAR_FILE} not readable" > /dev/stderr
  exit 1
fi

. ${VAR_FILE}

MAX=`expr ${#suffixes[*]} - 1`
for i in $(seq 0 ${MAX}); do
  echo ${suffixes[i]}
  sh runcollect.sh ${FLAGS} collect ${stats_strings[i]} ${SUFFIX_PREFIX}${suffixes[i]} ${VAR_FILE}
  sh merge_results.sh res/collect/${DATASET}_*_${ERR}_*_${SUFFIX_PREFIX}${suffixes[i]}.csv > res/merge/${DATASET}_summ_${SUFFIX_PREFIX}${suffixes[i]}.csv
  sh csv_to_latextable.sh res/merge/${DATASET}_summ_${SUFFIX_PREFIX}${suffixes[i]}.csv > res/tables/${DATASET}_summ_${SUFFIX_PREFIX}${suffixes[i]}.tex
done

