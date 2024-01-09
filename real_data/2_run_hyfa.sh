#!/bin/bash
cd ./sexde

while getopts "t:" option
do
    case "${option}"
    in
    t) tissue=${OPTARG};;
    esac
done

export PYTHONPATH=/lib/python3.8/site-packages/:$PYTHONPATH
export R_LIBS=/R/x86_64-pc-linux-gnu-library/3.5:$R_LIBS
out=./

mkdir -p ${out}/hyfa
mkdir -p ${out}/sexde/aux

# for tissue in $(cat ./data/queue/tissue.txt); do
# HYFA
/s/bin/python3.8 -m papermill ./HYFA/hyfa_impute.ipynb ${out}/hyfa/${tissue}.ipynb -p target_tissue ${tissue} -p out ${out}/hyfa

# Sex-differential gene expression
for name in lab_y lab_yhat unlab_yhat; do
    /s/pkg/linux64/R/3.5.1/bin/Rscript sexde.R \
    --data ./data \
    --outdir ${out}/sexde \
    --outname ${name} \
    --ge ${out}/hyfa/${tissue}.${name}.csv &
done

wait

