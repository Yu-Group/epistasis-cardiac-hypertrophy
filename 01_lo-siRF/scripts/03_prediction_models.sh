#!/bin/bash

# Fit prediction models to assess prediction check
thrs=("" "_binary_thr0.15" "_binary_thr0.2" "_binary_thr0.25")
dems=(0 1)
nsnps=1000
for thr in "${thrs[@]}"; do
	for dem in "${dems[@]}"; do
		Rscript 03_prediction_models.R "iLVM"${thr} ${nsnps} ${dem} >& "03_prediction_models_iLVM"${thr}"_"${dem}".out"
		python3 03_prediction_models.py "iLVM"${thr} ${nsnps} ${dem}
	done
done