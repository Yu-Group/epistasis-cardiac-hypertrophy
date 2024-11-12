#!/bin/bash

# Fit siRF and get importances for prioritization
thrs=("_binary_thr0.15" "_binary_thr0.2" "_binary_thr0.25")
nsnps=1000
for thr in "${thrs[@]}"; do
	Rscript sirf_local_stability_no_hypertension.R "iLVM"${thr} ${nsnps} >& "sirf_local_stability_no_hypertension_iLVM"${thr}".out"
done