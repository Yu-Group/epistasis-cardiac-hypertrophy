#!/bin/bash

# Fit siRF and get importances for prioritization
thrs=("_binary_thr0.15" "_binary_thr0.2" "_binary_thr0.25")
nsnps=1000
for thr in "${thrs[@]}"; do
	Rscript 04_irf_local_stability.R "iLVM"${thr} ${nsnps} >& "04_irf_local_stability_iLVM"${thr}".out"
done