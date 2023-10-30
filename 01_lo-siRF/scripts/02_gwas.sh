#!/bin/bash

# Do dimension reduction via GWAS
phenotypes=("iLVM_norm" "LVM_norm")
for phenotype in "${phenotypes[@]}"; do
	Rscript 02a_gwas_bolt_lmm.R ${phenotype} >& "02a_gwas_bolt_lmm_"${phenotype}".out"
	for chrom in {1..22}; do
		Rscript 02b_gwas_plink.R ${phenotype} ${chrom} >& "02b_gwas_plink_"${phenotype}"_"${chrom}".out"
	done
done