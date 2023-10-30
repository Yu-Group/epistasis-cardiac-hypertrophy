import sys
import os
from os.path import join as oj
sys.path.append(oj(os.path.dirname(__file__), '..', 'functions'))
from fit_python_models import fitModels

reg_models = ["kernel_ridge"]
class_models = ["svm"]

pheno_name = sys.argv[1]
nsnps = int(sys.argv[2])
include_dems = bool(int(sys.argv[3]))
if "binary" in pheno_name:
  models = class_models
else:
  models = reg_models

out = fitModels(pheno_name, nsnps, include_dems, models)
