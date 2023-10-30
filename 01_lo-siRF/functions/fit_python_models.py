import pandas as pd
import numpy as np
import os

from sklearn.svm import SVC
from sklearn.kernel_ridge import KernelRidge
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from autosklearn.regression import AutoSklearnRegressor

def fitModels(pheno_name, nsnps, include_dems, models):
    # get correct path
    if include_dems:
        save_path = "../results/" + pheno_name + "_gwas_filtered_nsnps" + str(nsnps) + "_dems"
    else:
        save_path = "../results/" + pheno_name + "_gwas_filtered_nsnps" + str(nsnps)
    
    # load in data
    Xtr = pd.read_csv(save_path + "/Xtrain.csv")
    Xts = pd.read_csv(save_path + "/Xtest.csv")
    Ytr = pd.read_csv(save_path + "/ytrain.csv")
    Yts = pd.read_csv(save_path + "/ytest.csv")
    
    ypred_dict = {}
    for model in models:
        print("Fitting " + model + "...")
        if model == "kernel_ridge":
            param_grid = {"kernelridge__alpha" : np.logspace(-4, 0, 5),
                          "kernelridge__gamma" : np.logspace(-6, -2, 5)}
            grid_search = GridSearchCV(make_pipeline(StandardScaler(), KernelRidge(kernel="rbf")), 
                                       param_grid)
            grid_search.fit(Xtr, Ytr)
            krr_best_model = grid_search.best_estimator_
            print(krr_best_model)
            print(pd.DataFrame.from_dict(grid_search.cv_results_))
            fitted_model = krr_best_model
            
        elif model == "mlp2":
            param_grid = {"mlpregressor__hidden_layer_sizes" : np.arange(50, 301, 50)}
            grid_search = GridSearchCV(make_pipeline(StandardScaler(), MLPRegressor(max_iter=500)), 
                                       param_grid)
            grid_search.fit(Xtr, Ytr.x.ravel())
            mlp_best_model = grid_search.best_estimator_
            print(mlp_best_model)
            print(pd.DataFrame.from_dict(grid_search.cv_results_))
            fitted_model = mlp_best_model
            
        elif model == "autosklearn":
            automl = AutoSklearnRegressor()
            automl.fit(Xtr, Ytr)
            print(automl.cv_results_)
            fitted_model = automl
            
        elif model == "svm":
            param_grid = {"svc__C" : np.logspace(-4, 4, 9)}
            grid_search = GridSearchCV(make_pipeline(StandardScaler(), SVC(kernel="rbf", probability=True)), 
                                       param_grid)
            grid_search.fit(Xtr, np.ravel(Ytr))
            svm_best_model = grid_search.best_estimator_
            print(svm_best_model)
            print(pd.DataFrame.from_dict(grid_search.cv_results_))
            fitted_model = svm_best_model
        
        if model == "svm":
            ypred_dict[model] = fitted_model.predict_proba(Xts)[:, 1].tolist()
        else:
            ypred_dict[model] = fitted_model.predict(Xts).tolist()
            
    ypred = pd.DataFrame.from_dict(ypred_dict)
    ypred_df = pd.concat([ypred[col].explode(col) for col in ypred.columns], axis = 1)
    ypred_df.to_csv(save_path + "/ypred.csv", index = False)
    
    # remove temporary files
    os.remove(save_path + "/Xtrain.csv")
    os.remove(save_path + "/Xtest.csv")
    os.remove(save_path + "/ytrain.csv")
    os.remove(save_path + "/ytest.csv")
    
    return
