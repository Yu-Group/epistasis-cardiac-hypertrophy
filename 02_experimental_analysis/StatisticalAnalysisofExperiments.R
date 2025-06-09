library(quantreg)
library(tidyverse)
library(doParallel)
library(quantreg)
library(WRS2)
library(rsample)
#-------------------------------------------------------------------------------
# Function to calculate and return p-values comparing scramble vs knocked down
# cell size (or other cell measurements)
# df: dataframe with cell measurements for each experiment. Should only include
# "good" cells after imaging filters
# batch: experimental batch to test (ex. 'Batch#2')
# flow: flowrate from experiment (ex. 'High')
# outlet: can consider Top (larger cells) vs Bottom (smaller cells) vs All
# cell_line: specify health (e273) vs disease (r403q)
# gene: specify gene knocked down (ccdc_igf1r, ccdc_ttn, ccdc, ttn, igf1r)
# quantile: which quantile to test, range from 0.5-0.9 
# nboot: no. of bootstrap samples for percentile test, 10000 in results
#-------------------------------------------------------------------------------

ScramblevsKD <- function(df,batch,flow,outlet,cell_line,gene,quantile,nboot){
  #----------------------------------------------------------------
  # filter good cells based on parameters
  #----------------------------------------------------------------
  mod_df <- df[df$batch==batch &
               df$flow==flow &
               df$outlet==outlet &
               df$cell_line == cell_line &
               df$gene ==gene,]
  
  #----------------------------------------------------------------
  # Run traditional Wilcoxon test for medians
  #----------------------------------------------------------------
  wilcox_test <- wilcox.test(x = mod_df$cell_size[mod_df$scramble=="Scramble"],
                             y =  mod_df$cell_size[mod_df$scramble=="KD"])
  p_val_median <- wilcox_test$p.value
  p_val_median <- c(batch,cell_line,gene,"median-trad",0.5,p_val_median)
  
  #----------------------------------------------------------------
  # Run percentile bootstrap quantile test
  #----------------------------------------------------------------
  quant_boot <- qcomhd(cell_size~scramble,q = quantile, nboot = nboot,data = mod_df)
  p_val_quant <- quant_boot$partable$p.value
  p_val_quant <- c(batch,cell_line,gene,"quantile-bootperc",quantile,p_val_quant)
  
  #----------------------------------------------------------------
  # stack and return
  #----------------------------------------------------------------
  p_val_all <- rbind(p_val_median,p_val_quant)
  colnames(p_val_all) <- c('Batch','Cell Line',"Gene",'Test','Quantile','p-value')
  return(p_val_all)
}

#----------------------------------------------------------------
# example
#----------------------------------------------------------------
load('df.Rdata') # load data
batch <- 'Batch#2'
flow <- 'High'
outlet <- 'All'
cell_line <- 'e273'
gene <- 'ccdc_igf1r'
quantile <- 0.5
nboot <- 10000
SvKD_p_value <- ScramblevsKD(df,batch,flow,outlet,cell_line,gene,quantile,nboot)


#-------------------------------------------------------------------------------
# Function to perform quantile regression testing to test non-additivity while accounting for batch effects
# df: dataframe with cell measurements across batches in different experiments
# (merged based on KD efficiency in manuscript results)
# note: the batchflow column should specify corresponding batches and their flow rates (ex. 'Batch#2-High')
# so they can be controlled for
# cell_line: specify healthy vs disease cell line to test (e273 vs r403q)
# interaction: the two genes knocked down together (ccdc_igf1r vs ccdc_ttn)
# quantile: the quantile level for testing in the regression (0.5-0.9)
# nboot: no. of bootstrap samples for percentile bootstrap regression t-test (10,000 in results)
#-------------------------------------------------------------------------------
NonAdditivity <- function(df,cell_line,interaction,quantile,nboot){
  mod_df <- df[df$cell_line == cell_line,]
  # list knockdown gene measurements to include in regression based on interaction
  interactions <- list()
  interactions[['ccdc_igf1r']] <- c('ccdc','igf1r','ccdc_igf1r')
  interactions[['ccdc_ttn']] <- c('ccdc','ttn','ccdc_ttn')
  genes <- interactions[[interaction]]
  
  #----------------------------------------------------------------
  # Adjust format with dummies for regression
  # gene1 - 0/1 column where it's 1 if the first gene was KD and 0 otherwise
  # gene2 - 0/1 column where it's 1 if the second gene was KD and 0 otherwise
  # Note that measurements from the scramble case have gene1=gene2=0
  # and measurements from interactions have gene1=gene2=1
  #----------------------------------------------------------------
  mod_df$gene1 <- 0
  mod_df$gene2 <- 0
  mod_df$gene1[mod_df$gene %in% genes[-2]] <- 1
  mod_df$gene2[mod_df$gene %in% genes[-1]] <- 1
  
  #----------------------------------------------------------------
  # Run tranditional quantile regression and collect t-test results
  #----------------------------------------------------------------
  quant_reg <- rq(formula = cell_size ~ gene1*gene2+batchflow,tau=quantile,data = mod_df)
  quant_reg_summary <- summary(quant_reg)
  p_val_quant_reg <- quant_reg_summary$coefficients[5,4]
  p_val_quant_reg <- c(interaction,"quantilereg-t_test",quantile,p_val_quant_reg)
  
  #----------------------------------------------------------------
  # Run quantile regression on bootstrap samples and collect test statistics
  #----------------------------------------------------------------
  nworkers <- detectCores()
  res_mclapply <- mclapply(1:nboot,function(trial){
    g1 <- mod_df[mod_df$gene %in% genes[1],]
    gb1 <- as.data.frame(bootstraps(g1, batchflow, times = 1)$splits[[1]])
    g2 <- mod_df[mod_df$gene %in% genes[2],]
    gb2 <- as.data.frame(bootstraps(g2, batchflow, times = 1)$splits[[1]])
    g12 <- mod_df[mod_df$gene %in% genes[3],]
    gb12 <- as.data.frame(bootstraps(g12, batchflow, times = 1)$splits[[1]])
    sc <- mod_df[!(mod_df$gene %in% genes),]
    scb <- as.data.frame(bootstraps(sc, batchflow, times = 1)$splits[[1]])
    
    bootsamp_df <- rbind(
      gb1, gb2, gb12, scb)

    quant_reg_bootsamp <- rq(formula = cell_size ~ gene1*gene2+batchflow,tau=quantile,data = bootsamp_df)
    broom::tidy(quant_reg_bootsamp)
  },mc.cores=nworkers)
  res_mclapply <- dplyr::bind_rows(res_mclapply)
  
  #----------------------------------------------------------------
  # Compute percentile bootstrap p-values
  #----------------------------------------------------------------
  b.orig <- quant_reg_summary$coefficients[which(rownames(quant_reg_summary$coefficients) == 'gene1:gene2'),1]
  t.orig <- quant_reg_summary$coefficients[which(rownames(quant_reg_summary$coefficients) == 'gene1:gene2'),3]
  b.boot <- res_mclapply$estimate[res_mclapply$term=='gene1:gene2']
  stderr.boot <- res_mclapply$std.error[res_mclapply$term=='gene1:gene2']
  t.boot <- (b.boot-b.orig)/stderr.boot
  p.upper <- mean(t.boot >= t.orig) 
  p.lower <- mean(t.boot <= t.orig) 
  p_val_quant_reg_percboot <- 2 * min(p.upper, p.lower)
  p_val_quant_reg_percboot <- c(interaction,"quantilereg-percboot-t_test",quantile,p_val_quant_reg_percboot)

  #----------------------------------------------------------------
  # stack and return
  #----------------------------------------------------------------
  p_val_all <- rbind(p_val_quant_reg,p_val_quant_reg_percboot)
  colnames(p_val_all) <- c('Interaction','Test','Quantile','p-value')
  return(p_val_all)
}


#----------------------------------------------------------------
# example
#----------------------------------------------------------------
load('df.Rdata') # load data
interaction <- 'ccdc_igf1r'
quantile <- 0.5
nboot <- 10000
NonAdd_p_value <- NonAdditivity(df,interaction,quantile,nboot)
