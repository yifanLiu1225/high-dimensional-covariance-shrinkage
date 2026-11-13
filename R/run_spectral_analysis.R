#!/usr/bin/env Rscript
# 在根目录运行 Rscript R/03_run_member_C.R --reps=500
args <- commandArgs(trailingOnly=TRUE)
getarg <- function(key,default) {
  x<-args[startsWith(args,paste0('--',key,'='))]
  if(length(x)) sub(paste0('^--',key,'='),'',x[1]) else default
}
reps <- as.integer(getarg('reps','500')); seed<-as.integer(getarg('seed','20261103'))
out<-getarg('out','.')
source('R/03_spectral_simulation.R',encoding='UTF-8')
source('R/03_C_figures.R',encoding='UTF-8')
raw<-run_spectral_C(reps,seed,out)
summarise_C(raw,out); plot_C(raw,out)
message('Simulation, summaries and figures finished.')
