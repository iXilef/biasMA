setwd("C:/Users/iXilef/OneDrive/Master/initial")
load("./data/ks_runs_final.RData")
load("./data/MWU_runs_final.RData")
library(topGO)

resultandseed <- function(result=NULL,seed=NULL)
{
  me <- list(
    classic = classic,
    elim = elim,
    weight=weight,
    lea=lea
  )
  
  # Set the name for the class
  class(me) <- append(class(me),"resultandseed")
  return(me)
}

orig_mwu <- get.data.parallel(data=mwu_runs_orig_data, what="MWU",sequ = seq(0, 0.35, 0.05),  iterations = 100,
                              annotfile = auxdata$original$annotation)

human_mwu <- get.data.parallel(data=mwu_runs_human_data, what="MWU",sequ = seq(0, 0.35, 0.05),  iterations = 100,
                               annotfile = auxdata$human$annotation)

rat_mwu <- get.data.parallel(data=mwu_runs_rat_data, what="MWU",sequ = seq(0, 0.35, 0.05),  iterations = 100,
                             annotfile = auxdata$rat$annotation)

# carefull changed get.data function for this

# orig_mwu_raw <- get.data.parallel(data=mwu_runs_orig_data, what="MWU",sequ = seq(0, 0.35, 0.05),  iterations = 100,
#                               annotfile = auxdata$original$annotation, correct="none")
# 
# human_mwu_raw <- get.data(data=mwu_runs_human_data, what="MWU",sequ = seq(0, 0.35, 0.05),  iterations = 100,
#                       annotfile = auxdata$human$annotation, correct="none")
# 
# rat_mwu_raw <- get.data.parallel(data=mwu_runs_rat_data, what="MWU",sequ = seq(0, 0.35, 0.05),  iterations = 100,
#                              annotfile = auxdata$rat$annotation, correct="none")


orig_results <- get.data.parallel(data=ks_runs_orig_data, sequ = seq(0, 0.35, 0.05), iterations = 100,
         annotfile = auxdata$original$annotation, GoDat=auxdata$original$GODat)

# orig_results_fdr <- get.data.parallel(data=ks_runs_orig_data, sequ = seq(0, 0.35, 0.05), iterations = 100,
#                           annotfile = auxdata$original$annotation, correct="fdr")

human_results <- get.data.parallel(data=ks_runs_human_data, sequ = seq(0, 0.35, 0.05), iterations = 100,
                          annotfile = auxdata$human$annotation, GoDat=auxdata$human$GOdat)

# human_results_fdr <- get.data.parallel(data=ks_runs_human_data, sequ = seq(0, 0.35, 0.05), iterations = 100,
#                           annotfile = auxdata$human$annotation, correct="fdr")

rat_results <- get.data.parallel(data=ks_runs_rat_data, sequ = seq(0, 0.35, 0.05), iterations = 100,
                          annotfile = auxdata$rat$annotation, GoDat=auxdata$rat$GODat)

# rat_results_fdr <- get.data.parallel(data=ks_runs_rat_data, sequ = seq(0, 0.35, 0.05), iterations = 100,
#                           annotfile = auxdata$rat$annotation, correct="fdr")



save(orig_results, orig_results_fdr, human_results, human_results_fdr, rat_results, rat_results_fdr, orig_mwu, human_mwu, rat_mwu, file="./data/results_final.RData")


save(orig_results, human_results, rat_results, file="./data/results_final_KS.RData")

save(orig_mwu, rat_mwu, human_mwu, file="./data/results_final_mwu.RData")

