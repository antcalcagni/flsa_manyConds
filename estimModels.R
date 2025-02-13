rm(list=ls())

setwd("/home/antonio/MEGA/Lavoro_sync/Conferences/2025_IES Bressanone/contribution_2/Rcode/")
source("FLSA.R"); source("utils.R")
library(Rdimtools)

if(!dir.exists("cache")){dir.create("cache")}; if(!dir.exists("results")){dir.create("results")}

K <- 10 #numer of CV splits
nocs <- seq(from=2,to=32)

cmds <- sapply(nocs,function(u){
  cmd <- paste0(
    "R CMD BATCH --no-restore --no-save '--args noc=", u, 
    " K=", K,
    "' estimModels_EXEC.R ", 
    " cache/out_noc_",u,".stdout"
  )
})

ncores <- 10
btcs <- split(cmds, ceiling(seq_along(cmds)/ncores)) 
for(jobs in btcs){ 
  cat("\n Running:");for(x in jobs){cat("\n\t",strsplit(x,split = "--args")[[1]][2])}
  
  batch_file <- tempfile("job_list_"); writeLines(jobs, batch_file)
  system(paste0("parallel -j ", ncores, " < ", batch_file)); unlink(batch_file)
}























