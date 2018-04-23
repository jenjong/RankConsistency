rm(list = ls())
setwd("~/task/RankConsistency")

rdata <- readLines(con = "real_cv.R", n= -1)
i = 1
for (i in 1:20)
{
  tmp <- rdata  
  tmp[[16]] <- paste0( "file_idx = ", i, sep='')
  cat(tmp, file = paste("./job1/task_", i,'.R', sep=''), sep = '\n')
}
