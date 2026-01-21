library(data.table)
library(MASS)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)
options(stringsAsFactors=F)


initial = args[1] #'PRS' or 'zero'
s3_output = args[2]
s4_output = args[3]


load(paste0(s3_output, '_1_testScoresum_initial_',initial,'.Rdata'))
prs_gt = test$score_gt

for (i in 2:22){
  load(load(paste0(s3_output, '_',chr,'_testScoresum_initial_',initial,'.Rdata')))
  prs_gt = prs_gt + test$score_gt
  
}
score = as.data.frame(prs_gt)
colnames(score)= c('prs_gt')
score$ID = rownames(score)
#scale the prs
score[, 1] <- scale(score[, 1])
write.table(score, paste0(s4_output, '_allchr_testScoresum_initial_',initial,'.txt'), quote=F, row.names=F, col.names=T, sep='\t')
