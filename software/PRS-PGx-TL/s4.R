library(data.table)
library(MASS)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)
options(stringsAsFactors=F)


initial = args[1] #'PRS' or 'zero'
best_criterion = args[2]
s3_output = args[3]
s4_output = args[4]


load(paste0(s3_output, '_1_testScoresum_initial_',initial,'_',best_criterion,'.Rdata'))
prs_g = test$score_g
prs_gt = test$score_gt

for (chr in 2:22){
  load(paste0(s3_output, '_',chr,'_testScoresum_initial_',initial,'_',best_criterion,'.Rdata'))
  prs_g = prs_g + test$score_g
  prs_gt = prs_gt + test$score_gt
  
}
score = as.data.frame(cbind(prs_g, prs_gt))
colnames(score)[1:2] = c('prs_g','prs_gt')
score$ID = rownames(score)
#scale the prs
score[, 1:2] <- scale(score[, 1:2])
write.table(score, paste0(s4_output, '_allchr_testScoresum_initial_',initial,'_',best_criterion,'.txt'), quote=F, row.names=F, col.names=T, sep='\t')
