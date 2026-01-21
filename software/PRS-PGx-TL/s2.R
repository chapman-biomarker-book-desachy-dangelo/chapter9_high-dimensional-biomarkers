library(data.table)
library(MASS)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)
options(stringsAsFactors=F)


X = args[1]
initial = args[2]
s1_output = args[3]
s2_output = args[4]
allchr = args[5] #yes: allchr in one file, no: 22 files for each chr 

### default arguments (optional)
covar <- gsub(x = args[grep(x = args, pattern = "covar=")], pattern = "covar=", replacement = "")
if (length(covar) == 0) {
    covar <- 'Tr'
} else {
    covars <- unlist(strsplit(covar, ","))
    covar <- paste0(c('Tr', covars), collapse = " + ")
}
print(covar)

##########  read required input
#X:  y and covar for the full-training set
#     need columns: ID, Y, covars (including Tr), 
#     


X = as.data.frame(fread(paste0(X),header=T))


########## sum score

########### test
collect <- function(j){
    if (allchr == 'yes') {
        load(paste0(s1_output,'_initial_',initial,'_innerCV_',j,'.Rdata'))
    } else {
        load(paste0(s1_output,'_1_initial_',initial,'_innerCV_',j,'.Rdata'))
    }
    prs_g = result$score_g
    prs_gt = result$score_gt
    prs_g_train = result$score_g_train
    prs_gt_train = result$score_gt_train
 
    tunning = result$tunning
    tunning = rbind(rep(0,ncol(tunning)), tunning)
    
    #add the prs across all chromosomes
    if (allchr == 'no') {
        for (chr in 2:22){
            load(paste0(s1_output,'_',chr,'_initial_',initial,'_innerCV_',j,'.Rdata'))
            prs_g = prs_g + result$score_g
            prs_gt = prs_gt + result$score_gt
            prs_g_train = prs_g_train + result$score_g_train
            prs_gt_train = prs_gt_train + result$score_gt_train
            
        }
    }


    df = X[X$ID %in% rownames(prs_g),]
    prs_g = prs_g[df$ID, ]
    prs_gt = prs_gt[df$ID, ]
    prs_gtt = prs_gt * df$Tr

    df_train = X[X$ID %in% rownames(prs_g_train),]
    prs_g_train = prs_g_train[df_train$ID, ]
    prs_gt_train = prs_gt_train[df_train$ID, ]
    prs_gtt_train = prs_gt_train * df_train$Tr


    model0=lm(paste0('Y ~ ',covar), data=df)
    df$res = model0$residuals
    model0_train=lm(paste0('Y ~ ',covar), data=df_train)
    df_train$res = model0_train$residuals

    tunning$R2 = tunning$R2_g = tunning$R2_gt = tunning$R2_train = tunning$R2_g_train = tunning$R2_gt_train = tunning$condR2_gt = tunning$condR2_gt_train = NA
    for (i in 1:nrow(tunning)){
        df$prs_g = prs_g[,i]
        df$prs_gtt = prs_gtt[,i]
        model1=lm(df$res ~ prs_g[,i] + prs_gtt[,i])
        tunning[i,'R2'] = summary(model1)$r.squared
        model2=lm(df$res ~ prs_g[,i])
        tunning[i,'R2_g'] = summary(model2)$r.squared
        model3=lm(df$res ~ prs_gtt[,i])
        tunning[i,'R2_gt'] = summary(model3)$r.squared
        model_g=lm(paste0('Y ~ prs_g + ',covar), data=df)
        df$res_g = model_g$residuals
        model4=lm(res_g ~ prs_gtt, data=df)
        tunning[i,'condR2_gt'] = summary(model4)$r.squared



        df_train$prs_g = prs_g_train[,i]
        df_train$prs_gtt = prs_gtt_train[,i]
        model1=lm(df_train$res ~ prs_g_train[,i] + prs_gtt_train[,i])
        tunning[i,'R2_train'] = summary(model1)$r.squared
        model2=lm(df_train$res ~ prs_g_train[,i])
        tunning[i,'R2_g_train'] = summary(model2)$r.squared
        model3=lm(df_train$res ~ prs_gtt_train[,i])
        tunning[i,'R2_gt_train'] = summary(model3)$r.squared
        model_g=lm(paste0('Y ~ prs_g + ',covar), data=df_train)
        df_train$res_g = model_g$residuals
        model4=lm(res_g ~ prs_gtt, data=df_train)
        tunning[i,'condR2_gt_train'] = summary(model4)$r.squared





    }
    tunning = tunning[,c("lambda","factor","steps","R2_g_train","R2_gt_train","condR2_gt_train","R2_train","R2_g","R2_gt","condR2_gt","R2")]
    return(tunning)

}

R2_1 = collect(1)
R2_2 = collect(2)
R2_3 = collect(3)
R2_4 = collect(4)

#max overall R2
R2 = as.data.frame(cbind(R2_1[,'R2'], R2_2[,'R2'],R2_3[,'R2'],R2_4[,'R2']))
mean_r2 = apply(R2, 1, mean)
#max cond R2
condR2 = as.data.frame(cbind(R2_1[,'condR2_gt'], R2_2[,'condR2_gt'],R2_3[,'condR2_gt'],R2_4[,'condR2_gt']))
mean_cond_r2 = apply(condR2, 1, mean)

compare = as.data.frame(cbind(R2_1[,1:3],mean_r2,mean_cond_r2))
best1 = compare[which.max(compare$mean_r2),]
print(best1)
#   lambda factor steps    mean_r2
#211    100     10    16 0.05026869

best2 = compare[which.max(compare$mean_cond_r2),]
print(best2)
#save the parameter 
write.table(best1, paste0(s2_output, '_initial_',initial,'_best_R2.txt'), quote=F, row.names=F, col.names=T, sep='\t')
write.table(best2, paste0(s2_output, '_initial_',initial,'_best_condR2.txt'), quote=F, row.names=F, col.names=T, sep='\t')
