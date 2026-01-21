library(data.table)
library(MASS)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)
options(stringsAsFactors=F)


sum_stats = args[1]
ped = args[2]
G = args[3]
bim = args[4]

initial = args[5] #'PRS' or 'zero'
num_snp = as.numeric(args[6]) #number of total snps with non-zero eff
s2_output = args[7]
s3_output = args[8]

### default arguments (optional)
covar <- gsub(x = args[grep(x = args, pattern = "covar=")], pattern = "covar=", replacement = "")
if (length(covar) == 0) {
  covar <- 'Tr'
} else {
  covars <- unlist(strsplit(covar, ","))
  covar <- paste0(c('Tr', covars), collapse = " + ")
}
print(covar)

#################### required input
#sum_stats: GWAS PRS weights 
#           need columns: chr, SNP, BP, Eff, Ref, Beta (need to be matched and flip)
#ped: y and covar for the full-training set
#     need columns: ID, Y, covars (including Tr), prs_g (from sum_stats)
#     
#G: genotype matrix ind x SNP 
#   need rownames as ID, colnames as SNP
#bim: allele info for G
#     need columns: SNP (match with sum_stats), Ref 

sum_stats = as.data.frame(fread(paste0(sum_stats),header=T))
ped = as.data.frame(fread(paste0(ped),header=T))
ped$prs_gt = ped$Tr * ped$prs_g
load(paste0(G))
bim = as.data.frame(fread(paste0(bim),header=T))
G2 = G #all samples
G = G2[rownames(G2) %in% ped$ID, ] #full-training samples

head(sum_stats)
head(ped)
G[1:10,1:10]
head(bim)

######### 0. match snp and flip +/-
bim2 = bim[,c('SNP','Ref')]
colnames(bim2)[2] = 'Ref_G'

bim_sum_stats=inner_join(sum_stats, bim2)
bim_sum_stats$Beta2=NA

#ref in prs matched with ref in training
flag1=which(bim_sum_stats$Ref==bim_sum_stats$Ref_G)
if (length(flag1)>0){  bim_sum_stats$Beta2[flag1]=bim_sum_stats$Beta[flag1]}
#eff in prs matched with ref in training
flag2=which(bim_sum_stats$Eff==bim_sum_stats$Ref_G)
if (length(flag2)>0){  bim_sum_stats$Beta2[flag2]=-bim_sum_stats$Beta[flag2]}
#table(bim_sum_stats$Beta2 == bim_sum_stats$Beta_matched)
bim_sum_stats=bim_sum_stats[is.na(bim_sum_stats$Beta2)==F,]


######### 1. split into 4 (3 for training, 1 for validation)
#skip this step 1

######### 2. adjust the phenotype 
adjust_y <- function(data){
    #data = ped
    #data = data[data$group != valid_gp, ]
    #adjust covariates and PRS_G!!!
    model0=lm(paste0('Y ~ prs_g + ',covar), data=data)
    data$res = model0$residuals
    model1=lm(res ~ prs_gt, data=data)

    y_res = data[,c('ID','res','Tr')]
    coeff = c(summary(model1)$coefficients['prs_gt','Estimate'])

    return(list(y_res=y_res,coeff = coeff))

}

#obj = adjust_y(ped, valid_gp)
######### 3. Gradient descent and keep all beta
GD <- function(obj, G, bim_sum_stats, step, nsnp, lambda, factor1){
    y_res = obj$y_res
    geno = G[y_res$ID, bim_sum_stats$SNP]
    geno_info=bim_sum_stats[,c('SNP','Ref_G','Beta2')]
    #impute NA as mean
    if (sum(is.na(geno))>0){
        for (j in 1:ncol(geno)){
            flag=which(is.na(geno[,j]))
            if (length(flag)>0){ geno[flag,j]=mean(geno[,j], na.rm=T)}
        }
    }
    geno_info$mean=colMeans(geno,na.rm=T)
    geno_info$sd=apply(geno, 2, sd)
    list1=which(geno_info$sd==0)
    if (length(list1)>0){
        geno_info=geno_info[-list1,]
        geno=geno[,-list1]
    }

    #add the geno x T (no NA in T)
    table(rownames(geno) == y_res$ID)
    GT = geno * y_res$Tr 
    geno_info$mean_gt = apply(GT, 2, mean)
    geno_info$sd_gt = apply(GT, 2, sd)
   
    table(colnames(GT) == colnames(G))
    colnames(GT) = paste0(colnames(GT), '_gt')

   
    #standardized 
    for (j in 1:ncol(geno)){
        geno[,j]=(geno[,j]-geno_info$mean[j])/geno_info$sd[j]
        GT[,j]=(GT[,j]-geno_info$mean_gt[j])/geno_info$sd_gt[j]
    }
    #X = [G GT]
    #X = as.matrix(cbind(geno, GT))
    X = as.matrix(GT)

    #check
    print(table(rownames(X) == y_res$ID))
    print(table(c(sub("_gt$", "", colnames(X))) == geno_info$SNP))
    
    GG=t(X)%*%X/(nrow(X)-1)
    gy=t(X)%*%(y_res$res)/(nrow(y_res)-1)
    if (nrow(geno_info)>0){
        #nsnp: number of total snps with non-zero eff
        
        if (initial == 'PRS'){betatemp_gt=geno_info$Beta2*geno_info$sd_gt*obj$coeff}
        if (initial == 'zero'){betatemp_gt=rep(0, nrow((geno_info)))}

        betatemp = betatemp_gt
        u0=gy-GG%*%betatemp
        beta.all=cbind(u0, betatemp)

        tunning = data.frame(lambda = NA, factor =NA, steps = NA)
        k=1
        betatemp=beta.all[,2]
        u0=beta.all[,1]
        while (k<=step){
            learningrate=1/nsnp*factor1
            if (learningrate>1){learningrate=1}
            print(learningrate)
            beta_old=betatemp
            betatemp = learningrate*u0+(1-lambda)*beta_old
            u0=u0-GG%*%(betatemp-beta_old)

            beta.all=cbind(beta.all,betatemp)
            k=k+1

            tunning = as.data.frame(rbind(tunning, c(lambda, factor1, k)))
        } 
  
        oc = data.frame(SNP = rownames(beta.all), Ref_G = c(geno_info$Ref_G),sd = c(geno_info$sd_gt), mean = c(geno_info$mean_gt))
        oc = as.data.frame(cbind(oc, beta.all))
        tunning = tunning[-1, ]

        #tunning = data.frame(factor = c(rep(1,step),rep(10,step),rep(50,step),rep(100,step)), steps = c(rep(1:step, 4)))
        return(list(oc=oc,tunning=tunning))

    }
}
#GD_obj = GD(obj, G, bim_sum_stats, 15)


######### 4. keep the best beta and save
best = read.table(paste0(s2_output, '_initial_',initial,'_best_R2.txt'),header=T)
best_lambda = best[1,'lambda']
best_factor = best[1,'factor']
best_step = best[1,'steps']
obj = adjust_y(ped)
GD_obj = GD(obj, G, bim_sum_stats, 200, num_snp, best_lambda, best_factor)

if (best_step == 0){
    beta = GD_obj$oc[,c(1:4, 6)]
} else {
    beta = GD_obj$oc[,c(1:4, 6+best_step-1)]
}
colnames(beta)[5] = 'best_coef'

oc_gt = beta

collect_beta = list(beta_gt=oc_gt)
save(collect_beta, file=paste0(s3_output, '_Beta_initial_',initial,'.Rdata'))


######### 5. apply the beta on the testing to get the scoresum (for this region) - optional
snps = sub("_gt$", "", rownames(oc_gt))
geno = G2[! rownames(G2) %in% ped$ID, snps] #keep testing samples not used in full-training set
if (nrow(geno) > 0) {
  #impute NA as mean
  if (sum(is.na(geno))>0){
    for (j in 1:ncol(geno)){
      flag=which(is.na(geno[,j]))
      if (length(flag)>0){ geno[flag,j]=mean(geno[,j], na.rm=T)}
    }
  }
  
  #get score sum for the testing set
  geno = as.matrix(geno)
  weight_gt = oc_gt$best_coef/oc_gt$sd
  score_gt = geno %*% weight_gt
  
  
  test = list(score_gt=score_gt)
  save(test, file=paste0(s3_output, '_testScoresum_initial_',initial,'.Rdata'))
  
  
}
