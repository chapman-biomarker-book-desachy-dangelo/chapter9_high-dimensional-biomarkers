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
output = args[7]

#sum_stats = "/home2/zhaiso/PROJECT/Book_PGxPRS/data/PRS-PGx-TL/sumstats.txt"
#ped = "/home2/zhaiso/PROJECT/Book_PGxPRS/data/PRS-PGx-TL/pheno.txt"
#G = "/home2/zhaiso/PROJECT/Book_PGxPRS/data/PRS-PGx-TL/geno.Rdata"
#bim = "/home2/zhaiso/PROJECT/Book_PGxPRS/data/PRS-PGx-TL/snp_info.txt"
#initial = "zero"
#num_snp = 500

### default arguments (optional)
covar <- "Tr"


#################### read required input
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
G = G[rownames(G) %in% ped$ID, ]

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
group_assignment <- rep(1:4, each = nrow(ped) %/% 4)
if (nrow(ped) %% 4 > 0){ group_assignment <- c(group_assignment, 1:(nrow(ped) %% 4))}
set.seed(123)  
# Randomly assign each individual to one of 4 groups
ped$group <- sample(group_assignment)
group_values <- as.numeric(names(table(ped$group)))
#valid_gp = group_values[1]

######### 2. adjust the phenotype 
adjust_y <- function(data, valid_gp){
    #data = ped
    #valid_gp=group_values[1]
    data = data[data$group != valid_gp, ]
    model0=lm(paste0('Y ~ ',covar), data=data)
    data$res = model0$residuals
    model1=lm(res ~ prs_g + prs_gt, data=data)

    y_res = data[,c('ID','res','Tr')]
    coeff = c(summary(model1)$coefficients['prs_g','Estimate'], summary(model1)$coefficients['prs_gt','Estimate'])

    return(list(y_res=y_res,coeff = coeff))

}

#obj = adjust_y(ped, valid_gp)
######### 3. Gradient descent and keep all beta
GD <- function(obj, G, bim_sum_stats, step = 15, nsnp){
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
    X = as.matrix(cbind(geno, GT))

    #check
    print(table(rownames(X) == y_res$ID))
    print(table(colnames(geno) == geno_info$SNP))
    
    GG=t(X)%*%X/(nrow(X)-1)
    gy=t(X)%*%(y_res$res)/(nrow(y_res)-1)
    if (nrow(geno_info)>0){
        #number of total snps with non-zero eff
        #nsnp = 24111 
        betatemp_g=geno_info$Beta2*geno_info$sd*obj$coeff[1]
        if (initial == 'PRS'){betatemp_gt=geno_info$Beta2*geno_info$sd_gt*obj$coeff[2]}
        if (initial == 'zero'){betatemp_gt=rep(0, length(betatemp_g))}
        betatemp = c(betatemp_g, betatemp_gt)
        u0=gy-GG%*%betatemp
        beta.all=cbind(u0, betatemp)

        tunning = data.frame(lambda = NA, factor =NA, steps = NA)
        for (lambda in c(0, 0.5, 0.99)){
            for (factor1 in c(1,10,50,70)){
                k=1
                betatemp=beta.all[,2]
                u0=beta.all[,1]
                while (k<=step){
                    learningrate=1/nsnp*factor1
                    if (learningrate>1){learningrate=1}
                    print(learningrate)
                    beta_old=betatemp
                    betatemp = learningrate*u0+(1-lambda)*beta_old
                    #betatemp = learningrate*u0+(1-learningrate*lambda)*beta_old
                    u0=u0-GG%*%(betatemp-beta_old)

                    beta.all=cbind(beta.all,betatemp)
                    k=k+1

                    tunning = as.data.frame(rbind(tunning, c(lambda, factor1, k)))
                } 
            }
        }

        oc = data.frame(SNP = rownames(beta.all), Ref_G = c(geno_info$Ref_G, geno_info$Ref_G),sd = c(geno_info$sd, geno_info$sd_gt), mean = c(geno_info$mean, geno_info$mean_gt))
        oc = as.data.frame(cbind(oc, beta.all))
        tunning = tunning[-1, ]

        #tunning = data.frame(factor = c(rep(1,step),rep(10,step),rep(50,step),rep(100,step)), steps = c(rep(1:step, 4)))
        return(list(oc=oc,tunning=tunning))

    }
}
#GD_obj = GD(obj, G, bim_sum_stats, 15)

######### 4. apply all the beta on the validation to get the scoresum (for this region)
get_score <- function(GD_obj, obj, G){
    #GD_obj = try
    oc = GD_obj$oc
    oc_g = oc[1:(nrow(oc)/2),]
    oc_gt = oc[(nrow(oc)/2+1):nrow(oc),]
    snps = rownames(oc_g)
    y_res = obj$y_res

    geno_valid = G[! rownames(G) %in% y_res$ID, snps]
    geno_train = G[rownames(G) %in% y_res$ID, snps]

    table(colnames(geno_valid) == rownames(oc_g))
    #get score sum for the validation set and the traning set
    geno_valid = as.matrix(geno_valid)
    geno_train = as.matrix(geno_train)
    weight_g = oc_g[,6]/oc_g$sd
    weight_gt = oc_gt[,6]/oc_gt$sd
    score_g = geno_valid %*% weight_g
    score_gt = geno_valid %*% weight_gt
    score_g_train = geno_train %*% weight_g
    score_gt_train = geno_train %*% weight_gt


    #get score sum for the training and validation set

     for (j in 7:ncol(oc_g)){
        weight_g = oc_g[,j]/oc_g$sd
        weight_gt = oc_gt[,j]/oc_gt$sd
        score_g = cbind(score_g, geno_valid %*% weight_g)
        score_gt = cbind(score_gt, geno_valid %*% weight_gt)
        score_g_train = cbind(score_g_train, geno_train %*% weight_g)
        score_gt_train = cbind(score_gt_train, geno_train %*% weight_gt)


     }
     return(list(score_g = score_g, score_gt=score_gt, score_g_train=score_g_train, score_gt_train=score_gt_train, tunning=GD_obj$tunning))

}

#result = get_score(try, obj, G)



######### run algorithm for each of the 4 folds

for (j in 1:length(group_values)){
    valid_gp = group_values[j]
    obj = adjust_y(ped, valid_gp)
    GD_obj = GD(obj, G, bim_sum_stats, 30, num_snp)
    result = get_score(GD_obj, obj, G)
    #save(result, file=paste0('/SFS/archive/data1/bardstms/shenjud/2024SummerIntern/real_data/TL_PGx_PRS/s1_output/penalized/lassosum_',CHR,'_outCV_',index,'_innerCV_',j,'.Rdata'))
    save(result, file=paste0(output,'_initial_',initial,'_innerCV_',j,'.Rdata'))
    print(j)
}

