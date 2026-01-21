require(mvtnorm)
require(lmtest)

# ##############################################################################
# GENERATE TRUE EFFECT SIZES
# ##############################################################################

generate.beta <- function(pcausal, nsnp, nblock, rhoY, rhoE, rhoC, h2, geno){
  corrmat <- matrix(c(1, rhoY, rhoY*rhoE, rhoC, rhoY*rhoC, rhoY*rhoC*rhoE,
                      rhoY, 1, rhoE, rhoY*rhoC, rhoC, rhoE*rhoC,
                      rhoY*rhoE, rhoE, 1, rhoY*rhoC*rhoE, rhoE*rhoC, rhoC,
                      rhoC, rhoY*rhoC, rhoY*rhoC*rhoE, 1, rhoY, rhoY*rhoE,
                      rhoY*rhoC, rhoC, rhoE*rhoC, rhoY, 1, rhoE,
                      rhoY*rhoC*rhoE, rhoE*rhoC, rhoC, rhoY*rhoE, rhoE, 1), ncol = 6, byrow = TRUE)
  
  coef <- double()
  sigma0 <- h2/(pcausal*nsnp)*corrmat
  mean0 <- rep(0,6)
  
  causal.range <- 51:100
  causal.index <- sort(sample(causal.range, nsnp*pcausal))
  causal.flag <- rep(0, nsnp); causal.flag[causal.index] <- 1
  
  for (i in 1:nsnp) {
    if(causal.flag[i] == 1){
      beta0 <- rmvnorm(1, mean = mean0, sigma = sigma0)
      beta0 <- as.vector(beta0)
    }
    if(causal.flag[i] == 0){
       beta0 <- rep(0,6)
    }
    coef <- rbind(coef, beta0)
  }
  coef <- as.data.frame(coef)
  colnames(coef) <- c("mu","b","a","phi","beta","alpha")
  rownames(coef) <- NULL
  coef <- cbind.data.frame(SNP=colnames(geno), coef)
  
  coef
}

# ##############################################################################
# GENERATE DATA
# ##############################################################################

generate.data <- function(geno, bim.geno, eff, se, h2.pgx, h2.dis){
  nsnp <- ncol(geno)
  nsubj <- nrow(geno)
  genox <- as.matrix(geno)
  
  # SIMULATE DISEASE GWAS SUMMARY STATISTICS
  mu.hat <- p.hat <- double()
  for (i in 1:nsnp) {
    x <- rnorm(1, mean = eff$mu[i], sd = se$SE[i]/2)
    mu.hat <- c(mu.hat, x)
    p.hat <- c(p.hat, 2*(1-pnorm(abs(x/se$SE[i]))))
  }
  dis.gwas.ss <- cbind.data.frame(SNP=colnames(geno), beta=mu.hat, se=se$SE, p = p.hat, N=se$N)
  dis.gwas.ss <- dis.gwas.ss %>% mutate(p=ifelse(p != 0, p, 1e-10)) %>% as.data.frame()
  dis.gwas.ss <- cbind.data.frame(CHR=bim.geno$CHR, POS=bim.geno$POS, RSID=bim.geno$RSID, A1=bim.geno$A1, A2=bim.geno$A2, dis.gwas.ss)
  
  # SIMULATE PGX DATA
  Tr <- rbinom(nsubj, 1, 0.5)
  sigma2 <- (1-h2.pgx)/h2.pgx*var(genox %*% eff$beta + (Tr*genox) %*% eff$alpha) - var(0.5*Tr)
  sigma <- sqrt(as.numeric(sigma2))
  e <- rnorm(nsubj, mean = 0, sd = sigma)
  y <- 0.5*Tr + genox %*% eff$beta + (Tr*genox) %*% eff$alpha + e
  y <- as.vector(as.numeric(y))
  dat.pgx <- cbind.data.frame(Y=y, Tr=Tr, geno)
  
  # SIMULATE TRAINING/VALIDATION/TESTING PGX DATA
  dat.train.pgx <- dat.pgx[1:3000,]
  dat.validation.pgx <- dat.pgx[3001:4000,]
  dat.test.pgx <- dat.pgx[4001:5000,]
  
  # SIMULATE PGX GWAS SUMMARY STATISTICS
  ss.pgx <- cal_ss(dat.train.pgx, eff)
  
  # SIMULATE DISEASE DATA
  sigma2 <- (1-h2.dis)/h2.dis*var(genox %*% eff$phi)
  sigma <- sqrt(as.numeric(sigma2))
  e <- rnorm(nsubj, mean = 0, sd = sigma)
  y <- genox %*% eff$phi + e
  y <- as.vector(as.numeric(y))
  dat.dis <- cbind.data.frame(Y=y, geno)
  
  # SIMULATE TRAINING/VALIDATION/TESTING PGX DATA
  dat.train.dis <- dat.dis[1:3000,]
  dat.validation.dis <- dat.dis[3001:4000,]
  dat.test.dis <- dat.dis[4001:5000,]
  
  re <- list(bim.geno = bim.geno, ss.dis = dis.gwas.ss, ss.pgx = ss.pgx, dat.train.pgx = dat.train.pgx, dat.validation.pgx = dat.validation.pgx, dat.test.pgx = dat.test.pgx, dat.train.dis = dat.train.dis, dat.validation.dis = dat.validation.dis, dat.test.dis = dat.test.dis, eff.ori=eff)
  re
}

cal_ss <- function(dat, eff){
  re <- double()
  for (k in 3:ncol(dat)) {
    df <- dat[,c(1,2,k)]
    colnames(df) <- c("Y","Tr","X")
    mod <- lm(Y ~ Tr*X, data = df)
    m <- summary(mod)
    modr <- lm(Y ~ Tr, data = df)
    fit <- waldtest(mod, modr)
    p0 <- fit$`Pr(>F)`[2]
    se1 <- m$coefficients[3,2]; se2 <- m$coefficients[4,2]
    p1 <- m$coefficients[3,4]; p2 <- m$coefficients[4,4]
    b1 <- rnorm(1, mean = eff$b[k-2], sd = se1/10)
    b2 <- rnorm(1, mean = eff$a[k-2], sd = se2/10)
    
    rex <- c(b1,se1,p1,b2,se2,p2,p0)
    re <- rbind(re, rex)
  }
  colnames(re) <- c("beta_est","beta_se","beta_p","alpha_est","alpha_se","alpha_p","2df_p")
  rownames(re) <- NULL
  re <- cbind.data.frame(SNP=colnames(dat)[3:ncol(dat)], re, N=nrow(dat))
  re
}









