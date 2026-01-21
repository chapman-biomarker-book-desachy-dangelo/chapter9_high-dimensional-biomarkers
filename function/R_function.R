optimal.grid.dis <- function(dat.validation, beta_grid, params){
  R2 <- double()
  for (i in 1:nrow(params)) {
    beta.est <- beta_grid[,i]
    if(sum(is.na(beta.est)) == 0){
      G.val <- dat.validation[,-1] %>% as.matrix()
      
      prs <- G.val %*% beta.est %>% as.vector()
      
      df <- cbind.data.frame(Y=dat.validation$Y, PRS=prs)
      fit <- lm(Y ~ PRS, data = df) %>% summary()
      
      R2 <- c(R2, fit$adj.r.squared)
    }
    if(sum(is.na(beta.est)) > 0){
      R2 <- c(R2, 0)
    }
  }
  names(R2) <- NULL
  
  df <- params[which.max(R2),]
  rownames(df) <- NULL
  
  df
}

optimal.grid.pgx <- function(dat.validation, beta_grid, params){
  R2 <- double()
  for (i in 1:nrow(params)) {
    beta.est <- beta_grid[,i]
    if(sum(is.na(beta.est)) == 0){
      G.val <- dat.validation[,-c(1,2)] %>% as.matrix()
      
      prs.g <- prs.gt <- G.val %*% beta.est %>% as.vector()
      
      df <- cbind.data.frame(Y=dat.validation$Y, Tr=dat.validation$Tr, PRSg=prs.g, PRSgt=prs.gt)
      fit <- lm(Y ~ Tr, data = df)
      resid <- fit$residuals
      df <- cbind.data.frame(df, resid=resid)
      fit <- lm(resid ~ PRSg + Tr:PRSgt, data = df) %>% summary()
      
      R2 <- c(R2, fit$adj.r.squared)
    }
    if(sum(is.na(beta.est)) > 0){
      R2 <- c(R2, 0)
    }
  }
  names(R2) <- NULL
  
  df <- params[which.max(R2),]
  rownames(df) <- NULL
  
  df
}

optimal.lambda <- function(dat.train, dat.validation, LAMBDA, method){
  R2 <- double()
  for (lambda in LAMBDA) {
    mod = PRS_PGx_Lasso(Y=dat.train$Y, Tr=dat.train$Tr, G=dat.train[,-c(1,2)], intercept = TRUE, lambda = lambda, method = method, alpha = 0.5)
    
    G.val <- dat.validation[,-c(1,2)] %>% as.matrix()
    
    prs.g <- G.val %*% mod$coef.G %>% as.vector()
    prs.gt <- G.val %*% mod$coef.G %>% as.vector()
    
    df <- cbind.data.frame(Y=dat.validation$Y, Tr=dat.validation$Tr, PRSg=prs.g, PRSgt=prs.gt)
    fit <- lm(Y ~ Tr, data = df)
    resid <- fit$residuals
    df <- cbind.data.frame(df, resid=resid)
    fit <- lm(resid ~ PRSg + Tr:PRSgt, data = df) %>% summary()
    
    R2 <- c(R2, fit$adj.r.squared)
  }
  names(R2) <- LAMBDA
  R2
}
