# Course-derived reference script
# English filename: 33_cibersort_xxdimmune_helper.R
# Original course path: see source_manifest.csv
# Role: deconvolution helper reference
# Adaptations applied for the skill package:
# - Filename flattened and translated to English.
# - Interactive working-directory selection replaced with SC_WORK_DIR/getwd() when present.
# - Example-specific object names, thresholds, metadata columns, and local filenames still require review.
# - Prefer scripts/course_adapted/ for runnable project workflows.
CoreAlg <- function(X, y){
  

  svn_itor <- 3
  
  res <- function(i){
    if(i==1){nus <- 0.25}
    if(i==2){nus <- 0.5}
    if(i==3){nus <- 0.75}
    model<-svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
    model
  }
  
  if(Sys.info()['sysname'] == 'Windows') out <- mclapply(1:svn_itor, res, mc.cores=1) else
    out <- mclapply(1:svn_itor, res, mc.cores=svn_itor)
  
  nusvm <- rep(0,svn_itor)
  corrv <- rep(0,svn_itor)
  

  t <- 1
  
  while(t <= svn_itor) {
    weights = t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights<0)]<-0
    w<-weights/sum(weights)
    u <- sweep(X,MARGIN=2,w,'*')
    k <- apply(u, 1, sum)
    nusvm[t] <- sqrt((mean((k - y)^2)))
    corrv[t] <- cor(k, y)
    t <- t + 1
  }
  

  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]]

  q <- t(model$coefs) %*% model$SV
  q[which(q<0)]<-0
  w <- (q/sum(q))
  
  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]
  
  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
  
}

doPerm <- function(perm, X, Y){
  message(paste0("进行doPerm循环,调用线程数： ",detectCores()))
  itor <- 1
  Ylist <- as.list(data.matrix(Y))
  dist<-list()
  registerDoParallel(round(detectCores()))
  dist<-foreach(itor = 1:perm) %dopar% {  

    library(limma)
    library(parallel)
    library(ggplot2)
    library(ggpubr)
    library(future.apply)
    library(car)
    library(ridge)
    library(e1071)
    library(preprocessCore)
    library(tcltk)
    library(limma)
    library(parallel)
    library(ggplot2)
    library(ggpubr)
    library(future.apply)
    library(car)
    library(ridge)
    library(e1071)
    library(preprocessCore)
    library(foreach)  
    library(doParallel)
    library(limma)
    library(parallel)
    library(ggplot2)
    library(ggpubr)
    library(future.apply)
    library(car)
    library(ridge)
    library(e1071)
    library(preprocessCore)
    library(tcltk)
    library(limma)
    library(parallel)
    library(ggplot2)
    library(ggpubr)
    library(future.apply)
    library(car)
    library(ridge)
    library(e1071)
    library(preprocessCore)
    library(foreach)  
    library(doParallel)
    library(data.table)
    CoreAlg <- function(X, y){

      svn_itor <- 3
      
      res <- function(i){
        if(i==1){nus <- 0.25}
        if(i==2){nus <- 0.5}
        if(i==3){nus <- 0.75}
        model<-svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
        model
      }
      
      if(Sys.info()['sysname'] == 'Windows') out <- mclapply(1:svn_itor, res, mc.cores=1) else
        out <- mclapply(1:svn_itor, res, mc.cores=svn_itor)
      
      nusvm <- rep(0,svn_itor)
      corrv <- rep(0,svn_itor)

      t <- 1
      
      while(t <= svn_itor) {
        weights = t(out[[t]]$coefs) %*% out[[t]]$SV
        weights[which(weights<0)]<-0
        w<-weights/sum(weights)
        u <- sweep(X,MARGIN=2,w,'*')
        k <- apply(u, 1, sum)
        nusvm[t] <- sqrt((mean((k - y)^2)))
        corrv[t] <- cor(k, y)
        t <- t + 1
      }

      rmses <- nusvm
      mn <- which.min(rmses)
      model <- out[[mn]]

      q <- t(model$coefs) %*% model$SV
      q[which(q<0)]<-0
      w <- (q/sum(q))
      
      mix_rmse <- rmses[mn]
      mix_r <- corrv[mn]
      
      newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
      
    }

    yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])

    yr <- (yr - mean(yr)) / sd(yr)
    
    result <- CoreAlg(X, yr)
    
    mix_r <- result$mix_r
    dist[[itor]]<-mix_r
    return(dist[[itor]])
  }
  stopImplicitCluster()
  dist<-do.call(rbind,dist)
  newList <- list("dist" = dist)
}

CIBERSORT <- function(sig_matrix, mixture_file, perm=0, QN=TRUE){
  library(e1071)
  library(parallel)
  library(preprocessCore)
  
  message("读取表达矩阵")
  X <- fread(sig_matrix,header=T,sep="\t",check.names=F,data.table=FALSE)
  rownames(X)=X[,1]
  X=X[,2:ncol(X)]
  Y <- fread(mixture_file, header=T, sep="\t",check.names=F,data.table=FALSE)
  rownames(Y)=Y[,1]
  Y=Y[,2:ncol(Y)]
  X <- na.omit(X)
  Y <- na.omit(Y)
  X <- data.matrix(X)
  Y <- data.matrix(Y)
  X<-as.matrix(X)
  Y<-as.matrix(Y)
  X <- X[order(rownames(X)),]
  Y <- Y[order(rownames(Y)),]
  
  P <- perm 

  if(max(Y) < 50) {Y <- 2^Y}
  
  if(QN == TRUE){
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }
  
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  Y <- Y[YintX,]
  XintY <- Xgns %in% row.names(Y)
  X <- X[XintY,]

  X <- (X - mean(X)) / sd(as.vector(X))
  
  if(P > 0) {nulldist <- sort(doPerm(P, X, Y)$dist)}

  header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")
  itor <- 1
  mixtures <- dim(Y)[2]
  pval <- 9999
  message(paste0("进行CIBERSORT循环,调用线程数： ",detectCores()))
  output=list()
    registerDoParallel(round(detectCores()))
    output<- foreach(itor = 1:mixtures) %dopar% {  
      library(limma)
      library(parallel)
      library(ggplot2)
      library(ggpubr)
      library(future.apply)
      library(car)
      library(ridge)
      library(e1071)
      library(preprocessCore)
      library(tcltk)
      library(limma)
      library(parallel)
      library(ggplot2)
      library(ggpubr)
      library(future.apply)
      library(car)
      library(ridge)
      library(e1071)
      library(preprocessCore)
      library(foreach)  
      library(doParallel)
      CoreAlg <- function(X, y){
        svn_itor <- 3
        res <- function(i){
          if(i==1){nus <- 0.25}
          if(i==2){nus <- 0.5}
          if(i==3){nus <- 0.75}
          model<-svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
          model
        }
        if(Sys.info()['sysname'] == 'Windows') out <- mclapply(1:svn_itor, res, mc.cores=1) else
          out <- mclapply(1:svn_itor, res, mc.cores=svn_itor)
        
        nusvm <- rep(0,svn_itor)
        corrv <- rep(0,svn_itor)

        t <- 1
        
        while(t <= svn_itor) {
          weights = t(out[[t]]$coefs) %*% out[[t]]$SV
          weights[which(weights<0)]<-0
          w<-weights/sum(weights)
          u <- sweep(X,MARGIN=2,w,'*')
          k <- apply(u, 1, sum)
          nusvm[t] <- sqrt((mean((k - y)^2)))
          corrv[t] <- cor(k, y)
          t <- t + 1
        }
        rmses <- nusvm
        mn <- which.min(rmses)
        model <- out[[mn]]
        q <- t(model$coefs) %*% model$SV
        q[which(q<0)]<-0
        w <- (q/sum(q))
        mix_rmse <- rmses[mn]
        mix_r <- corrv[mn]
        newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
        
      }
    y <- Y[,itor]
    y <- (y - mean(y)) / sd(y)
    result <- CoreAlg(X, y)
    w <- result$w
    mix_r <- result$mix_r
    mix_rmse <- result$mix_rmse
    if(P > 0) {pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))}
    out <- c(colnames(Y)[itor],w,pval,mix_r,mix_rmse)
    output[[itor]]<-out
    return(output[[itor]])
    }
    stopImplicitCluster()
    output<-do.call(cbind,output)
  output<-t(output)
  output<-rbind(header,output)
  write.table(output, file="CIBERSORT-Results.txt", sep="\t", row.names=F, col.names=F, quote=F)
  message("运行结束，结果导出")
}

