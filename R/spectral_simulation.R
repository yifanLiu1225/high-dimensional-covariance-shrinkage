# 从项目根目录 source；或用 Rscript R/03_run_member_C.R --reps=500。
source('R/03_C_models_estimators.R',encoding='UTF-8')
run_spectral_C <- function(reps=500L,seed=20261103L,out='.') {
  stopifnot(reps>=2L)
  dirs <- c('data/member_C/spectra','tables/member_C','figures/member_C')
  for(d in dirs) dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
  config <- list(reps=reps,seed=seed,regimes=regimes_C,structures=structures_C,
                 methods=methods_C,denominator='n-1',tol_relative=1e-8)
  saveRDS(config,file.path(out,'data/member_C/config.rds'))
  all_rows <- vector('list',20); ci <- 0L
  for(i in seq_len(nrow(regimes_C))) for(st in structures_C) {
    ci <- ci+1L; n <- regimes_C$n[i]; p <- regimes_C$p[i]
    Sigma <- make_covariance_C(p,st); L <- chol(Sigma); truth <- eigen(Sigma,symmetric=TRUE)
    # 每个条件有独立确定性种子，顺序运行或单条件重现都可复核。
    condition_seed <- seed+1000L*ci
    set.seed(condition_seed)
    spectra <- array(NA_real_,c(reps,p,length(methods_C)),dimnames=list(NULL,NULL,methods_C))
    rows <- vector('list',reps)
    for(r in seq_len(reps)) {
      X <- matrix(rnorm(n*p),n,p)%*%L
      S <- cov(X); mu <- mean(diag(S))
      es <- eigen(S,symmetric=TRUE,only.values=st!='spiked')
      ls <- estimate_intensity_C(S,n,'spherical'); ld <- estimate_intensity_C(S,n,'diagonal')
      lambdas <- c(0,.25,.5,1,ls,ld)
      targets <- c('none',rep('spherical',4),'diagonal')
      A_diag <- linear_shrinkage_C(S,make_shrinkage_target_C(S,'diagonal'),ld)
      ed <- eigen(A_diag,symmetric=TRUE,only.values=st!='spiked')
      metric <- vector('list',length(methods_C))
      for(j in seq_along(methods_C)) {
        lambda <- lambdas[j]
        if(j==1L) {A <- S; d <- es$values; U <- es$vectors}
        else if(j==6L) {A <- A_diag;d <- ed$values;U <- ed$vectors}
        else {
          A <- (1-lambda)*S+lambda*diag(mu,p)
          d <- (1-lambda)*es$values+lambda*mu
          U <- if(lambda<1-1e-12) es$vectors else NULL
        }
        if(st!='spiked') U <- NULL
        v <- spectral_metrics_C(A,Sigma,d,U,truth)
        spectra[r,,j] <- d
        metric[[j]] <- data.frame(regime=regimes_C$regime[i],n=n,p=p,c_n=p/n,c_df=p/(n-1),
          structure=st,rep=r,seed=condition_seed,method=methods_C[j],lambda=lambda,mu=mu,
          as.list(v),fingerprint_S11=S[1,1],check.names=FALSE)
      }
      rows[[r]] <- do.call(rbind,metric)
      # 小型完整样本夹具供独立 Python 程序复核，不保存全量 X。
      if(r==1L&&st=='spiked'&&regimes_C$regime[i]=='moderate') {
        for(nm in c('X','S','Sigma','A_diag')) write.csv(get(nm),
          file.path(out,'data/member_C',paste0('audit_',nm,'.csv')),row.names=FALSE)
      }
    }
    part <- do.call(rbind,rows); all_rows[[ci]] <- part
    key <- paste(regimes_C$regime[i],st,sep='_')
    saveRDS(list(spectra=spectra,Sigma=Sigma,truth=truth$values,config=config,
      regime=regimes_C$regime[i],n=n,p=p,condition_seed=condition_seed),
      file.path(out,'data/member_C/spectra',paste0(key,'.rds')),compress='gzip')
    message(sprintf('[%02d/20] %s: %d repeats complete',ci,key,reps))
  }
  raw <- do.call(rbind,all_rows); rownames(raw)<-NULL
  con <- gzfile(file.path(out,'tables/member_C/replicate_metrics.csv.gz'),'wt')
  write.csv(raw,con,row.names=FALSE,na='NA'); close(con)
  writeLines(c(capture.output(sessionInfo()),paste('RNGkind:',paste(RNGkind(),collapse=', '))),
             file.path(out,'data/member_C/sessionInfo.txt'))
  raw
}
summarise_C <- function(raw,out='.') {
  metrics <- c('lambda','dmin','dmax','kappa','zero_fraction','frobenius','spectral',
               'dmax_abs_error','sorted_spectrum_correlation','alignment1','alignment2','subspace_loss')
  groups <- split(raw,interaction(raw$regime,raw$structure,raw$method,drop=TRUE))
  stats <- lapply(groups,function(d) do.call(rbind,lapply(metrics,function(m) {
    x <- d[[m]]; finite <- x[is.finite(x)]; k <- length(finite); B <- nrow(d)
    avg <- if(any(is.infinite(x))) Inf else if(k) mean(finite) else NA_real_
    ss <- if(any(is.infinite(x))) NA_real_ else if(k>1L) sd(finite) else NA_real_
    se <- ss/sqrt(k); delta <- if(k>1L) qt(.975,k-1)*se else NA_real_
    data.frame(d[1,c('regime','n','p','c_n','c_df','structure','method')],metric=m,
      reps=B,n_finite=k,n_missing=sum(is.na(x)),n_infinite=sum(is.infinite(x)),
      mean=avg,sd=ss,mcse=se,ci_lower=avg-delta,ci_upper=avg+delta,
      median=if(k) median(x,na.rm=TRUE) else if(any(is.infinite(x))) Inf else NA_real_)
  })))
  summary <- do.call(rbind,stats); rownames(summary)<-NULL
  write.csv(summary,file.path(out,'tables/member_C/summary_long.csv'),row.names=FALSE)
  # 按同一份样本做成对差值，避免把不同方法视为独立样本。
  baseline <- raw[raw$method=='sample',]
  keys <- c('regime','structure','rep')
  paired <- merge(raw,baseline[,c(keys,'frobenius','spectral')],by=keys,suffixes=c('','_sample'))
  paired$delta_frobenius <- paired$frobenius-paired$frobenius_sample
  paired$delta_spectral <- paired$spectral-paired$spectral_sample
  ps <- do.call(rbind,lapply(split(paired,interaction(paired$regime,paired$structure,paired$method)),function(d) {
    do.call(rbind,lapply(c('delta_frobenius','delta_spectral'),function(m) {
      x<-d[[m]]; se<-sd(x)/sqrt(length(x)); t<-qt(.975,length(x)-1)
      data.frame(d[1,c('regime','structure','method')],metric=m,mean=mean(x),sd=sd(x),mcse=se,
        ci_lower=mean(x)-t*se,ci_upper=mean(x)+t*se)
    }))
  }))
  write.csv(ps,file.path(out,'tables/member_C/paired_differences.csv'),row.names=FALSE)
  sing <- do.call(rbind,lapply(groups,function(d) {
    fail <- sum(d$rank<d$p); B <- nrow(d); z <- qnorm(.975); ph <- fail/B
    center <- (ph+z*z/(2*B))/(1+z*z/B)
    half <- z*sqrt(ph*(1-ph)/B+z*z/(4*B*B))/(1+z*z/B)
    data.frame(d[1,c('regime','structure','method')],failures=fail,reps=B,rate=ph,
      wilson_lower=center-half,wilson_upper=center+half)
  }))
  write.csv(sing,file.path(out,'tables/member_C/singularity_rates.csv'),row.names=FALSE)
  edges <- transform(regimes_C,c_n=p/n,c_df=p/(n-1))
  edges$mp_lower <- (1-sqrt(edges$c_df))^2; edges$mp_upper <- (1+sqrt(edges$c_df))^2
  edges$mp_atom <- pmax(0,1-1/edges$c_df)
  edges$exact_zero_count <- pmax(0,edges$p-(edges$n-1))
  edges$kappa_proxy <- ifelse(edges$c_df<1,edges$mp_upper/edges$mp_lower,Inf)
  write.csv(edges,file.path(out,'tables/member_C/mp_reference.csv'),row.names=FALSE)
  invisible(summary)
}
