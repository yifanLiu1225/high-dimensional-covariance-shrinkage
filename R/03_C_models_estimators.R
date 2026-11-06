# C 模块：零第三方依赖的模型、估计量和指标。
# 与上传的 A/B 快照参数一致。B 的三个函数逐式移植，见 provenance/README.md。
regimes_C <- data.frame(regime=c('low','moderate','high','near_singular','singular'),
                       n=c(500L,200L,200L,150L,100L),p=c(10L,50L,150L,140L,200L))
structures_C <- c('identity','ar1','block','spiked')
methods_C <- c('sample','sphere_025','sphere_050','sphere_100','sphere_auto','diag_auto')

make_covariance_C <- function(p, structure, rho=0.6) {
  stopifnot(p>=2L, structure %in% structures_C)
  if (structure=='identity') return(diag(p))
  if (structure=='ar1') return(rho^abs(outer(seq_len(p),seq_len(p),'-')))
  if (structure=='block') {
    groups <- rep(seq_len(ceiling(p/10)),each=10,length.out=p)
    S <- matrix(0.1,p,p); S[outer(groups,groups,'==')] <- 0.7
    diag(S) <- 1; return(S)
  }
  # 保留调用者的 RNG 状态。协方差方向固定，但不重置模拟样本流。
  existed <- exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE)
  if (existed) old <- get('.Random.seed',envir=.GlobalEnv)
  on.exit(if(existed) assign('.Random.seed',old,envir=.GlobalEnv) else
    if(exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE)) rm('.Random.seed',envir=.GlobalEnv))
  set.seed(1000+p)
  Q <- qr.Q(qr(matrix(rnorm(p*p),p)))
  Q %*% diag(c(6,3,rep(1,p-2))) %*% t(Q)
}
make_shrinkage_target_C <- function(S,target) {
  if(target=='spherical') return(diag(mean(diag(S)),ncol(S)))
  if(target=='diagonal') return(diag(diag(S),ncol(S)))
  stop('Unsupported target')
}
linear_shrinkage_C <- function(S,T,lambda) (1-lambda)*S+lambda*T
estimate_intensity_C <- function(S,n,target) {
  stopifnot(n>2L,target %in% c('spherical','diagonal'))
  m <- n-1; den <- (m-1)*(m+2)
  V <- ((m-2)*S^2+m*outer(diag(S),diag(S)))/den
  if(target=='spherical') {
    trace_square <- (m^2*sum(S^2)-m*sum(diag(S))^2)/den
    numerator <- sum(V)-2*trace_square/(m*ncol(S))
  } else numerator <- sum(V)-sum(diag(V))
  distance <- sum((S-make_shrinkage_target_C(S,target))^2)
  if(distance<=.Machine$double.eps*sum(S^2)) return(1)
  min(1,max(0,numerator/distance))
}
mp_edges_C <- function(c) c(lower=(1-sqrt(c))^2,upper=(1+sqrt(c))^2)
mp_density_C <- function(x,c) {
  e <- mp_edges_C(c); z <- numeric(length(x)); ok <- x>e[1]&x<e[2]&x>0
  z[ok] <- sqrt((e[2]-x[ok])*(x[ok]-e[1]))/(2*pi*c*x[ok]); z
}
mp_atom_C <- function(c) max(0,1-1/c)
# finite-c proxy uses df=n-1 because S is centered and unbiased.
pca_reference_C <- function(alpha,c) {
  if(alpha<=1+sqrt(c)) return(c(sample_limit=(1+sqrt(c))^2,alignment=0))
  c(sample_limit=alpha*(1+c/(alpha-1)),alignment=(1-c/(alpha-1)^2)/(1+c/(alpha-1)))
}
spectral_metrics_C <- function(A,Sigma,eig=NULL,U=NULL,truth=NULL) {
  p <- ncol(A)
  if(is.null(eig)) eig <- eigen(A,symmetric=TRUE,only.values=TRUE)$values
  tol <- max(1,abs(eig[1]))*1e-8
  # 未把微小负值截掉后冒充正定。原始最小值与容差分开记录。
  rank <- sum(eig>tol)
  E <- A-Sigma; err_eig <- eigen(E,symmetric=TRUE,only.values=TRUE)$values
  true_eig <- if(is.null(truth)) eigen(Sigma,symmetric=TRUE,only.values=TRUE)$values else truth$values
  sorted_cor <- if(sd(true_eig)<1e-12||sd(eig)<1e-12) NA_real_ else cor(eig,true_eig)
  # Pearson 相关：按特征值降序配对，区别于“先排序再做 Spearman”的恒为1指标。
  out <- c(dmax=eig[1],dmin=eig[p],rank=rank,zero_fraction=(p-rank)/p,
           kappa=if(rank==p) eig[1]/eig[p] else Inf,
           frobenius=sum(E^2),spectral=max(abs(err_eig)),
           dmax_abs_error=abs(eig[1]-true_eig[1]),
           sorted_spectrum_correlation=sorted_cor,
           alignment1=NA_real_,alignment2=NA_real_,subspace_loss=NA_real_)
  if(!is.null(U)&&!is.null(truth)) {
    out['alignment1'] <- sum(U[,1]*truth$vectors[,1])^2
    out['alignment2'] <- sum(U[,2]*truth$vectors[,2])^2
    out['subspace_loss'] <- max(0,4-2*sum(crossprod(U[,1:2],truth$vectors[,1:2])^2))
  }
  out
}
