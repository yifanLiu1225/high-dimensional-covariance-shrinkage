# D 组：考察协方差收缩对马氏距离估计和 Hotelling 均值检验的影响。
# RStudio：打开本文件，点击代码编辑区，再按 Ctrl+Alt+R。
# Windows 控制台：Rscript --vanilla 04_run_member_D.R
# 正式设置是每种条件各 500 次校准和评估。请等待控制台显示“完成”。

# R 在少数 Windows 终端会以非 UTF-8 区域设置启动。
if (.Platform$OS.type == 'windows' && !isTRUE(l10n_info()[['UTF-8']])) {
  invisible(suppressWarnings(Sys.setlocale('LC_CTYPE', '.UTF-8')))
}

# 1. 总体模型与估计器
member_D_design <- function() {
  data.frame(regime = c('low', 'moderate', 'high', 'near_singular', 'singular'),
             n = c(500, 200, 200, 150, 100), p = c(10, 50, 150, 140, 200))
}

member_D_iris_cor <- function() {
  x <- as.matrix(iris[, 1:4])
  for (j in seq_len(ncol(x))) x[, j] <- x[, j] - ave(x[, j], iris$Species)
  cor(x)
}

member_D_sigma <- function(p, structure) {
  index <- seq_len(p)
  if (structure == 'identity') return(diag(p))
  if (structure == 'ar1') return(0.6 ^ abs(outer(index, index, '-')))
  if (structure == 'block') {
    groups <- (index - 1) %/% 10
    S <- matrix(0.1, p, p)
    S[outer(groups, groups, '==')] <- 0.7
    diag(S) <- 1
    return(S)
  }
  if (structure == 'spiked') {
    # B、C 用同一种子固定主方向。和 C 一样，在退出时恢复抽样随机数状态。
    had_seed <- exists('.Random.seed', envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get('.Random.seed', envir = .GlobalEnv)
    on.exit(if (had_seed) assign('.Random.seed', old_seed, envir = .GlobalEnv) else
      rm('.Random.seed', envir = .GlobalEnv))
    set.seed(1000 + p)
    Q <- qr.Q(qr(matrix(rnorm(p * p), p)))
    return(Q %*% diag(c(6, 3, rep(1, p - 2))) %*% t(Q))
  }
  if (structure == 'semi_real') {
    stopifnot(p >= 4)
    # 保留标准化 iris 物种内相关结构。其余维度是与前4维独立的标准正态噪声。
    S <- diag(p)
    S[1:4, 1:4] <- member_D_iris_cor()
    return(S)
  }
  stop('Unknown covariance structure: ', structure)
}

member_D_strengths <- function(S, m) {
  # m 为有效协方差自由度。普通模拟 m=n-1；iris 合并三个物种的残差，m=90-3。
  p <- ncol(S)
  variances <- diag(S)
  scale_squared <- sum(S^2)
  denominator <- (m - 1) * (m + 2)
  V <- ((m - 2) * S^2 + m * outer(variances, variances)) / denominator
  trace_square <- (m^2 * scale_squared - m * sum(variances)^2) / denominator
  sphere <- diag(mean(variances), p)
  diagonal <- diag(variances, p)

  # 目标已与 S 重合时取强度 1。阈值随数据单位一起缩放，沿用 B、C 的判定。
  intensity <- function(numerator, target) {
    distance <- sum((S - target)^2)
    if (distance <= .Machine$double.eps * scale_squared) return(1)
    min(1, max(0, numerator / distance))
  }
  c(sphere_auto = intensity(sum(V) - 2 * trace_square / (m * p), sphere),
    diag_auto = intensity(sum(V) - sum(diag(V)), diagonal))
}

member_D_estimates <- function(S, m, strengths = member_D_strengths(S, m)) {
  p <- ncol(S)
  sphere <- diag(mean(diag(S)), p)
  diagonal <- diag(diag(S), p)
  list(sample = S,
       sphere_auto = (1 - strengths[1]) * S + strengths[1] * sphere,
       diag_auto = (1 - strengths[2]) * S + strengths[2] * diagonal)
}

member_D_wilson <- function(k,N) {
  if(!N)return(c(NA_real_,NA_real_))
  z<-qnorm(.975);ph<-k/N;a<-(ph+z*z/(2*N))/(1+z*z/N)
  b<-z*sqrt(ph*(1-ph)/N+z*z/(4*N*N))/(1+z*z/N)
  pmax(0,pmin(1,c(a-b,a+b)))
}

# 2. 模拟及鸢尾花留出实验
run_member_D <- function(B=500L,Bcal=500L,root='.',seed=20261107L) {
  stopifnot(B>=20,Bcal>=20)
  td<-file.path(root,'tables/member_D');dd<-file.path(root,'data/member_D')
  dir.create(td,recursive=TRUE,showWarnings=FALSE);dir.create(dd,recursive=TRUE,showWarnings=FALSE)
  settings<-member_D_design();structures<-c('identity','ar1','block','spiked','semi_real')
  methods<-c('sample','sphere_auto','diag_auto');raws<-list();tests<-list();criticals<-list();k<-0L
  for(i in seq_len(nrow(settings)))for(st in structures){
    k<-k+1L;set.seed(seed+1000L*k)
    n<-settings$n[i];p<-settings$p[i];regime<-settings$regime[i];m<-n-1
    Sigma<-member_D_sigma(p,st);C<-chol(Sigma);P<-chol2inv(C)
    delta<-sqrt(10/n)*as.numeric(t(C)[,1])
    cat(sprintf('%02d/25 %s %s n=%d p=%d\n',k,st,regime,n,p));flush.console()
    simulate<-function(evaluate=FALSE){
      # 正态样本n-1个独立对比，与通常的中心化S具有完全相同的Wishart分布。
      X<-matrix(rnorm(m*p),m,p)%*%C;S<-crossprod(X)/m
      strengths <- member_D_strengths(S, m)
      estimates <- member_D_estimates(S, m, strengths)
      z<-as.numeric(rnorm(p)%*%C) # sqrt(n)*(样本均值-mu0)，独立于S
      condition_before_noise <- NA_real_
      if (evaluate) {
        v <- as.numeric(rnorm(p) %*% C)
        truth <- as.numeric(crossprod(v, P %*% v))
        if (st == 'semi_real') {
          # 同一样本先取前4维，再与加入独立噪声后的矩阵比较。
          baseline <- S[1:4, 1:4]
          condition_before_noise <- norm(baseline, 'O') * norm(chol2inv(chol(baseline)), 'O')
        }
      }
      ans<-lapply(methods,function(method){
        A<-estimates[[method]]
        # Cholesky判断和求解，无广义逆；样本p>=n确定不可逆。
        ch<-if(method=='sample' && p>=n)NULL else tryCatch(chol(A),error=function(e)NULL)
        valid<-!is.null(ch)
        qfun<-function(w)if(valid)sum(backsolve(ch,w,transpose=TRUE)^2) else NA_real_
        null<-qfun(z)
        if(!evaluate)return(null)
        dq <- qfun(v)
        inverse <- if (valid) chol2inv(ch) else NULL
        condition_after_noise <- if (st != 'semi_real') NA_real_ else
          if (valid) norm(A, 'O') * norm(inverse, 'O') else Inf
        data.frame(method,inverse_failure=!valid,frobenius=sum((A-Sigma)^2),
          inverse_error=if(valid)sum((inverse-P)^2)/sum(P^2) else NA_real_,
          condition_before_noise=condition_before_noise, condition_after_noise=condition_after_noise,
          distance_ratio=dq/truth,distance_rel_error=abs(dq-truth)/truth,
          true_distance=truth,distance=dq,t_null=null,t_alt=qfun(z+sqrt(n)*delta),
          lambda=if(method=='sample')0 else unname(strengths[method]))
      })
      if(evaluate)do.call(rbind,ans) else setNames(unlist(ans),methods)
    }
    cal<-replicate(Bcal,simulate(FALSE))
    cuts<-apply(cal,1,function(z)if(all(is.na(z)))NA_real_ else unname(quantile(z,.95,na.rm=TRUE,type=1)))
    a<-do.call(rbind,lapply(seq_len(B),function(b){x<-simulate(TRUE);x$rep<-b;x}))
    a$structure<-st;a$regime<-regime;a$n<-n;a$p<-p;a$condition_seed<-seed+1000L*k
    raws[[k]]<-a
    criticals[[k]]<-data.frame(structure=st,regime,n,p,method=methods,threshold=cuts,Bcal)
    classic<-if(p<n)p*(n-1)/(n-p)*qf(.95,p,n-p) else NA_real_
    for(mt in methods)for(rule in c('classical_F','oracle_MC'))for(state in c('null','alt')){
      b<-a[a$method==mt,];cut<-if(rule=='classical_F')classic else cuts[mt]
      v<-b[[paste0('t_',state)]];good<-is.finite(v)&is.finite(cut);N<-sum(good)
      count<-if(N)sum(v[good]>cut) else NA_integer_;ci<-if(N)member_D_wilson(count,N) else c(NA,NA)
      tests[[length(tests)+1L]]<-data.frame(structure=st,regime,n,p,method=mt,rule,state,threshold=cut,N,
        rejected=count,rate=if(N)count/N else NA_real_,lo=ci[1],hi=ci[2])
    }
  }
  raw<-do.call(rbind,raws);test<-do.call(rbind,tests)
  metrics<-c('frobenius','inverse_error','distance_ratio','distance_rel_error','lambda')
  summarize<-function(x){
    N<-sum(is.finite(x));v<-x[is.finite(x)]
    if(!N)return(c(N=0,mean=NA,sd=NA,mcse=NA,lo=NA,hi=NA))
    se<-sd(v)/sqrt(N);ci<-qt(.975,N-1)*se
    c(N=N,mean=mean(v),sd=sd(v),mcse=se,lo=mean(v)-ci,hi=mean(v)+ci)
  }
  groups<-split(raw,interaction(raw$structure,raw$regime,raw$method,drop=TRUE))
  summary<-do.call(rbind,lapply(groups,function(a)do.call(rbind,lapply(metrics,function(me){
    ss<-as.data.frame(as.list(summarize(a[[me]])))
    cbind(a[1,c('structure','regime','n','p','method')],metric=me,ss)
  }))))
  paired<-list()
  for(st in structures)for(reg in settings$regime)for(mt in methods[-1])for(me in c('frobenius','inverse_error','distance_rel_error')){
    a<-raw[raw$structure==st&raw$regime==reg&raw$method==mt,]
    b<-raw[raw$structure==st&raw$regime==reg&raw$method=='sample',]
    stopifnot(identical(a$rep,b$rep))
    paired[[length(paired)+1L]]<-cbind(data.frame(structure=st,regime=reg,method=mt,metric=paste0('delta_',me)),as.data.frame(as.list(summarize(a[[me]]-b[[me]]))))
  }
  failure<-aggregate(inverse_failure~structure+regime+n+p+method,raw,mean)
  write.csv(summary,file.path(td,'summary_long.csv'),row.names=FALSE)
  write.csv(test,file.path(td,'rejection_rates.csv'),row.names=FALSE)
  write.csv(do.call(rbind,paired),file.path(td,'paired_differences.csv'),row.names=FALSE)
  write.csv(failure,file.path(td,'inverse_failure_rates.csv'),row.names=FALSE)
  # 原始记录保留每次的前后配对，汇总表同时给出有限样本数及中位数。
  semi <- raw[raw$structure == 'semi_real', ]
  condition_summary <- do.call(rbind, lapply(settings$regime, function(regime_name) {
    group <- semi[semi$regime == regime_name, ]
    baseline <- group$condition_before_noise[group$method == 'sample']
    conditions <- c(list(before_noise = baseline),
                    setNames(lapply(methods, function(method_name) {
                      group$condition_after_noise[group$method == method_name]
                    }), methods))
    do.call(rbind, lapply(names(conditions), function(name) {
      x <- conditions[[name]]
      data.frame(regime = regime_name, n = group$n[1],
                 p = if (name == 'before_noise') 4 else group$p[1],
                 method = name, N = length(x), n_finite = sum(is.finite(x)),
                 median_condition_1 = median(x))
    }))
  }))
  write.csv(condition_summary, file.path(td, 'semi_real_condition.csv'), row.names=FALSE)
  write.csv(do.call(rbind,criticals),file.path(td,'oracle_critical_values.csv'),row.names=FALSE)
  con<-gzfile(file.path(dd,'replicate_metrics.csv.gz'),'wt');write.csv(raw,con,row.names=FALSE);close(con)
  write.csv(member_D_iris_cor(),file.path(td,'iris_within_species_correlation.csv'))
  capture.output(sessionInfo(),file=file.path(dd,'sessionInfo.txt'))
  writeLines(c(paste('seed',seed),paste('evaluation repetitions',B),paste('independent calibration repetitions',Bcal),
    'Gaussian samples; alpha=.05; alternative n*delta^T Sigma^-1 delta=10.',
    'Semi-real population: iris within-species correlation block plus independent standard normal noise.',
    'Condition comparison: exact 1-norm condition numbers on the same semi-real sample, p=4 and p>=10.',
    'oracle_MC uses true Sigma, only an ideal simulation benchmark.',
    'Wilson intervals condition on the estimated threshold, excluding calibration uncertainty.'),file.path(dd,'settings.txt'))
  invisible(list(raw=raw,summary=summary,tests=test))
}

# 真实低维控制组：iris物种内残差的留出高斯预测评分，绝不虚构真实Σ。
run_member_D_iris <- function(root='.',B=500L){
  set.seed(20261201);records<-list()
  for(b in seq_len(B)){
    train<-unlist(lapply(split(seq_len(nrow(iris)),iris$Species),function(ii)sample(ii,30)))
    test<-setdiff(seq_len(nrow(iris)),train)
    x<-as.matrix(iris[train,1:4]);y<-as.matrix(iris[test,1:4])
    for(sp in levels(iris$Species)){
      ii<-iris$Species[train]==sp;jj<-iris$Species[test]==sp
      mu<-colMeans(x[ii,,drop=FALSE]);x[ii,]<-sweep(x[ii,,drop=FALSE],2,mu);y[jj,]<-sweep(y[jj,,drop=FALSE],2,mu)
    }
    scales<-sqrt(colSums(x*x)/(length(train)-3));x<-sweep(x,2,scales,'/');y<-sweep(y,2,scales,'/')
    S<-crossprod(x)/(length(train)-3)
    estimates<-member_D_estimates(S,length(train)-3)
    for(mt in names(estimates)){
      C<-chol(estimates[[mt]]);scores<-2*sum(log(diag(C)))+rowSums(t(backsolve(C,t(y),transpose=TRUE))^2)
      records[[length(records)+1L]]<-data.frame(rep=b,method=mt,score=mean(scores))
    }
  }
  raw<-do.call(rbind,records)
  write.csv(raw,file.path(root,'tables/member_D/iris_holdout_scores.csv'),row.names=FALSE)
  invisible(raw)
}

# 3. 距离误差、检验表现和应用对照
plot_member_D <- function(root='.') {
  if(!requireNamespace('ggplot2',quietly=TRUE))stop('绘图需要ggplot2。请在R中安装install.packages("ggplot2")')
  library(ggplot2)
  td<-file.path(root,'tables/member_D');fd<-file.path(root,'figures/member_D')
  dir.create(fd,recursive=TRUE,showWarnings=FALSE)
  s<-read.csv(file.path(td,'summary_long.csv'));r<-read.csv(gzfile(file.path(root,'data/member_D/replicate_metrics.csv.gz')))
  t<-read.csv(file.path(td,'rejection_rates.csv'));iris_scores<-read.csv(file.path(td,'iris_holdout_scores.csv'))
  regimes<-c('low','moderate','high','near_singular','singular')
  labels<-c('低维','中等','高维','近奇异','奇异')
  prep<-function(d){
    d$regime<-factor(d$regime,regimes,labels)
    d$method<-factor(d$method,c('sample','sphere_auto','diag_auto'),c('样本协方差','球形自动收缩','对角自动收缩'))
    d$structure<-factor(d$structure,c('identity','ar1','block','spiked','semi_real'),c('Identity','AR(1)','Block','Spiked','半真实'))
    d
  }
  s<-prep(s);r<-prep(r);t<-prep(t)
  pal<-c('样本协方差'='#444E5C','球形自动收缩'='#087F8C','对角自动收缩'='#D65A83')
  previous_theme <- theme_get()
  on.exit(theme_set(previous_theme), add=TRUE)
  B <- length(unique(r$rep))
  theme_set(theme_minimal(base_size=13,base_family='Microsoft YaHei')+
    theme(panel.grid.minor=element_blank(),legend.position='bottom',legend.title=element_blank(),
          plot.title=element_text(face='bold',size=19,color='#17384B'),plot.subtitle=element_text(size=12,color='#627383'),
          strip.text=element_text(face='bold',size=14),axis.title=element_text(size=12),plot.margin=margin(16,18,12,14)))
  savefig<-function(g,name,w=12,h=7){
    ggsave(file.path(fd,paste0(name,'.png')),g,width=w,height=h,dpi=190)
    ggsave(file.path(fd,paste0(name,'.pdf')),g,width=w,height=h,device=cairo_pdf)
  }
  # 图 1 比较四种理论总体下的距离平方误差，观察不同样本量与维数组合下各估计法的表现。
  # 纵轴为相对绝对误差，越低表示估计越接近真实距离平方。
  a<-subset(s,structure!='半真实'&metric=='distance_rel_error'&is.finite(mean))
  g<-ggplot(a,aes(regime,mean,color=method,group=method))+
    geom_line(linewidth=.65)+geom_point(size=2)+geom_errorbar(aes(ymin=pmax(lo,1e-5),ymax=hi),width=.12)+
    facet_wrap(~structure,nrow=2)+scale_color_manual(values=pal)+scale_y_log10()+
    labs(title='普通与收缩协方差的马氏距离平方误差',subtitle=sprintf('独立新观测，均值已知；点为%d次均值，误差线为95%% Monte Carlo区间', B),x=NULL,y='平均相对绝对误差（对数刻度）',caption='样本协方差在奇异组无逆矩阵，因此没有对应点；连线仅辅助阅读不同设计点。')
  savefig(g,'01_mahalanobis_error')
  # 图 2 用高维组的箱线图检查距离平方是否被系统放大，以及每次抽样的波动。
  # 比值 1 表示估计值等于真实值。
  a<-subset(r,structure!='半真实'&regime==labels[3]&is.finite(distance_ratio))
  g<-ggplot(a,aes(method,distance_ratio,fill=method))+geom_boxplot(outlier.size=.5,width=.55)+
    geom_hline(yintercept=1,linetype=2,color='#B76835')+facet_wrap(~structure,nrow=1)+
    scale_fill_manual(values=pal)+scale_y_log10()+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    labs(title='高维组的距离膨胀与抽样波动',subtitle=sprintf('n=200，p=150；每种结构%d次；虚线1表示估计距离等于真实距离', B),x=NULL,y='估计距离平方 / 真实距离平方（对数刻度）')
  savefig(g,'02_distance_distribution',12,5.5)
  # 图 3 检查原假设成立时的拒绝率是否接近 5%，分别使用经典 F 阈值和理想模拟阈值。
  a<-subset(t,structure!='半真实'&state=='null'&is.finite(rate));a$rule<-factor(a$rule,c('classical_F','oracle_MC'),c('直接使用经典F临界值','理想模拟校准（已知真实Σ）'))
  pos<-position_dodge(width=.45)
  g<-ggplot(a,aes(regime,rate,color=method))+geom_hline(yintercept=.05,linetype=2,color='#B76835')+
    geom_errorbar(aes(ymin=lo,ymax=hi),position=pos,width=.13)+geom_point(position=pos,size=1.8)+
    facet_grid(rule~structure)+scale_color_manual(values=pal)+scale_y_continuous(labels=function(x)paste0(100*x,'%'))+
    labs(title='检验拒绝率：可逆性改善以后还要检查标尺',subtitle='H0为真，名义水平5%；误差线为95% Wilson区间',x=NULL,y='原假设拒绝率',caption='经典F仅适用于普通S且n>p，收缩套F只作反例。\n理想模拟校准使用真实Σ，不是未知总体下的现成检验。')+
    theme(axis.text.x=element_text(size=7),strip.text=element_text(size=11))
  savefig(g,'03_null_rejection_rates',14,7.5)
  # 图 4 在指定均值偏移下比较检出率。阈值来自独立的原假设校准样本，并按方法分别估计。
  a<-subset(t,structure!='半真实'&state=='alt'&rule=='oracle_MC'&is.finite(rate))
  g<-ggplot(a,aes(regime,rate,color=method,group=method))+geom_line(linewidth=.6)+geom_point(size=2)+
    geom_errorbar(aes(ymin=lo,ymax=hi),width=.12)+facet_wrap(~structure,nrow=1)+scale_color_manual(values=pal)+
    scale_y_continuous(limits=c(0,1),labels=function(x)paste0(100*x,'%'))+
    labs(title='功效取决于结构与均值偏移方向',subtitle='理想校准下，一个预先指定的方向；非中心参数n δᵀΣ⁻¹δ=10',x=NULL,y='备择下拒绝率',caption='仅报告这一备择设定，不由此宣称收缩对所有方向都更有检验力。')+theme(axis.text.x=element_text(size=8))
  savefig(g,'04_power',13,5.5)
  # 图 5 保留鸢尾花的四维相关结构并加入独立噪声，比较矩阵损失和距离平方误差。
  a<-subset(s,structure=='半真实'&metric%in%c('frobenius','distance_rel_error')&is.finite(mean))
  a$metric<-factor(a$metric,c('frobenius','distance_rel_error'),c('平方Frobenius损失','马氏距离平方的相对绝对误差'))
  g<-ggplot(a,aes(regime,mean,color=method,group=method))+geom_line()+geom_point(size=2.5)+
    facet_wrap(~metric,scales='free_y')+scale_y_log10()+scale_color_manual(values=pal)+
    labs(title='半真实实验：保留真实相关性，加入独立噪声',subtitle=sprintf('前4维采用iris物种内相关矩阵，新增维度为独立标准正态噪声；%d次评估', B),x=NULL,y='Monte Carlo均值（对数刻度）',caption='前4维及新增变量均从指定总体模拟；实际测量只用于估计前4维相关结构。')
  savefig(g,'05_semi_real',12,6)
  # 图 6 比较真实低维鸢尾花留出数据上的高斯预测评分，观察收缩的收益。评分越小越好。
  iris_scores$method<-factor(iris_scores$method,c('sample','sphere_auto','diag_auto'),names(pal))
  g<-ggplot(iris_scores,aes(method,score,fill=method))+geom_boxplot(width=.5,outlier.size=.6)+
    scale_fill_manual(values=pal)+labs(title='真实低维控制组：iris留出预测评分',subtitle=sprintf('每次每物种30条训练、20条测试；训练集估计物种均值和尺度，重复%d次随机划分', length(unique(iris_scores$rep))),x=NULL,y='留出高斯评分（越小越好）',caption='评分 = log|Σ̂| + 平均马氏距离平方，省略共同常数。重复划分复用同一数据集，不视为500个独立真实样本。')
  savefig(g,'06_iris_control',10,6)
  # 图 7 用理论曲线说明小特征值求逆后会产生很大的权重，以及球形收缩如何限制该权重。
  d<-seq(.005,1,length.out=300);a<-data.frame(eigenvalue=rep(d,2),weight=c(1/d,1/(.75*d+.25)),method=rep(c('样本逆特征值','固定λ=0.25球形收缩'),each=length(d)))
  g<-ggplot(a,aes(eigenvalue,weight,color=method))+geom_line(linewidth=1.1)+
    scale_color_manual(values=c('#087F8C','#444E5C'))+labs(title='小特征值怎样被求逆放大',subtitle='理论示意：平均方差τ=1，固定收缩强度λ=0.25',x='样本特征值 d',y='该方向在逆矩阵中的权重')+
    coord_cartesian(ylim=c(0,25))
  savefig(g,'07_inverse_weight',10,5.5)
  # 图 8 用高维 AR(1) 的原假设统计量直方图对照两类阈值，解释收缩后照用经典阈值时为何很少拒绝。
  null_draws <- subset(r, structure == 'AR(1)' & regime == '高维' & is.finite(t_null))
  cuts <- subset(t, structure == 'AR(1)' & regime == '高维' & state == 'null')
  cuts$threshold_source <- factor(cuts$rule, levels = c('classical_F', 'oracle_MC'),
                                  labels = c('经典 F 临界值', '本方法的理想校准值'))
  g <- ggplot(null_draws, aes(t_null)) +
    geom_histogram(bins = 52, fill = '#477888', color = 'white', linewidth = .18) +
    geom_vline(data = cuts, aes(xintercept = threshold, linetype = threshold_source,
                                color = threshold_source), linewidth = .85) +
    facet_grid(method ~ .) +
    scale_linetype_manual(values = c('solid', 'dashed')) +
    scale_color_manual(values = c('#333D44', '#AF3E20')) +
    labs(title = '同一条经典临界线，在三种统计量分布中的位置',
         subtitle = sprintf('AR(1)，n = 200，p = 150。每行来自 %d 次独立评估，直方条为落入各区间的次数。', length(unique(null_draws$rep))),
         x = '原假设下的二次型统计量', y = '模拟次数',
         caption = '右尾超过临界线即拒绝。虚线取自另一组已知真实协方差的模拟。') +
    theme(strip.text.y = element_text(angle = 0), legend.box = 'horizontal')
  savefig(g, '08_hotelling_null_distribution', 11, 6.5)
  invisible(TRUE)
}

# 4. 自动定位结果目录，运行正式模拟并出图
member_D_output_parent <- function() {
  # RStudio 的“运行全部”按编辑区的文件位置保存；纯控制台粘贴时使用当前目录。
  path <- character()
  if (interactive() && requireNamespace('rstudioapi', quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    path <- tryCatch(rstudioapi::getSourceEditorContext()$path,
                     error = function(e) character())
  }
  if (!length(path) || !nzchar(path)) {
    cli <- grep('^--file=', commandArgs(), value = TRUE)
    path <- sub('^--file=', '', cli)
  }
  # 保留团队通过整份脚本调用的兼容方式。
  if (!length(path) || identical(path, '-')) {
    paths <- unlist(lapply(sys.frames(), function(frame) frame$ofile), use.names = FALSE)
    if (length(paths)) path <- tail(paths, 1)
  }
  if (!length(path) || !nzchar(path) || !file.exists(path)) return(getwd())
  parent <- dirname(normalizePath(path, winslash = '/', mustWork = TRUE))
  if (basename(parent) == 'R') dirname(parent) else parent
}

member_D_run_all <- function() {
  if (!requireNamespace('ggplot2', quietly = TRUE)) {
    stop('首次运行需要 ggplot2。请在 RStudio 的 Console 输入 install.packages("ggplot2")，安装后再运行整个文件。',
         call. = FALSE)
  }
  if (!capabilities('cairo')) {
    stop('当前 R 缺少 Cairo 绘图支持。请在 RStudio 中选用支持 Cairo 的 R。', call. = FALSE)
  }
  reps <- suppressWarnings(as.integer(Sys.getenv('D_REPS', '500')))
  cal_reps <- suppressWarnings(as.integer(Sys.getenv('D_CAL_REPS', '500')))
  if (anyNA(c(reps, cal_reps)) || min(reps, cal_reps) < 20L) {
    stop('D_REPS 和 D_CAL_REPS 都需要是至少 20 的整数。', call. = FALSE)
  }

  parent <- member_D_output_parent()
  output <- file.path(parent, 'D_复现结果')
  if (file.exists(output)) output <- tempfile('D_复现结果_', tmpdir = parent)
  staging <- tempfile('.D_复现中_', tmpdir = parent)
  if (!dir.create(staging)) stop('无法创建结果目录，请把脚本放在可写文件夹。', call. = FALSE)

  message('开始计算。每组评估 ', reps, ' 次，另有 ', cal_reps, ' 次校准。')
  run_member_D(B = reps, Bcal = cal_reps, root = staging)
  run_member_D_iris(root = staging, B = reps)
  message('模拟结束，正在绘制距离误差、检验表现和应用对照。')
  plot_member_D(root = staging)

  if (!file.rename(staging, output)) {
    stop('结果已生成，请到这个目录查看：', staging, call. = FALSE)
  }
  figures <- normalizePath(file.path(output, 'figures/member_D'),
                           winslash = '/', mustWork = TRUE)
  message('完成。PNG 和 PDF 图表目录：', figures)
  invisible(figures)
}

member_D_run_all()
