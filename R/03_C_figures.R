# 全部图来自已保存的正式重复实验；坐标轴为英文，报告提供中文图注。
plot_C <- function(raw,out='.') {
  figdir <- file.path(out,'figures/member_C')
  cols <- c(sample='#64748b',sphere_025='#84a39a',sphere_050='#e69f00',sphere_100='#cc6677',sphere_auto='#0072b2',diag_auto='#009e73')
  labels <- c(sample='Sample',sphere_025='Sphere .25',sphere_050='Sphere .50',sphere_100='Sphere 1.0',sphere_auto='Sphere auto',diag_auto='Diagonal auto')
  draw <- function(name,w,h,fun) {
    png(file.path(figdir,paste0(name,'.png')),width=w,height=h,res=150)
    par(mar=c(4,4,3,1),mgp=c(2.4,.7,0),las=1,cex=.85)
    fun(); dev.off()
    pdf(file.path(figdir,paste0(name,'.pdf')),width=w/150,height=h/150,useDingbats=FALSE)
    par(mar=c(4,4,3,1),mgp=c(2.4,.7,0),las=1,cex=.85)
    fun(); dev.off()
  }
  draw('01_mp_identity',1800,1150,function() {
    par(mfrow=c(2,3),mar=c(4,4,5,1))
    for(i in 1:5) {
      d <- readRDS(file.path(out,'data/member_C/spectra',paste0(regimes_C$regime[i],'_identity.rds')))
      vals <- as.vector(d$spectra[,,'sample']); c <- d$p/(d$n-1); e <- mp_edges_C(c)
      tol <- 1e-8*max(1,max(vals)); positive <- vals[vals>tol]
      # 密度按全体 p 个特征值归一化，而非对正特征值重新归一化。
      h <- hist(positive,breaks=60,plot=FALSE)
      h$density <- h$counts/(length(vals)*diff(h$breaks))
      xx <- seq(e[1]+1e-7,e[2]-1e-7,length.out=1800); yy <- mp_density_C(xx,c)
      plot(h,freq=FALSE,col='#dbe7ef',border='white',xlim=c(0,max(e[2],max(positive))*1.05),
           ylim=c(0,max(h$density,yy)*1.1),main='',xlab='Eigenvalue',ylab='Density (mass / all p)')
      title(main=sprintf('%s: n=%d, p=%d',d$regime,d$n,d$p),line=3,cex.main=1)
      lines(xx,yy,col='#0072b2',lwd=2); abline(v=e,lty=2,col='#e69f00')
      mtext(sprintf('c_df=%.3f; zero mass: %.3f (MP %.3f)',c,mean(vals<=tol),mp_atom_C(c)),side=3,line=1.25,cex=.7)
    }
    plot.new(); legend('center',legend=c('Pooled empirical continuous part','MP density (Gaussian identity)','MP continuous support edges','Zero atom shown as text, not density'),
      lty=c(1,1,2,NA),lwd=c(7,2,1,NA),col=c('#dbe7ef','#0072b2','#e69f00','black'),bty='n',cex=.9)
  })
  draw('02_spectra_four_structures',1800,1150,function() {
    par(mfrow=c(2,2),mar=c(4,5,3,1))
    for(st in structures_C) {
      d <- readRDS(file.path(out,'data/member_C/spectra',paste0('high_',st,'.rds')))
      means <- apply(d$spectra,c(2,3),mean)
      use <- c('sample','sphere_auto','diag_auto','sphere_100')
      matplot(1:d$p,means[,use],type='l',lty=1,lwd=2,col=cols[use],log='y',
        xlab='Eigenvalue rank (descending)',ylab='Eigenvalue (log scale)',main=paste(st,'(n=200, p=150)'))
      lines(1:d$p,d$truth,lwd=2,lty=2,col='black')
      legend('topright',legend=c('Population',labels[use]),col=c('black',cols[use]),lty=c(2,rep(1,4)),bty='n',cex=.7)
    }
  })
  draw('03_condition_numbers',1800,1250,function() {
    par(mfrow=c(2,2),mar=c(5,4,3,1))
    use <- c('sample','sphere_auto','diag_auto')
    for(st in structures_C) {
      d <- raw[raw$structure==st & raw$method%in%use,]
      groups <- lapply(seq_len(5),function(i) lapply(use,function(m) {
        x<-d$kappa[d$regime==regimes_C$regime[i]&d$method==m]; log10(x[is.finite(x)])
      })); groups <- unlist(groups,recursive=FALSE)
      names(groups)<-rep(c('S','SA','DA'),5)
      boxplot(groups,col=rep(cols[use],5),outline=FALSE,ylab='log10(condition number), finite only',main=st,xaxt='n')
      axis(1,at=seq(2,14,3),labels=c('low','mod','high','near','sing'))
      # 不把 Inf 当有限值画，明确标记最后一组样本协方差不可逆。
      u <- par('usr'); text(13,u[4]-.08*diff(u[3:4]),sprintf('S: Inf\n%d/%d',sum(d$method=='sample'&d$regime=='singular'),sum(d$method=='sample'&d$regime=='singular')),cex=.65,col='#64748b')
      legend('topleft',legend=labels[use],fill=cols[use],bty='n',cex=.7)
    }
  })
  draw('04_extreme_eigenvalues',1800,1600,function() {
    par(mfrow=c(4,2),mar=c(4,5,3,1),cex.axis=.85)
    use<-c('sample','sphere_auto','diag_auto')
    for(st in structures_C) for(metric in c('dmin','dmax')) {
      d<-raw[raw$structure==st&raw$method%in%use,]
      plot(0,0,type='n',xlim=c(1,5),ylim=range(d[[metric]]),xaxt='n',xlab='Regime',ylab=metric,main=st)
      axis(1,1:5,labels=c('L','M','H','N','S'))
      for(m in use) {
        avg<-se<-numeric(5)
        for(i in 1:5) {x<-d[[metric]][d$method==m&d$regime==regimes_C$regime[i]];avg[i]<-mean(x);se[i]<-sd(x)/sqrt(length(x))}
        lines(1:5,avg,col=cols[m],lwd=2,type='b',pch=16)
        segments(1:5,avg-qt(.975,length(x)-1)*se,1:5,avg+qt(.975,length(x)-1)*se,col=cols[m]); segments((1:5)-.035,avg-qt(.975,length(x)-1)*se,(1:5)+.035,avg-qt(.975,length(x)-1)*se,col=cols[m]); segments((1:5)-.035,avg+qt(.975,length(x)-1)*se,(1:5)+.035,avg+qt(.975,length(x)-1)*se,col=cols[m])
      }
      if(st=='identity') legend('topleft',legend=labels[use],col=cols[use],lty=1,bty='n',cex=.6)
    }
  })
  draw('05_risk_and_spectral_error',1800,1600,function() {
    par(mfrow=c(4,2),mar=c(4,5,3,1),cex.axis=.85)
    use<-c('sample','sphere_100','sphere_auto','diag_auto')
    for(st in structures_C) for(metric in c('frobenius','spectral')) {
      d<-raw[raw$structure==st&raw$method%in%use,]
      avg<-sapply(use,function(m) sapply(regimes_C$regime,function(r) mean(d[[metric]][d$method==m&d$regime==r])))
      matplot(1:5,avg,type='b',log='y',col=cols[use],lty=1,pch=16,xaxt='n',xlab='Regime',ylab=paste('Mean',metric),main=st)
      axis(1,1:5,labels=c('L','M','H','N','S'))
      if(st=='identity') legend('topleft',legend=labels[use],col=cols[use],lty=1,bty='n',cex=.58)
    }
  })
  draw('06_spiked_pca',1800,850,function() {
    par(mfrow=c(1,3),mar=c(4,4,3,1))
    d<-raw[raw$structure=='spiked',]
    for(metric in c('alignment1','alignment2','subspace_loss')) {
      boxplot(lapply(regimes_C$regime,function(r) d[[metric]][d$method=='sample'&d$regime==r]),
        names=c('low','mod','high','near','sing'),ylim=if(metric!='subspace_loss')c(0,1)else c(0,4),
        col='#dbe7ef',main=metric,ylab=if(metric=='subspace_loss')'Squared projector error'else'Squared cosine',xlab='Regime',outline=FALSE)
      if(metric!='subspace_loss') {
        alpha<-if(metric=='alignment1')6 else 3
        # Paul (2007) reference only plotted for c_df<1, matching stated theorem scope.
        refs<-sapply(1:4,function(i) pca_reference_C(alpha,regimes_C$p[i]/(regimes_C$n[i]-1))['alignment'])
        points(1:4,refs,pch=18,col='#e69f00',cex=1.4)
      }
      mtext('Sphere lambda<1 has the same directions; lambda=1 undefined',side=3,line=.2,cex=.58)
    }
  })
  draw('07_shrinkage_tradeoff',1500,650,function() {
    par(mfrow=c(1,2),mar=c(5,4,3,1))
    use<-methods_C[c(1,2,3,5,4)]
    for(metric in c('kappa','frobenius')) {
      d<-raw[raw$structure=='spiked'&raw$regime=='high',]
      boxplot(lapply(use,function(m)d[[metric]][d$method==m]),names=labels[use],col=cols[use],log='y',
        outline=FALSE,main=paste('Spiked, n=200 p=150:',metric),ylab=metric,cex.axis=.7)
    }
  })
  invisible(TRUE)
}
