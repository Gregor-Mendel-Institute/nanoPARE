# arg = commandArgs(trailingOnly=TRUE)
# metaplot_file = arg[1]
TSS_scores = list()
TES_scores = list()
for(samplename in c('fb1_1to1','fb2_1to1','fb3_1to1','1ng_fb_1','1ng_fb_2','1ng_fb_3','100pg_fb_1','100pg_fb_2','100pg_fb_3','10pg_fb_1','10pg_fb_2','10pg_fb_3')){
  metaplot_file = paste('L:/members/Schon/Writing/nanoPARE/data_generated/EndGraph/',samplename,'/TSS_metaplot_10.tab',sep='')
  readtypes = c()
  if(grepl('TSS',metaplot_file,fixed=T)){readtypes = append(readtypes,'TSS')}
  if(grepl('TES',metaplot_file,fixed=T)){readtypes = append(readtypes,'TES')}
  bookend_colors=c(BODY="#808285",TSS="#1C75BC",TES="#BE1E2D")
  
  auc=function(x)sum(x[which(x>0)])
  calculate_scalefactor=function(value_vector,maxpoint,direction='both'){
    if(direction=='both'){
      optimize(function(x){abs(auc(value_vector+x)-1)},lower = -maxpoint,upper = maxpoint)$minimum
    }else{
      if(direction==TRUE){
        optimize(function(x){abs(auc(value_vector+x)-1)},lower = 0,upper = maxpoint)$minimum
      }else{
        if(direction==FALSE){
          optimize(function(x){abs(auc(value_vector+x)-1)},lower = -maxpoint,upper = 0)$minimum
        }
      }
    }
  }
  
  metaplot_table=read.table(metaplot_file,header = F,stringsAsFactors = F)
  colnames(metaplot_table)=paste(rep(readtypes,each=2),rep(c('start','end'),length(readtypes)),sep='_')
  
  midpoint=floor(nrow(metaplot_table)/2)
  overall_maxpoint=max(abs(metaplot_table))
  
  pdf(paste(paste(readtypes,samplename,sep='_'),'_metaplot.pdf',sep=''),useDingbats = F)
  if(length(readtypes)==2){
  par(mfrow=c(2,2),mar=c(4,4,3,2))
  }else{
  par(mfrow=c(1,2),mar=c(4,4,3,2))
  }
  minuscolor = bookend_colors['BODY']
  for(readtype in readtypes){
  for(side in c('start','end')){
    cat(paste(samplename,side,auc(score),sep='\t'),'\n')
  pluscolor = bookend_colors[readtype]
  score=metaplot_table[,paste(readtype,side,sep='_')]
  if(side=='start'){
    plot(score,type='l',lwd=1,col='black',xaxt='n',xlab = "Position relative to TSS",ylab=paste(readtype,' reads - body reads',sep=''),frame.plot = F,ylim=c(-.025,.15))#ylim=c(-overall_maxpoint,overall_maxpoint))
    TSS_scores[[samplename]]=score
    polygon(x=c(midpoint+1,nrow(metaplot_table),nrow(metaplot_table),midpoint+1),y=c(1,1,-1,-1),col=adjustcolor('gold',alpha.f = .1),border = NA)
  }else{
    plot(score,type='l',lwd=1,col='black',xaxt='n',xlab = "Position relative to TES",ylab=paste(readtype,' reads - body reads',sep=''),frame.plot = F,ylim=c(-.025,.15))#,ylim=c(-overall_maxpoint,overall_maxpoint))
    TES_scores[[samplename]]=score
    polygon(x=c(1,midpoint,midpoint,1),y=c(1,1,-1,-1),col=adjustcolor('gold',alpha.f = .1),border = NA)
  }
  score_above=score
  score_above[score_above<0]=0
  score_below=score
  score_below[score_below>0]=0
  axis(1,c(1,midpoint,length(score)),labels = c(-midpoint,0,midpoint))
  polygon(x=c(1:length(score),length(score):1),y=c(score_below,rep(0,length(score))),col = minuscolor,border = NA)
  polygon(x=c(1:length(score),length(score):1),y=c(score_above,rep(0,length(score))),col = pluscolor,border = NA)
  abline(h=0,lty=1)
  abline(v=midpoint,lty=3)
  }
  }
  trash <- dev.off()
}


