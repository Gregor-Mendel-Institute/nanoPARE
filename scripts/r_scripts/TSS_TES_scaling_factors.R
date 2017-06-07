arg = commandArgs(trailingOnly=TRUE)
if(length(arg)==1){
  metaplot_file = arg[1]
}else{
  metaplot_file = 'TSS_TES_metaplot.tab'
}

readtypes = c()
if(grepl(arg[1],'TSS',fixed=T)){readtypes = append(readtypes,'TSS')}
if(grepl(arg[1],'TES',fixed=T)){readtypes = append(readtypes,'TES')}

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
if(ncol(metaplot_table)==4){
	colnames(metaplot_table)=c('TSS_start','TSS_end','TES_start','TES_end')
}else{
	colnames(metaplot_table)=c(paste(readtypes[1],'start',sep='_'),paste(readtypes[1],'end',sep='_'))
}
midpoint=floor(nrow(metaplot_table)/2)
overall_maxpoint=max(abs(metaplot_table))

if('TSS'%in%readtypes){
	opt_direction=1>auc(metaplot_table[,"TSS_start"])
	maxpoint=max(abs(metaplot_table[,"TSS_start"]))
	TSS_optimum = calculate_scalefactor(metaplot_table[,"TSS_start"],maxpoint,opt_direction)
	TSS_scaling_factor = auc(-metaplot_table[,"TSS_start"]-TSS_optimum) / auc(-metaplot_table[,"TSS_start"])
	TSS_bandwidth      = min(which(sapply(1:midpoint,function(x)auc(metaplot_table[(midpoint-x):(midpoint+x),"TSS_start"]+TSS_optimum))>=0.6826))
	cat(paste(c(TSS_scaling_factor,' ',TSS_bandwidth,' '),collapse=''))
}
if('TES'%in%readtypes){
	opt_direction=1>auc(metaplot_table[,"TES_end"])
	maxpoint=max(abs(metaplot_table[,"TES_end"]))
	TES_optimum = calculate_scalefactor(metaplot_table[,"TES_end"],maxpoint,opt_direction)
	TES_scaling_factor = auc(-metaplot_table[,"TES_end"]-TES_optimum) / auc(-metaplot_table[,"TES_end"])
	TES_bandwidth      = min(which(sapply(1:midpoint,function(x)auc(metaplot_table[(midpoint-x):(midpoint+x),"TES_end"]+TES_optimum))>=0.6826))
	cat(paste(c(TES_scaling_factor,' ',TES_bandwidth,' '),collapse=''))
}
