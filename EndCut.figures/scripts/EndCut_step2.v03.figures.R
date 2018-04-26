#EndCut_step2.v03.figures.R
#code used to generate Fig. 4 and associated figures in Schon/Kellner et al. (unpublished)
#input: *detectedSites.txt files
#output: various graphs as represented in Fig. 4 of Schon/Kellner et al. and associated supplemental files

biocLite("survcomp")
library(stringr)
library("pheatmap")
library("vioplot")

#set filepaths
setwd('/Users/michael.nodine/Desktop/R_Desktop/nanoPARE/180228/')
fileRoot = '/Users/michael.nodine/Desktop/R_Desktop/nanoPARE/180228/'
dataRoot = '/Users/michael.nodine/Desktop/R_Desktop/nanoPARE/180228/'
printRoot = '/Users/michael.nodine/Desktop/R_Desktop/nanoPARE/180228/graphs/'

################################################FUNCTIONS################################################
get.sites <- function(dir,sample,req.p.val) {
  #retrieves pre-compiled sites from tables
  sites = data.frame()
  sites = read.delim(paste(dataRoot,'results/',dir,'/',sample,'.detected.sites.tsv',sep=''),sep='\t', row.names=NULL, header=TRUE,stringsAsFactors=FALSE,strip.white=TRUE)
  
  #make another interaction column using sRNA and AGI (e.g. some AGIs have identical common names like SPL13, ARF, etc.)
  mir.AGI = c()
  
  for (row in rownames(sites)) {
    mir.AGI.sub = paste(sites[row,'sRNA'],':',strsplit(sites[row,'target.transcript'],"[.]")[[1]][1],sep='')
    mir.AGI = c(mir.AGI,mir.AGI.sub)
  }
  
  rownames(sites) = make.unique(mir.AGI)
  
  #select for significant ones
  sites.sig = subset(sites, adj.p.val < req.p.val & fold.change > 1.0 & slice.site.rpm > 0.1)
  return(list(sites,sites.sig))
}

get.dists.individual <- function(dir,sample,transcript.v,sRNA.v,allen.min,allen.max,flank.min,flank.max) {
  
  test = read.delim(paste(dataRoot,'results/distances/',dir,'/',sample,'.pred.sites.100nt.overlap.tsv',sep=''),
                    sep='\t', row.names=NULL, header=TRUE,stringsAsFactors=FALSE,strip.white=TRUE)
  
  #normalize values in test
  test.bg = read.delim(paste(dataRoot,'bedgraphs/',sample,'.W.transcript.capmasked.bedgraph',sep=''),
                       sep='\t', row.names=NULL, header=FALSE,stringsAsFactors=FALSE,strip.white=TRUE)
  
  total.reads = sum(test.bg$V4)
  
  test.norm = data.frame(test, 'Reads.norm'=test$Reads/(total.reads/1000000))
  
  test.sub1 = subset(test.norm, sRNA %in% sRNA.v & Transcript %in% transcript.v & Allen.Score >= allen.min & 
                       Allen.Score < allen.max & Diff >= flank.min & Diff <= flank.max, select=c('Diff','Reads.norm'))
  
  #tally up reads for each position
  test.dist.df = data.frame()
  
  for (i in flank.min:flank.max) {
    test.sub = subset(test.sub1, Diff==i, select=Reads.norm)
    i.sum = as.numeric(colSums(test.sub))
    test.dist.df = rbind(test.dist.df, data.frame('i.sum'=i.sum))
  }
  
  rownames(test.dist.df) = c(flank.min:flank.max)
  
  #subset for those with >0 reads
  test.dist.df.sub = subset(test.dist.df, i.sum > 0)
  
  return(list(test.dist.df,test.dist.df.sub))
}

print.example <- function(dir,transcript,sRNA,common,flank.min,flank.max) {
  sample.v=c('fb1_1to1','fb2_1to1','fb3_1to1','xrn1_fb_1','xrn1_fb_2','xrn1_fb_3','xrn4_M_1','xrn4_M_2','xrn4_M_3')
  allen.min=0
  allen.max=6.1
  
  dists.list.exp = list()
  
  for (sample in sample.v) {
    dists = get.dists.individual(dir,sample,c(transcript),c(sRNA),allen.min,allen.max,flank.min,flank.max)[[2]]
    dists.list.exp[[sample]] = dists
  }
  
  #plot just averages + SE to simplify
  get.avg.prop <- function(df.list) {
    df.merge = merge(df.list[[1]],df.list[[2]], by='row.names', all=TRUE)
    rownames(df.merge) = df.merge$Row.names
    df.merge.1 = subset(df.merge, select=c(i.sum.x,i.sum.y))
    df.merge.2 = merge(df.merge.1, df.list[[3]], by='row.names', all=TRUE)
    rownames(df.merge.2) = df.merge.2$Row.names
    df.merge.3 = subset(df.merge.2, select=c(i.sum.x,i.sum.y,i.sum))
    #replace NAs with zeros
    df.merge.3[is.na(df.merge.3)] <- 0
    df.avg = rowMeans(df.merge.3)
    df.avg.df = data.frame('avg.prop'=df.avg,'sd.prop'=apply(df.merge.3, 1, sd))
    return(df.avg.df)
  }
  
  fb.avg.df = get.avg.prop(list(dists.list.exp[[1]],dists.list.exp[[2]],dists.list.exp[[3]]))
  xrn1.avg.df = get.avg.prop(list(dists.list.exp[[4]],dists.list.exp[[5]],dists.list.exp[[6]]))
  xrn4.avg.df = get.avg.prop(list(dists.list.exp[[7]],dists.list.exp[[8]],dists.list.exp[[9]]))
  
  add.se <- function(avg.df) {
    for (row in rownames(avg.df)) {
      arrows(as.numeric(row), (avg.df[row,'avg.prop'] - avg.df[row,'sd.prop']), as.numeric(row), (avg.df[row,'avg.prop'] + avg.df[row,'sd.prop']),
             length=0.02, angle=90, code=3, lwd=1)
    }
  }
  
  plotFile = paste(printRoot,'/figures/distances/',dir,'/',sRNA,'.',common,'.distance.rpm.pdf',sep='')
  pdf(plotFile, useDingbats=FALSE, width=5, height=4)
  plot(as.numeric(rownames(fb.avg.df)),fb.avg.df$avg.prop, xlim=c(flank.min,flank.max), xlab="Distance from predicted cleavage site", ylab="Mean RPM", 
       las=1, col='#2171b5', ylim=c(0,max(c(fb.avg.df$avg.prop + fb.avg.df$sd.prop,xrn1.avg.df$avg.prop + xrn1.avg.df$sd.prop,
                                            xrn4.avg.df$avg.prop + xrn4.avg.df$sd.prop), na.rm=TRUE)),pch=19, main=paste(sRNA,':',common,sep=''))
  abline(v=0,lty=2)
  add.se(fb.avg.df)
  add.se(xrn1.avg.df)
  add.se(xrn4.avg.df)
  points(rownames(fb.avg.df),fb.avg.df$avg.prop,col='#2171b5',pch=19)
  points(rownames(xrn1.avg.df),xrn1.avg.df$avg.prop,col='#cb181d',pch=19)
  points(rownames(xrn4.avg.df),xrn4.avg.df$avg.prop,col='#238b45',pch=19)
  legend(x='topleft', legend=c('Col-0 (-Xrn1)','Col-0 (+Xrn1)','xrn4'), col=c('#2171b5','#cb181d','#238b45'), pch=c(19,19,19), pt.cex=1.2, bty='n')
  dev.off()

}

get.sRNAs.anno <- function(sRNA.type) {
  df = read.delim(paste(dataRoot,'Ath_annotations/anno.',sRNA.type,'.tsv',sep=''),sep='\t', row.names=1, header=TRUE,stringsAsFactors=FALSE,strip.white=TRUE)
  return(df)
}

get.union.2 <- function(df.list) {
  #returns vector of number of 1) union of all unique interactions, 2) interactions with FDRs < 0.05 in 2/3 reps and 3) FDRs < 0.1 in 2/3 reps
  
  unique.interactions = unique(c(rownames(df.list[[1]]),rownames(df.list[[2]]),rownames(df.list[[3]])))
  
  fdr_05.v = c()
  fdr_1.v = c()
  
  for (ui in unique.interactions) {
    pval1 = 1
    pval2 = 1
    pval3 = 1
    t05 = 0
    t1 = 0
    if (ui %in% rownames(df.list[[1]]) == TRUE) {pval1 = df.list[[1]][ui,'adj.p.val']}
    if (ui %in% rownames(df.list[[2]]) == TRUE) {pval2 = df.list[[2]][ui,'adj.p.val']}
    if (ui %in% rownames(df.list[[3]]) == TRUE) {pval3 = df.list[[3]][ui,'adj.p.val']}
    if (pval1 < 0.05) {t05 = t05 + 1}
    if (pval2 < 0.05) {t05 = t05 + 1}
    if (pval3 < 0.05) {t05 = t05 + 1}
    if (pval1 < 0.1) {t1 = t1 + 1}
    if (pval2 < 0.1) {t1 = t1 + 1}
    if (pval3 < 0.1) {t1 = t1 + 1}
    if (t05 >= 2) {fdr_05.v = c(fdr_05.v,ui)}
    if (t1 >= 2) {fdr_1.v = c(fdr_1.v,ui)}
  }
  return(list(unique.interactions,fdr_05.v,fdr_1.v))
}

get.sites <- function(dir,sample,req.p.val) {

  #retrieves pre-compiled sites from tables
  sites = data.frame()
  sites = read.delim(paste(dataRoot,'results/',dir,'/',sample,'.detected.sites.tsv',sep=''),sep='\t', row.names=NULL, header=TRUE,stringsAsFactors=FALSE,strip.white=TRUE)
  
  #make another interaction column using sRNA and AGI (e.g. some AGIs have identical common names like SPL13, ARF, etc.)
  mir.AGI = c()
  
  for (row in rownames(sites)) {
    mir.AGI.sub = paste(sites[row,'sRNA'],':',strsplit(sites[row,'target.transcript'],"[.]")[[1]][1],sep='')
    mir.AGI = c(mir.AGI,mir.AGI.sub)
  }
  
  rownames(sites) = make.unique(mir.AGI)
  
  #select for significant ones
  sites.sig = subset(sites, adj.p.val < req.p.val & fold.change > 1.0 & slice.site.rpm > 0.1)
  return(list(sites,sites.sig))
}

get.dists.select <- function(dir,sample,transcript.v,sRNA.v,allen.min,allen.max,flank.min,flank.max) {
  
  test = read.delim(paste(dataRoot,'results/distances/',dir,'/',sample,'.pred.sites.100nt.overlap.tsv',sep=''),
                    sep='\t', row.names=NULL, header=TRUE,stringsAsFactors=FALSE,strip.white=TRUE)
  
  test.sub1 = subset(test, sRNA %in% sRNA.v & Transcript %in% transcript.v & Allen.Score >= allen.min & 
                       Allen.Score <= allen.max & Diff >= flank.min & Diff <= flank.max, select=c('Diff','Reads'))
  
  #tally up reads for each position
  test.dist.df = data.frame()
  
  for (i in flank.min:flank.max) {
    test.sub = subset(test.sub1, Diff==i, select=Reads)
    i.sum=as.numeric(colMeans(test.sub))
    test.dist.df = rbind(test.dist.df, data.frame('i.sum'=i.sum))
  }
  
  rownames(test.dist.df) = c(flank.min:flank.max)
  
  #normalize values in test
  test.bg = read.delim(paste(dataRoot,'bedgraphs/',sample,'.W.transcript.capmasked.bedgraph',sep=''),
                       sep='\t', row.names=NULL, header=FALSE,stringsAsFactors=FALSE,strip.white=TRUE)
  
  total.reads = sum(test.bg$V4)
  
  #replace NAs with zeros
  test.dist.df[is.na(test.dist.df)] <- 0
  
  test.norm = data.frame(test.dist.df, 'Reads.norm'=test.dist.df$i.sum/(total.reads/1000000))
  
  return(test.norm)
}

get.all.site.proportions.transcript <- function(dir,sample,transcript.v,sRNA.v,allen.min,allen.max,flank.min,flank.max) {

  test = read.delim(paste(dataRoot,'results/distances/',dir,'/',sample,'.pred.sites.100nt.overlap.tsv',sep=''),
                    sep='\t', row.names=NULL, header=TRUE,stringsAsFactors=FALSE,strip.white=TRUE)
  
  test.sub1 = subset(test, sRNA %in% sRNA.v & Transcript %in% transcript.v & Allen.Score >= allen.min & 
                       Allen.Score <= allen.max & Diff >= flank.min & Diff <= flank.max, select=c('Transcript','sRNA','Diff','Reads'))
  
  #get total number of reads per uncapped transcript
  test.bg = read.delim(paste(dataRoot,'bedgraphs/',sample,'.W.transcript.unmasked.bedgraph',sep=''),
                       sep='\t', row.names=NULL, header=FALSE,stringsAsFactors=FALSE,strip.white=TRUE)
  
  #tally up reads for each position
  test.prop.df=data.frame(matrix(ncol=101,nrow=0))
  names(test.prop.df) = as.character(c(flank.min:flank.max))
  prop.name=c()
  
  for (i in 1:length(sRNA.v)) {
    transcript.reads = as.numeric(sum(as.numeric(subset(test.bg, V1==transcript.v[i])$V4)))
    prop.name=c(prop.name,paste(sRNA.v[i],':',transcript.v[i],sep=""))
    prop.v = c()    
    for (j in flank.min:flank.max) {
      prop = 0
      prop = as.numeric(colSums(subset(test.sub1, sRNA==sRNA.v[i] & Transcript==transcript.v[i] & Diff==j, select=Reads))/transcript.reads)
      prop.v = c(prop.v,prop)
    }
    test.prop.df = rbind(test.prop.df, prop.v)
  }
  
  names(test.prop.df) = as.character(c(flank.min:flank.max))
  rownames(test.prop.df) = prop.name  
  #replace NAs with zeros
  test.prop.df[is.na(test.prop.df)] <- 0
  
  return(test.prop.df)
}

get.all.site.proportions.transcript.pred <- function(dir,sample,transcript.v,sRNA.v,allen.min,allen.max,flank.min,flank.max) {

  test = read.delim(paste(dataRoot,'results/distances/',dir,'/',sample,'.pred.sites.100nt.overlap.tsv',sep=''),
                    sep='\t', row.names=NULL, header=TRUE,stringsAsFactors=FALSE,strip.white=TRUE)
  
  test.sub1 = subset(test, Transcript %in% transcript.v & sRNA %in% sRNA.v & Allen.Score >= allen.min & Allen.Score <= allen.max & Diff >= flank.min & Diff <= flank.max, 
                     select=c('Transcript','sRNA','Diff','Reads'))
  
  #get total number of reads per uncapped transcript
  test.bg = read.delim(paste(dataRoot,'bedgraphs/',sample,'.W.transcript.unmasked.bedgraph',sep=''), sep='\t', row.names=NULL, header=FALSE,stringsAsFactors=FALSE,strip.white=TRUE)
  test.bg.sub = subset(test.bg, V1 %in% transcript.v)
  
  #tally up number of reads per uncapped transcript
  transcript.v.unique = unique(test.bg.sub$V1)
  
  transcript.total.v=c()
  
  for (transcript in transcript.v.unique) {
    transcript.total=sum(subset(test.bg.sub, V1==transcript)$V4)
    transcript.total.v=c(transcript.total.v,transcript.total)
  }
  
  transcript.bg.df = data.frame("transcript.total"=transcript.total.v)
  rownames(transcript.bg.df) = transcript.v.unique
  
  #tally up reads for each position
  test.prop.df=data.frame()
  
  for (i in flank.min:flank.max) {
    i.prop=c()
    test.sub = subset(test.sub1, Diff==i, select=c(Transcript,Reads))
    for (row in rownames(test.sub)) {
      i.prop=c(i.prop,test.sub[row,"Reads"]/transcript.bg.df[test.sub[row,"Transcript"],"transcript.total"])
    }
    test.prop.df = rbind(test.prop.df, data.frame('prop'=mean(i.prop)))
  }
  
  rownames(test.prop.df) = c(flank.min:flank.max)
  
  #replace NAs with zeros
  test.prop.df[is.na(test.prop.df)] <- 0
  
  #transpose and change names
  test.prop.df = t(test.prop.df)
  rownames(test.prop.df) = "Predictions"
  
  return(test.prop.df)
}

#########Below shows how data were plotted Fig. 4 and associated supplementary figures#############
##first, retrieve all detected sRNA/target cleavage sites

req.p.val = 0.05

##miRNAs
###bias corrected and masked cap features (internal)
fb_1.c.w.anno.mir = get.sites('detected.sites.capmasked','fb_1.anno.mir',req.p.val)
fb_2.c.w.anno.mir = get.sites('detected.sites.capmasked','fb_2.anno.mir',req.p.val)
fb_3.c.w.anno.mir = get.sites('detected.sites.capmasked','fb_3.anno.mir',req.p.val)
xrn1_1.c.w.anno.mir = get.sites('detected.sites.capmasked','xrn1_1.anno.mir',req.p.val)
xrn1_2.c.w.anno.mir = get.sites('detected.sites.capmasked','xrn1_2.anno.mir',req.p.val)
xrn1_3.c.w.anno.mir = get.sites('detected.sites.capmasked','xrn1_3.anno.mir',req.p.val)
xrn4_1.c.w.anno.mir = get.sites('detected.sites.capmasked','xrn4_1.anno.mir',req.p.val)
xrn4_2.c.w.anno.mir = get.sites('detected.sites.capmasked','xrn4_2.anno.mir',req.p.val)
xrn4_3.c.w.anno.mir = get.sites('detected.sites.capmasked','xrn4_3.anno.mir',req.p.val)
d234_fb4.c.w.anno.mir = get.sites('detected.sites.capmasked','d234_fb4.anno.mir',req.p.val)
d234_fb5.c.w.anno.mir = get.sites('detected.sites.capmasked','d234_fb5.anno.mir',req.p.val)
d234_fb6.c.w.anno.mir = get.sites('detected.sites.capmasked','d234_fb6.anno.mir',req.p.val)
###just bias corrected (public)
GSM278334.u.w.anno.mir = get.sites('detected.sites.unmasked','GSM278334.anno.mir',req.p.val)
GSM278335.u.w.anno.mir = get.sites('detected.sites.unmasked','GSM278335.anno.mir',req.p.val)
GSM280226.u.w.anno.mir = get.sites('detected.sites.unmasked','GSM280226.anno.mir',req.p.val)
GSM280227.u.w.anno.mir = get.sites('detected.sites.unmasked','GSM280227.anno.mir',req.p.val)
SRX003076.u.w.anno.mir = get.sites('detected.sites.unmasked','SRX003076.anno.mir',req.p.val)
SRX003078.u.w.anno.mir = get.sites('detected.sites.unmasked','SRX003078.anno.mir',req.p.val)
SRX283759.u.w.anno.mir = get.sites('detected.sites.unmasked','SRX283759.anno.mir',req.p.val)
SRX283760.u.w.anno.mir = get.sites('detected.sites.unmasked','SRX283760.anno.mir',req.p.val)
SRX1559792.u.w.anno.mir = get.sites('detected.sites.unmasked','SRX1559792.anno.mir',req.p.val)
SRX1980216.u.w.anno.mir = get.sites('detected.sites.unmasked','SRX1980216.anno.mir',req.p.val)
SRX1980220.u.w.anno.mir = get.sites('detected.sites.unmasked','SRX1980220.anno.mir',req.p.val)
SRX1429178.u.w.anno.mir = get.sites('detected.sites.unmasked','SRX1429178.anno.mir',req.p.val)
SRX1429179.u.w.anno.mir = get.sites('detected.sites.unmasked','SRX1429179.anno.mir',req.p.val)
SRX384362.u.w.anno.mir = get.sites('detected.sites.unmasked','SRX384362.anno.mir',req.p.val)

##tasiRNAs
###bias corrected and masked cap features (internal)
fb_1.c.w.anno.tas = get.sites('detected.sites.capmasked','fb_1.anno.tas',req.p.val)
fb_2.c.w.anno.tas = get.sites('detected.sites.capmasked','fb_2.anno.tas',req.p.val)
fb_3.c.w.anno.tas = get.sites('detected.sites.capmasked','fb_3.anno.tas',req.p.val)
xrn1_1.c.w.anno.tas = get.sites('detected.sites.capmasked','xrn1_1.anno.tas',req.p.val)
xrn1_2.c.w.anno.tas = get.sites('detected.sites.capmasked','xrn1_2.anno.tas',req.p.val)
xrn1_3.c.w.anno.tas = get.sites('detected.sites.capmasked','xrn1_3.anno.tas',req.p.val)
xrn4_1.c.w.anno.tas = get.sites('detected.sites.capmasked','xrn4_1.anno.tas',req.p.val)
xrn4_2.c.w.anno.tas = get.sites('detected.sites.capmasked','xrn4_2.anno.tas',req.p.val)
xrn4_3.c.w.anno.tas = get.sites('detected.sites.capmasked','xrn4_3.anno.tas',req.p.val)
d234_fb4.c.w.anno.tas = get.sites('detected.sites.capmasked','d234_fb4.anno.tas',req.p.val)
d234_fb5.c.w.anno.tas = get.sites('detected.sites.capmasked','d234_fb5.anno.tas',req.p.val)
d234_fb6.c.w.anno.tas = get.sites('detected.sites.capmasked','d234_fb6.anno.tas',req.p.val)
###just bias corrected (public)
GSM278334.u.w.anno.tas = get.sites('detected.sites.unmasked','GSM278334.anno.tas',req.p.val)
GSM278335.u.w.anno.tas = get.sites('detected.sites.unmasked','GSM278335.anno.tas',req.p.val)
GSM280226.u.w.anno.tas = get.sites('detected.sites.unmasked','GSM280226.anno.tas',req.p.val)
GSM280227.u.w.anno.tas = get.sites('detected.sites.unmasked','GSM280227.anno.tas',req.p.val)
SRX003076.u.w.anno.tas = get.sites('detected.sites.unmasked','SRX003076.anno.tas',req.p.val)
SRX003078.u.w.anno.tas = get.sites('detected.sites.unmasked','SRX003078.anno.tas',req.p.val)
SRX283759.u.w.anno.tas = get.sites('detected.sites.unmasked','SRX283759.anno.tas',req.p.val)
SRX283760.u.w.anno.tas = get.sites('detected.sites.unmasked','SRX283760.anno.tas',req.p.val)
SRX1559792.u.w.anno.tas = get.sites('detected.sites.unmasked','SRX1559792.anno.tas',req.p.val)
SRX1980216.u.w.anno.tas = get.sites('detected.sites.unmasked','SRX1980216.anno.tas',req.p.val)
SRX1980220.u.w.anno.tas = get.sites('detected.sites.unmasked','SRX1980220.anno.tas',req.p.val)
SRX1429178.u.w.anno.tas = get.sites('detected.sites.unmasked','SRX1429178.anno.tas',req.p.val)
SRX1429179.u.w.anno.tas = get.sites('detected.sites.unmasked','SRX1429179.anno.tas',req.p.val)
SRX384362.u.w.anno.tas = get.sites('detected.sites.unmasked','SRX384362.anno.tas',req.p.val)


###############################################Fig. 4A#######################################################################

##miR173-5p targets
for (name in c("miR173-5p:AT2G27400","miR173-5p:AT2G39675","miR173-5p:AT2G39681")) {
  if (name %in% rownames(fb_1.c.w.anno.mir[[2]]) == TRUE) {
    print.example("transcript_bedgraph_capmasked",fb_1.c.w.anno.mir[[2]][name,"target.transcript"],fb_1.c.w.anno.mir[[2]][name,"sRNA"],
                  fb_1.c.w.anno.mir[[2]][name,"common.name"],-50,50)
  }else if (name %in% rownames(fb_2.c.w.anno.mir[[2]]) == TRUE) {
    print.example("transcript_bedgraph_capmasked",fb_2.c.w.anno.mir[[2]][name,"target.transcript"],fb_2.c.w.anno.mir[[2]][name,"sRNA"],
                  fb_1.c.w.anno.mir[[2]][name,"common.name"],-50,50)
  }
  
}

###############################################Fig. 4B#######################################################################

#Generate sRNA.v and transcript.v for all sites that were called significant in at least one Col-0 (-Xrn1) sample
fb.c.w.anno.mir.union = get.union.2(list(fb_1.c.w.anno.mir[[2]],fb_2.c.w.anno.mir[[2]],fb_3.c.w.anno.mir[[2]]))[[1]]

sRNA.v=c()
transcript.v=c()

for (name in sort(fb.c.w.anno.mir.union)) {
  if (name %in% rownames(fb_1.c.w.anno.mir[[2]]) == TRUE) {
    sRNA.v=c(sRNA.v,fb_1.c.w.anno.mir[[2]][name,"sRNA"])
    transcript.v=c(transcript.v,fb_1.c.w.anno.mir[[2]][name,"target.transcript"])
  }else if (name %in% rownames(fb_2.c.w.anno.mir[[2]]) == TRUE) {
    sRNA.v=c(sRNA.v,fb_2.c.w.anno.mir[[2]][name,"sRNA"])
    transcript.v=c(transcript.v,fb_2.c.w.anno.mir[[2]][name,"target.transcript"])
  }else if (name %in% rownames(fb_3.c.w.anno.mir[[2]]) == TRUE) {
    sRNA.v=c(sRNA.v,fb_3.c.w.anno.mir[[2]][name,"sRNA"])
    transcript.v=c(transcript.v,fb_3.c.w.anno.mir[[2]][name,"target.transcript"])
  }
}

dir='transcript_bedgraph_capmasked'
sample.v=c('fb1_1to1','fb2_1to1','fb3_1to1','xrn1_fb_1','xrn1_fb_2','xrn1_fb_3','xrn4_M_1','xrn4_M_2','xrn4_M_3')
allen.min=0
allen.max=6.1
flank.min = -50
flank.max = 50

dists.list.exp.all = list()

sample.v=c('fb1_1to1','fb2_1to1','fb3_1to1','xrn1_fb_1','xrn1_fb_2','xrn1_fb_3','xrn4_M_1','xrn4_M_2','xrn4_M_3')

for (sample in sample.v) {
  dists = get.all.site.proportions.transcript(dir,sample,transcript.v,sRNA.v,allen.min,allen.max,flank.min,flank.max)
  dists.list.exp.all[[sample]] = dists
} 

anno.miRNAs = get.sRNAs.anno('miRNA')
 
#get predicted
predictions = subset(read.delim('/Users/michael.nodine/Desktop/R_Desktop/nanoPARE/171021/GSTAr/anno.mir.tas.fa.GSTAr',sep='\t', row.names=NULL, header=TRUE, skip=7, stringsAsFactors=FALSE,strip.white=TRUE), AllenScore <= 6.0 & Query %in% rownames(anno.miRNAs), select=c(Transcript,Query))
 
dists.list.exp.pred.all = list()
 
for (sample in sample.v) {
  dists = get.all.site.proportions.transcript.pred(dir,sample,predictions$Transcript,predictions$Query,allen.min,allen.max,flank.min,flank.max)
  dists.list.exp.pred.all[[sample]] = dists
} 

fb.avg=colMeans(rbind(dists.list.exp.all[[1]],dists.list.exp.all[[2]],dists.list.exp.all[[3]])) * 100
fb.sd=apply(rbind(colMeans(dists.list.exp.all[[1]]),colMeans(dists.list.exp.all[[2]]),colMeans(dists.list.exp.all[[3]])), 2, sd) * 100
fb.se=fb.sd/sqrt(3)
xrn1.avg=colMeans(rbind(dists.list.exp.all[[4]],dists.list.exp.all[[5]],dists.list.exp.all[[4]])) * 100
xrn1.sd=apply(rbind(colMeans(dists.list.exp.all[[4]]),colMeans(dists.list.exp.all[[5]]),colMeans(dists.list.exp.all[[6]])), 2, sd) * 100
xrn1.se = xrn1.sd/sqrt(3)
xrn4.avg=colMeans(rbind(dists.list.exp.all[[7]],dists.list.exp.all[[8]],dists.list.exp.all[[9]])) * 100
xrn4.sd=apply(rbind(colMeans(dists.list.exp.all[[7]]),colMeans(dists.list.exp.all[[8]]),colMeans(dists.list.exp.all[[9]])), 2, sd) * 100
xrn4.se = xrn4.sd/sqrt(3)

fb.pred.avg=colMeans(rbind(dists.list.exp.pred.all[[1]],dists.list.exp.pred.all[[2]],dists.list.exp.pred.all[[3]])) * 100
xrn1.pred.avg=colMeans(rbind(dists.list.exp.pred.all[[4]],dists.list.exp.pred.all[[5]],dists.list.exp.pred.all[[6]])) * 100
xrn4.pred.avg=colMeans(rbind(dists.list.exp.pred.all[[7]],dists.list.exp.pred.all[[8]],dists.list.exp.pred.all[[9]])) * 100

add.se.transcript <- function(avg.df,se.df,mp) {
  for (i in c(1:length(avg.df))) {
    if (as.numeric(avg.df[i]) - as.numeric(se.df[i]) > 0) {
      arrows(mp[i], as.numeric(avg.df[i] - se.df[i]), mp[i], as.numeric(avg.df[i] + se.df[i]),
             length=0.015, angle=90, code=3, lwd=0.5)      
    }
  }
}

#determine p-values for on-target cleavage between conditions
##generate vectors of number of read 5' ends mapping to expected clevage position
fb.cleavage.site.v = as.numeric(rowMeans(rbind(dists.list.exp.all[[1]]["0"],dists.list.exp.all[[2]]["0"],dists.list.exp.all[[3]]["0"])))
xrn1.cleavage.site.v = as.numeric(rowMeans(rbind(dists.list.exp.all[[4]]["0"],dists.list.exp.all[[5]]["0"],dists.list.exp.all[[6]]["0"])))
xrn4.cleavage.site.v = as.numeric(rowMeans(rbind(dists.list.exp.all[[7]]["0"],dists.list.exp.all[[8]]["0"],dists.list.exp.all[[9]]["0"])))

wt.vs.xrn1.pval=ks.test(fb.cleavage.site.v,xrn1.cleavage.site.v,alternative="less")$p.value
wt.vs.xrn4.pval=ks.test(fb.cleavage.site.v,xrn4.cleavage.site.v,alternative="greater")$p.value

#determine fold-changes (for text)
mean(fb.cleavage.site.v)/mean(xrn1.cleavage.site.v)   #9.73648
mean(xrn4.cleavage.site.v)/mean(fb.cleavage.site.v)   #1.841652

#generate label key
fb.xrn1.sym=""
fb.xrn4.sym=""

if (wt.vs.xrn1.pval < 0.001) {
  fb.xrn1.sym="***"
}else if (wt.vs.xrn1.pval < 0.01) {
  fb.xrn1.sym="**"
}else if (wt.vs.xrn1.pval < 0.05) {
  fb.xrn1.sym="*"
}

if (wt.vs.xrn4.pval < 0.001) {
  fb.xrn4.sym="***"
}else if (wt.vs.xrn4.pval < 0.01) {
  fb.xrn4.sym="**"
}else if (wt.vs.xrn4.pval < 0.05) {
  fb.xrn4.sym="*"
}

col.v = c('#2171b5','#cb181d','#238b45')

plotFile = paste(printRoot,'/figures/fb.c.w.union.distance.vs.transcript.proportions.pdf',sep='')
pdf(plotFile, useDingbats=FALSE, width=2.753145, height=4.306425)
par(mfrow=c(3,1))
par(mar=c(0,4,2,2))
ymax=max(c(fb.avg+fb.se,xrn1.avg+xrn1.se,xrn4.avg+xrn4.se)) + 0.5
line.width=1
##fb and predictions
mp=barplot(fb.avg, xlab=NULL, ylab="", main="", beside=FALSE, col="#2171b5", las=1, ylim=c(0,ymax), xaxt="n", border=NA)
#abline(h=0,lwd=1)
axis(1, at=c(mp[1],mp[26],mp[51],mp[76],mp[101]), labels=c("","","","",""), padj=-0.25)
#add se bars
add.se.transcript(fb.avg,fb.se,mp)
#add trend line for predictions
par(new=FALSE)
lines(mp,as.numeric(fb.pred.avg), col="black", lty=1, lwd=1)
#add legend
text(x=mp[90],y=12, labels=c("Col-0\n(-Xrn1)"), pos=1)
##xrn1 and predictions
par(mar=c(1,4,1,2))
mp=barplot(xrn1.avg, xlab=NULL, ylab="Mean read proportion", main="", beside=FALSE, col="#cb181d", las=1, ylim=c(0,ymax), xaxt="n", border=NA)
axis(1, at=c(mp[1],mp[26],mp[51],mp[76],mp[101]), labels=c("","","","",""), padj=-0.25)
#add se bars
add.se.transcript(xrn1.avg,xrn1.se,mp)
#add trend line for predictions
par(new=FALSE)
lines(mp,as.numeric(xrn1.pred.avg), col="black", lty=1, lwd=1)
#add legend
text(x=mp[90],y=12, labels=c("Col-0\n(+Xrn1)"), pos=1) #, adj=1)
#add asterisks denoting significance
text(mp[51],xrn1.avg["0"] + xrn1.se["0"] - 0.8,labels=fb.xrn1.sym,pos=3,cex=1)
##xrn4 and predictions
par(mar=c(2,4,0,2))
mp=barplot(xrn4.avg, xlab=NULL, ylab="", main="", beside=FALSE, col="#238b45", las=1, ylim=c(0,ymax), xaxt="n", border=NA)
#add se bars
add.se.transcript(xrn4.avg,xrn4.se,mp)
#add axis
axis(1, at=c(mp[1],mp[26],mp[51],mp[76],mp[101]), labels=c(-50,-25,0,25,50), padj=-0.25)
#add trend line for predictions
par(new=FALSE)
lines(mp,as.numeric(xrn4.pred.avg), col="black", lty=1, lwd=1)
#add legend
text(x=mp[90],y=12, labels=c("xrn4"), pos=1) #, adj=1)
#add asterisks denoting significance
text(mp[51],xrn4.avg["0"] + xrn4.se["0"] - 0.8,labels=fb.xrn4.sym,pos=3,cex=1, xpd="TRUE")
dev.off()

#########################Supplemental Figure. NanoPARE read percentages along predicted targets (related to Fig. 4B)#########################
#plot 5' end proportions for predictions using various Allen scores (all [0-6], low [0-2], medium [2-4] and high [4-6])
#Generate sRNA.v and transcript.v for all sites that were called significant in at least one Col-0 (-Xrn1) sample
anno.miRNAs = get.sRNAs.anno('miRNA')
anno.tasiRNAs = get.sRNAs.anno('tasiRNA')

dir='transcript_bedgraph_capmasked'
sample.v=c('fb1_1to1','fb2_1to1','fb3_1to1','xrn1_fb_1','xrn1_fb_2','xrn1_fb_3','xrn4_M_1','xrn4_M_2','xrn4_M_3')

#get predicted
predictions = subset(read.delim('/Users/michael.nodine/Desktop/R_Desktop/nanoPARE/171021/GSTAr/anno.mir.tas.fa.GSTAr',sep='\t', row.names=NULL, header=TRUE, skip=7, stringsAsFactors=FALSE,strip.white=TRUE), AllenScore <= 6.0 & Query %in% rownames(anno.miRNAs), select=c(Transcript,Query))

##all
dists.list.exp.pred.all = list()

for (sample in sample.v) {
  dists = get.all.site.proportions.transcript.pred(dir,sample,predictions$Transcript,predictions$Query,0,6,-50,50)
  dists.list.exp.pred.all[[sample]] = dists
} 

##high
dists.list.exp.pred.high = list()

for (sample in sample.v) {
  dists = get.all.site.proportions.transcript.pred(dir,sample,predictions$Transcript,predictions$Query,4.01,6,-50,50)
  dists.list.exp.pred.high[[sample]] = dists
} 

##medium
dists.list.exp.pred.medium = list()

for (sample in sample.v) {
  dists = get.all.site.proportions.transcript.pred(dir,sample,predictions$Transcript,predictions$Query,2.01,4,-50,50)
  dists.list.exp.pred.medium[[sample]] = dists
} 

##low
dists.list.exp.pred.low = list()

for (sample in sample.v) {
  dists = get.all.site.proportions.transcript.pred(dir,sample,predictions$Transcript,predictions$Query,0,2,-50,50)
  dists.list.exp.pred.low[[sample]] = dists
} 


#get means and se for each pred set
##all
fb.pred.all.avg=colMeans(rbind(dists.list.exp.pred.all[[1]],dists.list.exp.pred.all[[2]],dists.list.exp.pred.all[[3]])) * 100
fb.pred.all.sd=apply(rbind(colMeans(dists.list.exp.pred.all[[1]]),colMeans(dists.list.exp.pred.all[[2]]),colMeans(dists.list.exp.pred.all[[3]])), 2, sd) * 100
fb.pred.all.se=fb.pred.all.sd/sqrt(3)
xrn1.pred.all.avg=colMeans(rbind(dists.list.exp.pred.all[[4]],dists.list.exp.pred.all[[5]],dists.list.exp.pred.all[[4]])) * 100
xrn1.pred.all.sd=apply(rbind(colMeans(dists.list.exp.pred.all[[4]]),colMeans(dists.list.exp.pred.all[[5]]),colMeans(dists.list.exp.pred.all[[6]])), 2, sd) * 100
xrn1.pred.all.se = xrn1.pred.all.sd/sqrt(3)
xrn4.pred.all.avg=colMeans(rbind(dists.list.exp.pred.all[[7]],dists.list.exp.pred.all[[8]],dists.list.exp.pred.all[[9]])) * 100
xrn4.pred.all.sd=apply(rbind(colMeans(dists.list.exp.pred.all[[7]]),colMeans(dists.list.exp.pred.all[[8]]),colMeans(dists.list.exp.pred.all[[9]])), 2, sd) * 100
xrn4.pred.all.se = xrn4.pred.all.sd/sqrt(3)
##high
fb.pred.high.avg=colMeans(rbind(dists.list.exp.pred.high[[1]],dists.list.exp.pred.high[[2]],dists.list.exp.pred.high[[3]])) * 100
fb.pred.high.sd=apply(rbind(colMeans(dists.list.exp.pred.high[[1]]),colMeans(dists.list.exp.pred.high[[2]]),colMeans(dists.list.exp.pred.high[[3]])), 2, sd) * 100
fb.pred.high.se=fb.pred.high.sd/sqrt(3)
xrn1.pred.high.avg=colMeans(rbind(dists.list.exp.pred.high[[4]],dists.list.exp.pred.high[[5]],dists.list.exp.pred.high[[4]])) * 100
xrn1.pred.high.sd=apply(rbind(colMeans(dists.list.exp.pred.high[[4]]),colMeans(dists.list.exp.pred.high[[5]]),colMeans(dists.list.exp.pred.high[[6]])), 2, sd) * 100
xrn1.pred.high.se = xrn1.pred.high.sd/sqrt(3)
xrn4.pred.high.avg=colMeans(rbind(dists.list.exp.pred.high[[7]],dists.list.exp.pred.high[[8]],dists.list.exp.pred.high[[9]])) * 100
xrn4.pred.high.sd=apply(rbind(colMeans(dists.list.exp.pred.high[[7]]),colMeans(dists.list.exp.pred.high[[8]]),colMeans(dists.list.exp.pred.high[[9]])), 2, sd) * 100
xrn4.pred.high.se = xrn4.pred.high.sd/sqrt(3)
##medium
fb.pred.medium.avg=colMeans(rbind(dists.list.exp.pred.medium[[1]],dists.list.exp.pred.medium[[2]],dists.list.exp.pred.medium[[3]])) * 100
fb.pred.medium.sd=apply(rbind(colMeans(dists.list.exp.pred.medium[[1]]),colMeans(dists.list.exp.pred.medium[[2]]),colMeans(dists.list.exp.pred.medium[[3]])), 2, sd) * 100
fb.pred.medium.se=fb.pred.medium.sd/sqrt(3)
xrn1.pred.medium.avg=colMeans(rbind(dists.list.exp.pred.medium[[4]],dists.list.exp.pred.medium[[5]],dists.list.exp.pred.medium[[4]])) * 100
xrn1.pred.medium.sd=apply(rbind(colMeans(dists.list.exp.pred.medium[[4]]),colMeans(dists.list.exp.pred.medium[[5]]),colMeans(dists.list.exp.pred.medium[[6]])), 2, sd) * 100
xrn1.pred.medium.se = xrn1.pred.medium.sd/sqrt(3)
xrn4.pred.medium.avg=colMeans(rbind(dists.list.exp.pred.medium[[7]],dists.list.exp.pred.medium[[8]],dists.list.exp.pred.medium[[9]])) * 100
xrn4.pred.medium.sd=apply(rbind(colMeans(dists.list.exp.pred.medium[[7]]),colMeans(dists.list.exp.pred.medium[[8]]),colMeans(dists.list.exp.pred.medium[[9]])), 2, sd) * 100
xrn4.pred.medium.se = xrn4.pred.medium.sd/sqrt(3)
##low
fb.pred.low.avg=colMeans(rbind(dists.list.exp.pred.low[[1]],dists.list.exp.pred.low[[2]],dists.list.exp.pred.low[[3]])) * 100
fb.pred.low.sd=apply(rbind(colMeans(dists.list.exp.pred.low[[1]]),colMeans(dists.list.exp.pred.low[[2]]),colMeans(dists.list.exp.pred.low[[3]])), 2, sd) * 100
fb.pred.low.se=fb.pred.low.sd/sqrt(3)
xrn1.pred.low.avg=colMeans(rbind(dists.list.exp.pred.low[[4]],dists.list.exp.pred.low[[5]],dists.list.exp.pred.low[[4]])) * 100
xrn1.pred.low.sd=apply(rbind(colMeans(dists.list.exp.pred.low[[4]]),colMeans(dists.list.exp.pred.low[[5]]),colMeans(dists.list.exp.pred.low[[6]])), 2, sd) * 100
xrn1.pred.low.se = xrn1.pred.low.sd/sqrt(3)
xrn4.pred.low.avg=colMeans(rbind(dists.list.exp.pred.low[[7]],dists.list.exp.pred.low[[8]],dists.list.exp.pred.low[[9]])) * 100
xrn4.pred.low.sd=apply(rbind(colMeans(dists.list.exp.pred.low[[7]]),colMeans(dists.list.exp.pred.low[[8]]),colMeans(dists.list.exp.pred.low[[9]])), 2, sd) * 100
xrn4.pred.low.se = xrn4.pred.low.sd/sqrt(3)

add.se.transcript <- function(avg.df,se.df,mp) {
  for (i in c(1:length(avg.df))) {
    if (as.numeric(avg.df[i]) - as.numeric(se.df[i]) > 0) {
      arrows(mp[i], as.numeric(avg.df[i] - se.df[i]), mp[i], as.numeric(avg.df[i] + se.df[i]),
             length=0.015, angle=90, code=3, lwd=0.5)      
    }
  }
}

#generate graph
col.v = c('#2171b5','#cb181d','#238b45')

plotFile = paste(printRoot,'/figures/supplemental/fb.c.w.union.distance.vs.transcript.proportions.predictions.pdf',sep='')
pdf(plotFile, useDingbats=FALSE, width=10, height=5)
par(mfcol=c(3,4))
ymax=max(c(fb.pred.all.avg+fb.pred.all.se,xrn4.pred.all.avg+xrn4.pred.all.se,
           fb.pred.high.avg+fb.pred.high.se,xrn4.pred.high.avg+xrn4.pred.high.se,
           fb.pred.medium.avg+fb.pred.medium.se,xrn4.pred.medium.avg+xrn4.pred.medium.se,
           fb.pred.low.avg+fb.pred.low.se,xrn4.pred.low.avg+xrn4.pred.low.se)) + 0.5

line.width=1

print.sub <- function(df.avg,df.se,print.col,df.name,mtext.name) {

  mp=barplot(df.avg, xlab=NULL, ylab="Mean read proportion", main="", beside=FALSE, col=print.col, las=1, ylim=c(0,ymax), 
             xaxt="n", yaxt="n", border=NA,xpd=FALSE)
  #add se bars
  add.se.transcript(df.avg,df.se,mp)
  #add legend
  text(x=mp[90],y=ymax-1, labels=c(df.name), pos=1) #, adj=1)
  #add axis
  axis(1, at=c(mp[1],mp[26],mp[51],mp[76],mp[101]), labels=c(-50,-25,0,25,50), padj=-0.25)
  axis(2, labels=TRUE, las=1)
  if (mtext.name != "") {
    mtext(mtext.name,3,0)
  }
  
  
}

##all
par(mar=c(0,4,2,2))
print.sub(fb.pred.all.avg,fb.pred.all.se,"#2171b5","Col-0\n(-Xrn1)","0-6 Allen scores")
par(mar=c(1,4,1,2))
print.sub(xrn1.pred.all.avg,xrn1.pred.all.se,"#cb181d","Col-0\n(+Xrn1)","")
par(mar=c(2,4,0,2))
print.sub(xrn4.pred.all.avg,xrn4.pred.all.se,"#238b45","xrn4","")
##high
par(mar=c(0,4,2,2))
print.sub(fb.pred.high.avg,fb.pred.high.se,"#2171b5","Col-0\n(-Xrn1)","4-6 Allen scores")
par(mar=c(1,4,1,2))
print.sub(xrn1.pred.high.avg,xrn1.pred.high.se,"#cb181d","Col-0\n(+Xrn1)","")
par(mar=c(2,4,0,2))
print.sub(xrn4.pred.high.avg,xrn4.pred.high.se,"#238b45","xrn4","")
##medium
par(mar=c(0,4,2,2))
print.sub(fb.pred.medium.avg,fb.pred.medium.se,"#2171b5","Col-0\n(-Xrn1)","2-4 Allen scores")
par(mar=c(1,4,1,2))
print.sub(xrn1.pred.medium.avg,xrn1.pred.medium.se,"#cb181d","Col-0\n(+Xrn1)","")
par(mar=c(2,4,0,2))
print.sub(xrn4.pred.medium.avg,xrn4.pred.medium.se,"#238b45","xrn4","")
##low
par(mar=c(0,4,2,2))
print.sub(fb.pred.low.avg,fb.pred.low.se,"#2171b5","Col-0\n(-Xrn1)","0-2 Allen scores")
par(mar=c(1,4,1,2))
print.sub(xrn1.pred.low.avg,xrn1.pred.low.se,"#cb181d","Col-0\n(+Xrn1)","")
par(mar=c(2,4,0,2))
print.sub(xrn4.pred.low.avg,xrn4.pred.low.se,"#238b45","xrn4","")
dev.off()

#####################################Figs. 4C/D#####################################
##this was output from EndCut_step2.v03.git.R

#####################################Figs. 4E/F#####################################
##get sites detected in at least one biorep from each condition (unions)
#miRNAs
fb.c.w.anno.mir.union = get.union.2(list(fb_1.c.w.anno.mir[[2]],fb_2.c.w.anno.mir[[2]],fb_3.c.w.anno.mir[[2]]))[[1]]
xrn1.c.w.anno.mir.union = get.union.2(list(xrn1_1.c.w.anno.mir[[2]],xrn1_2.c.w.anno.mir[[2]],xrn1_3.c.w.anno.mir[[2]]))[[1]]
xrn4.c.w.anno.mir.union = get.union.2(list(xrn4_1.c.w.anno.mir[[2]],xrn4_2.c.w.anno.mir[[2]],xrn4_3.c.w.anno.mir[[2]]))[[1]]
d234.c.w.anno.mir.union = get.union.2(list(d234_fb4.c.w.anno.mir[[2]],d234_fb5.c.w.anno.mir[[2]],d234_fb6.c.w.anno.mir[[2]]))[[1]]
#tasiRNAs
fb.c.w.anno.tas.union = get.union.2(list(fb_1.c.w.anno.tas[[2]],fb_2.c.w.anno.tas[[2]],fb_3.c.w.anno.tas[[2]]))[[1]]
xrn1.c.w.anno.tas.union = get.union.2(list(xrn1_1.c.w.anno.tas[[2]],xrn1_2.c.w.anno.tas[[2]],xrn1_3.c.w.anno.tas[[2]]))[[1]]
xrn4.c.w.anno.tas.union = get.union.2(list(xrn4_1.c.w.anno.tas[[2]],xrn4_2.c.w.anno.tas[[2]],xrn4_3.c.w.anno.tas[[2]]))[[1]]
d234.c.w.anno.tas.union = get.union.2(list(d234_fb4.c.w.anno.tas[[2]],d234_fb5.c.w.anno.tas[[2]],d234_fb6.c.w.anno.tas[[2]]))[[1]]
##get sites detected in at least 2/3 bioreps from each condition (high conf. set)
##miRNAs
fb.c.w.anno.mir.set = get.union.2(list(fb_1.c.w.anno.mir[[2]],fb_2.c.w.anno.mir[[2]],fb_3.c.w.anno.mir[[2]]))[[2]]
xrn1.c.w.anno.mir.set = get.union.2(list(xrn1_1.c.w.anno.mir[[2]],xrn1_2.c.w.anno.mir[[2]],xrn1_3.c.w.anno.mir[[2]]))[[2]]
xrn4.c.w.anno.mir.set = get.union.2(list(xrn4_1.c.w.anno.mir[[2]],xrn4_2.c.w.anno.mir[[2]],xrn4_3.c.w.anno.mir[[2]]))[[2]]
d234.c.w.anno.mir.set = get.union.2(list(d234_fb4.c.w.anno.mir[[2]],d234_fb5.c.w.anno.mir[[2]],d234_fb6.c.w.anno.mir[[2]]))[[2]]
##tasiRNAs
fb.c.w.anno.tas.set = get.union.2(list(fb_1.c.w.anno.tas[[2]],fb_2.c.w.anno.tas[[2]],fb_3.c.w.anno.tas[[2]]))[[2]]
xrn1.c.w.anno.tas.set = get.union.2(list(xrn1_1.c.w.anno.tas[[2]],xrn1_2.c.w.anno.tas[[2]],xrn1_3.c.w.anno.tas[[2]]))[[2]]
xrn4.c.w.anno.tas.set = get.union.2(list(xrn4_1.c.w.anno.tas[[2]],xrn4_2.c.w.anno.tas[[2]],xrn4_3.c.w.anno.tas[[2]]))[[2]]
d234.c.w.anno.tas.set = get.union.2(list(d234_fb4.c.w.anno.tas[[2]],d234_fb5.c.w.anno.tas[[2]],d234_fb6.c.w.anno.tas[[2]]))[[2]]

#Fig. 4E
plotFile=paste(printRoot,'/figures/fig.4e.anno.miRNA.c.w.strip.pdf',sep='')
pdf(plotFile, useDingbats=FALSE, width=4, height=5)
ymax=90
stripchart(list(c(length(fb_1.c.w.anno.mir[[2]][,1]),length(fb_2.c.w.anno.mir[[2]][,1]),length(fb_3.c.w.anno.mir[[2]][,1])),
                c(length(xrn1_1.c.w.anno.mir[[2]][,1]),length(xrn1_2.c.w.anno.mir[[2]][,1]),length(xrn1_3.c.w.anno.mir[[2]][,1])),
                c(length(xrn4_1.c.w.anno.mir[[2]][,1]),length(xrn4_2.c.w.anno.mir[[2]][,1]),length(xrn4_3.c.w.anno.mir[[2]][,1])),
                c(length(d234_fb4.c.w.anno.mir[[2]][,1]),length(d234_fb5.c.w.anno.mir[[2]][,1]),length(d234_fb6.c.w.anno.mir[[2]][,1]))),
           method='jitter', jitter=0.2, main=NULL, vertical=TRUE, las=1, pch=1, cex=1.2, col='black',
           group.names=c('Col-0\n-Xrn1','Col-0\n+Xrn1','xrn4','dcl234'), xlab=NULL, ylab='Number of miRNA target sites (p.adj < 0.05)',ylim=c(0,ymax), cex.axis=1.0)
abline(h=0, lty=2)
#high-confidence set (2/3)
stripchart(list(c(length(fb.c.w.anno.mir.set)),c(length(xrn1.c.w.anno.mir.set)),c(length(xrn4.c.w.anno.mir.set)),c(length(d234.c.w.anno.mir.set))),
           method='jitter', jitter=0.0, vertical=TRUE, pch=19, las=1, col='#b2182b', cex=1.2, xlab=NULL, ylim=c(0,ymax), cex.axis=1.0, add=TRUE)
#union
stripchart(list(c(length(fb.c.w.anno.mir.union)),c(length(xrn1.c.w.anno.mir.union)),c(length(xrn4.c.w.anno.mir.union)),c(length(d234.c.w.anno.mir.union))),
           method='jitter', jitter=0.0, vertical=TRUE, pch=17, las=1, col='#b2182b', cex=1.2, xlab=NULL, ylim=c(0,ymax), cex.axis=1.0, add=TRUE)
legend(x='topleft', legend=c('Bioreps','Union','High Conf.'), col=c('black','#b2182b','#b2182b'), pch=c(1,17,19), pt.cex=1.2, bty='n')
dev.off()

#Fig. 4F
plotFile=paste(printRoot,'/figures/fig.4f.anno.tasiRNA.c.w.strip.pdf',sep='')
pdf(plotFile, useDingbats=FALSE, width=4, height=5)
ymax=18
stripchart(list(c(length(fb_1.c.w.anno.tas[[2]][,1]),length(fb_2.c.w.anno.tas[[2]][,1]),length(fb_3.c.w.anno.tas[[2]][,1])),
                c(length(xrn1_1.c.w.anno.tas[[2]][,1]),length(xrn1_2.c.w.anno.tas[[2]][,1]),length(xrn1_3.c.w.anno.tas[[2]][,1])),
                c(length(xrn4_1.c.w.anno.tas[[2]][,1]),length(xrn4_2.c.w.anno.tas[[2]][,1]),length(xrn4_3.c.w.anno.tas[[2]][,1])),
                c(length(d234_fb4.c.w.anno.tas[[2]][,1]),length(d234_fb5.c.w.anno.tas[[2]][,1]),length(d234_fb6.c.w.anno.tas[[2]][,1]))),
           method='jitter', jitter=0.2, main=NULL, vertical=TRUE, las=1, pch=1, cex=1.2, col='black',
           group.names=c('Col-0\n-Xrn1','Col-0\n+Xrn1','xrn4','dcl234'), xlab=NULL, ylab='Number of tasiRNA target sites (p.adj < 0.05)',ylim=c(0,ymax), cex.axis=1.0)
abline(h=0, lty=2)
#high-confidence set (2/3)
stripchart(list(c(length(fb.c.w.anno.tas.set)),c(length(xrn1.c.w.anno.tas.set)),c(length(xrn4.c.w.anno.tas.set)),c(length(d234.c.w.anno.tas.set))),
           method='jitter', jitter=0.0, vertical=TRUE, pch=19, las=1, col='#2166ac', cex=1.2, xlab=NULL, ylim=c(0,ymax), cex.axis=1.0, add=TRUE)
#union
stripchart(list(c(length(fb.c.w.anno.tas.union)),c(length(xrn1.c.w.anno.tas.union)),c(length(xrn4.c.w.anno.tas.union)),c(length(d234.c.w.anno.tas.union))),
           method='jitter', jitter=0.0, vertical=TRUE, pch=17, las=1, col='#2166ac', cex=1.2, xlab=NULL, ylim=c(0,ymax), cex.axis=1.0, add=TRUE)
legend(x='topleft', legend=c('Bioreps','Union','High Conf.'), col=c('black','#2166ac','#2166ac'), pch=c(1,17,19), pt.cex=1.2, bty='n')
dev.off()

######FOR TEXT#######
##miRNAs
#determine fold-changes and p-values
#Col-0 (-Xrn1) vs. Col-0 (+Xrn1)
mean(c(35,28,35))/mean(c(3,5,6))                                  #7
t.test(c(35,28,35),c(3,5,6), alternative = "greater")$p.value     #0.001515833
#xrn4 vs. Col-0 (-Xrn1)
mean(c(47,66,45))/mean(c(35,28,35))                               #1.612245
t.test(c(35,28,35),c(47,66,45), alternative = "less")$p.value     #0.04152325

##tasiRNAs
#determine fold-changes and p-values
#Col-0 (-Xrn1) vs. Col-0 (+Xrn1)
mean(c(4,2,1))/mean(c(0,1,0))                               #7
t.test(c(4,2,1),c(0,1,0), alternative = "greater")$p.value  #0.0697175
#xrn4 vs. Col-0 (-Xrn1)
mean(c(8,9,6))/mean(c(4,2,1))                               #3.285714
t.test(c(4,2,1),c(8,9,6), alternative = "less")$p.value     #0.006443344


###################################################Fig. 4G###############################################
##miRNAs
#get list of dfs for all interactions detected
df.inList.total = list(fb_1.c.w.anno.mir[[1]],fb_2.c.w.anno.mir[[1]],fb_3.c.w.anno.mir[[1]],xrn1_1.c.w.anno.mir[[1]],xrn1_2.c.w.anno.mir[[1]],xrn1_3.c.w.anno.mir[[1]],
                       xrn4_1.c.w.anno.mir[[1]],xrn4_2.c.w.anno.mir[[1]],xrn4_3.c.w.anno.mir[[1]],d234_fb4.c.w.anno.mir[[1]],d234_fb5.c.w.anno.mir[[1]],d234_fb6.c.w.anno.mir[[1]])

#get vector of unique names of high confidence interactions
high.conf.v = unique(c(fb.c.w.anno.mir.set,xrn1.c.w.anno.mir.set,xrn4.c.w.anno.mir.set,d234.c.w.anno.mir.set))

#generate vector of names
df.names = c('fb_1','fb_2','fb_3','xrn1_1','xrn1_2','xrn1_3','xrn4_1','xrn4_2','xrn4_3','d234_fb4','d234_fb5','d234_fb6')

gaps.col = c(3,6,9) 

rpm.max = 5

inList.bioreps.sub = list()

df.rpm = data.frame(row.names=high.conf.v)

#build interaction.v for later assignment to rownames
interaction.v = c()
for (si in high.conf.v) {
  interaction=""
  for (df in df.inList.total) {
    if (si %in% rownames(df) == TRUE) {
      interaction=df[si,"interaction"]
    }
  }
  interaction.v=c(interaction.v,interaction)
}

for (df in df.inList.total) {
  rpm.v = c()
  for (si in high.conf.v) {
    if (si %in% rownames(df)) {
      rpm = as.numeric(df[si,'slice.site.rpm'])
    }else {
      rpm = 0
    }
    rpm.v = c(rpm.v,rpm)
  }
  df.rpm  = cbind(df.rpm,rpm.v)
}

names(df.rpm) = df.names
rownames(df.rpm) = make.unique(interaction.v)

df.rpm[(df.rpm > rpm.max)] <- rpm.max  

#build annotation df with verified sites and those not present in at least one biorep from cap.plus union
df.anno = data.frame(row.names=high.conf.v, stringsAsFactors = FALSE)


for (row in rownames(df.anno)) {
  verified = 'no'
  anno.row = c(verified)
  for (df in df.inList.total) {
    if (row %in% rownames(df) == TRUE) {
      verified = as.character(df[row,'verified.'])
    }
  }
  anno.row = c(verified)
  df.anno=rbind(df.anno, anno.row, stringsAsFactors=FALSE)
}

rownames(df.anno) = make.unique(interaction.v)
names(df.anno) = c("5' RACE")

plotFile = paste(printRoot,'/figures/fig.4g.miRNA.rpm.heatmap.pdf',sep='')
pdf(plotFile)
pheatmap(log((df.rpm * 10) + 1,10), cluster_rows=T, cluster_cols=F, legend=T, main='miRNA slice site RPTM + 1 (log10)', las=2, treeheight_row = 0, fontsize_row=8, fontsize_col=8, 
         show_rownames=TRUE, na.rm=F, scale="none", gaps_col=gaps.col, annotation_row=df.anno, cellwidth=10, cellheight=8)
dev.off()


##tasiRNAs
#get list of dfs for all interactions detected
df.inList.total = list(fb_1.c.w.anno.tas[[1]],fb_2.c.w.anno.tas[[1]],fb_3.c.w.anno.tas[[1]],xrn1_1.c.w.anno.tas[[1]],xrn1_2.c.w.anno.tas[[1]],xrn1_3.c.w.anno.tas[[1]],
                       xrn4_1.c.w.anno.tas[[1]],xrn4_2.c.w.anno.tas[[1]],xrn4_3.c.w.anno.tas[[1]],d234_fb4.c.w.anno.tas[[1]],d234_fb5.c.w.anno.tas[[1]],d234_fb6.c.w.anno.mir[[1]])

#get vector of unique names of high confidence interactions
high.conf.v = unique(c(fb.c.w.anno.tas.set,xrn1.c.w.anno.tas.set,xrn4.c.w.anno.tas.set,d234.c.w.anno.tas.set))

#generate vector of names
df.names = c('fb_1','fb_2','fb_3','xrn1_1','xrn1_2','xrn1_3','xrn4_1','xrn4_2','xrn4_3','d234_fb4','d234_fb5','d234_fb6')

gaps.col = c(3,6,9) 

rpm.max = 5

inList.bioreps.sub = list()

df.rpm = data.frame(row.names=high.conf.v)

#build interaction.v for later assignment to rownames
interaction.v = c()
for (si in high.conf.v) {
  interaction=""
  for (df in df.inList.total) {
    if (si %in% rownames(df) == TRUE) {
      interaction=df[si,"interaction"]
    }
  }
  interaction.v=c(interaction.v,interaction)
}

for (df in df.inList.total) {
  rpm.v = c()
  for (si in high.conf.v) {
    if (si %in% rownames(df)) {
      rpm = as.numeric(df[si,'slice.site.rpm'])
    }else {
      rpm = 0
    }
    rpm.v = c(rpm.v,rpm)
  }
  df.rpm  = cbind(df.rpm,rpm.v)
}

names(df.rpm) = df.names
rownames(df.rpm) = make.unique(interaction.v)

df.rpm[(df.rpm > rpm.max)] <- rpm.max  

#build annotation df with verified sites and those not present in at least one biorep from cap.plus union
df.anno = data.frame(row.names=high.conf.v, stringsAsFactors = FALSE)
#names(df.anno) = c('common','verified')

for (row in rownames(df.anno)) {
  verified = 'no'
  anno.row = c(verified)
  for (df in df.inList.total) {
    if (row %in% rownames(df) == TRUE) {
      verified = as.character(df[row,'verified.'])
    }
  }
  anno.row = c(verified)
  df.anno=rbind(df.anno, anno.row, stringsAsFactors=FALSE)
}

rownames(df.anno) = make.unique(interaction.v)
names(df.anno) = c("5' RACE")

plotFile = paste(printRoot,'/figures/fig.4g.tasiRNA.rpm.heatmap.pdf',sep='')
pdf(plotFile)
pheatmap(log((df.rpm * 10) + 1,10), cluster_rows=T, cluster_cols=F, legend=T, main='tasiRNA slice site RPTM + 1 (log10)', las=2, treeheight_row = 0, fontsize_row=8, fontsize_col=8, 
         show_rownames=TRUE, na.rm=F, scale="none", gaps_col=gaps.col, annotation_row=df.anno, cellwidth=10, cellheight=8)
dev.off()


#####################################Figs. 4H#####################################
#get unions for public sets
##miRNAs
addo.u.w.anno.mir.union = unique(c(rownames(GSM278334.u.w.anno.mir[[2]]),rownames(GSM278335.u.w.anno.mir[[2]])))
willman.u.w.anno.mir.union = unique(c(rownames(SRX283759.u.w.anno.mir[[2]]),rownames(SRX283760.u.w.anno.mir[[2]])))
hou.u.w.anno.mir.union = unique(c(rownames(SRX1559792.u.w.anno.mir[[2]]),rownames(SRX1980216.u.w.anno.mir[[2]]),rownames(SRX1980220.u.w.anno.mir[[2]])))
yu.u.w.anno.mir.union = unique(c(rownames(SRX1429178.u.w.anno.mir[[2]]),rownames(SRX1429179.u.w.anno.mir[[2]])))

##tasiRNAs
addo.u.w.anno.tas.union = unique(c(rownames(GSM278334.u.w.anno.tas[[2]]),rownames(GSM278335.u.w.anno.tas[[2]])))
willman.u.w.anno.tas.union = unique(c(rownames(SRX283759.u.w.anno.tas[[2]]),rownames(SRX283760.u.w.anno.tas[[2]])))
hou.u.w.anno.tas.union = unique(c(rownames(SRX1559792.u.w.anno.tas[[2]]),rownames(SRX1980216.u.w.anno.tas[[2]]),rownames(SRX1980220.u.w.anno.tas[[2]])))
yu.u.w.anno.tas.union = unique(c(rownames(SRX1429178.u.w.anno.tas[[2]]),rownames(SRX1429179.u.w.anno.tas[[2]])))

plotFile=paste(printRoot,'/figures/fig.4h.public.compare.strip.pdf',sep='')
pdf(plotFile, useDingbats=FALSE, width=8, height=5)
par(mar=c(7,4,4,2))
ymax=180
##miRNAs
stripchart(list(c(length(fb_1.c.w.anno.mir[[2]][,1]),length(fb_2.c.w.anno.mir[[2]][,1]),length(fb_3.c.w.anno.mir[[2]][,1])),
                c(length(xrn4_1.c.w.anno.mir[[2]][,1]),length(xrn4_2.c.w.anno.mir[[2]][,1]),length(xrn4_3.c.w.anno.mir[[2]][,1])),
                c(length(GSM278334.u.w.anno.mir[[2]][,1]),length(GSM278335.u.w.anno.mir[[2]][,1])),
                length(GSM280226.u.w.anno.mir[[2]][,1]),
                length(GSM280227.u.w.anno.mir[[2]][,1]),
                length(SRX003076.u.w.anno.mir[[2]][,1]),
                length(SRX003078.u.w.anno.mir[[2]][,1]),
                c(length(SRX283759.u.w.anno.mir[[2]][,1]),length(SRX283760.u.w.anno.mir[[2]][,1])),
                c(length(SRX1559792.u.w.anno.mir[[2]][,1]),length(SRX1980216.u.w.anno.mir[[2]][,1]),length(SRX1980220.u.w.anno.mir[[2]][,1])),
                c(length(SRX1429178.u.w.anno.mir[[2]][,1]),length(SRX1429179.u.w.anno.mir[[2]][,1])),
                length(SRX384362.u.w.anno.mir[[2]][,1])),
           method='jitter', jitter=0.2, main=NULL, vertical=TRUE, las=2, pch=19, cex=1.2, col=c('#b2182b'),
           group.names=c('Col-0','xrn4','Col-0\nAddo-Quaye','Col-0\nGerman','xrn4\nGerman et al.',
                         'Col-0\nGregory','xrn4\nGregory','Col-0\nWillman et al.','Col-0\nHou et al.',
                         'Col-0 +CHX\nYu et al.','Col-0\nCreasy et al.'), 
           xlab=NULL, ylab='Number of sRNA target sites (p.adj < 0.05)',ylim=c(0,ymax), cex.axis=1.0)
#union
stripchart(list(c(length(fb.c.w.anno.mir.union)),c(length(xrn4.c.w.anno.mir.union)),c(length(addo.u.w.anno.mir.union)),c(),c(),c(),c(),
                c(length(willman.u.w.anno.mir.union)),c(length(hou.u.w.anno.mir.union)),c(length(yu.u.w.anno.mir.union))),
           method='jitter', jitter=0.0, vertical=TRUE, pch=17, las=1, col=c('#b2182b'), cex=1.2, 
           xlab=NULL, ylim=c(0,ymax), cex.axis=1.0, add=TRUE)
##tasiRNAs
stripchart(list(c(length(fb_1.c.w.anno.tas[[2]][,1]),length(fb_2.c.w.anno.tas[[2]][,1]),length(fb_3.c.w.anno.tas[[2]][,1])),
                c(length(xrn4_1.c.w.anno.tas[[2]][,1]),length(xrn4_2.c.w.anno.tas[[2]][,1]),length(xrn4_3.c.w.anno.tas[[2]][,1])),
                c(length(GSM278334.u.w.anno.tas[[2]][,1]),length(GSM278335.u.w.anno.tas[[2]][,1])),
                length(GSM280226.u.w.anno.tas[[2]][,1]),
                length(GSM280227.u.w.anno.tas[[2]][,1]),
                length(SRX003076.u.w.anno.tas[[2]][,1]),
                length(SRX003078.u.w.anno.tas[[2]][,1]),
                c(length(SRX283759.u.w.anno.tas[[2]][,1]),length(SRX283760.u.w.anno.tas[[2]][,1])),
                c(length(SRX1559792.u.w.anno.tas[[2]][,1]),length(SRX1980216.u.w.anno.tas[[2]][,1]),length(SRX1980220.u.w.anno.tas[[2]][,1])),
                c(length(SRX1429178.u.w.anno.tas[[2]][,1]),length(SRX1429179.u.w.anno.tas[[2]][,1])),
                length(SRX384362.u.w.anno.tas[[2]][,1])),
           method='jitter', jitter=0.2, main=NULL, vertical=TRUE, las=2, pch=19, cex=1.2, col=c('#2166ac'),
           group.names=c('Col-0','xrn4','Col-0\nAddo-Quaye','Col-0\nGerman','xrn4\nGerman et al.',
                         'Col-0\nGregory','xrn4\nGregory','Col-0\nWillman et al.','Col-0\nHou et al.',
                         'Col-0 +CHX\nYu et al.','Col-0\nCreasy et al.'), 
           xlab=NULL, ylab='Number of sRNA target sites (p.adj < 0.05)',ylim=c(0,ymax), cex.axis=1.0, add=TRUE)
#union
stripchart(list(c(length(fb.c.w.anno.tas.union)),c(length(xrn4.c.w.anno.tas.union)),c(length(addo.u.w.anno.tas.union)),c(),c(),c(),c(),
                c(length(willman.u.w.anno.tas.union)),c(length(hou.u.w.anno.tas.union)),c(length(yu.u.w.anno.tas.union))),
           method='jitter', jitter=0.0, vertical=TRUE, pch=17, las=1, col=c('#2166ac'), cex=1.2, 
           xlab=NULL, ylim=c(0,ymax), cex.axis=1.0, add=TRUE)
legend(x='topleft', legend=c('Bioreps (miRNAs)','Bioreps (tasiRNAs)','Union (miRNAs)','Union (tasiRNAs)'), 
       col=c('#b2182b','#2166ac','#b2182b','#2166ac'), pch=c(19,19,17,17), pt.cex=1.2, bty='n')

dev.off()

#######FOR TEXT######
public.mir.unique = unique(c(rownames(GSM278334.u.w.anno.mir[[2]]),rownames(GSM278335.u.w.anno.mir[[2]]),rownames(GSM280226.u.w.anno.mir[[2]]),
                             rownames(GSM280227.u.w.anno.mir[[2]]),rownames(SRX003076.u.w.anno.mir[[2]]),rownames(SRX003078.u.w.anno.mir[[2]]),
                             rownames(SRX283759.u.w.anno.mir[[2]]),rownames(SRX283760.u.w.anno.mir[[2]]),rownames(SRX1559792.u.w.anno.mir[[2]]),
                             rownames(SRX1980216.u.w.anno.mir[[2]]),rownames(SRX1980220.u.w.anno.mir[[2]]),rownames(SRX1429178.u.w.anno.mir[[2]]),
                             rownames(SRX1429179.u.w.anno.mir[[2]]),rownames(SRX384362.u.w.anno.mir[[2]])))
#length(public.mir.unique) = 263

public.tas.unique = unique(c(rownames(GSM278334.u.w.anno.tas[[2]]),rownames(GSM278335.u.w.anno.tas[[2]]),rownames(GSM280226.u.w.anno.tas[[2]]),
                             rownames(GSM280227.u.w.anno.tas[[2]]),rownames(SRX003076.u.w.anno.tas[[2]]),rownames(SRX003078.u.w.anno.tas[[2]]),
                             rownames(SRX283759.u.w.anno.tas[[2]]),rownames(SRX283760.u.w.anno.tas[[2]]),rownames(SRX1559792.u.w.anno.tas[[2]]),
                             rownames(SRX1980216.u.w.anno.tas[[2]]),rownames(SRX1980220.u.w.anno.tas[[2]]),rownames(SRX1429178.u.w.anno.tas[[2]]),
                             rownames(SRX1429179.u.w.anno.tas[[2]]),rownames(SRX384362.u.w.anno.tas[[2]])))
#length(public.tas.unique) = 22


#total number of unique miRNA and tasiRNA targets identified in public datasets as well as all datasets combined
##capmasked
all.mir.unique = unique(c(fb.c.w.anno.mir.union,xrn1.c.w.anno.mir.union,xrn4.c.w.anno.mir.union,d234.c.w.anno.mir.union,
                        rownames(GSM278334.u.w.anno.mir[[2]]),rownames(GSM278335.u.w.anno.mir[[2]]),rownames(GSM280226.u.w.anno.mir[[2]]),
                        rownames(GSM280227.u.w.anno.mir[[2]]),rownames(SRX003076.u.w.anno.mir[[2]]),rownames(SRX003078.u.w.anno.mir[[2]]),
                        rownames(SRX283759.u.w.anno.mir[[2]]),rownames(SRX283760.u.w.anno.mir[[2]]),rownames(SRX1559792.u.w.anno.mir[[2]]),
                        rownames(SRX1980216.u.w.anno.mir[[2]]),rownames(SRX1980220.u.w.anno.mir[[2]]),rownames(SRX1429178.u.w.anno.mir[[2]]),
                        rownames(SRX1429179.u.w.anno.mir[[2]]),rownames(SRX384362.u.w.anno.mir[[2]])))
#length(all.mir.unique) = 296

all.tas.unique = unique(c(fb.c.w.anno.tas.union,xrn1.c.w.anno.tas.union,xrn4.c.w.anno.tas.union,d234.c.w.anno.tas.union,
                          rownames(GSM278334.u.w.anno.tas[[2]]),rownames(GSM278335.u.w.anno.tas[[2]]),rownames(GSM280226.u.w.anno.tas[[2]]),
                          rownames(GSM280227.u.w.anno.tas[[2]]),rownames(SRX003076.u.w.anno.tas[[2]]),rownames(SRX003078.u.w.anno.tas[[2]]),
                          rownames(SRX283759.u.w.anno.tas[[2]]),rownames(SRX283760.u.w.anno.tas[[2]]),rownames(SRX1559792.u.w.anno.tas[[2]]),
                          rownames(SRX1980216.u.w.anno.tas[[2]]),rownames(SRX1980220.u.w.anno.tas[[2]]),rownames(SRX1429178.u.w.anno.tas[[2]]),
                          rownames(SRX1429179.u.w.anno.tas[[2]]),rownames(SRX384362.u.w.anno.tas[[2]])))
#length(all.tas.unique) = 29

np.mir.unique = unique(c(fb.c.w.anno.mir.union,xrn1.c.w.anno.mir.union,xrn4.c.w.anno.mir.union,d234.c.w.anno.mir.union))
#length(np.mir.unique) = 126

np.tas.unique = unique(c(fb.c.w.anno.tas.union,xrn1.c.w.anno.tas.union,xrn4.c.w.anno.tas.union,d234.c.w.anno.tas.union))
#length(np.tas.unique) = 18
