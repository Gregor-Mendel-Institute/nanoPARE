#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#EndCut_step2.v03.git.R
#takes output from EndCut_step1.py and for each sRNA to be tested:
# 1) Retrieves corresponding shuffled entries and computes corresponding ecdfs to generate background models for 
#     a) fold-change between slice-site and flanking 20-nt and 50-nt regions and 
#     b) Allen Scores
# 2) Computes p-values that sites for sRNAs being tested have higher fold-changes and lower Allen Scores than 
#    expected based on the corresponding background models
#input:
#     args: [1] sample name (e.g. fb1_1to1), [2] directory where *detectedSites.txt files are located (e.g. transcript_bedgraph_capmasked/W), 
#           [3] file suffix (e.g. W.transcript.capmasked), [4] filepath where files are located (e.g. '/Users/michael.nodine/Desktop/R_Desktop/nanoPARE/180228/'),
#output:
#     1. Graphs illustrating cumulative fractions of fold-changes and Allen scores of shuffled and non-shuffled sRNAs being tested
#     2. Tab-delimited files of various statistics, including p-values, for all detected sites

if (length(args) != 4) {
  stop("Four arguments are required", call.=FALSE)
} else if (length(args)==4) {
  name = args[1]
  dir = args[2]
  suff = args[3]
  dataRoot = args[4]
}

print(paste("sample name:",name))
print(paste("directory:",dir))
print(paste("file ending:",suff))
print(paste("Working directory:",dataRoot))

setwd(dataRoot)
dir.create(file.path(dataRoot,"graphs"))
printRoot = paste(dataRoot,"graphs",sep="")
dir.create(file.path(printRoot,dir,"cdfs/"))
dir.create(file.path(dataRoot,"results"))

library(stringr)
library(survcomp)

################################################FUNCTIONS################################################
get.sRNAs.anno <- function(sRNA.type) {
  df = read.delim(paste(dataRoot,'Ath_annotations/anno.',sRNA.type,'.tsv',sep=''),sep='\t', row.names=1, header=TRUE,stringsAsFactors=FALSE,strip.white=TRUE)
  return(df)
}

get.common <- function() {
  anno = read.delim(paste(dataRoot,'Ath_annotations/gene_aliases_20130130.txt',sep=''),sep='\t', row.names=NULL, header=TRUE,stringsAsFactors=FALSE,strip.white=TRUE)
  names.new = make.names(anno$locus_name, unique=TRUE)
  rownames(anno) = names.new
  return(anno)
}

get.5pTargets <- function(sRNA.type) {
  targets = read.delim(paste(dataRoot,'Ath_annotations/sRNA_targets_5pRACE.tsv',sep=''),sep='\t', row.names=NULL, header=TRUE,stringsAsFactors=FALSE,strip.white=TRUE)
  targets.sub=subset(targets, select=c(sRNA,targetAGI))
  tar.fam.v=c()
  if (sRNA.type=='miRNA') {
    tar.fam=substr(as.character(targets.sub$sRNA),1,6)  
  }
  if (sRNA.type=='tasiRNA') {
    tar.fam=substr(as.character(targets.sub$sRNA),1,4)    
  }
  targets.sub.2=data.frame(targets.sub,'fam'=tar.fam)
  return(targets.sub.2)
}

get.detected.sites <- function(dir,sample,suff,window,sRNA.df,shuff.num) {

  if (substr(rownames(sRNA.df)[1],1,3) == 'TAS') {sRNA.type='tasiRNA'}
  if (substr(rownames(sRNA.df)[1],1,3) == 'miR') {sRNA.type='miRNA'}
  
  #get df of fc and Allen frequencies
  test=data.frame()
  test = read.delim(paste(dataRoot,'results/',dir,'/',sample,'/anno.mir.tas.fa.GSTAr_',sample,'.',suff,'_win_',window,'_detectedSites.txt',sep=''),sep='\t', row.names=NULL, header=TRUE,stringsAsFactors=FALSE,strip.white=TRUE)
  ##subset for specific sRNA isoforms
  test = subset(test, sRNA %in% rownames(sRNA.df))
  
  #sort by Allen, rpm and sRNA
  test = test[order(test$AllenScore,test$slice.site.rpm,test$sRNA, decreasing=c(FALSE,TRUE,FALSE)),]
  
  #build slice site identifier, calculate new fold-change (with pseudocount) and generate AGI column
  if (length(test$Transcript) > 0) {
    slice.id = paste(test$Transcript,'_',test$Slice.site,sep='')
    AGI.v = str_split_fixed(as.character(test$Transcript), "[.]",length(test$Transcript))[,1]
  }else {
    slice.id = c()
    AGI.v = c()
  }
  
  fc.new =  (test$slice.site.rpm + 1)/(test$flankMax.rpm + 1) #c()
  test = data.frame(test,'fc.new'=fc.new, 'slice.id'=slice.id, 'AGI'=AGI.v)
  
  #remove self:self hits
  ##get df of sRNAs and their precursors
  target.AGI.df = read.delim(paste(dataRoot,'Ath_annotations/miRBase21_and_TAIR10_miRNA_tasiRNA_precursors.txt',sep=''),sep='\t', row.names=1, header=TRUE,stringsAsFactors=FALSE,strip.white=TRUE)
  
  test.new = data.frame()
  
  for (row in rownames(test)) {
    if (as.character(test[row,'AGI']) != target.AGI.df[test[row,'sRNA'],'precursor']) {
      test.new=rbind(test.new,test[row,])
    }
  }
  
  #remove redundant slice.sites +/- 5 nt; if >1, then take one with lower Allen score...if same, then higher rpm...if same, then first alphabetical sRNA
  ##get unique transcripts with detectable sites
  transcripts.uni.v = unique(test.new$Transcript)
  
  test.new.2 = data.frame()

  for (transcript in transcripts.uni.v) {
    #get entries with same transcript
    transcript.sub = subset(test.new, Transcript==transcript)
    #add first entry (top RPM) to df
    transcript.sub.1 = transcript.sub[1,]
    test.new.2=rbind(test.new.2,transcript.sub.1)
    #look for additional entries that are more than 5 nt away
    for (row in rownames(transcript.sub)) {
      test.new.2.sites = as.numeric(subset(test.new.2, Transcript==transcript)$Slice.site)
      if (min(abs(as.numeric(transcript.sub[row,'Slice.site']) - test.new.2.sites)) > 5) {
        test.new.2=rbind(test.new.2,transcript.sub[row,])
      }
    }
  }

  test = test.new.2
  
  #get corresponding dataframes from shuffled sets; generate single dataframe from all shuffled
  ##note: not necessary to remove self:self hits (don't exist), also don't remove redundant slice sites because
  ##randomized sRNAs are highly unlikely to overlap and this becomes practically difficult due to combining all shuff entries
  shuff = data.frame()
  shuff = as.data.frame(setNames(replicate(length(names(test)),numeric(0), simplify = F), names(names(test))))
  for (i in 0:(shuff.num-1)) {
    shuff.sub = read.delim(paste(dataRoot,'results/',dir,'/',sample,'/anno.mir.tas.',i,'.shuffled.fa.GSTAr_',sample,'.',suff,'_win_',window,'_detectedSites.txt',sep=''),sep='\t', row.names=NULL, header=TRUE,stringsAsFactors=FALSE,strip.white=TRUE)
    shuff.sub = subset(shuff.sub, sRNA %in% rownames(sRNA.df))
    shuff = rbind(shuff,shuff.sub)
  }
  
  fc.new = (shuff$slice.site.rpm + 1)/(shuff$flankMax.rpm + 1)
  shuff = data.frame(shuff, 'fc.new'=fc.new) #, 'slice.id'=slice.id)
  
  shuff[is.na(shuff)] <- 0
  
  test = test[order(test$sRNA),]
  
  #for test, add common names and whether they have been verified by 5pRACE
  common = get.common()
  common.v = c()
  
  raceTars = get.5pTargets(sRNA.type)
  
  tar.v = c()
  
  for (name in rownames(test)) {
    agi=as.character(test[name,'AGI'])
    if (sRNA.type=='miRNA') {sRNA.fam=substr(test[name,'sRNA'],1,6)}
    if (sRNA.type=='tasiRNA')  {sRNA.fam=substr(test[name,'sRNA'],1,4)}
    if (agi %in% rownames(common)) {common.sub = common[agi,'symbol']}else {common.sub = 'NA'}
    common.v = c(common.v, common.sub)
    #note: if multiple sRNAs have been verified then just take one for simplicity
    #for miRNAs, test whether first six characters match (miRNA family); for tasiRNAs, test whether full match
    raceTars.sub = subset(raceTars, agi==targetAGI & sRNA.fam==fam)
    if (length(raceTars.sub[,1] > 0)) {tar.sub='yes'}else {tar.sub = 'no'}
    
    tar.v = c(tar.v, tar.sub)
  }
  
  test = data.frame(test, 'common.name'=common.v, 'validated.by.5pRACE'=tar.v)
  
  return(list(test,shuff))
}

get.target.sites <- function(dir,sites,sample.name) {

  test = sites[[1]]
  shuff = sites[[2]]
  
  #for each sRNA, extract fc and allen vals from corresponding shuffled sets,
  #combine and generate ecdf, compute p-vals and plot graphs
  ##get all unique sRNAs
  sRNA.v = unique(test$sRNA)

  ##for each sRNA, extract dfs and add to lists
  sRNA.shuff.list = list()
  sRNA.test.list = list()
  
  for (sRNA.u in sRNA.v) {
    sRNA.test = subset(test, sRNA == sRNA.u)
    sRNA.test.list[[sRNA.u]] = sRNA.test
    sRNA.shuff = subset(shuff, sRNA == sRNA.u)
    sRNA.shuff.list[[sRNA.u]] = sRNA.shuff
  }
  
  ##for each sRNA, generate fc and Allen ecdfs with shuffled set and compute p-values for each detected sRNA:target in tests
  extra.names = c('fc.p.val','allen.p.val')
  final = as.data.frame(setNames(replicate(length(c(names(test),extra.names)),numeric(0), simplify = F), c(names(test),extra.names)))
  
  for (sRNA.u in sRNA.v) {
    ind.test = sRNA.test.list[[sRNA.u]]
    ind.shuff = sRNA.shuff.list[[sRNA.u]]
    #compute p.vals
    fc.p.val = 1 - ecdf(ind.shuff$fc.new)(ind.test$fc.new)
    allen.p.val = ecdf(ind.shuff$AllenScore)(ind.test$AllenScore)
    
    final = rbind(final,data.frame(subset(test,sRNA==sRNA.u),'fc.p.val'=fc.p.val,'allen.p.val'=allen.p.val))
    
    plotFile=paste(printRoot,'/',dir,'/cdfs/',sRNA.u,'.',sample.name,'.fc.ecdf.pdf',sep='')
    pdf(plotFile, useDingbats=FALSE, width=5, height=8)
    par(mfrow = c(2,1))
    plot(NULL, main=paste('Detected', sRNA.u, 'sites in',sample.name), xlab='Fold-changes', ylab='Proportion', las=1, xlim=c(0,max(c(ind.shuff$fc.new,ind.test$fc.new),na.rm=T)), ylim=c(0,1.0))
    lines(ecdf(ind.shuff$fc.new), pch=1)
    color='coral'
    points(ind.test$fc.new,ecdf(ind.shuff$fc.new)(ind.test$fc.new), pch=19, col=color)
    legend(x='bottomright', legend=c(paste('Test (n = ',length(rownames(ind.test)),')',sep=''),paste('Shuffled (n = ',length(rownames(ind.shuff)),')',sep='')), col=c(color,'black'), bty='n', pch=c(19,1))
    text(ind.test$fc.new,ecdf(ind.shuff$fc.new)(ind.test$fc.new),labels=paste(ind.test$Transcript,'/',ind.test$common.name,sep=''),col=color,pos=1,cex=0.6)
    
    plot(NULL, main=NULL, xlab='Allen Scores', ylab='Proportion', las=1, xlim=c(0,max(c(ind.shuff$AllenScore,ind.test$AllenScore),na.rm=T)), ylim=c(0,1.0))
    lines(ecdf(ind.shuff$AllenScore), pch=1)
    color='slateblue'
    points(ind.test$AllenScore,ecdf(ind.shuff$AllenScore)(ind.test$AllenScore), pch=19, col=color)
    text(ind.test$AllenScore,ecdf(ind.shuff$AllenScore)(ind.test$AllenScore),labels=paste(ind.test$Transcript,'/',ind.test$common.name,sep=''),col=color,pos=3,cex=0.6)
    legend(x='topleft', legend=c(paste('Test (n = ',length(rownames(ind.test)),')',sep=''),paste('Shuffled (n = ',length(rownames(ind.shuff)),')',sep='')), col=c(color,'black'), bty='n', pch=c(19,1) )
    dev.off()
  }
  
  #combine p-values Fishers method
  fishers.v = c()
  
  for (row in rownames(final)) {
    fishers.combined.p.val = survcomp::combine.test(c(final[row,'fc.p.val'],final[row,'allen.p.val']), method='fisher')
    fishers.v = c(fishers.v,fishers.combined.p.val)    
  }
  
  combined.p.val.adj = p.adjust(fishers.v, method='fdr')
  
  final.plus = data.frame(final,'fishers.combined.p.val'=fishers.v,'fishers.combined.p.val.adj'=combined.p.val.adj)
  
  final.plus = final.plus[order(final.plus$fishers.combined.p.val.adj),]
  
  sig.hits = subset(final.plus, fc.new > 1.0 & fishers.combined.p.val.adj < 0.05 & slice.site.rpm > 0.1)
  
  return(list(final.plus,sig.hits))
}

get.tars <- function(dir,sample,suff,sRNAs) {

  detected.sites.50 = get.detected.sites(dir,sample,suff,50,sRNAs,1000)
  detected.sites.20 = get.detected.sites(dir,sample,suff,20,sRNAs,1000)
  
  detected.tars.50 = get.target.sites(dir,detected.sites.50,paste(sample,'.50',sep=''))
  detected.tars.20 = get.target.sites(dir,detected.sites.20,paste(sample,'.20',sep=''))
  
  detected.sites = data.frame()
  target.sites = data.frame()
  
  if (length(detected.tars.50[[1]][,1]) > 0 & length(detected.tars.20[[1]][,1]) > 0) {
    detected.sites = merge(detected.tars.20[[1]],detected.tars.50[[1]], by='slice.id', all=TRUE)
  }else if (length(detected.tars.50[[1]][,1]) > 0) {
    detected.sites = detected.tars.50[[1]]
  }else if (length(detected.tars.20[[1]][,1]) > 0) {
    detected.sites = detected.tars.20[[1]]
  }
  
  if (length(detected.tars.50[[2]][,1]) > 0 & length(detected.tars.20[[2]][,1]) > 0) {
    target.sites = merge(detected.tars.20[[2]],detected.tars.50[[2]], by='slice.id', all=TRUE)
  }else if (length(detected.tars.50[[2]][,1]) > 0) {
    target.sites = detected.tars.50[[2]]
  }else if (length(detected.tars.20[[2]][,1]) > 0) {
    target.sites = detected.tars.20[[2]]
  }
  
  #reformat target.sites to unique sRNA:common (or AGI)
  name.new = c()
  
  if (length(detected.tars.50[[2]][,1]) > 0 & length(detected.tars.20[[2]][,1]) > 0) {
    for (name in rownames(target.sites)) {
      if (is.na(target.sites[name,'sRNA.x']) == FALSE) {sRNA=as.character(target.sites[name,'sRNA.x'])}else {sRNA=as.character(target.sites[name,'sRNA.y'])}
      if (is.na(target.sites[name,'common.name.x']) == FALSE & target.sites[name,'common.name.x'] != 'NA') {
        target=as.character(target.sites[name,'common.name.x'])
      } else if (is.na(target.sites[name,'common.name.y']) == FALSE & target.sites[name,'common.name.y'] != 'NA') {target=as.character(target.sites[name,'common.name.y'])
      } else if (is.na(target.sites[name,'Transcript.x']) == FALSE) {target=as.character(target.sites[name,'Transcript.x'])
      } else {target=as.character(target.sites[name,'Transcript.y'])}
      
      name.new.sub = paste(sRNA,':',target,sep='')
      name.new = c(name.new,name.new.sub)
    }
  }
  
  if (length(detected.tars.50[[2]][,1]) == 0 | length(detected.tars.20[[2]][,1]) == 0) {
    for (name in rownames(target.sites)) {
      sRNA=as.character(target.sites[name,'sRNA'])
      target=as.character(target.sites[name,'common.name'])
      
      name.new.sub = paste(sRNA,':',target,sep='')
      name.new = c(name.new,name.new.sub)
    }
  }
  
  target.sites= data.frame(target.sites,'interaction'=name.new)
  
  name.new = c()
  
  for (name in rownames(detected.sites)) {
    if (is.na(detected.sites[name,'sRNA.x']) == FALSE) {sRNA=as.character(detected.sites[name,'sRNA.x'])}else {sRNA=as.character(detected.sites[name,'sRNA.y'])}
    if (is.na(detected.sites[name,'common.name.x']) == FALSE & detected.sites[name,'common.name.x'] != 'NA') {
      target=as.character(detected.sites[name,'common.name.x'])
    } else if (is.na(detected.sites[name,'common.name.y']) == FALSE & detected.sites[name,'common.name.y'] != 'NA') {target=as.character(detected.sites[name,'common.name.y'])
    } else if (is.na(detected.sites[name,'Transcript.x']) == FALSE) {target=as.character(detected.sites[name,'Transcript.x'])
    } else {target=as.character(detected.sites[name,'Transcript.y'])}
    
    name.new.sub = paste(sRNA,':',target,sep='')
    name.new = c(name.new,name.new.sub)
  }
  
  detected.sites= data.frame(detected.sites,'interaction'=name.new)
  
  return(list(detected.sites,target.sites))
}

print.tables <- function(df.total,outFile) {

  df.total.new = data.frame()
  
  #compute minimum combined adj. p-values for each interactions
  p.val.min.v = c()
  sRNA.v = c()
  slice.v = c()
  allen.v = c()
  rpm.v = c()
  fc.v = c()
  verify.v = c()
  common.v = c()
  
  for (row in rownames(df.total)) {
    verify.sub = min(as.character(df.total[row,'validated.by.5pRACE.x']),as.character(df.total[row,'validated.by.5pRACE.y']), na.rm=TRUE)
    common.sub = min(as.character(df.total[row,'common.name.x']),as.character(df.total[row,'common.name.y'], na.rm=TRUE))
    p.val.min.sub = min(df.total[row,'fishers.combined.p.val.adj.x'],df.total[row,'fishers.combined.p.val.adj.y'], na.rm=TRUE)
    if (p.val.min.sub == df.total[row,'fishers.combined.p.val.adj.x']) {
      sRNA.sub = df.total[row,'sRNA.x']
      allen.sub = df.total[row,'AllenScore.x']
      rpm.sub = df.total[row,'slice.site.rpm.x']
      fc.sub= df.total[row,'fc.new.x']
      slice.sub= df.total[row,'Slice.site.x']
    }else {
      sRNA.sub = df.total[row,'sRNA.y']
      allen.sub = df.total[row,'AllenScore.y']
      rpm.sub = df.total[row,'slice.site.rpm.y']
      fc.sub= df.total[row,'fc.new.y']
      slice.sub= df.total[row,'Slice.site.y']
    }
    p.val.min.v = c(p.val.min.v,p.val.min.sub)
    sRNA.v = c(sRNA.v,sRNA.sub)
    allen.v = c(allen.v,allen.sub)
    rpm.v = c(rpm.v,rpm.sub)
    fc.v = c(fc.v,fc.sub)
    slice.v = c(slice.v,slice.sub)
    verify.v = c(verify.v,verify.sub)
    common.v = c(common.v,common.sub)
  }
  
  df.total.new = data.frame('target.transcript'=df.total$Transcript.x, 'sRNA'=sRNA.v, 'slice.site'=slice.v, 'Allen.score'=allen.v, 'slice.site.rpm'=rpm.v, 'fold.change'=fc.v,
                            'verified?'=verify.v, 'common.name'=common.v, 'interaction'=as.character(df.total$interaction), 'adj.p.val'=p.val.min.v)
  
  df.total.new = df.total.new[order(df.total.new$adj.p.val),]
  
  write.table(df.total.new, file=paste(dataRoot,'/results/',outFile,'.detected.sites.tsv',sep=''), quote=FALSE, row.names=FALSE, col.names=TRUE, sep='\t')
}
#########################################################################################################

#get annotated miRNAs and tasiRNAs to test for cleavage sites
print("Getting annotated miRNAs and tasiRNAs...")
anno.miRNAs = get.sRNAs.anno('miRNA')
anno.tasiRNAs = get.sRNAs.anno('tasiRNA')

print("Testing for tasiRNA targets...")
sample.tas = get.tars(dir,name,suff,anno.tasiRNAs)
print("Testing for miRNA targets...")
sample.mir = get.tars(dir,name,suff,anno.miRNAs)

#write to tsv files
print("Writing to outfiles...")
##tasiRNAs
print.tables(sample.tas[[1]],paste(dir,'/',name,'/',name,'.anno.tas.allen',sep=''))
##miRNAs
print.tables(sample.mir[[1]],paste(dir,'/',name,'/',name,'.anno.mir.allen',sep=''))
