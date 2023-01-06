### plotting Translational efficiencies comparing different clusters of MZ genes
library(ggplot2)
library(dplyr)
library(cowplot)

theme_ameres <- function (type) { 
  
  types_plot = c("boxplot","barplot","violin")
  if(type %in% types_plot ==T)
  {
    
    if(type == "boxplot"){
      theme(legend.title=element_blank(),axis.text.x = element_text(margin=margin(10,15,10,15,"pt"),size = 15),axis.text.y = element_text(margin=margin(5,15,10,5,"pt"),size = 15), axis.ticks.length = unit(-0.25 , "cm"),legend.position="bottom",axis.line = element_line(colour = "black", lineend = "round"), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18))   
    }
    
    if(type == "barplot"){
      theme(legend.title=element_blank(),axis.text.x = element_text(margin=margin(10,15,10,15,"pt"),size = 15, hjust = 1),axis.text.y = element_text(margin=margin(5,15,10,5,"pt"),size = 15), axis.ticks.length = unit(-0.25 , "cm"),legend.position="bottom",axis.line = element_line(colour = "black", lineend = "round"), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18))   
    }
    
    
  }
  
}

masterTable_data = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure5/data//ClassifiedGenes_masterTable.txt",stringsAsFactors = F, header = T,sep = "\t")
masterTable_data =  masterTable_data %>% dplyr::filter(description != 'MZ from Z') %>% dplyr::filter(description != 'MZ from M')
#### rate of increase of MZ genes... 
MZgenes = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure4/data/clusteredBasedOnTCaccumulation.bed",sep="\t",stringsAsFactors = F)
MZgenes = MZgenes %>% dplyr::mutate(gene=V10,class=V12) %>% dplyr::select('gene','V11','class')
MZgenes = MZgenes %>% dplyr::mutate(external_gene_name = gene,MZrate = class) %>% dplyr::select('external_gene_name','MZrate','V11')


masterTable_data_fractionU = masterTable_data %>% dplyr::select(dplyr::contains("TE_"))
metaData = masterTable_data[,c(1:3)]
masterTable_data_fractionU  = cbind.data.frame(masterTable_data_fractionU,metaData)
masterTable_data_fractionU =  plyr::join(masterTable_data_fractionU,MZgenes,by='external_gene_name')            
masterTable_data_fractionU$MZrate = factor(masterTable_data_fractionU$MZrate)
timepoints = colnames(masterTable_data_fractionU)[1:5]
times = c('2_4cell','256_cell','1K_cell','Dome','Shield')



pdf("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure5/plots/TranslationEffeciency_Allgenes.pdf", height = 4, width=4)

for(i in 1:length(timepoints)){
  Mgenes = masterTable_data_fractionU  %>%
    dplyr::mutate(fractionU=get(timepoints[i])) %>% dplyr::filter(class== "M") %>% dplyr::mutate(V11='M')%>% dplyr::select(c('MZrate','fractionU','V11'))
  Mstablegenes = masterTable_data_fractionU  %>%
    dplyr::mutate(fractionU=get(timepoints[i])) %>% dplyr::filter(class== "M-stable") %>% dplyr::mutate(V11='M-stable')%>% dplyr::select(c('MZrate','fractionU','V11'))
  
  
  b= masterTable_data_fractionU %>%
    dplyr::mutate(fractionU=get(timepoints[i])) %>% dplyr::select(c('MZrate','fractionU','V11'))
  
  b = rbind.data.frame(Mgenes,b)
  b= b[complete.cases(b$V11) ,]
  
  b_split = split(b,b$V11)
  pval_ks= ks.test(b_split$`1`$fractionU,b_split$`4`$fractionU)$p.value
  
  b_split= lapply(b_split,function(x) x[order(x[,'fractionU'],decreasing = F),])
  b_split = lapply(b_split,function(x) x %>% dplyr::mutate(fractionU_cumulative  = cumsum(fractionU)))
  b_split=lapply(b_split,function(x) x%>% dplyr::mutate(cumsum_frac = fractionU_cumulative/max(fractionU_cumulative)))
  
  b_split_2hpf = do.call(rbind,b_split)
  b_split_2hpf = b_split_2hpf %>% dplyr::group_by(V11) %>% dplyr::mutate(n=n())  %>% dplyr::mutate(category_sample = paste0(V11,"(n=",n,")"))
  b_split_2hpf$V11  = as.factor(b_split_2hpf$V11)
  p =  ggplot(b_split_2hpf, aes(log10(fractionU), color=V11,group=V11)) +
    stat_ecdf() + theme_cowplot()
  p = p+ theme_ameres(type = 'barplot')
  print(p)
}

b = masterTable_data_fractionU  %>% dplyr::mutate(fractionU = TE_Shield/TE_256) %>% dplyr::filter(class=="M" | class=="MZ" | class == "M-stable") %>% dplyr::select(c('class','V11','fractionU','external_gene_name'))
Mgenes = b %>% dplyr::filter(class=="M")%>% dplyr::mutate(V11=ifelse(class=='M','M',V11))
Mstable = b %>% dplyr::filter(class=="M-stable")%>% dplyr::mutate(V11=ifelse(class=='M-stable','M-stable',V11))
other_b =  b %>% dplyr::filter(class=="MZ")
b = rbind.data.frame(other_b,Mgenes,Mstable)
b= b[complete.cases(b$fractionU) ,] 


      b_split = split(b,b$V11)
      pval_ks_4= ks.test(b_split$`M-stable`$fractionU,b_split$`4`$fractionU)$p.value
      pval_ks_1= ks.test(b_split$`1`$fractionU,b_split$`M-stable`$fractionU)$p.value
      pval_ks_2= ks.test(b_split$`2`$fractionU,b_split$`M-stable`$fractionU)$p.value
      pval_ks_3= ks.test(b_split$`3`$fractionU,b_split$`M-stable`$fractionU)$p.value
      pval_ks_Mstable = ks.test(b_split$`M-stable`$fractionU,b_split$M$fractionU)$p.value
      
      b_split= lapply(b_split,function(x) x[order(x[,'fractionU'],decreasing = F),])
      b_split = lapply(b_split,function(x) x %>% dplyr::mutate(fractionU_cumulative  = cumsum(fractionU)))
      b_split=lapply(b_split,function(x) x%>% dplyr::mutate(cumsum_frac = fractionU_cumulative/max(fractionU_cumulative)))
      
      b_split_2hpf = do.call(rbind,b_split)
      b_split_2hpf = b_split_2hpf %>% dplyr::group_by(V11) %>% dplyr::mutate(n=n())  %>% dplyr::mutate(category_sample = paste0(V11,"(n=",n,")"))
      b_split_2hpf$V11 = as.factor(b_split_2hpf$V11)
      
      p =  ggplot(b_split_2hpf, aes(log2(fractionU), color=category_sample,group=V11)) +
        stat_ecdf()  + theme_cowplot()
      p = p + annotate( "text", x = 1, y = 0.5, label = pval_ks_1)+  annotate( "text", x = 1, y = 0.4, label = pval_ks_2)+
        annotate( "text", x = 1, y = 0.3, label = pval_ks_3) +  annotate( "text", x = 1, y = 0.2, label = pval_ks_4)+
        annotate( "text", x = 1, y = 0.1, label = pval_ks_Mstable)+
        theme_ameres(type = 'barplot') + ylab('Cumulative fraction') + xlab ('Fraction uridylation')
      print(p)
      
      p = ggplot(b_split_2hpf,aes(y=log10(fractionU),x=category_sample)) + geom_boxplot(outlier.shape = NA) + ggpubr::stat_compare_means(ref.group = 'M(n=1262)') + theme_cowplot()
      p + theme_ameres(type = 'barplot')
      for(i in 1:5){
        masterTable_data_fractionU[,i][ which(masterTable_data_fractionU[,i]=="Inf")] <-0
      }
      
      masterTable_data_fractionU = masterTable_data_fractionU %>% dplyr::filter(class== 'MZ' | class=="M" | class=="M-stable")
      M =  masterTable_data_fractionU %>% dplyr::filter(class=="M") %>% dplyr::mutate(V11=ifelse(class=="M","M",V11)) 
      
      Mstable =  masterTable_data_fractionU %>% dplyr::filter(class=="M-stable") %>% dplyr::mutate(V11=ifelse(class=="M-stable","M-stable",V11)) 
      othergenes = masterTable_data_fractionU %>% dplyr::filter(class=="MZ")
      masterTable_data_fractionU = rbind.data.frame(M,Mstable,othergenes)
      
            masterTable_data_fractionU_split = split(masterTable_data_fractionU,masterTable_data_fractionU$V11,T)
      
      
        increaseOvertime = reshape::melt(lapply(masterTable_data_fractionU_split,function(y) apply(y[,c(1:5)],2,function(x) mean(x,na.rm=T,trim=0.01))))
        increaseOvertime$time = times
        q = ggpubr::ggline(increaseOvertime,x='time',y='value',color='L1',ylab='Fraction uridylation') + theme_ameres(type = 'barplot') 
        print(q)
        
      
     
dev.off()  
pvals = cbind.data.frame(c('cluster1','cluster2','cluster3','cluster4','M-unastable'),c(pval_ks_1,pval_ks_2,pval_ks_3,pval_ks_4,pval_ks_Mstable))
colnames(pvals) = c('sample','pval_kstest')
write.table(pvals, "/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure5/data/pvals_mstable_translationalEfficiency.txt",
            sep="\t",quote = F, row.names = F)
TEdata_write = do.call(rbind.data.frame,masterTable_data_fractionU_split)
write.table(TEdata_write,"/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure5/data/teData.txt",sep="\t",quote = F,row.names = F)
write.table(increaseOvertime,"/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure5/data/meanFractionTE.txt",
            sep="\t",quote = F,row.names = F)
###### CAI value


masterTable_data_CAI = masterTable_data %>% dplyr::select(dplyr::contains("CAI"))
metaData = masterTable_data[,c(1:3)]
masterTable_data_CAI  = cbind.data.frame(masterTable_data_CAI,metaData)
masterTable_data_CAI = masterTable_data_CAI[complete.cases(masterTable_data_CAI),]
masterTable_data_CAI =  plyr::join(masterTable_data_CAI,MZgenes,by='external_gene_name')            
masterTable_data_CAI = masterTable_data_CAI[complete.cases(masterTable_data_CAI ),]
masterTable_data_CAI$MZrate = factor(masterTable_data_CAI$MZrate)


######## combining CAI with length of 3' UTRs... 


###### divide into 4 categories based on the CAI and length of UTR (same as what was done in Mishima et al )
masterTable_data = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures//figure5/data//ClassifiedGenes_masterTable.txt",stringsAsFactors = F, header = T,sep = "\t")

classifiedGenes = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/packages/data//countingWindows_classified_conversions.txt",
                             sep="\t",stringsAsFactors = F,header = T)
classifiedGenes$V4 = classifiedGenes$gene
CWconnection = read.table("//Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure4/data//finalEnds_intronsRemoved.bed",sep="\t",stringsAsFactors = F)
CWconnection =  plyr::join(CWconnection,classifiedGenes,by='V4')
CWconnection = CWconnection%>% dplyr::mutate(length = V3-V2)
CWconnection = CWconnection %>% dplyr::group_by(V4) %>% dplyr::mutate(meanLength = mean(V3-V2))
CWconnection = CWconnection[!duplicated(CWconnection$V4),] %>% dplyr::mutate(external_gene_name = V4) %>% dplyr::ungroup() %>%
  dplyr::select(c(external_gene_name, meanLength))
masterTable_data = masterTable_data %>% dplyr::select(c(external_gene_name,meanCAI,meanLength,class))
masterTable_data =  plyr::join(masterTable_data,MZgenes)

cutoff_CAI = median(masterTable_data$meanCAI,na.rm = T)
cutoff_Length = median(masterTable_data$meanLength,na.rm = T)



#################### chekcing for enrichment of genes................ #############################################

###### there are 4 comparisons for slow and 4 for fast genes
#### A - high CAI and  low UTR length - q4
#### B - high CAI anfd high UTR length - q2
#### C - low CAI and low UTR -q3
#### D - low CAI and high UTR - q1

masterTable_data$MZrate[is.na(masterTable_data$MZrate) ] <-'Diff'
masterTable_data = masterTable_data %>% dplyr::mutate(V11 = ifelse(is.na(V11),class,V11))
masterTable_data = masterTable_data[complete.cases(masterTable_data),]
classes_diff = c(1,2,3,4,"M","M-stable")
#### A 
q4_enrichment = c()
q4_depletion = c()
q4_fraction = c()

for(i in 1:length(classes_diff)){
  fast_A = nrow(masterTable_data %>% dplyr::filter(V11 == classes_diff[i] ) %>% dplyr::filter(meanCAI >= cutoff_CAI & meanLength <= cutoff_Length))
  allFast = nrow(masterTable_data %>% dplyr::filter(V11 == classes_diff[i] )) - (fast_A)
  
  nonFast_A = nrow(masterTable_data %>% dplyr::filter(V11 != classes_diff[i] ) %>% dplyr::filter(meanCAI >= cutoff_CAI & meanLength <= cutoff_Length))
  allNonFast =  nrow(masterTable_data %>% dplyr::filter(V11 != classes_diff[i] )) - nonFast_A
  
  q4_fraction = c(q4_fraction,fast_A/(fast_A+allFast))
  q4_depletion = c(q4_depletion,fisher.test(x = matrix(c(fast_A,nonFast_A,allFast,allNonFast),nrow = 2,byrow = T), alternative = 'less')$p.value ) ### depletion
  q4_enrichment = c(q4_enrichment,fisher.test(x = matrix(c(fast_A,nonFast_A,allFast,allNonFast),nrow = 2,byrow = T), alternative = 'greater')$p.value  )### enrichment
  
}

### B (q2)
q2_enrichment = c()
q2_depletion = c()
q2_fraction = c()

for(i in 1:length(classes_diff)){
        fast_A = nrow(masterTable_data %>% dplyr::filter(V11 ==classes_diff[i]  ) %>% dplyr::filter(meanCAI > cutoff_CAI & meanLength > cutoff_Length))
        allFast = nrow(masterTable_data %>% dplyr::filter(V11 == classes_diff[i] )) - (fast_A)
        
        nonFast_A = nrow(masterTable_data %>% dplyr::filter(V11 != classes_diff[i] ) %>% dplyr::filter(meanCAI > cutoff_CAI & meanLength > cutoff_Length))
        allNonFast =  nrow(masterTable_data %>% dplyr::filter(V11 != classes_diff[i] )) - nonFast_A
        q2_fraction = c(q2_fraction,fast_A/(fast_A+allFast))
        
        q2_depletion = c(q2_depletion,fisher.test(x = matrix(c(fast_A,nonFast_A,allFast,allNonFast),nrow = 2,byrow = T), alternative = 'less')$p.value)
        q2_enrichment = c(q2_enrichment,fisher.test(x = matrix(c(fast_A,nonFast_A,allFast,allNonFast),nrow = 2,byrow = T), alternative = 'greater')$p.value)

}

#### C q3 

q3_enrichment = c()
q3_depletion = c()
q3_fraction = c()

for(i in 1:length(classes_diff)){
  
        fast_A = nrow(masterTable_data %>% dplyr::filter(V11 == classes_diff[i] ) %>% dplyr::filter(meanCAI <= cutoff_CAI & meanLength <= cutoff_Length))
        allFast = nrow(masterTable_data %>% dplyr::filter(V11 == classes_diff[i] )) - (fast_A)
        
        nonFast_A = nrow(masterTable_data %>% dplyr::filter(V11 != classes_diff[i] ) %>% dplyr::filter(meanCAI <= cutoff_CAI & meanLength <= cutoff_Length))
        allNonFast =  nrow(masterTable_data %>% dplyr::filter(V11 != classes_diff[i] )) - nonFast_A
        q3_fraction = c(q3_fraction,fast_A/(fast_A+allFast))
        
        q3_depletion = c(q3_depletion, fisher.test(x = matrix(c(fast_A,nonFast_A,allFast,allNonFast),nrow = 2,byrow = T), alternative = 'less')$p.value)
        q3_enrichment = c(q3_enrichment,fisher.test(x = matrix(c(fast_A,nonFast_A,allFast,allNonFast),nrow = 2,byrow = T), alternative = 'greater')$p.value)

}


##### D q1 


q1_enrichment = c()
q1_depletion = c()
q1_fraction = c()

for(i in 1:length(classes_diff)){
        fast_A = nrow(masterTable_data %>% dplyr::filter(V11 == classes_diff[i] ) %>% dplyr::filter(meanCAI < cutoff_CAI & meanLength > cutoff_Length))
        allFast = nrow(masterTable_data %>% dplyr::filter(V11 == classes_diff[i]  )) - (fast_A)
        
        nonFast_A = nrow(masterTable_data %>% dplyr::filter(V11 != classes_diff[i]  ) %>% dplyr::filter(meanCAI < cutoff_CAI & meanLength > cutoff_Length))
        allNonFast =  nrow(masterTable_data %>% dplyr::filter(V11 != classes_diff[i]  )) - nonFast_A
        q1_fraction = c(q1_fraction,fast_A/(fast_A+allFast))
        
        q1_depletion = c(q1_depletion,fisher.test(x = matrix(c(fast_A,nonFast_A,allFast,allNonFast),nrow = 2,byrow = T), alternative = 'less')$p.value)
        q1_enrichment = c(q1_enrichment,fisher.test(x = matrix(c(fast_A,nonFast_A,allFast,allNonFast),nrow = 2,byrow = T), alternative = 'greater')$p.value)
}

fractions_samples = data.frame(q1 = q1_fraction,q2 = q2_fraction, q3 = q3_fraction, q4 = q4_fraction)
row.names(fractions_samples) = classes_diff

depletions_samples = data.frame(q1 = q1_depletion,q2 = q2_depletion, q3 = q3_depletion,q4 = q4_depletion)
row.names(depletions_samples) = classes_diff

enrichments_samples = data.frame(q1 = q1_enrichment,q2 = q2_enrichment, q3 = q3_enrichment,q4 = q4_enrichment)
row.names(enrichments_samples) = classes_diff


fast = masterTable_data %>% dplyr::filter(V11 == 4)  
fast = fast[complete.cases(fast),] %>% dplyr::mutate(meanLength = log10(meanLength))


####### plotting the same plits but also adding miR430targets...
pdf("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure4/plots//CAI_MZgenes_includingMiRNA.pdf", height = 4, width=4)

        
        allTargets = read.table('/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure4/data//mir430_definedTargets.txt',
                                sep='\t',stringsAsFactors = F)
        allTargets = allTargets %>% dplyr::filter(UTRsites==T & downReg == T)
        
        p =   ggplot2::ggplot(data=fast,aes(x=meanCAI,y=meanLength)) + ggrastr::geom_point_rast(alpha=0.2) + geom_hline(yintercept =log10(cutoff_Length) ,linetype = 'dashed',color='black') + 
          geom_vline(xintercept =cutoff_CAI ,linetype = 'dashed',color='black') +theme_cowplot() + xlab('CAI') + ylab('3 UTR length')
        p = p + theme_ameres(type = 'barplot') + ggtitle(paste0('4 genes',nrow(fast)))
        p + ggrastr::geom_point_rast(data = fast[fast$external_gene_name %in%  allTargets$gene,],aes(x=meanCAI,y=meanLength),colour='purple')
        
       
        fast = masterTable_data %>% dplyr::filter(V11 == 3)  
        fast = fast[complete.cases(fast),] %>% dplyr::mutate(meanLength = log10(meanLength))
        
        p =   ggplot2::ggplot(data=fast,aes(x=meanCAI,y=meanLength)) + ggrastr::geom_point_rast(alpha=0.2) + geom_hline(yintercept =log10(cutoff_Length) ,linetype = 'dashed',color='black') + 
          geom_vline(xintercept =cutoff_CAI ,linetype = 'dashed',color='black') +theme_cowplot() + xlab('CAI') + ylab('3 UTR length')
        p = p + theme_ameres(type = 'barplot') + ggtitle(paste0('3 genes',nrow(fast)))
        p + ggrastr::geom_point_rast(data = fast[fast$external_gene_name %in%  allTargets$gene,],aes(x=meanCAI,y=meanLength),colour='purple')
        
         
        
        fast = masterTable_data %>% dplyr::filter(V11 == 2)  
        fast = fast[complete.cases(fast),] %>% dplyr::mutate(meanLength = log10(meanLength))
        
        p =   ggplot2::ggplot(data=fast,aes(x=meanCAI,y=meanLength)) + ggrastr::geom_point_rast(alpha=0.2) + geom_hline(yintercept =log10(cutoff_Length) ,linetype = 'dashed',color='black') + 
          geom_vline(xintercept =cutoff_CAI ,linetype = 'dashed',color='black') +theme_cowplot() + xlab('CAI') + ylab('3 UTR length')
        p = p + theme_ameres(type = 'barplot') + ggtitle(paste0('2 genes',nrow(fast)))
        p + ggrastr::geom_point_rast(data = fast[fast$external_gene_name %in%  allTargets$gene,],aes(x=meanCAI,y=meanLength),colour='purple')
        
        
        slow = masterTable_data %>% dplyr::filter(V11 == 1)  
        slow = slow[complete.cases(slow),] %>% dplyr::mutate(meanLength = log10(meanLength))
        
        p =   ggplot2::ggplot(data=slow,aes(x=meanCAI,y=meanLength)) + ggrastr::geom_point_rast(alpha=0.2) + geom_hline(yintercept =log10(cutoff_Length) ,linetype = 'dashed',color='black') + 
          geom_vline(xintercept =cutoff_CAI ,linetype = 'dashed',color='black') +theme_cowplot() + xlab('CAI') + ylab('3 UTR length')
        p = p + theme_ameres(type = 'barplot') + ggtitle(paste0('slow genes',nrow(slow)))
        p+ ggrastr::geom_point_rast(data = slow[slow$external_gene_name %in%  allTargets$gene,],aes(x=meanCAI,y=meanLength),colour='purple')
        
         M = masterTable_data %>% dplyr::filter(V11 == 'M')  
         M = M[complete.cases(M),] %>% dplyr::mutate(meanLength = log10(meanLength))

        p =   ggplot2::ggplot(data=M,aes(x=meanCAI,y=meanLength)) + ggrastr::geom_point_rast(alpha=0.2) + geom_hline(yintercept =log10(cutoff_Length) ,linetype = 'dashed',color='black') +
          geom_vline(xintercept =cutoff_CAI ,linetype = 'dashed',color='black') +theme_cowplot() + xlab('CAI') + ylab('3 UTR length')
        p = p + theme_ameres(type = 'barplot') + ggtitle(paste0('M genes',nrow(M)))
        p+ ggrastr::geom_point_rast(data = M[M$external_gene_name %in%  allTargets$gene,],aes(x=meanCAI,y=meanLength),colour='purple')


        Mstable = masterTable_data %>% dplyr::filter(V11 == 'M-stable')
        Mstable= Mstable[complete.cases(Mstable),] %>% dplyr::mutate(meanLength = log10(meanLength))

        p =   ggplot2::ggplot(data=Mstable,aes(x=meanCAI,y=meanLength)) + ggrastr::geom_point_rast(alpha=0.2) + geom_hline(yintercept =log10(cutoff_Length) ,linetype = 'dashed',color='black') +
          geom_vline(xintercept =cutoff_CAI ,linetype = 'dashed',color='black') +theme_cowplot() + xlab('CAI') + ylab('3 UTR length')
        p = p + theme_ameres(type = 'barplot') + ggtitle(paste0('Mstable genes',nrow(Mstable)))
        p+ ggrastr::geom_point_rast(data = Mstable[Mstable$external_gene_name %in%  allTargets$gene,],aes(x=meanCAI,y=meanLength),colour='purple')

dev.off()