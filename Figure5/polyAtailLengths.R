#### plotting polyA tail lengths of different classes of MZ genes

library(ggplot2)
library(dplyr)

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

masterTable_data = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures//figure5/data//ClassifiedGenes_masterTable.txt",
                              stringsAsFactors = F, header = T,sep = "\t")
masterTable_data =  masterTable_data %>% dplyr::filter(description != 'MZ from Z') %>% dplyr::filter(description != 'MZ from M')
#### rate of increase of MZ genes... 
MZgenes = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure4/data/clusteredBasedOnTCaccumulation.bed",sep="\t",stringsAsFactors = F)
MZgenes = MZgenes %>% dplyr::mutate(gene=V10,class=V12) %>% dplyr::select('gene','V11','class')
MZgenes = MZgenes %>% dplyr::mutate(external_gene_name = gene,MZrate = class) %>% dplyr::select('external_gene_name','MZrate','V11')

masterTable_data_fractionU = masterTable_data %>% dplyr::select(dplyr::contains("polyA"))
metaData = masterTable_data[,c(1:3)]
masterTable_data_fractionU  = cbind.data.frame(masterTable_data_fractionU,metaData)
masterTable_data_fractionU =  plyr::join(masterTable_data_fractionU,MZgenes,by='external_gene_name')            
masterTable_data_fractionU$MZrate = factor(masterTable_data_fractionU$MZrate)
timepoints = colnames(masterTable_data_fractionU)[1:6]
times = c(0,1,2,4,6,8)
dir.create('/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure5//plots/')
pdf("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure5//plots/polyAtails_ClustersMZgenes.pdf", height = 4, width=4)

for(i in 1:length(timepoints)){
  Mgenes = masterTable_data_fractionU  %>%
    dplyr::mutate(fractionU=get(timepoints[i])) %>% dplyr::filter(class== "M") %>% dplyr::mutate(V11='M')%>% dplyr::select(c('MZrate','fractionU','V11'))
  Mstablegenes = masterTable_data_fractionU  %>%
    dplyr::mutate(fractionU=get(timepoints[i])) %>% dplyr::filter(class== "M-stable") %>% dplyr::mutate(V11='M-stable')%>% dplyr::select(c('MZrate','fractionU','V11'))
  
  
  b= masterTable_data_fractionU %>%
    dplyr::mutate(fractionU=get(timepoints[i])) %>% dplyr::select(c('MZrate','fractionU','V11'))
  
  b = rbind.data.frame(Mgenes,Mstablegenes,b)
  b= b[complete.cases(b$V11) ,]
  
  b_split = split(b,b$V11)
  #names(b_split) = paste0('Cluster ',c(1:4))
  pval_ks= ks.test(b_split$`1`$fractionU,b_split$`4`$fractionU)$p.value
  
  b_split= lapply(b_split,function(x) x[order(x[,'fractionU'],decreasing = F),])
  b_split = lapply(b_split,function(x) x %>% dplyr::mutate(fractionU_cumulative  = cumsum(fractionU)))
  b_split=lapply(b_split,function(x) x%>% dplyr::mutate(cumsum_frac = fractionU_cumulative/max(fractionU_cumulative)))
  
  b_split_2hpf = do.call(rbind,b_split)
  b_split_2hpf = b_split_2hpf %>% dplyr::group_by(V11) %>% dplyr::mutate(n=n())  %>% dplyr::mutate(category_sample = paste0(V11,"(n=",n,")"))
  b_split_2hpf$V11  = as.factor(b_split_2hpf$V11)
  p =  ggplot(b_split_2hpf, aes(log10(fractionU), color=V11,group=V11)) +
    stat_ecdf() + theme_cowplot()
  p = p + annotate("text", x = 1, y = 0.5, label = pval_ks) + theme_ameres(type = 'barplot')+ coord_fixed()
  #print(p)
  
  
  
}

b = masterTable_data_fractionU  %>% dplyr::mutate(fractionU = X6hpf_polyA/X2hpf_polyA) %>% dplyr::filter(class=="M" | class=="MZ" | class == "M-stable") %>% dplyr::select(c('class','V11','fractionU','external_gene_name'))
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
  theme_ameres(type = 'barplot') + ylab('Cumulative fraction') + xlab ('poly(A) tail length (log2)')
print(p)

p = ggplot(b_split_2hpf,aes(y=log2(fractionU),x=category_sample)) + geom_boxplot(outlier.shape = NA) + ggpubr::stat_compare_means(ref.group = 'M(n=1013)') + theme_cowplot()
p + theme_ameres(type = 'barplot')


masterTable_data_fractionU = masterTable_data_fractionU %>% dplyr::filter(class== 'MZ' | class=="M" | class=="M-stable")

M =  masterTable_data_fractionU %>% dplyr::filter(class=="M") %>% dplyr::mutate(V11=ifelse(class=="M","M",V11)) 
Mstable =  masterTable_data_fractionU %>% dplyr::filter(class=="M-stable") %>% dplyr::mutate(V11=ifelse(class=="M-stable","M-stable",V11)) 
othergenes = masterTable_data_fractionU %>% dplyr::filter(class=="MZ")
masterTable_data_fractionU = rbind.data.frame(M,Mstable,othergenes)
masterTable_data_fractionU_split = split(masterTable_data_fractionU,masterTable_data_fractionU$V11,T)
increaseOvertime = reshape::melt(lapply(masterTable_data_fractionU_split,function(y) apply(y[,c(1:6)],2,function(x) mean(x,na.rm=T))))
increaseOvertime$time = times
q = ggpubr::ggline(increaseOvertime,x='time',y='value',color='L1',ylab='poly(A) tail length') + theme_ameres(type = 'barplot') 
print(q)

dev.off()
polyAlengths = do.call(rbind.data.frame,masterTable_data_fractionU)
write.table(polyAlengths,"/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure5/data/polyAlengths.txt",
            sep="\t",quote = F,row.names = F)
write.table(increaseOvertime,"/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure5/data/meanpolyAtailLength.txt",
            sep="\t",quote = F,row.names = F)
pvals = cbind.data.frame(c('cluster1','cluster2','cluster3','cluster4','M-unastable'),c(pval_ks_1,pval_ks_2,pval_ks_3,pval_ks_4,pval_ks_Mstable))
colnames(pvals) = c('sample','pval_kstest')
write.table(pvals, "/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure5/data/pvals_mstable_polyAtails.txt",
            sep="\t",quote = F, row.names = F)
