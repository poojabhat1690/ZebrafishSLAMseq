library(reshape)
theme_ameres <- function (type) { 
  
  types_plot = c("boxplot","barplot","violin")
  if(type %in% types_plot ==T)
  {
    
    if(type == "boxplot"){
      theme(legend.title=element_blank(),axis.text.x = ggplot2::element_text(margin=ggplot2::margin(10,15,10,15,"pt"),size = 15),axis.text.y = element_text(margin=ggplot2::margin(5,15,10,5,"pt"),size = 15), axis.ticks.length = unit(-0.25 , "cm"),legend.position="bottom",axis.line = element_line(colour = "black", lineend = "round"), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18))   
    }
    
    if(type == "barplot"){
      theme(legend.title=element_blank(),axis.text.x = ggplot2::element_text(margin=ggplot2::margin(10,15,10,15,"pt"),size = 15),axis.text.y = element_text(margin=ggplot2::margin(5,15,10,5,"pt"),size = 15), axis.ticks.length = unit(-0.25 , "cm"),legend.position="bottom",axis.line = element_line(colour = "black", lineend = "round"), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18))   
    }
    
    
  }
  
}


heynEtal = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/packages/data/figure2/heynEtal_s1.txt",
                      sep="\t",stringsAsFactors = F, header = T)
ensembl_geneNames_dr7 = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/packages/data/figure2/allGenes_dr7.txt",
                                   sep="\t",stringsAsFactors = F, header = T)
ensembl_geneNames_dr7 = ensembl_geneNames_dr7 %>% dplyr::distinct(ensembl_gene_id,.keep_all=T)
heynEtal = plyr::join(heynEtal,ensembl_geneNames_dr7)
heynEtal = heynEtal %>% dplyr::filter(gene_biotype=='protein_coding')

dr11 = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/packages/data/figure2//allGenes_dr11.txt",
                  sep="\t",stringsAsFactors = F, header = T)

`%!in%` = Negate(`%in%`)
heynEtal_dr11 = heynEtal[heynEtal$external_gene_name %in% dr11$external_gene_name,]
heynEtal_dr11 = heynEtal_dr11[-which(heynEtal_dr11$external_gene_name == ""),]
heynEtal_dr11 = heynEtal_dr11[-which(heynEtal_dr11$Gene_name  == ""),]


classifiedGenes = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/packages/data/countingWindows_classified_conversions.txt",sep="\t",stringsAsFactors = F, header = T)
classifiedGenes = classifiedGenes %>% dplyr::mutate(external_gene_name = gene) %>% dplyr::select(c('external_gene_name','class'))

RPMs_all  = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/packages/data/RPM_allCws.txt",
                       sep="\t",stringsAsFactors = F, header = T)
RPMs_all = RPMs_all %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Inj_R3_TP9), sum, na.rm = TRUE)

splitReplicates = function(dataFrameToSplit,condition,metadata_add){
  dataFrameToSplit_condition = dataFrameToSplit[,grep(condition,colnames(dataFrameToSplit ))]
  # dataFrameToSplit_condition_R1 = dataFrameToSplit_condition[,grep("R1",colnames(dataFrameToSplit_condition ))]
  dataFrameToSplit_condition_R2 = dataFrameToSplit_condition[,grep("R2",colnames(dataFrameToSplit_condition ))]
  dataFrameToSplit_condition_R3 = dataFrameToSplit_condition[,grep("R3",colnames(dataFrameToSplit_condition ))]
  mean_repl = (dataFrameToSplit_condition_R2+dataFrameToSplit_condition_R3)/2
  #dataFrameToSplit_condition_R1 = cbind.data.frame(dataFrameToSplit_condition_R1,metadata_add)
  dataFrameToSplit_condition_R2 = cbind.data.frame(dataFrameToSplit_condition_R2,metadata_add)
  dataFrameToSplit_condition_R3 = cbind.data.frame(dataFrameToSplit_condition_R3,metadata_add)
  mean_repl = cbind.data.frame(mean_repl,metadata_add)
  splitReplicates = list(dataFrameToSplit_condition_R2,dataFrameToSplit_condition_R3,mean_repl)
  names(splitReplicates) = c("R2","R3","mean")
  return(splitReplicates)
}

RPMs_mean = splitReplicates(RPMs_all,condition = "Inj",metadata_add = RPMs_all[,1])[[3]]
heynEtal_dr11 = heynEtal_dr11[heynEtal_dr11$external_gene_name %in% RPMs_mean$name,]

heynEtal_dr11 = plyr::join(heynEtal_dr11,classifiedGenes)

noclassifiedInSLAMseq = heynEtal_dr11[!complete.cases(heynEtal_dr11),] ### is this because of low coverage?

RPMs_mean$external_gene_name = RPMs_mean$name
noclassifiedInSLAMseq = plyr::join(noclassifiedInSLAMseq,RPMs_mean)

RPMS_nonDetected = noclassifiedInSLAMseq %>% dplyr::select(dplyr::contains("Inj_"))
colnames(RPMS_nonDetected) = c('0.75','2','2.5','3','3.5','4','4.5','5','5.5')
RPMS_nonDetected =  reshape::melt(RPMS_nonDetected)
RPMS_nonDetected$variable = as.character(RPMS_nonDetected$variable)
RPMS_nonDetected$variable  = as.numeric(RPMS_nonDetected$variable)
ggpubr::ggboxplot(RPMS_nonDetected,x='variable',y='value',numeric.x.axis=T)
library(ggplot2)
p= ggplot(RPMS_nonDetected,aes(x=variable,y=value, group=variable)) + geom_boxplot()  + cowplot::theme_cowplot()
p

#### harvey et al... 
classifiedGenes = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/packages/data/countingWindows_classified_conversions.txt",sep="\t",stringsAsFactors = F, header = T)
classifiedGenes = classifiedGenes %>% dplyr::mutate(external_gene_name = gene) %>% dplyr::select(c('external_gene_name','class'))

harveyEtal = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/packages/data/figure2/harveyEtal.txt",
                        sep="\t", stringsAsFactors = F, header = T)
harveyEtal_Z = harveyEtal %>% dplyr::filter(class == "Z")
harveyEtal_MZ = harveyEtal %>% dplyr::filter(class == "MZ")
harveyEtal_Z = data.frame(ensembl_gene_id = unlist(strsplit(harveyEtal_Z$Gene.ID,",",T)),class = "Z") %>% dplyr::distinct(ensembl_gene_id,.keep_all=T)
harveyEtal_MZ = data.frame(ensembl_gene_id = unlist(strsplit(harveyEtal_MZ$Gene.ID,",",T)),class = "MZ")%>% dplyr::distinct(ensembl_gene_id,.keep_all=T)
totalHarvey = rbind.data.frame(harveyEtal_Z,harveyEtal_MZ)

dr11 = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/packages/data/figure2/allGenes_dr11.txt",
                  sep="\t",stringsAsFactors = F, header = T)

totalHarvey = plyr::join(totalHarvey,dr11) %>% dplyr::distinct(ensembl_gene_id,.keep_all=T)



RPMs_all  = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/packages/data/RPM_allCws.txt",
                       sep="\t",stringsAsFactors = F, header = T)
RPMs_all = RPMs_all %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Inj_R3_TP9), sum, na.rm = TRUE)

splitReplicates = function(dataFrameToSplit,condition,metadata_add){
  dataFrameToSplit_condition = dataFrameToSplit[,grep(condition,colnames(dataFrameToSplit ))]
  # dataFrameToSplit_condition_R1 = dataFrameToSplit_condition[,grep("R1",colnames(dataFrameToSplit_condition ))]
  dataFrameToSplit_condition_R2 = dataFrameToSplit_condition[,grep("R2",colnames(dataFrameToSplit_condition ))]
  dataFrameToSplit_condition_R3 = dataFrameToSplit_condition[,grep("R3",colnames(dataFrameToSplit_condition ))]
  mean_repl = (dataFrameToSplit_condition_R2+dataFrameToSplit_condition_R3)/2
  #dataFrameToSplit_condition_R1 = cbind.data.frame(dataFrameToSplit_condition_R1,metadata_add)
  dataFrameToSplit_condition_R2 = cbind.data.frame(dataFrameToSplit_condition_R2,metadata_add)
  dataFrameToSplit_condition_R3 = cbind.data.frame(dataFrameToSplit_condition_R3,metadata_add)
  mean_repl = cbind.data.frame(mean_repl,metadata_add)
  splitReplicates = list(dataFrameToSplit_condition_R2,dataFrameToSplit_condition_R3,mean_repl)
  names(splitReplicates) = c("R2","R3","mean")
  return(splitReplicates)
}

RPMs_mean = splitReplicates(RPMs_all,condition = "Inj",metadata_add = RPMs_all[,1])[[3]]
heynEtal_dr11 = heynEtal_dr11[heynEtal_dr11$external_gene_name %in% RPMs_mean$name,]

heynEtal_dr11 = plyr::join(heynEtal_dr11,classifiedGenes)

noclassifiedInSLAMseq = heynEtal_dr11[!complete.cases(heynEtal_dr11),] ### is this because of low coverage?

RPMs_mean$external_gene_name = RPMs_mean$name
  RPMs_mean =  plyr::join(RPMs_mean,dr11)

totalHarvey = totalHarvey[totalHarvey$ensembl_gene_id %in% RPMs_mean$ensembl_gene_id,]

totalHarvey =  plyr::join(totalHarvey,classifiedGenes,by='external_gene_name')
table(totalHarvey[,8])
notIdInSLAMseq = totalHarvey[(is.na(totalHarvey[,8])),]
notIdInSLAMseq= plyr::join(notIdInSLAMseq,RPMs_mean,by='ensembl_gene_id')
notIdInSLAMseq = notIdInSLAMseq[!duplicated(notIdInSLAMseq$ensembl_gene_id),]

classifiedGenes = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/packages/data/countingWindows_classified_conversions.txt",sep="\t",stringsAsFactors = F, header = T)
totalHarvey_M = totalHarvey[which(totalHarvey[,8] == "M" ),]
totalHarvey_M_stable = totalHarvey[which(totalHarvey[,8] == "M-stable" ),]

misClassified_SLAMseq_M = classifiedGenes[classifiedGenes$gene %in% totalHarvey_M$external_gene_name,]
misClassified_SLAMseq_M = misClassified_SLAMseq_M %>% dplyr::select(dplyr::contains("Inj"))
colnames(misClassified_SLAMseq_M) = c(0.75,2,2.5,3,3.5,4,4.5,5,5.5)
misClassified_SLAMseq_M = melt(misClassified_SLAMseq_M)


misClassified_SLAMseq_Mstable = classifiedGenes[classifiedGenes$gene %in% totalHarvey_M_stable$external_gene_name,]
misClassified_SLAMseq_Mstable = misClassified_SLAMseq_Mstable %>% dplyr::select(dplyr::contains("Inj"))
colnames(misClassified_SLAMseq_Mstable) = c(0.75,2,2.5,3,3.5,4,4.5,5,5.5)
misClassified_SLAMseq_Mstable = melt(misClassified_SLAMseq_Mstable)
misClassified_SLAMseq_M$variable = as.numeric(as.character(misClassified_SLAMseq_M$variable))
misClassified_SLAMseq_Mstable$variable = as.numeric(as.character(misClassified_SLAMseq_Mstable$variable))

pdf("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/packages/plots/figure2/tcConversionRates_misclassifiedInHarveyEtAl.pdf",height = 3, width=6)
      p= ggplot(misClassified_SLAMseq_M,aes(x=as.numeric(variable),y=value, group=variable)) + geom_boxplot()  + cowplot::theme_cowplot()
      p = p + ylim(c(0,0.15)) + theme_ameres(type = 'barplot') + ylab('TC conversion rate')+xlab('time (hpf)') + ggtitle(paste0('mis-classified as maternal=',nrow(misClassified_SLAMseq_M)/9))
      p + scale_x_continuous(breaks=c(0.75,2,2.5,3,3.5,4,4.5,5,5.5),labels=c(0.75,2,2.5,3,3.5,4,4.5,5,5.5))
      
      p= ggplot(misClassified_SLAMseq_Mstable,aes(x=variable,y=value, group=variable)) + geom_boxplot(outlier.shape=NA)  + cowplot::theme_cowplot()
      p = p + ylim(c(0,0.15))+ theme_ameres(type = 'barplot') + ylab('TC conversion rate')+xlab('time (hpf)')+ ggtitle(paste0('mis-classified as maternal=',nrow(misClassified_SLAMseq_Mstable)/9))
      p + scale_x_continuous(breaks=c(0.75,2,2.5,3,3.5,4,4.5,5,5.5),labels=c(0.75,2,2.5,3,3.5,4,4.5,5,5.5))
      
dev.off()


#### also comparing with lee et al 2013, where they used inhibition using ActD to identify zygotic transcripts.... 


classifiedGenes = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/packages/data/countingWindows_classified_conversions.txt",sep="\t",stringsAsFactors = F, header = T)
classifiedGenes = classifiedGenes %>% dplyr::mutate(external_gene_name = gene) %>% dplyr::select(c('external_gene_name','class'))


########### processing lee et al... 

leeEtal = read.delim("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/packages/data/figure2/leeEtal_complete.txt",
                     sep="\t", stringsAsFactors = F, header = T,skip = 23)
t4 = leeEtal %>% dplyr::filter(WT4_txed == 'E' | WT4_txed =='I' | WT4_txed == "EI")
t6 = leeEtal %>% dplyr::filter(WT6_txed == 'E' | WT6_txed =='I' | WT6_txed == "EI")
u1u2 = leeEtal %>% dplyr::filter(U1U2_txed == 'E' | U1U2_txed =='I' | U1U2_txed == "EI")
leeEtal = rbind.data.frame(t4,t6, u1u2)
leeEtal = leeEtal[!duplicated(leeEtal$Symbol),]

leeEtal = leeEtal %>% dplyr::mutate(external_gene_name = Symbol, ensembl_gene_id = Gene_id) %>%  
  dplyr::select(c('external_gene_name','ensembl_gene_id','Maternal_contr'))
dr11 = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/packages/data/figure2/allGenes_dr11.txt",
                  sep="\t",stringsAsFactors = F, header = T)
leeEtal  = plyr::join(leeEtal, dr11,by='external_gene_name') 

leeEtal= leeEtal[!duplicated(leeEtal$external_gene_name),]
leeEtal = leeEtal[complete.cases(leeEtal),] #### using this only to get the genes which we can detect in dr11
        
leeEtal = plyr::join(leeEtal, classifiedGenes)

notClassifiedSLAMseq = leeEtal[!complete.cases(leeEtal$class),]

RPMs_all  = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/packages/data/RPM_allCws.txt",
                       sep="\t",stringsAsFactors = F, header = T)
RPMs_all = RPMs_all %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Inj_R3_TP9), sum, na.rm = TRUE)


RPMs_mean = splitReplicates(RPMs_all,condition = "Inj",metadata_add = RPMs_all[,1])[[3]]


RPMs_mean$external_gene_name = RPMs_mean$name
RPMs_mean =  plyr::join(RPMs_mean,dr11)
notClassifiedSLAMseq = RPMs_mean[RPMs_mean$external_gene_name %in% notClassifiedSLAMseq$external_gene_name ,]
notClassifiedSLAMseq = notClassifiedSLAMseq[!duplicated(notClassifiedSLAMseq$external_gene_name),]


notClassifiedSLAMseq_RPMS = notClassifiedSLAMseq %>% dplyr::select(dplyr::contains("Inj_"))
colnames(notClassifiedSLAMseq_RPMS) = c('0.75','2','2.5','3','3.5','4','4.5','5','5.5')
notClassifiedSLAMseq_RPMS =  reshape::melt(notClassifiedSLAMseq_RPMS)
notClassifiedSLAMseq_RPMS$variable = as.character(notClassifiedSLAMseq_RPMS$variable)
notClassifiedSLAMseq_RPMS$variable  = as.numeric(notClassifiedSLAMseq_RPMS$variable)
notClassifiedSLAMseq_RPMS$value = log10(notClassifiedSLAMseq_RPMS$value+1)
ggpubr::ggboxplot(notClassifiedSLAMseq_RPMS,x='variable',y='value',numeric.x.axis=T)
library(ggplot2)

      
      p= ggplot(notClassifiedSLAMseq_RPMS,aes(x=variable,y=value, group=variable)) + geom_boxplot(outlier.shape = NA)  + cowplot::theme_cowplot()
      p + theme_ameres(type = 'barplot') + ylab("RPM")+ scale_x_continuous(breaks=c(0.75,2,2.5,3,3.5,4,4.5,5,5.5),labels=c(0.75,2,2.5,3,3.5,4,4.5,5,5.5))



##### ribominus RNAseq... 

#TPM_RNAseq = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/externalDatat/fromAndi_ribo0VsPolyA/TPM_STAR_dr11.txt",sep="\t",stringsAsFactors = F, header = T)
      TPM_RNAseq = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/packages/data/figure2/TPM_STAR_dr11.txt",sep="\t",stringsAsFactors = F, header = T)
      
      TPM_RNAseq$ensembl_gene_id = row.names(TPM_RNAseq)

notClassifiedSLAMseq =  plyr::join(notClassifiedSLAMseq,TPM_RNAseq)
notClassifiedSLAMseq = notClassifiedSLAMseq %>% dplyr::mutate(polyAT4 = (X49986_wt.2.T4.polyA + X49982_wt.1.T4.polyA + X49990_wt.3.T4.polyA)/3, ribominusT4 =  (X50018_wt.3.T4.RiboCop +X50014_wt.2.T4.RiboCop + X50010_wt.1.T4.RiboCop)/3)
notClassifiedSLAMseq = notClassifiedSLAMseq %>% dplyr::mutate(polyAT4 = log10(polyAT4+1), ribominusT4 = log10(ribominusT4+1))%>% 
  dplyr::select(c('polyAT4','ribominusT4'))


`%!in%` = Negate(`%in%`)


TPM_RNAseq_other = TPM_RNAseq[TPM_RNAseq$ensembl_gene_id %!in%  notClassifiedSLAMseq$ensembl_gene_id,]
TPM_RNAseq_other = TPM_RNAseq_other %>% dplyr::mutate(polyAT4 = (X49986_wt.2.T4.polyA + X49982_wt.1.T4.polyA + X49990_wt.3.T4.polyA)/3, ribominusT4 =  (X50018_wt.3.T4.RiboCop +X50014_wt.2.T4.RiboCop + X50010_wt.1.T4.RiboCop)/3)
TPM_RNAseq_other = TPM_RNAseq_other %>% dplyr::mutate(polyAT4 = log10(polyAT4+1), ribominusT4 = log10(ribominusT4+1))%>% 
  dplyr::select(c('polyAT4','ribominusT4'))


notClassifiedSLAMseq_melt =  reshape::melt(notClassifiedSLAMseq) %>% dplyr::mutate(type='notClassified')
TPM_RNAseq_other_melt = reshape::melt(TPM_RNAseq_other) %>% dplyr::mutate(type = 'otherTPM')
totalRPMs = rbind.data.frame(notClassifiedSLAMseq_melt,TPM_RNAseq_other_melt)

pdf("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/packages/plots/figure2/RPMundetectedByLeeEtal.pdf",height = 6, width=4)

  ggpubr::ggboxplot(totalRPMs,x='type',y='value',facet.by = 'variable', outlier.shape=NA, ylab = 'log10 (TPM)',ylim = c(0,3.5),bxp.errorbar = TRUE, fill='type', palette = 'Set1') + ggpubr::stat_compare_means(method = "wilcox.test",
   label.y = 3.3) + stat_summary(fun.data = function(x) data.frame(y=3, label = paste("Mean=",mean(x))), geom="text")  + theme_ameres(type = 'barplot')

dev.off()


##### compare what we identify to these... 

heynEtal = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/packages/data/figure2/heynEtal_s1.txt",
                      sep="\t",stringsAsFactors = F, header = T)
ensembl_geneNames_dr7 = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/packages/data/figure2//allGenes_dr7.txt",
                                   sep="\t",stringsAsFactors = F, header = T)
ensembl_geneNames_dr7 = ensembl_geneNames_dr7 %>% dplyr::distinct(ensembl_gene_id,.keep_all=T)
heynEtal = plyr::join(heynEtal,ensembl_geneNames_dr7)
heynEtal = heynEtal %>% dplyr::select(c('ensembl_gene_id')) 

harveyEtal = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/packages/data/figure2/harveyEtal.txt",
                        sep="\t", stringsAsFactors = F, header = T)
harveyEtal_Z = harveyEtal %>% dplyr::filter(class == "Z")
harveyEtal_MZ = harveyEtal %>% dplyr::filter(class == "MZ")
harveyEtal_Z = data.frame(ensembl_gene_id = unlist(strsplit(harveyEtal_Z$Gene.ID,",",T)),class = "Z") %>% dplyr::distinct(ensembl_gene_id,.keep_all=T)
harveyEtal_MZ = data.frame(ensembl_gene_id = unlist(strsplit(harveyEtal_MZ$Gene.ID,",",T)),class = "MZ")%>% dplyr::distinct(ensembl_gene_id,.keep_all=T)
totalHarvey = rbind.data.frame(harveyEtal_Z,harveyEtal_MZ) %>% dplyr::select('ensembl_gene_id')

leeEtal = read.delim("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/packages/data/figure2/leeEtal_complete.txt",
                     sep="\t", stringsAsFactors = F, header = T,skip = 23)
t4 = leeEtal %>% dplyr::filter(WT4_txed == 'E' | WT4_txed =='I' | WT4_txed == "EI")
t6 = leeEtal %>% dplyr::filter(WT6_txed == 'E' | WT6_txed =='I' | WT6_txed == "EI")
u1u2 = leeEtal %>% dplyr::filter(U1U2_txed == 'E' | U1U2_txed =='I' | U1U2_txed == "EI")
leeEtal = rbind.data.frame(t4,t6, u1u2)
leeEtal = leeEtal[!duplicated(leeEtal$Symbol),]

leeEtal = leeEtal %>% dplyr::mutate(external_gene_name = Symbol, ensembl_gene_id = Gene_id) %>%  
  dplyr::select(c('ensembl_gene_id'))

leeEtal_write = leeEtal
leeEtal_write$sample = 'Lee et.al 2013'

heynEtal_write = heynEtal
heynEtal_write$sample = 'Heyn et.al 2013'

totalHarvey_write = totalHarvey
totalHarvey_write$sample = 'Harvey et al 2014 '
total_write = rbind.data.frame(leeEtal_write,heynEtal_write,totalHarvey_write)
write.table(total_write,file = "/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/tables/otherStudies.txt",
            sep="\t",row.names = F,quote = F)


totalOtherStudies = rbind.data.frame(leeEtal,heynEtal, totalHarvey)
totalOtherStudies = totalOtherStudies[!duplicated(totalOtherStudies),]
totalOtherStudies = data.frame(ensembl_gene_id = totalOtherStudies)

##### checking how many of these are in dr11... 
dr11 = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/packages/data/figure2//allGenes_dr11.txt",
                  sep="\t",stringsAsFactors = F, header = T)
totalOtherStudies = plyr::join(totalOtherStudies, dr11)
 totalOtherStudies  = totalOtherStudies[!duplicated(totalOtherStudies$ensembl_gene_id),]
 totalOtherStudies  = totalOtherStudies[complete.cases(totalOtherStudies),]
 
 
 classifiedGenes = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/packages/data/countingWindows_classified_conversions.txt",sep="\t",stringsAsFactors = F, header = T)
 classifiedGenes = classifiedGenes %>% dplyr::mutate(external_gene_name = gene) %>% dplyr::select(c('external_gene_name','class'))

 totalOtherStudies = totalOtherStudies %>% dplyr::select(c('external_gene_name')) %>% dplyr::mutate(curriculum="Z", semester='other studies')
 classifiedGenes = classifiedGenes %>% dplyr::mutate(curriculum = class, semester='this study') %>% dplyr::select(c('external_gene_name','curriculum','semester')) 
 
 total = rbind.data.frame(totalOtherStudies,classifiedGenes)
 library(ggalluvial)
 total = as.data.frame(total)
 total$curriculum <- as.factor(total$curriculum)
 total = total[!duplicated(total),]
 
 pdf("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/packages/plots/figure2//allStudies_comparisonToSLAMseq.pdf",height=3,width=3)
 
       p = ggplot(total,
              aes(x = semester, stratum = curriculum, alluvium = external_gene_name,
                  fill = curriculum, label = curriculum)) +
         scale_fill_brewer(type = "qual", palette = "Set1") +
         geom_flow(stat = "alluvium", lode.guidance = "rightleft",
                   color = "darkgray") +
         geom_stratum() +
         theme(legend.position = "bottom") +
          geom_label(stat = "stratum") + xlab('Study') + cowplot::theme_cowplot()
       p+  theme_ameres(type = 'barplot') +
         geom_stratum(alpha = 0) +
         geom_text(stat = "stratum")
       
 dev.off()
 
 
 `%!in%` = Negate(`%in%`)
 
onlySLAMseq = classifiedGenes[ classifiedGenes$external_gene_name %!in% totalOtherStudies$external_gene_name,]
onlySLAMseq_Z = onlySLAMseq %>% dplyr::filter(curriculum == "Z")
onlySLAMseq_MZ = onlySLAMseq %>% dplyr::filter(curriculum == "MZ")

conversionRates  = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/packages/data/countingWindows_classified_conversions.txt",sep="\t",stringsAsFactors = F, header = T)
conversionRates_MZonlySLAMseq = conversionRates[conversionRates$gene %in% onlySLAMseq_MZ$external_gene_name,]
conversionRates_MZonlySLAMseq = conversionRates_MZonlySLAMseq %>% dplyr::select(dplyr::contains("Inj"))
colnames(conversionRates_MZonlySLAMseq) =  c('0.75','2','2.5','3','3.5','4','4.5','5','5.5')
conversionRates_MZonlySLAMseq = reshape::melt(conversionRates_MZonlySLAMseq)
conversionRates_MZonlySLAMseq$variable = as.numeric(as.character(conversionRates_MZonlySLAMseq$variable))
pdf("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/packages/plots/figure2/TCconversions_MZonlySLAMseq.pdf",height = 3, width=6)

  p= ggplot(conversionRates_MZonlySLAMseq,aes(x=variable,y=value, group=variable)) + geom_boxplot(outlier.shape = NA) +ylim(c(0,0.15)) + cowplot::theme_cowplot() + ggtitle(paste("n=",nrow(conversionRates_MZonlySLAMseq)/9))
  p = p + theme_ameres(type = 'barplot') + ylab("TC conversion rate")+ scale_x_continuous(breaks=c(0.75,2,2.5,3,3.5,4,4.5,5,5.5),labels=c(0.75,2,2.5,3,3.5,4,4.5,5,5.5)) + theme_ameres(type = 'barplot')+ stat_boxplot(geom = "errorbar")
print(p)
dev.off()


conversionRates_ZonlySLAMseq = conversionRates[conversionRates$gene %in% onlySLAMseq_Z$external_gene_name,]
conversionRates_ZonlySLAMseq = conversionRates_ZonlySLAMseq %>% dplyr::select(dplyr::contains("Inj"))
colnames(conversionRates_ZonlySLAMseq) =  c('0.75','2','2.5','3','3.5','4','4.5','5','5.5')
conversionRates_ZonlySLAMseq = reshape::melt(conversionRates_ZonlySLAMseq)
conversionRates_ZonlySLAMseq$variable = as.numeric(as.character(conversionRates_ZonlySLAMseq$variable))
pdf("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/packages/plots/figure2/TCconversions_ZonlySLAMseq.pdf",height = 3, width=6)

  p= ggplot(conversionRates_ZonlySLAMseq,aes(x=variable,y=value, group=variable)) + geom_boxplot(outlier.shape = NA) +ylim(c(0,0.15)) + cowplot::theme_cowplot() + ggtitle(paste("n=",nrow(conversionRates_ZonlySLAMseq)/9))
  p = p + theme_ameres(type = 'barplot') + ylab("TC conversion rate")+ scale_x_continuous(breaks=c(0.75,2,2.5,3,3.5,4,4.5,5,5.5),labels=c(0.75,2,2.5,3,3.5,4,4.5,5,5.5)) + theme_ameres(type = 'barplot')+ stat_boxplot(geom = "errorbar")
  print(p)
dev.off()