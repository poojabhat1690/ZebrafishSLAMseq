#### script used to plot the relationship between CAI and miR430 targets

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

miRNA430Targets = read.table("//Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure4/data//miRNA430Tarets_utrs.txt",
                             sep="\t",stringsAsFactors = F, header = T)


CWconnection = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure4/data//finalSplitAnnotations.bed",sep="\t",stringsAsFactors = F)
CWconnection=  CWconnection %>% dplyr::mutate(a_Gene_ID = V7) %>% dplyr::select(c('a_Gene_ID','V4'))
CWconnection_genes = data.frame(V4=unique(CWconnection$V4),starts=T)


miRNA430Targets= plyr::join(miRNA430Targets,CWconnection)
miRNA430Targets=  miRNA430Targets[!duplicated(miRNA430Targets),]

classifiedGenes = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/packages/data/countingWindows_classified_conversions.txt",
                             sep="\t",stringsAsFactors = F,header = T)
classifiedGenes_info = classifiedGenes %>% dplyr::select(c('gene','class','description'))

miRNA430Targets$gene = miRNA430Targets$V4
miRNA430Targets = plyr::join(classifiedGenes_info,miRNA430Targets,by='gene')

miRNA430Targets = miRNA430Targets %>% dplyr::select(c('gene','class','description','Site_type')) %>%   dplyr::group_by(gene) %>% dplyr::mutate(allSites_UTR = paste(Site_type,collapse = ","))
miRNA430Targets= miRNA430Targets[!duplicated(miRNA430Targets$gene),]



#### getting the ORF targets... 

orfTargets = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure4/data/orfTargets.txt",
                        sep="\t",stringsAsFactors = F, header = T)


orfTargets = orfTargets %>% dplyr::select(c('a_Gene_ID','miRNA_family_ID','Site_type')) %>% dplyr::group_by(a_Gene_ID) %>%
  dplyr::mutate(allSites_ORF=paste(Site_type,collapse = ","))
orfTargets = orfTargets[!duplicated(orfTargets$a_Gene_ID),]

orfTargets$ensembl_transcript_id = orfTargets$a_Gene_ID
dr11_connection = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure4/data//allGenes_dr11.txt",
                             sep="\t",stringsAsFactors = F, header = T)

dr11_connection= dr11_connection %>% dplyr::select(c('ensembl_transcript_id','external_gene_name','ensembl_gene_id'))
orfTargets = as.data.frame(orfTargets)
orfTargets  = plyr::join(orfTargets,dr11_connection)
orfTargets = orfTargets[!duplicated(orfTargets$external_gene_name),]
orfTargets$gene = orfTargets$external_gene_name


classifiedGenes_info=  plyr::join(classifiedGenes_info,orfTargets)
classifiedGenes_info= classifiedGenes_info[!is.na(classifiedGenes_info$miRNA_family_ID),]

classifiedGenes_info_orf=  classifiedGenes_info %>% dplyr::select(c('gene','allSites_ORF'))
classifiedGenes_info_targets = plyr::join(as.data.frame(miRNA430Targets),as.data.frame(classifiedGenes_info_orf))
classifiedGenes_info_targets = as.data.frame(classifiedGenes_info_targets)
mir430Containing = classifiedGenes_info_targets[!is.na(classifiedGenes_info_targets$Site_type),]

###### adding to this the informatio

######## the direct targets will be upregulated in the asence of mir430... for this taking data from mishima tomari...


#### find the genes upregulated in in the presence of a miR430 oligo


countMat = read.table("//Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure4/data//counts_MishimaTomari_6hpf.txt",
                      sep="\t",stringsAsFactors = F, header = T)
counts = countMat[,c(7:12)]
colnames(counts) = c('unt1','unt2','unt3','mir430_1','mir430_2','mir430_3')
row.names(counts)= countMat$Geneid
coldata = data.frame(condition=c('untreated','untreated','untreated','treated','treated','treated'),
                     type=rep('single-read',6))
rownames(coldata) = c('unt1','unt2','unt3','mir430_1','mir430_2','mir430_3')
library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
plotMA
res$ensembl_gene_id = rownames(res)
res =  plyr::join(as.data.frame(res),dr11_connection)
res = res[!duplicated(res$external_gene_name),]
res$gene = res$external_gene_name
res$external_gene_name = NULL
classifiedGenes_info_targets = plyr::join(classifiedGenes_info_targets,res,by='gene')

classifiedGenes_info_targets = classifiedGenes_info_targets %>% dplyr::mutate(UTRsites = ifelse(allSites_UTR == "NA",F,T))  %>%
  dplyr::mutate(ORFsites = ifelse(allSites_ORF == "NA",F,T)) %>% dplyr::mutate(downReg = ifelse(padj<0.01 & log2FoldChange < -0.58, T, F ))

classifiedGenes_info_targets = classifiedGenes_info_targets %>%
  dplyr::select(c('gene','class','description','UTRsites','ORFsites','downReg'))


MZgenes = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure4/data//clusteredBasedOnTCaccumulation.bed",sep="\t",stringsAsFactors = F)
MZgenes$gene = MZgenes$V10
MZgenes = MZgenes %>% dplyr::select(c(gene,V11))
masterTable_data = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure5/data//ClassifiedGenes_masterTable.txt",stringsAsFactors = F, header = T,sep = "\t")
colnames(MZgenes) = c('gene','V3')

Mgenes = masterTable_data %>% dplyr::filter(class == 'M') %>% dplyr::mutate(MZrate = 'M') %>%
  dplyr::select('external_gene_name','MZrate') %>% dplyr::mutate(gene = external_gene_name,V3 = MZrate) %>% dplyr::select(gene,V3)
Zgenes = masterTable_data %>% dplyr::filter(class == 'Z') %>% dplyr::mutate(MZrate = 'Z')%>%
  dplyr::select('external_gene_name','MZrate') %>% dplyr::mutate(gene = external_gene_name,V3 = MZrate)%>% dplyr::select(gene,V3)
Mstablegenes = masterTable_data %>% dplyr::filter(class == 'M-stable') %>% dplyr::mutate(MZrate = 'Mstable') %>% 
  dplyr::select('external_gene_name','MZrate') %>% dplyr::mutate(gene = external_gene_name,V3 = MZrate) %>% dplyr::select(gene,V3)

MZgenes = rbind.data.frame(MZgenes,Mgenes,Zgenes,Mstablegenes)

allTargets = plyr::join(classifiedGenes_info_targets,MZgenes)
allTargets$V3[is.na(allTargets$V3)]<-0

CWconnection_genes$gene = as.character(CWconnection_genes$V4)
CWconnection_genes = CWconnection_genes %>% dplyr::select(c(gene,starts))
allTargets =  plyr::join(allTargets,CWconnection_genes,by='gene')
allTargets = allTargets[!is.na(allTargets$starts),]
masterTable_data$gene = masterTable_data$external_gene_name
allTargets = allTargets %>% dplyr::select(-c(class ,description ))

allTargets = plyr::join(allTargets,masterTable_data,by='gene')


################################## 

masterTable_data = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure5/data//ClassifiedGenes_masterTable.txt",stringsAsFactors = F, header = T,sep = "\t")

classifiedGenes = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/packages/data//countingWindows_classified_conversions.txt",
                             sep="\t",stringsAsFactors = F,header = T)
classifiedGenes$V4 = classifiedGenes$gene
CWconnection = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure4/data/finalEnds_intronsRemoved.bed",sep="\t",stringsAsFactors = F)
CWconnection =  plyr::join(CWconnection,classifiedGenes,by='V4')
CWconnection = CWconnection%>% dplyr::mutate(length = V3-V2)
CWconnection = CWconnection %>% dplyr::group_by(V4) %>% dplyr::mutate(meanLength = mean(V3-V2))
CWconnection = CWconnection[!duplicated(CWconnection$V4),] %>% dplyr::mutate(external_gene_name = V4) %>% dplyr::ungroup() %>%
  dplyr::select(c(external_gene_name, meanLength))
masterTable_data = masterTable_data %>% dplyr::select(c(external_gene_name,meanCAI,meanLength))
masterTable_data = masterTable_data %>% dplyr::mutate(gene=external_gene_name) %>% dplyr::select(c(gene,meanLength))

allTargets  = plyr::join(allTargets,masterTable_data)
#allTargets = allTargets[which(allTargets$meanLength<5000),]
cutoffCAI = median(allTargets$meanCAI,na.rm=T)
cutoffLenth = median(allTargets$meanLength,na.rm = T)
`%!in%` = Negate(`%in%`)
write.table(allTargets,'/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/MZdynamics/miRNAanalysis/mir430_definedTargets.txt',sep = '\t',quote = F)

classes_genes = unique(allTargets$V3)
q1_enriched = c()
q1_depleted = c()
###### q1

for(i in 1:length(classes_genes)){
  
  
  targetSet = allTargets %>% dplyr::filter(V3 == classes_genes[i] & meanCAI<cutoffCAI & meanLength>cutoffLenth)
  nonTargetSet = allTargets[allTargets$gene  %!in% targetSet$gene,]
  
  mirs_inTarget = nrow(targetSet %>% dplyr::filter(UTRsites==T & downReg == T))
  nonMirs_inTarget =  nrow(targetSet) - mirs_inTarget
  
  mirs_inNONTarget = nrow(nonTargetSet %>% dplyr::filter(UTRsites==T & downReg == T))
  nonMirs_inNOnTarget =  nrow(nonTargetSet) - mirs_inNONTarget
  
  
  q1_depleted = c(q1_depleted,fisher.test(matrix(c(mirs_inTarget,mirs_inNONTarget,nonMirs_inTarget,nonMirs_inNOnTarget),nrow = 2, byrow = T), 
                                          alternative = 'less')$p.value)
  q1_enriched = c(q1_enriched,fisher.test(matrix(c(mirs_inTarget,mirs_inNONTarget,nonMirs_inTarget,nonMirs_inNOnTarget),nrow = 2, byrow = T), 
                                          alternative = 'greater')$p.value)
  
}


############################## q2 slow genes
q2_enriched = c()
q2_depleted = c()
for(i in 1:length(classes_genes)){
  targetSet = allTargets %>% dplyr::filter(V3 == classes_genes[i] & meanCAI>cutoffCAI & meanLength>cutoffLenth)
  nonTargetSet = allTargets[allTargets$gene  %!in% targetSet$gene,]
  
  mirs_inTarget = nrow(targetSet %>% dplyr::filter(UTRsites==T & downReg == T))
  nonMirs_inTarget =  nrow(targetSet) - mirs_inTarget
  
  mirs_inNONTarget = nrow(nonTargetSet %>% dplyr::filter(UTRsites==T & downReg == T))
  nonMirs_inNOnTarget =  nrow(nonTargetSet) - mirs_inNONTarget
  
  
  q2_depleted = c(q2_depleted,fisher.test(matrix(c(mirs_inTarget,mirs_inNONTarget,nonMirs_inTarget,nonMirs_inNOnTarget),
                                                 nrow = 2, byrow = T), alternative = 'less')$p.value)
  q2_enriched = c(q2_enriched,fisher.test(matrix(c(mirs_inTarget,mirs_inNONTarget,nonMirs_inTarget,nonMirs_inNOnTarget),
                                                 nrow = 2, byrow = T), alternative = 'greater')$p.value)
  
}
############################## q3 slow genes

q3_enriched = c()
q3_depleted = c()
for(i in 1:length(classes_genes)){
  targetSet = allTargets %>% dplyr::filter(V3 == classes_genes[i] & meanCAI<cutoffCAI & meanLength<cutoffLenth)
  nonTargetSet = allTargets[allTargets$gene  %!in% targetSet$gene,]
  
  mirs_inTarget = nrow(targetSet %>% dplyr::filter(UTRsites==T & downReg == T))
  nonMirs_inTarget =  nrow(targetSet) - mirs_inTarget
  
  mirs_inNONTarget = nrow(nonTargetSet %>% dplyr::filter(UTRsites==T & downReg == T))
  nonMirs_inNOnTarget =  nrow(nonTargetSet) - mirs_inNONTarget
  
  
  q3_depleted = c(q3_depleted,fisher.test(matrix(c(mirs_inTarget,mirs_inNONTarget,nonMirs_inTarget,nonMirs_inNOnTarget),nrow = 2, byrow = T), alternative = 'less')$p.value)
  q3_enriched = c(q3_enriched,fisher.test(matrix(c(mirs_inTarget,mirs_inNONTarget,nonMirs_inTarget,nonMirs_inNOnTarget),nrow = 2, byrow = T), alternative = 'greater')$p.value)
  
}
##### q4 slow genes

q4_enriched = c()
q4_depleted = c()
for(i in 1:length(classes_genes)){
  
  targetSet = allTargets %>% dplyr::filter(V3 == classes_genes[i] & meanCAI>cutoffCAI & meanLength<cutoffLenth)
  nonTargetSet = allTargets[allTargets$gene  %!in% targetSet$gene,]
  
  mirs_inTarget = nrow(targetSet %>% dplyr::filter(UTRsites==T & downReg == T))
  nonMirs_inTarget =  nrow(targetSet) - mirs_inTarget
  
  mirs_inNONTarget = nrow(nonTargetSet %>% dplyr::filter(UTRsites==T & downReg == T))
  nonMirs_inNOnTarget =  nrow(nonTargetSet) - mirs_inNONTarget
  
  
  q4_depleted = c(q4_depleted,fisher.test(matrix(c(mirs_inTarget,mirs_inNONTarget,nonMirs_inTarget,nonMirs_inNOnTarget),nrow = 2, byrow = T), alternative = 'less')$p.value)
  q4_enriched = c(q4_enriched,fisher.test(matrix(c(mirs_inTarget,mirs_inNONTarget,nonMirs_inTarget,nonMirs_inNOnTarget),nrow = 2, byrow = T), alternative = 'greater')$p.value)
  
}
#######################################
allEnriched = cbind.data.frame(q1_enriched,q2_enriched,q3_enriched,q4_enriched)
rownames(allEnriched) = classes_genes

allDepleted = cbind.data.frame(q1_depleted,q2_depleted,q3_depleted,q4_depleted)
rownames(allDepleted) = classes_genes
