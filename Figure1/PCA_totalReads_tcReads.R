
library(DESeq2)

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
allReads = read.table('/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_may2019/dataTables_expandedCounting/totalCounts_allCws.txt',
                      sep="\t",stringsAsFactors = F, header = T)
allReads =  allReads %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Untreated_TP9), sum, na.rm = TRUE)
allReads_Inj =  allReads %>% dplyr::select(-dplyr::contains('Inc'))%>% dplyr::select(-name)

sampleInfo = data.frame(sample = colnames(allReads_Inj))

sampleInfo = data.frame(                        condition = unlist(lapply(strsplit(as.character(sampleInfo$sample),"_",T),function(x) x[1])))
sampleInfo$replicate = c(rep('R2',each=9),rep('R3',9),rep('Untreated',3))
sampleInfo$timepoint = c(paste0('TP',c(1:9)),paste0('TP',c(1:9)),'TP1','TP6','TP9')
dds = DESeqDataSetFromMatrix(countData = as.matrix(allReads_Inj),colData =sampleInfo,design = condition~timepoint )
vsd_minusInc = vst(dds, blind=FALSE)
p_rnaseq = DESeq2::plotPCA(vsd_minusInc, intgroup=c("replicate"),returnData=T)
p_rnaseq$time = c(paste0('TP',c(1:9)),paste0('TP',c(1:9)),'TP1','TP6','TP9')
percentVar_rnaseq <- round(100 * attr(p_rnaseq, "percentVar")) 
pca_rnaseq = ggpubr::ggscatter(p_rnaseq,x='PC1',y='PC2',color = 'replicate',palette = 'Dark2' ,
                               label = 'time',repel = T,size=2.5,ylab=paste0('PC2'," ",percentVar_rnaseq[2]," %"),
                               xlab = paste0('PC1'," ",percentVar_rnaseq[1]," %"),title = paste0('RNAseq data (n= ',
                                                                                                 nrow(allReads),")") )+ 
  theme_ameres(type = 'barplot') 

pdf('/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/plots/PCAplots.pdf', height=4, width=5)

print(pca_rnaseq)


##### also now getting TC reads only...

#### Deseq2 dosnt allow decimal values.. so just taking the total vals without error subtraction

tcReads = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_may2019/dataTables_expandedCounting/numberOfreadsWithTC.txt",
                     sep="\t", stringsAsFactors = F, header = T)
tcReads =  tcReads %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Untreated_TP9), sum, na.rm = TRUE)
tcReads = tcReads %>% dplyr::select(-dplyr::contains('Inc'))
tcReads = tcReads %>% dplyr::select(Inj_R2_TP1:Untreated_TP9)


dds = DESeqDataSetFromMatrix(countData = tcReads,colData =sampleInfo,design = condition~timepoint )
vsd_minusInc = vst(dds, blind=FALSE)
p_rnaseq = DESeq2::plotPCA(vsd_minusInc, intgroup=c("replicate"),returnData=T)
p_rnaseq$time = c(paste0('TP',c(1:9)),paste0('TP',c(1:9)),'TP1','TP6','TP9')
percentVar_rnaseq <- round(100 * attr(p_rnaseq, "percentVar")) 
pca_rnaseq = ggpubr::ggscatter(p_rnaseq,x='PC1',y='PC2',color = 'replicate',palette = 'Dark2' ,
                               label = 'time',repel = T,size=2.5,ylab=paste0('PC2'," ",percentVar_rnaseq[2]," %"),
                               xlab = paste0('PC1'," ",percentVar_rnaseq[1]," %"),title = paste0('RNAseq data (n= ',
                                                                                                 nrow(tcReads),")") )+ 
  theme_ameres(type = 'barplot') 
print(pca_rnaseq)

dev.off()