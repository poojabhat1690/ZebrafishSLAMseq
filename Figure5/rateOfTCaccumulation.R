###### script clustering gene expression from MZ contstant genes, generating heatmaps and checking for candidates of cytoplasmic polyadenylation.

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
library(dplyr)
`%!in%` = Negate(`%in%`)
library(ggplot2)


totalReads = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures//figure2/data//totalCounts_allCws.txt",
                        sep="\t", header = T, stringsAsFactors = F)
tcReads = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure2/data//numberOfreadsWithTC.txt",
                     sep="\t", stringsAsFactors = F, header = T)



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

totalReads_Inj = totalReads %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Inj_R3_TP9), sum, na.rm = TRUE)
tcReads_Inj = tcReads %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Inj_R3_TP9), sum, na.rm = TRUE)

tcReads_Inj[,2:19]= tcReads_Inj[,2:19]/0.82
totalReads_Inj_mean = splitReplicates(dataFrameToSplit = totalReads_Inj, condition = "Inj",metadata_add = totalReads_Inj$name)[[3]]
tcReads_Inj_mean = splitReplicates(dataFrameToSplit = tcReads_Inj, condition = "Inj",metadata_add = tcReads_Inj$name)[[3]]



classifiedGenes = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/packages/data//countingWindows_classified_conversions.txt",
                             sep="\t",stringsAsFactors = F,header = T)
classifiedGenes_MZ = classifiedGenes %>% dplyr::filter(description=="MZ pure")


totalReads_Inj_mean = totalReads_Inj_mean[totalReads_Inj_mean$metadata_add %in% classifiedGenes_MZ$gene,]
tcReads_Inj_mean = tcReads_Inj_mean[tcReads_Inj_mean$metadata_add %in% classifiedGenes_MZ$gene,]

fractionDisplaced = tcReads_Inj_mean[,c(1:9)]/totalReads_Inj_mean[,c(1:9)]
fractionDisplaced$gene = tcReads_Inj_mean$metadata_add


a = cluster::clara(fractionDisplaced[,c(1:9)],k = 4)
fractionDisplaced$cluster = a$clustering
fractionDisplaced = fractionDisplaced %>% dplyr::mutate(category = ifelse(cluster !=1, 'fast','slow'))



##### now also plotting the RPM of genes in these clusters

RPM_reads = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure2/data//RPM_allCws.txt",
                       sep="\t",stringsAsFactors = F, header = T)
RPM_reads =  RPM_reads %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Inj_R3_TP9), sum, na.rm = TRUE)
RPM_mean = splitReplicates(dataFrameToSplit = RPM_reads, condition = "Inj",metadata_add = RPM_reads$name)[[3]]
RPM_mean$gene = RPM_mean$metadata_add
fractionDisplaced_tmp = fractionDisplaced
colnames(fractionDisplaced) = c(paste0('fraction_',colnames(fractionDisplaced)[1:9]),'gene','cluster','category')
fractionDisplaced = plyr::join(fractionDisplaced,RPM_mean,by='gene')
fractionDisplaced$log10RPM = log10(fractionDisplaced$Inj_R2_TP1)
fractionDisplaced = fractionDisplaced %>% dplyr::group_by(category) %>% dplyr::mutate(n=n()) %>% dplyr::mutate(category_sample = paste0(category,"(n=",n,")"))
fractionDisplaced = fractionDisplaced[-grep('mt-',fractionDisplaced$metadata_add),]

write.table(fractionDisplaced, "/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure4/data//clusteredBasedOnTCaccumulation_fractionTC.bed",sep="\t",quote = F,row.names = F,col.names = F)

#write.table(fractionDisplaced,"//groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/MZdynamics/differentClassesOfTranscripts/clusteredBasedOnTCaccumulation.bed",sep="\t",quote = F,row.names = F,col.names = F)


b = t(apply(fractionDisplaced[,c("Inj_R2_TP1", "Inj_R2_TP2", "Inj_R2_TP3", "Inj_R2_TP4", "Inj_R2_TP5", "Inj_R2_TP6", "Inj_R2_TP7", "Inj_R2_TP8", "Inj_R2_TP9" )],1,function(x) x/max(x)))

cluster1_fraction = fractionDisplaced %>% dplyr::filter(cluster==1)
cluster2_fraction = fractionDisplaced %>% dplyr::filter(cluster==2) %>% dplyr::mutate(cluster=3)
cluster3_fraction = fractionDisplaced %>% dplyr::filter(cluster==3) %>% dplyr::mutate(cluster=2)
cluster4_fraction = fractionDisplaced %>% dplyr::filter(cluster==4) 
fractionDisplaced = rbind.data.frame(cluster1_fraction,cluster2_fraction,cluster3_fraction,cluster4_fraction)

write.table(fractionDisplaced,"/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure4/data//clusteredBasedOnTCaccumulation.txt",sep="\t",quote = F,row.names = F,col.names = F)

fractionDisplaced_split = split(fractionDisplaced,fractionDisplaced$cluster,T)
fractionDisplaced_split_copy = fractionDisplaced_split
numberPerCluster = reshape::melt(lapply(fractionDisplaced_split,nrow))
numberPerCluster$numberGenes = numberPerCluster$value
numberPerCluster = numberPerCluster %>% dplyr::select('numberGenes','L1')

pdf("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure4/plots//clusters_accumulation.pdf",
    height = 4, width = 5)
        fractionDisplaced_allGenes = fractionDisplaced
        fractionDisplaced_allGenes$cluster = 0
        fractionDisplaced = rbind.data.frame(fractionDisplaced,fractionDisplaced_allGenes)
        ggpubr::ggviolin(fractionDisplaced,x='cluster',y='log10RPM', add='boxplot', ylab = "log10 (RPM 0.75hpf)",fill = 'cluster', palette = "Set1") +
          theme_ameres(type = 'barplot') + ggpubr::stat_compare_means(ref.group = '0',method.args = list(alternative ='less'))
        
        ggpubr::ggbarplot(numberPerCluster,x='L1',y='value',label = 'value',fill='black',xlab='cluster',ylab='Number of genes')
        fractionDisplaced_split = reshape::melt(lapply(fractionDisplaced_split,function(x) colMeans(x[,c(1:9)])))
        fractionDisplaced_split$time = c(0.75,2,2.5,3,3.5,4,4.5,5,5.5)
        
        plyr::join(fractionDisplaced_split,numberPerCluster)
        ggpubr::ggline(fractionDisplaced_split,x='time',y='value',color = 'L1',ylab="fraction of TC reads", palette = 'Set1',numeric.x.axis = F) + 
          theme_ameres(type = 'barplot')   + ylim(c(0,1))
        

dev.off()



fractionDisplaced_split_copy = lapply(fractionDisplaced_split_copy,function(x) x[order(x$Inj_R2_TP9,decreasing = T),])
fractionDisplaced_split_copy_highest = lapply(fractionDisplaced_split_copy,function(x) x)
fractionDisplaced_split_copy$`4`$fraction_Inj_R2_TP9[which(fractionDisplaced_split_copy$`4`[,c(9)]>1)] <-1

pdf("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure4//plots//clusters_accumulation_topgenes_RAWVALUES.pdf", height = 4, width = 3)
    
    NMF::aheatmap((fractionDisplaced_split_copy$`1`[,c(5:9)]),Rowv = NA,Colv = NA,annRow = F,annCol = F,breaks = seq(0,1,0.025),color= RColorBrewer::brewer.pal(n = 9,name = 'Reds'))
    NMF::aheatmap((fractionDisplaced_split_copy$`2`[,c(5:9)]),Rowv = NA,Colv = NA,annRow = F,annCol = F,breaks = seq(0,1,0.025),color= RColorBrewer::brewer.pal(n = 9,name = 'Reds'))
    NMF::aheatmap((fractionDisplaced_split_copy$`3`[,c(5:9)]),Rowv = NA,Colv = NA,annRow = F,annCol = F,breaks = seq(0,1,0.025),color= RColorBrewer::brewer.pal(n = 9,name = 'Reds'))
    NMF::aheatmap((fractionDisplaced_split_copy$`4`[,c(5:9)]),Rowv = NA,Colv = NA,annRow = F,annCol = F,breaks = seq(0,1,0.025),color= RColorBrewer::brewer.pal(n = 9,name = 'Reds'))
    
    #NMF::aheatmap(log10(fractionDisplaced_split_copy$`3`[,c(5:9)]+1),Rowv = NA,Colv = NA,annRow = F,annCol = F,breaks = seq(0,0.4,0.01),color= RColorBrewer::brewer.pal(n = 9,name = 'Reds'))
    
    # NMF::aheatmap(log10(fractionDisplaced_split_copy$`1`[,c(13:21)]),Rowv = NA,Colv = NA,annRow = F,annCol = F,breaks = seq(0,1,0.01))
    # NMF::aheatmap(log10(fractionDisplaced_split_copy$`2`[,c(13:21)]),Rowv = NA,Colv = NA,annRow = F,annCol = F)
    # NMF::aheatmap(log10(fractionDisplaced_split_copy$`3`[,c(13:21)]),Rowv = NA,Colv = NA,annRow = F,annCol = F)
    # NMF::aheatmap(log10(fractionDisplaced_split_copy$`4`[,c(13:21)]),Rowv = NA,Colv = NA,annRow = F,annCol = F)


dev.off()


#### i also wan to check if these genes undergo CPA

CPA = read.table('/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure5/data/winata_CPA.txt',
                 sep="\t",stringsAsFactors = F, header = T)
nonCPA =read.table('/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure5/data/winata_nonCPA.txt',sep="\t",
                   stringsAsFactors = F, header = T)
allGenes_winataEtall = rbind.data.frame(CPA,nonCPA)


fractionDisplaced = fractionDisplaced %>% dplyr::ungroup()
cluster4 = fractionDisplaced[fractionDisplaced$cluster==4,]%>% dplyr::select(metadata_add,cluster)
cluster1 = fractionDisplaced[fractionDisplaced$cluster==1,]%>% dplyr::select(metadata_add,cluster)
cluster2 = fractionDisplaced[fractionDisplaced$cluster==2,]%>% dplyr::select(metadata_add,cluster)
cluster3 = fractionDisplaced[fractionDisplaced$cluster==3,]%>% dplyr::select(metadata_add,cluster)
Mgenes = classifiedGenes %>% dplyr::filter(class=='M') %>% dplyr::mutate(cluster=5,metadata_add=gene)%>% dplyr::select(metadata_add,cluster)
Zgenes = classifiedGenes %>% dplyr::filter(class=='Z') %>% dplyr::mutate(cluster=6,metadata_add=gene)%>% dplyr::select(metadata_add,cluster)
MstableGenes = classifiedGenes %>% dplyr::filter(class=='M-stable') %>% dplyr::mutate(cluster=7,metadata_add=gene)%>% dplyr::select(metadata_add,cluster)
MZother =  classifiedGenes %>% dplyr::filter(description=='MZ from M' | description == 'MZ from Z') %>% dplyr::mutate(cluster=8,metadata_add=gene)%>% dplyr::select(metadata_add,cluster)


allGenes = rbind.data.frame(cluster1,cluster2,cluster3,cluster4,Mgenes,Zgenes,MstableGenes,MZother)
allGenes= allGenes[allGenes$metadata_add %in% allGenes_winataEtall$gene.name,] 
allGenes = allGenes %>% dplyr::mutate(gene.name=metadata_add)
allGenes_winataEtall = plyr::join(allGenes_winataEtall,allGenes)
allGenes_winataEtall = allGenes_winataEtall[complete.cases(allGenes_winataEtall),]


pval = c()
oddsRatio = c()
inClusterCPA_nums = c()
allInluster = c()
for(i in 1:8){
  inCluster_CPA = nrow(allGenes_winataEtall %>% dplyr::filter(cluster==i & class == 'CPA'))
  NOTinCluster_CPA = nrow(allGenes_winataEtall %>% dplyr::filter(cluster!=i & class == 'CPA'))
  inCluster_NOTCPA = nrow(allGenes_winataEtall %>% dplyr::filter(cluster==i & class != 'CPA'))
  NOTinCluster_NOTCPA = nrow(allGenes_winataEtall %>% dplyr::filter(cluster!=i & class != 'CPA'))
  allInluster = c(allInluster,nrow(allGenes_winataEtall %>% dplyr::filter(cluster==i)))
  
  pval = c(pval,fisher.test(matrix(c(inCluster_CPA,inCluster_NOTCPA,NOTinCluster_CPA,NOTinCluster_NOTCPA),byrow = T,nrow = 2))$p.value)
  oddsRatio = c(oddsRatio,fisher.test(matrix(c(inCluster_CPA,inCluster_NOTCPA,NOTinCluster_CPA,NOTinCluster_NOTCPA),byrow = T,nrow = 2))$estimate)
  inClusterCPA_nums = c(inClusterCPA_nums, inCluster_CPA)
  
}


rations_pvals = data.frame(pval = unlist(pval),oddsRatio=unlist(oddsRatio),numberInCluster=unlist(inClusterCPA_nums),cluster=c(1,2,3,4,'M','Z','M-stable','otherMZ'))
rations_pvals$pval = p.adjust(rations_pvals$pval)
rations_pvals$oddsRatio= log10(rations_pvals$oddsRatio) 
rations_pvals = rations_pvals %>% dplyr::mutate(significance=ifelse(pval<0.01,'sig','nonSig'))

rations_pvals = rations_pvals %>% dplyr::filter(cluster != "Z") %>% dplyr::filter(cluster != "otherMZ")
pdf('/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure4/plots//fishers_cpa.pdf', height=3, width=5)

  p = ggpubr::ggbarplot(rations_pvals,x='cluster',y='oddsRatio',fill='significance',ylab='Odds Ratio (log10)',palette = 'Set1',label = rations_pvals$numberInCluster) +
 geom_hline(yintercept = 0) + theme_ameres(type = 'barplot')
  print(p)
dev.off()

### maternal genes for whichthere isa large increase poly(A) 

maternalGenes_highIncreasePolyA  = b %>% dplyr::filter(class=="M" & fractionU>2) 
allCPA = allGenes_winataEtall %>% dplyr::filter(class=='CPA')
allNONCPA = allGenes_winataEtall %>% dplyr::filter(class!='CPA')

cpa_nonCPA = c(length(intersect(maternalGenes_highIncreasePolyA$external_gene_name,allCPA$gene.name)),
length(intersect(maternalGenes_highIncreasePolyA$external_gene_name,allNONCPA$gene.name)))
cpa_nonCPA = data.frame(cpa_nonCPA,type=c('cpa','nonCPA'))
ggpubr::ggbarplot(cpa_nonCPA,x='type',y='cpa_nonCPA')
sum(maternalGenes_highIncreasePolyA$external_gene_name %in% allGenes_winataEtall$gene.name)
