##### Classification of gene expression dynamics based on gene expression and T>C conversion rates

library(cluster)
library(factoextra)
library(dplyr)
`%!in%` = Negate(`%in%`)
options(scipen=999)

##### installing ggraster from github
#install.packages('devtools')
#devtools::install_github('VPetukhov/ggrastr')
library(ggrastr)
library(ggplot2)
library(reshape)
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


#### classification of gene expression patterns by classification of Gene expression patterns and based on the presence of conversions... 
#### i want to read in the metadata from the annotation used... 

errorRates = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_may2019//errorRates_predicted_observed.txt",header = T)
RPMs_all = read.table("///Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_may2019/dataTables_expandedCounting/RPM_allCws.txt",sep="\t",stringsAsFactors = F,header=T)
conversionRates_all = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017//analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_may2019//dataTables_expandedCounting//perGenes_conversion.txt",sep="\t",stringsAsFactors = F,header=T)

conversionRates_all[is.na(conversionRates_all)] <- 0
# conversionRates_all = as.matrix(conversionRates_all)
# conversionRates_all[which(conversionRates_all == "Inf")]<-0
conversionRates_all = as.data.frame(conversionRates_all)
conversionRates_all_injection = conversionRates_all %>% dplyr::select(-starts_with("Inc"))
conversionRates_all_injection = conversionRates_all_injection %>% dplyr::select(starts_with("Inj"))
#### i need to check now if this is above or below the background... if below the background, i want to set it to 0


for(i in 1:ncol(conversionRates_all_injection)){
  conversionRates_all_injection[,i] = conversionRates_all_injection[,i] - errorRates$predicted[i] 
}
conversionRates_all_injection = as.matrix(conversionRates_all_injection)
conversionRates_all_injection[which(conversionRates_all_injection <0 )] <- 0
conversionRates_all_injection[grep("Inf",conversionRates_all_injection)]<-0
conversionRates_all_injection = as.data.frame(conversionRates_all_injection)

RPMs_all = RPMs_all %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Inj_R3_TP9), sum, na.rm = TRUE)
conversionRates_all_injection$name = conversionRates_all$name

##### i now want to split this by replicate and create the mean
conversionRates_all_injection_split  = splitReplicates(dataFrameToSplit = conversionRates_all_injection,condition = "Inj",metadata_add = conversionRates_all_injection$name)
RPMs_all_split = splitReplicates(dataFrameToSplit = RPMs_all,condition = "Inj",metadata_add = RPMs_all[,1])


library(textshape)
#### getting the mean and normalizing this.. 
RPMs_all_split_mean = column_to_rownames(RPMs_all_split$mean,10)
conversionRates_all_injection_mean = column_to_rownames(conversionRates_all_injection_split$mean,10)
RPMs_all_split_mean = RPMs_all_split_mean[apply(RPMs_all_split_mean,1,max) > 5,]
RPMs_all_split_mean$gene = rownames(RPMs_all_split_mean)
write.table(RPMs_all_split_mean,"/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_may2019/dataTables_expandedCounting/RPM_countingWindows_above5rpm.txt",sep="\t",quote = F,row.names = F)
RPMs_all_split_mean$gene = NULL
RPMs_all_split_mean_normalized = as.data.frame(t(apply(RPMs_all_split_mean,1,function(x) x/max(x))))
clusterGenes = cluster::clara(RPMs_all_split_mean_normalized,3)
RPMs_all_split_mean$cluster = clusterGenes$clustering
names_RPMs = rownames(RPMs_all_split_mean)
RPMs_all_split_mean = RPMs_all_split_mean %>% dplyr::mutate(pre = Inj_R2_TP1+Inj_R2_TP2, post = Inj_R2_TP8 + Inj_R2_TP9) %>% 
  dplyr::mutate(preVsPost = pre/post , log10PreVsPost = log10(pre/post)) 
rownames(RPMs_all_split_mean) = names_RPMs
ggpubr::ggviolin(RPMs_all_split_mean,y='log10PreVsPost',x='cluster', fill='black')


#### cluster 1 = MZ, cluster 2 = M, cluster 3 = Zygotic
RPMs_all_split_mean_normalized$cluster = clusterGenes$clustering
RPMs_all_split_mean_normalized_split  = split(RPMs_all_split_mean_normalized,RPMs_all_split_mean_normalized$cluster,T)


### i just want to plot these clusters now... 
RPMs_all_split_mean_normalized_split = lapply(RPMs_all_split_mean_normalized_split,function(x) x[,-10])
RPMs_all_split_mean_normalized_split_melt = lapply(RPMs_all_split_mean_normalized_split,function(x) melt(x))
for(i in 1:3){
  RPMs_all_split_mean_normalized_split_melt[[i]]$cluster = i
}
RPMs_all_split_mean_normalized_split_melt = do.call(rbind.data.frame,RPMs_all_split_mean_normalized_split_melt)
ggpubr::ggboxplot(RPMs_all_split_mean_normalized_split_melt,x='variable',y='value',facet.by = 'cluster' )

ggpubr::ggviolin(RPMs_all_split_mean_normalized_split_melt,x='variable',y='value',facet.by = 'cluster' ,add = 'boxplot',fill='red',alpha = 0.4)



###### now checking if these have any TC conversions

initialMZ = rownames(RPMs_all_split_mean_normalized_split$`1`)
initialM = rownames(RPMs_all_split_mean_normalized_split$`2`)
initialZ = rownames(RPMs_all_split_mean_normalized_split$`3`)

#### checking conversion rates in purely M transcripts. 
conversionRates_all_injection_split$mean = conversionRates_all_injection_split$mean %>% dplyr::mutate(name = metadata_add) %>% dplyr::select(-metadata_add)
conversionRates_initalM =  conversionRates_all_injection_split$mean[conversionRates_all_injection_split$mean$name  %in% initialM  ,]%>% column_to_rownames(loc = 10)
clustering_MconversionRates = clara(conversionRates_initalM,3)$clustering
conversionRates_initalM$clusters = clustering_MconversionRates

clusterAndPlot = function(dataf,nClus){
  dataf = dataf %>% dplyr::select(c(Inj_R2_TP1:Inj_R2_TP9))
  dataf_normalized = as.data.frame(t(apply(dataf,1,function(x) x/max(x))))
  dataf_normalized[is.na(dataf_normalized)]<-0
  dataf_normalized$cluster = clara(dataf_normalized,nClus)$clustering
  dataf$cluster = dataf_normalized$cluster
  dataf_split  = split(dataf,dataf$cluster,T)
  
  
  ### i just want to plot these clusters now... 
  dataf_split = lapply(dataf_split,function(x) x[,-10])
  dataf_split_melt = lapply(dataf_split,function(x) melt(x))
  for(i in 1:nClus){
    dataf_split_melt[[i]]$cluster = i
  }
  dataf_split_melt = do.call(rbind.data.frame,dataf_split_melt)
  p =  ggpubr::ggboxplot(dataf_split_melt,x='variable',y='value',facet.by = 'cluster' )
  returnList = vector("list",2)
  returnList[[1]] = dataf
  returnList[[2]] = p
  return(returnList)
}



################ only maternal genes ####################


clusterConversions_M = clusterAndPlot(dataf = conversionRates_initalM,nClus = 3)
### only cluster 1,2 seems to have no conversions...the other clusters have some conversions. - so cluster 1 is a purely maternal cluster

purelyMaternal = clusterConversions_M[[1]] %>% tibble::rownames_to_column('gene')   %>% dplyr::filter(cluster !=3) %>% 
  dplyr::mutate(class = "M", description = "decreasing gene expression - no TC conversions") %>%  tibble::column_to_rownames('gene')


#### the othere genes are MZ
MZ_fromM =  clusterConversions_M[[1]] %>% tibble::rownames_to_column('gene')  %>% dplyr::filter(cluster ==3) %>% 
  dplyr::mutate(class = "MZ", description = "decreasing gene expression -  TC conversions present") %>%  tibble::column_to_rownames('gene')


###### initial zygotic genes      

conversionRates_initalZ =  conversionRates_all_injection_split$mean[conversionRates_all_injection_split$mean$name  %in% initialZ  ,]%>% column_to_rownames(loc = 10)
### all these have increasing gene expression - but if they are expressed at the first 2 time points with >1 RPM -> MZ  

RPM_initialZ = RPMs_all_split_mean %>% tibble::rownames_to_column('gene') %>% dplyr::filter(cluster == 3) %>%  tibble::column_to_rownames('gene')
pureZ = RPM_initialZ %>% tibble::rownames_to_column('gene') %>% dplyr::filter(Inj_R2_TP1<1 | Inj_R2_TP2 <1) %>% 
  dplyr::mutate(class="Z",description = "increasing gene expression - absent in initial stages") %>%  tibble::column_to_rownames('gene')
pureZ_names = rownames(pureZ)

MZ_fromZ =  RPM_initialZ %>% tibble::rownames_to_column('gene') 
MZ_fromZ = MZ_fromZ[MZ_fromZ$gene %!in% pureZ_names,]
MZ_fromZ = MZ_fromZ %>%  dplyr::mutate(class = "MZ", description = "increasing gene expression - present from 1st two stages ") %>% tibble::column_to_rownames('gene') 
MZ_fromZ_names = rownames(MZ_fromZ)


pureZ = conversionRates_initalZ[rownames(conversionRates_initalZ) %in% pureZ_names,]  %>% tibble::rownames_to_column('gene') %>%
  dplyr::mutate(class="Z",description = "increasing gene expression - absent in initial stages")  %>% tibble::column_to_rownames('gene') 
MZ_fromZ = conversionRates_initalZ[rownames(conversionRates_initalZ) %in% MZ_fromZ_names,] %>%  dplyr::mutate(class="MZ",description ="increasing gene expression - present from 1st two stages ") 
rownames(MZ_fromZ) = MZ_fromZ_names
####### the initial MZ genes can either be true MZ or can be maternal long lived genes

conversionRates_initalMZ =  conversionRates_all_injection_split$mean[conversionRates_all_injection_split$mean$name  %in% initialMZ  ,] %>% column_to_rownames(loc = 10)
conversionRates_initalMZ_cluster = clusterAndPlot(dataf = conversionRates_initalMZ,nClus = 3)

###### looks like cluster 1,2 has no TC conversions - these are maternal transcripts that are just degrated at later time points. 

Mstable = conversionRates_initalMZ_cluster[[1]] %>% tibble::rownames_to_column('gene') %>% dplyr::filter(cluster !=3) %>% dplyr::mutate(class = "M-stable", description = "constant gene expression - noTC ") %>% tibble::column_to_rownames('gene') 
MZ = conversionRates_initalMZ_cluster[[1]] %>% tibble::rownames_to_column('gene') %>% dplyr::filter(cluster ==3) %>%  dplyr::mutate(class = "MZ", description = "constant gene expression - increase in TC ") %>% tibble::column_to_rownames('gene') 



##### total genes
purelyMaternal =  purelyMaternal %>% dplyr::select(c(Inj_R2_TP1,Inj_R2_TP2,Inj_R2_TP3,Inj_R2_TP4,Inj_R2_TP5,Inj_R2_TP6,Inj_R2_TP7,Inj_R2_TP8,Inj_R2_TP9,class,description)) %>% 
  tibble::rownames_to_column('gene') %>%
  dplyr::mutate(description = "M pure")
MZ_fromM = MZ_fromM %>% dplyr::select(c(Inj_R2_TP1,Inj_R2_TP2,Inj_R2_TP3,Inj_R2_TP4,Inj_R2_TP5,Inj_R2_TP6,Inj_R2_TP7,Inj_R2_TP8,Inj_R2_TP9,class,description)) %>% 
  tibble::rownames_to_column('gene') %>% dplyr::mutate(description = "MZ from M")
MZ_fromZ = MZ_fromZ  %>% dplyr::select(c(Inj_R2_TP1,Inj_R2_TP2,Inj_R2_TP3,Inj_R2_TP4,Inj_R2_TP5,Inj_R2_TP6,Inj_R2_TP7,Inj_R2_TP8,Inj_R2_TP9,class,description)) %>%
  tibble::rownames_to_column('gene') %>% dplyr::mutate(description = "MZ from Z")
pureZ = pureZ  %>% dplyr::select(c(Inj_R2_TP1,Inj_R2_TP2,Inj_R2_TP3,Inj_R2_TP4,Inj_R2_TP5,Inj_R2_TP6,Inj_R2_TP7,Inj_R2_TP8,Inj_R2_TP9,class,description))  %>% 
  tibble::rownames_to_column('gene') %>% dplyr::mutate(description = "Z pure")
Mstable = Mstable  %>% dplyr::select(c(Inj_R2_TP1,Inj_R2_TP2,Inj_R2_TP3,Inj_R2_TP4,Inj_R2_TP5,Inj_R2_TP6,Inj_R2_TP7,Inj_R2_TP8,Inj_R2_TP9,class,description)) %>% 
  tibble::rownames_to_column('gene') %>% dplyr::mutate(description = "M stable")
MZ = MZ  %>% dplyr::select(c(Inj_R2_TP1,Inj_R2_TP2,Inj_R2_TP3,Inj_R2_TP4,Inj_R2_TP5,Inj_R2_TP6,Inj_R2_TP7,Inj_R2_TP8,Inj_R2_TP9,class,description)) %>% 
  tibble::rownames_to_column('gene') %>% dplyr::mutate(description = "MZ pure")

totalSamples  = rbind(purelyMaternal, MZ_fromM,MZ_fromZ,pureZ,Mstable,MZ)

write.table(totalSamples,"/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_may2019//dataTables_expandedCounting/countingWindows_classified_conversions.txt",sep="\t",quote = F,row.names = F)

library(PKNCA)
####### plotting 

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/classificationOfTranscripts/classification_genes.pdf",
    height = 4,width=6)
    

    maternalGenes  = totalSamples %>% dplyr::filter(class == "M"| class=="M-stable") %>% dplyr::select(Inj_R2_TP1:Inj_R2_TP9)
    maternalGenes_mean = data.frame(meanSamples = colMeans(maternalGenes),time = c(0.75,2,2.5,3,3.5,4,4.5,5,5.5),type = 'Maternal' )
    
    maternalGenes = melt(maternalGenes)
    p = ggpubr::ggboxplot(maternalGenes,x='variable',y='value',title = paste("Maternal genes n = (",nrow(maternalGenes)/9,")"), ylim=c(0,0.25),
                          ylab = 'Fraction TC', xlab = 'timepoint (hpf)',fill='red',outlier.shape = NA)
    p + theme_ameres(type = 'barplot')
    
    
    ZygoticGenes  = totalSamples %>% dplyr::filter(class == "Z") %>% dplyr::select(Inj_R2_TP1:Inj_R2_TP9)
    ZygoticGenes_mean = data.frame(meanSamples = colMeans(ZygoticGenes),time = c(0.75,2,2.5,3,3.5,4,4.5,5,5.5),type = 'Zygotic' )
    
    ZygoticGenes = melt(ZygoticGenes)
    p = ggpubr::ggboxplot(ZygoticGenes,x='variable',y='value',title = paste("Zygotic genes n = (",nrow(ZygoticGenes)/9,")"), ylim=c(0,0.25),
                          ylab = 'Fraction TC', xlab = 'timepoint (hpf)',fill='red',outlier.shape = NA)
    p + theme_ameres(type = 'barplot')
    
    
    MZGenes  = totalSamples %>% dplyr::filter(class == "MZ") %>% dplyr::select(Inj_R2_TP1:Inj_R2_TP9)
    MZGenes_mean = data.frame(meanSamples = colMeans(MZGenes),time = c(0.75,2,2.5,3,3.5,4,4.5,5,5.5),type = 'MZ' )
    
    MZGenes = melt(MZGenes)
    p = ggpubr::ggboxplot(MZGenes,x='variable',y='value',title = paste("M-Z genes genes n = (",nrow(MZGenes)/9,")"),
                          ylim=c(0,0.25),ylab = 'Fraction TC', xlab = 'timepoint (hpf)',fill='red',outlier.shape = NA)
    p + theme_ameres(type = 'barplot')

    MstableGenes  = totalSamples %>% dplyr::filter(class == "M-stable") %>% dplyr::select(Inj_R2_TP1:Inj_R2_TP9)
    MstableGenes_mean = data.frame(meanSamples = colMeans(MstableGenes),time = c(0.75,2,2.5,3,3.5,4,4.5,5,5.5),type = 'M-stable' )
    
    MstableGenes = melt(MstableGenes)
    p = ggpubr::ggboxplot(MstableGenes,x='variable',y='value',title = paste("M-stable genes genes n = (",nrow(MstableGenes)/9,")"), ylim=c(0,0.25),
                          ylab = 'Fraction TC', xlab = 'timepoint (hpf)',fill='red',outlier.shape = NA)
    p + theme_ameres(type = 'barplot')

    maternalGenes = maternalGenes %>% dplyr::mutate(type = 'Maternal')
    ZygoticGenes = ZygoticGenes %>% dplyr::mutate(type = 'Zygotic')
    MstableGenes = MstableGenes %>% dplyr::mutate(type = 'Maternal stable')
    MZGenes = MZGenes    %>% dplyr::mutate(type = 'MZ')
    
total = rbind.data.frame(maternalGenes_mean,ZygoticGenes_mean,MstableGenes_mean,MZGenes_mean)
total_means = total %>% dplyr::mutate(log10 = log10(meanSamples)) 
ggpubr::ggline(total_means,x = 'time',y='log10',group='type',col='type',palette = 'Set1',ylab = 'log10(fraction TC)')




   ### now plotting the RPMs... 

    
    maternalGenes  = totalSamples %>% dplyr::filter( class == "M-stable" | class == "M") %>% dplyr::select('gene')
    RPMs_all_M = RPMs_all_split_mean_normalized[rownames(RPMs_all_split_mean_normalized) %in% maternalGenes$gene,] %>% dplyr::select(-'cluster')
    RPMs_all_M = melt(RPMs_all_M)
    p = ggpubr::ggboxplot(RPMs_all_M,x='variable',y='value',title = paste("Maternal genes n = (",nrow(RPMs_all_M)/9,")"), ylim=c(0,1),
                          ylab = 'Normalized RPM ', xlab = 'timepoint (hpf)',fill='blue',outlier.shape = NA)
    p + theme_ameres(type = 'barplot')
    
    
    zygoticGenes  = totalSamples %>% dplyr::filter(class == "Z") %>% dplyr::select('gene')
    RPMs_all_Z = RPMs_all_split_mean_normalized[rownames(RPMs_all_split_mean_normalized) %in% zygoticGenes$gene,] %>% dplyr::select(-'cluster')
    RPMs_all_Z = melt(RPMs_all_Z)
    p = ggpubr::ggboxplot(RPMs_all_Z,x='variable',y='value',title = paste("Zygotic genes n = (",nrow(RPMs_all_Z)/9,")"), ylim=c(0,1),
                          ylab = 'Normalized RPM', xlab = 'timepoint (hpf)',fill='blue',outlier.shape = NA)
    p + theme_ameres(type = 'barplot')


    
    maternalZygoticGenes =  totalSamples %>% dplyr::filter(class == "MZ") %>% dplyr::select('gene')
    RPMs_all_MZ = RPMs_all_split_mean_normalized[rownames(RPMs_all_split_mean_normalized) %in% maternalZygoticGenes$gene,] %>% dplyr::select(-'cluster')
    RPMs_all_MZ = melt(RPMs_all_MZ)
    p = ggpubr::ggboxplot(RPMs_all_MZ,x='variable',y='value',title = paste("MZ genes n = (",nrow(RPMs_all_MZ)/9,")"), ylim=c(0,1),
                          ylab = 'Normalized RPM', xlab = 'timepoint (hpf)',fill='blue',outlier.shape = NA)
    p + theme_ameres(type = 'barplot')
    

    
    maternalStableGenes =  totalSamples %>% dplyr::filter(class == "M-stable") %>% dplyr::select('gene')
    RPMs_all_Mstable = RPMs_all_split_mean_normalized[rownames(RPMs_all_split_mean_normalized) %in% maternalStableGenes$gene,] %>% dplyr::select(-'cluster')
    RPMs_all_Mstable = melt(RPMs_all_Mstable)
    p = ggpubr::ggboxplot(RPMs_all_Mstable,x='variable',y='value',title = paste("Mstable genes n = (",nrow(RPMs_all_Mstable)/9,")"), ylim=c(0,1),
                          ylab = 'Normalized RPM', xlab = 'timepoint (hpf)',fill='blue',outlier.shape = NA)
    p + theme_ameres(type = 'barplot')

dev.off()

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/classificationOfTranscripts/classification_genes_TClog.pdf",height = 4,width=4)
    
ggpubr::ggline(total_means,x = 'time',y='log10',group='type',col='type',palette = 'Set1',ylab = 'log10(fraction TC)') + theme_ameres(type = 'barplot')

dev.off()

