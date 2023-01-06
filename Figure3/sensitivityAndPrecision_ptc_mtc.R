######## identification of transcripts that are expressed before the main wave of ZGA.. 
#### the idea is to use 3 parameters to find these transcripts... 
#### should be greater than the background conversion rate
#### shotuld have percentage TC conversion rate greater than a set of purely zygotic genes. 
#### should have multiple TC conversions..

library(ggalluvial)
library(dplyr)
library(ggbeeswarm)

theme_ameres <- function (type) { 
  
  types_plot = c("boxplot","barplot","violin")
  if(type %in% types_plot ==T)
  {
    
    if(type == "boxplot"){
      theme(legend.title=element_blank(),axis.text.x = element_text(margin=ggplot2::margin(10,15,10,15,"pt"),size = 15),axis.text.y = element_text(margin=ggplot2::margin(5,15,10,5,"pt"),size = 15), axis.ticks.length = unit(-0.25 , "cm"),legend.position="bottom",axis.line = element_line(colour = "black", lineend = "round"), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18))   
    }
    
    if(type == "barplot"){
      theme(legend.title=element_blank(),axis.text.x = element_text(margin=ggplot2::margin(10,15,10,15,"pt"),size = 15, hjust = 1),axis.text.y = element_text(margin=ggplot2::margin(5,15,10,5,"pt"),size = 15), axis.ticks.length = unit(-0.25 , "cm"),legend.position="bottom",axis.line = element_line(colour = "black", lineend = "round"), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18))   
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


#### i need to get the background conversions per replicate and subtract from the conversion rate. 
errorRates = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline//paper_ZebrafishSLAMseq/other/rawfigures/figure2/data//errorRates_predicted_observed.txt",header = T)
conversionRates_all = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure2/data/perGenes_conversion.txt",sep="\t",stringsAsFactors = F,header=T)

conversionRates_all[is.na(conversionRates_all)] <- 0
conversionRates_all = as.data.frame(conversionRates_all)
conversionRates_all_injection = conversionRates_all %>% dplyr::select( c(dplyr::contains("Inj"),'Untreated_TP1','Untreated_TP6','Untreated_TP9'))

for(i in 1:ncol(conversionRates_all_injection)){
  conversionRates_all_injection[,i] = conversionRates_all_injection[,i] - errorRates$predicted[i] 
}

conversionRates_all_injection = as.matrix(conversionRates_all_injection)
conversionRates_all_injection[which(conversionRates_all_injection <0 )] <- 0
conversionRates_all_injection[grep("Inf",conversionRates_all_injection)]<-0
conversionRates_all_injection = as.data.frame(conversionRates_all_injection)
conversionRates_all_injection$name = conversionRates_all$name

##### splitting into replicates                    
conversionRates_all_injection_split = splitReplicates(dataFrameToSplit = conversionRates_all_injection,condition = "Inj",
                                                      metadata_add = conversionRates_all_injection$name)

conversionRates_untreated = conversionRates_all_injection[,c('Untreated_TP1','Untreated_TP6','Untreated_TP9')]


######## now i want to get the percentage TC conversions... i.e (#TC converted * 100 )/ (#TC covered* number of TC containgin reads)
##### i.e conversion rate / total number of reads. 

#### reading in the total number of reads... 

readCounts = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures//figure2/data/totalCounts_allCws.txt",
                        stringsAsFactors = F,sep="\t",header = T)
readCounts = readCounts %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Inj_R3_TP9), sum, na.rm = TRUE) ### per gene read count
readCounts_split = splitReplicates(dataFrameToSplit = readCounts,condition = "Inj",
                                   metadata_add = readCounts$name)


### reads with single TC

readsWithTC = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure2/data//numberOfreadsWithTC.txt",
                         stringsAsFactors = F, sep="\t",header = T)                      

classifiedGenes = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/packages/data/countingWindows_classified_conversions.txt",
                             sep="\t",stringsAsFactors = F, header = T)
classifiedGenes = classifiedGenes %>% dplyr::select(c(gene,class))
readsWithTC = readsWithTC %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Untreated_TP9), sum, na.rm = TRUE) ### per gene read count
reads_untreated = readsWithTC[,c('Untreated_TP1','Untreated_TP6','Untreated_TP9')]
readsWithTC_split = splitReplicates(dataFrameToSplit = readsWithTC,condition = "Inj",
                                    metadata_add = readsWithTC$name)
#### creating the perentage TC score by multiplying this conversion TC with the read counts. 

percentageTCscore = conversionRates_all_injection_split

for(i in 1:length(readCounts_split)){
  percentageTCscore[[i]][,c(1:9)] = (conversionRates_all_injection_split[[i]][,c(1:9)] * 100) / readsWithTC_split[[i]][,c(1:9)]
}  

ptc_untreated = (conversionRates_untreated * 100)/reads_untreated
ptc_untreated$metadata_add= percentageTCscore$R2$metadata_add

exampleGenes = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure1/data/info_exampleMandZ.txt",
                          sep="\t",stringsAsFactors = F, header = T)
exampleGenes_Z = exampleGenes %>% dplyr::filter(type=="Zygotic")
exampleGenes_M = exampleGenes %>% dplyr::filter(type=="Maternal")
exampleGenes_Z$metadata_add = exampleGenes_Z$external_gene_name
exampleGenes_M$metadata_add = exampleGenes_M$external_gene_name


percentageTCscore_R1 = percentageTCscore$R2
percentageTCscore_R1[is.na(percentageTCscore_R1)]<-0

zygoticGenes_ptc_r1 = plyr::join(exampleGenes_Z,percentageTCscore_R1,by='metadata_add')
maternalGenes_ptc_r1 = plyr::join(exampleGenes_M,percentageTCscore_R1,by='metadata_add')

percentageTCscore_R2 = percentageTCscore$R3
percentageTCscore_R2[is.na(percentageTCscore_R2)]<-0


zygoticGenes_ptc_r2 = plyr::join(exampleGenes_Z,percentageTCscore_R2,by='metadata_add')
maternalGenes_ptc_r2 = plyr::join(exampleGenes_M,percentageTCscore_R2,by='metadata_add')

ptc_untreated[is.na(ptc_untreated)]<-0
zygoticGenes_ptc_untreated = plyr::join(exampleGenes_Z,ptc_untreated,by='metadata_add')
maternalGenes_ptc_untreated = plyr::join(exampleGenes_M,ptc_untreated,by='metadata_add')


seqs_samples = seq(2,0,by = -0.01)
TPR_seqs_r1 = c()
FPR_seqs_r1 = c()

TPR_seqs_r2 = c()
FPR_seqs_r2 = c()

FPR_seqs_unt_t1 = c()
TPR_seqs_unt_t1=c()

FPR_seqs_unt_t6 = c()
TPR_seqs_unt_t6=c()

FPR_seqs_unt_t9 = c()
TPR_seqs_unt_t9=c()

for(i in 1:length(seqs_samples)){
  
  TP_r1 = nrow(zygoticGenes_ptc_r1[zygoticGenes_ptc_r1$Inj_R2_TP9 >= seqs_samples[i],])
  FN_r1 = nrow(zygoticGenes_ptc_r1[zygoticGenes_ptc_r1$Inj_R2_TP9 < seqs_samples[i],])
  TPR_r1 = TP_r1/(TP_r1+FN_r1)
  TPR_seqs_r1 = c(TPR_seqs_r1,TPR_r1)
  
  FP_r1 =  nrow(maternalGenes_ptc_r1[maternalGenes_ptc_r1$Inj_R2_TP9 >= seqs_samples[i],])
  TN_r1=nrow(maternalGenes_ptc_r1[maternalGenes_ptc_r1$Inj_R2_TP9 < seqs_samples[i],])
  FPR_r1 = FP_r1/(FP_r1+TN_r1)
  FPR_seqs_r1 = c(FPR_seqs_r1,FPR_r1)
  
  TP_r2 = nrow(zygoticGenes_ptc_r2[zygoticGenes_ptc_r2$Inj_R3_TP9 >= seqs_samples[i],])
  FN_r2 = nrow(zygoticGenes_ptc_r2[zygoticGenes_ptc_r2$Inj_R3_TP9 < seqs_samples[i],])
  TPR_r2 = TP_r2/(TP_r2+FN_r2)
  TPR_seqs_r2 = c(TPR_seqs_r2,TPR_r2)
  
  FP_r2 =  nrow(maternalGenes_ptc_r2[maternalGenes_ptc_r2$Inj_R3_TP9 >= seqs_samples[i],])
  TN_r2=nrow(maternalGenes_ptc_r2[maternalGenes_ptc_r2$Inj_R3_TP9 < seqs_samples[i],])
  FPR_r2 = FP_r2/(FP_r2+TN_r2)
  FPR_seqs_r2 = c(FPR_seqs_r2,FPR_r2)
  
  ### untreated time points.
  
  ##### untreated tp1
  
  TP_unt = nrow(zygoticGenes_ptc_untreated[zygoticGenes_ptc_untreated$Untreated_TP1  >= seqs_samples[i],])
  FN_unt = nrow(zygoticGenes_ptc_untreated[zygoticGenes_ptc_untreated$Untreated_TP1 < seqs_samples[i],])
  TPR_unt = TP_unt/(TP_unt+FN_unt)
  TPR_seqs_unt_t1 = c(TPR_seqs_unt_t1,TPR_unt)
  
  
  FP_unt =  nrow(maternalGenes_ptc_untreated[maternalGenes_ptc_untreated$Untreated_TP1 >= seqs_samples[i],])
  TN_unt=nrow(maternalGenes_ptc_untreated[maternalGenes_ptc_untreated$Untreated_TP1 < seqs_samples[i],])
  FPR_unt = FP_unt/(FP_unt+TN_unt)
  FPR_seqs_unt_t1 = c(FPR_seqs_unt_t1,FPR_unt)

  #### untreated tp6
  TP_unt = nrow(zygoticGenes_ptc_untreated[zygoticGenes_ptc_untreated$Untreated_TP6  >= seqs_samples[i],])
  FN_unt = nrow(zygoticGenes_ptc_untreated[zygoticGenes_ptc_untreated$Untreated_TP6 < seqs_samples[i],])
  TPR_unt = TP_unt/(TP_unt+FN_unt)
  TPR_seqs_unt_t6 = c(TPR_seqs_unt_t6,TPR_unt)
  
  
  FP_unt =  nrow(maternalGenes_ptc_untreated[maternalGenes_ptc_untreated$Untreated_TP6 >= seqs_samples[i],])
  TN_unt=nrow(maternalGenes_ptc_untreated[maternalGenes_ptc_untreated$Untreated_TP6 < seqs_samples[i],])
  FPR_unt = FP_unt/(FP_unt+TN_unt)
  FPR_seqs_unt_t6 = c(FPR_seqs_unt_t6,FPR_unt)
  
  #### untreated tp9
  TP_unt = nrow(zygoticGenes_ptc_untreated[zygoticGenes_ptc_untreated$Untreated_TP9  >= seqs_samples[i],])
  FN_unt = nrow(zygoticGenes_ptc_untreated[zygoticGenes_ptc_untreated$Untreated_TP9 < seqs_samples[i],])
  TPR_unt = TP_unt/(TP_unt+FN_unt)
  TPR_seqs_unt_t9 = c(TPR_seqs_unt_t9,TPR_unt)
  
  
  FP_unt =  nrow(maternalGenes_ptc_untreated[maternalGenes_ptc_untreated$Untreated_TP9 >= seqs_samples[i],])
  TN_unt=nrow(maternalGenes_ptc_untreated[maternalGenes_ptc_untreated$Untreated_TP9 < seqs_samples[i],])
  FPR_unt = FP_unt/(FP_unt+TN_unt)
  FPR_seqs_unt_t9 = c(FPR_seqs_unt_t9,FPR_unt)
  
}

roc_r1 = data.frame(FP=FPR_seqs_r1,TP=TPR_seqs_r1,replicate='R1')
roc_r2 = data.frame(FP=FPR_seqs_r2,TP=TPR_seqs_r2,replicate='R2')
roc_unt_t1 = data.frame(FP=FPR_seqs_unt_t1,TP=TPR_seqs_unt_t1,replicate='Untreated_t1')
roc_unt_t6 = data.frame(FP=FPR_seqs_unt_t6,TP=TPR_seqs_unt_t6,replicate='Untreated_t6')
roc_unt_t9 = data.frame(FP=FPR_seqs_unt_t9,TP=TPR_seqs_unt_t9,replicate='Untreated_t9')

roc_ptc= rbind.data.frame(roc_r1,roc_r2,roc_unt_t1,roc_unt_t6,roc_unt_t9)
roc_ptc$threshold = seqs_samples

roc_ptc = ggplot(roc_ptc,aes(x=FP,y=TP,group=replicate,color=replicate))+ geom_line() + geom_path() + theme_cowplot()
roc_ptc = roc_ptc + theme_ameres(type = 'barplot')

roc_ptc_values= rbind.data.frame(roc_r1,roc_r2,roc_unt_t1,roc_unt_t6,roc_unt_t9)


roc_ptc_values$threshold = seqs_samples




roc_ptc_values_split= split(roc_ptc_values,roc_ptc_values$replicate,T)

roc_ptc_values_split = lapply(roc_ptc_values_split,function(x) x[order(x$TP,decreasing = T),] %>% dplyr::filter(FP==0))
ptcThresholds = lapply(roc_ptc_values_split,function(x) x[1,])
ptcThresholds = do.call(rbind.data.frame,ptcThresholds)


#################

pdf("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure2/plots/fp_tp_Ptc_Mtc.pdf", height = 3, width=6)
      
      
      TP_FP_pTC_R1 = ggplot(roc_ptc_values %>% dplyr::filter(replicate=='R1'),aes(x=threshold,y=TP)) +
        geom_line() + geom_line(data=roc_ptc_values %>% dplyr::filter(replicate=='R1'),aes(x=threshold,y=FP),color='red') + theme_cowplot()
      TP_FP_pTC_R1 = TP_FP_pTC_R1 + theme_ameres(type = 'barplot') + ggtitle('pTC R1')+ geom_vline(xintercept = ptcThresholds[1,'threshold'],linetype='dashed')
      
      TP_FP_pTC_R2 = ggplot(roc_ptc_values %>% dplyr::filter(replicate=='R2'),aes(x=threshold,y=TP)) +
        geom_line() + geom_line(data=roc_ptc_values %>% dplyr::filter(replicate=='R2'),aes(x=threshold,y=FP),color='red')+ theme_cowplot()
      TP_FP_pTC_R2 = TP_FP_pTC_R2 + theme_ameres(type = 'barplot')+ ggtitle('pTC R2')+ geom_vline(xintercept = ptcThresholds[2,'threshold'],linetype='dashed')
      
      TP_FP_pTC_R1+TP_FP_pTC_R2

      
      
      untreatedPtc_threshold_TP1 = reshape2::melt(ptc_untreated,id.vars=c('metadata_add')) %>% dplyr::filter(value >0 ) %>% dplyr::filter(variable == "Untreated_TP1")
      untreatedPtc_threshold_TP6 = reshape2::melt(ptc_untreated,id.vars=c('metadata_add')) %>% dplyr::filter(value >0 ) %>% dplyr::filter(variable == "Untreated_TP6")
      untreatedPtc_threshold_TP9 = reshape2::melt(ptc_untreated,id.vars=c('metadata_add')) %>% dplyr::filter(value >0 ) %>% dplyr::filter(variable == "Untreated_TP9")
      
      total_ptc_untreated = rbind.data.frame(untreatedPtc_threshold_TP1,untreatedPtc_threshold_TP6,untreatedPtc_threshold_TP9)
      total_ptc_untreated$ptc = log10(total_ptc_untreated$value)
      
      interval_tp1 = log10(t.test(untreatedPtc_threshold_TP1$value,conf.level = 0.95)$conf.int[2])
      interval_tp6 = log10(t.test(untreatedPtc_threshold_TP6$value,conf.level = 0.95)$conf.int[2])
      interval_tp9 = log10(t.test(untreatedPtc_threshold_TP9$value,conf.level = 0.95)$conf.int[2])
      
      ptc_untreated= ggpubr::ggdensity(total_ptc_untreated,x='ptc',color = 'variable', fill='variable') + theme_ameres(type = 'barplot') + 
        geom_vline(xintercept = interval_tp6)+ geom_vline(xintercept = interval_tp1) + geom_vline(xintercept = interval_tp9)+
        geom_vline(xintercept = log10(ptcThresholds[1,'threshold']),linetype='dashed') + geom_vline(xintercept = log10(ptcThresholds[2,'threshold']),linetype='dashed') 

      print(ptc_untreated)
      
dev.off()




########################## gretting the genes which are above the defined ptc thresholds ##################

percentageTCscore = lapply(percentageTCscore,function(x) x[x$metadata_add %in% classifiedGenes$gene ,])
percentageTCscore$R2  = percentageTCscore$R2 %>% dplyr::mutate_at(vars(dplyr::matches("Inj")),.funs = funs(pass= if_else(.> ptcThresholds$threshold[1], T,F )))
percentageTCscore$R3  = percentageTCscore$R3 %>% dplyr::mutate_at(vars(dplyr::matches("Inj")),.funs = funs(pass= if_else(.> ptcThresholds$threshold[2], T,F )))
genes_names = percentageTCscore$R2$metadata_add


t1_R1 = percentageTCscore$R2 %>% dplyr::filter(Inj_R2_TP1_pass == T) %>% dplyr::mutate(rep='R1',ptc=Inj_R2_TP1) %>% dplyr::select(c('metadata_add','rep','ptc'))
t1_R2 = percentageTCscore$R3 %>% dplyr::filter(Inj_R3_TP1_pass == T) %>% dplyr::mutate(rep='R2',ptc=Inj_R3_TP1) %>% dplyr::select(c('metadata_add','rep','ptc'))
total_t1 = rbind.data.frame(t1_R1,t1_R2) %>% dplyr::mutate(time = 0.75)


t2_R1 = percentageTCscore$R2 %>% dplyr::filter(Inj_R2_TP2_pass == T) %>% dplyr::mutate(rep='R1',ptc=Inj_R2_TP2) %>% dplyr::select(c('metadata_add','rep','ptc'))
t2_R2 = percentageTCscore$R3 %>% dplyr::filter(Inj_R3_TP2_pass == T) %>% dplyr::mutate(rep='R2',ptc=Inj_R3_TP2) %>% dplyr::select(c('metadata_add','rep','ptc'))
total_t2 = rbind.data.frame(t2_R1,t2_R2) %>% dplyr::mutate(time = 2)

t3_R1 = percentageTCscore$R2 %>% dplyr::filter(Inj_R2_TP3_pass == T) %>% dplyr::mutate(rep='R1',ptc=Inj_R2_TP3) %>% dplyr::select(c('metadata_add','rep','ptc'))
t3_R2 = percentageTCscore$R3 %>% dplyr::filter(Inj_R3_TP3_pass == T) %>% dplyr::mutate(rep='R2',ptc=Inj_R3_TP3) %>% dplyr::select(c('metadata_add','rep','ptc'))
total_t3 = rbind.data.frame(t3_R1,t3_R2) %>% dplyr::mutate(time = 2.5)


t4_R1 = percentageTCscore$R2 %>% dplyr::filter(Inj_R2_TP4_pass == T) %>% dplyr::mutate(rep='R1',ptc=Inj_R2_TP4) %>% dplyr::select(c('metadata_add','rep','ptc'))
t4_R2 = percentageTCscore$R3 %>% dplyr::filter(Inj_R3_TP4_pass == T) %>% dplyr::mutate(rep='R2',ptc=Inj_R3_TP4) %>% dplyr::select(c('metadata_add','rep','ptc'))
total_t4 = rbind.data.frame(t4_R1,t4_R2) %>% dplyr::mutate(time = 3)


t5_R1 = percentageTCscore$R2 %>% dplyr::filter(Inj_R2_TP5_pass == T) %>% dplyr::mutate(rep='R1',ptc=Inj_R2_TP5) %>% dplyr::select(c('metadata_add','rep','ptc'))
t5_R2 = percentageTCscore$R3 %>% dplyr::filter(Inj_R3_TP5_pass == T) %>% dplyr::mutate(rep='R2',ptc=Inj_R3_TP5) %>% dplyr::select(c('metadata_add','rep','ptc'))
total_t5 = rbind.data.frame(t5_R1,t5_R2) %>% dplyr::mutate(time = 3.5)


t6_R1 = percentageTCscore$R2 %>% dplyr::filter(Inj_R2_TP6_pass == T) %>% dplyr::mutate(rep='R1',ptc=Inj_R2_TP6) %>% dplyr::select(c('metadata_add','rep','ptc'))
t6_R2 = percentageTCscore$R3 %>% dplyr::filter(Inj_R3_TP6_pass == T) %>% dplyr::mutate(rep='R2',ptc=Inj_R3_TP6) %>% dplyr::select(c('metadata_add','rep','ptc'))
total_t6 = rbind.data.frame(t6_R1,t6_R2) %>% dplyr::mutate(time = 4)


t7_R1 = percentageTCscore$R2 %>% dplyr::filter(Inj_R2_TP7_pass == T) %>% dplyr::mutate(rep='R1',ptc=Inj_R2_TP7) %>% dplyr::select(c('metadata_add','rep','ptc'))
t7_R2 = percentageTCscore$R3 %>% dplyr::filter(Inj_R3_TP7_pass == T) %>% dplyr::mutate(rep='R2',ptc=Inj_R3_TP7) %>% dplyr::select(c('metadata_add','rep','ptc'))
total_t7 = rbind.data.frame(t7_R1,t7_R2) %>% dplyr::mutate(time = 4.5)

t8_R1 = percentageTCscore$R2 %>% dplyr::filter(Inj_R2_TP8_pass == T) %>% dplyr::mutate(rep='R1',ptc=Inj_R2_TP8) %>% dplyr::select(c('metadata_add','rep','ptc'))
t8_R2 = percentageTCscore$R3 %>% dplyr::filter(Inj_R3_TP8_pass == T) %>% dplyr::mutate(rep='R2',ptc=Inj_R3_TP8) %>% dplyr::select(c('metadata_add','rep','ptc'))
total_t8 = rbind.data.frame(t8_R1,t8_R2) %>% dplyr::mutate(time = 5)

t9_R1 = percentageTCscore$R2 %>% dplyr::filter(Inj_R2_TP9_pass == T) %>% dplyr::mutate(rep='R1',ptc=Inj_R2_TP9) %>% dplyr::select(c('metadata_add','rep','ptc'))
t9_R2 = percentageTCscore$R3 %>% dplyr::filter(Inj_R3_TP9_pass == T) %>% dplyr::mutate(rep='R2',ptc=Inj_R3_TP9) %>% dplyr::select(c('metadata_add','rep','ptc'))
total_t9 = rbind.data.frame(t9_R1,t9_R2) %>% dplyr::mutate(time = 5.5)

total_allGenes = rbind.data.frame(total_t1,total_t2,total_t3,total_t4,total_t5,total_t6,total_t7,total_t8,total_t9)
write.table(total_allGenes,'/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure2/data/genesAtDifferentTimepoints.txt',
            sep = '\t',quote = F,row.names = F,col.names = F)

totalList = percentageTCscore
save(  totalList,file = "/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure2/data/allGenes_percentageTC_multiTC.Rdata")


### once I known the genes above the threshold, I want to know which of these have conversion rates at later time points.. 
      ### this is done based on conversion rates
            
            errorRates = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure2/data/errorRates_predicted_observed.txt",header = T)
            conversionRates_all = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure2/data//perGenes_conversion.txt",sep="\t",stringsAsFactors = F,header=T)
            
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
            
            conversionRates_all_injection$name = conversionRates_all$name

            
            
            plotHeatmap = function(inTP1,conversionRates_all_injection){
              conversions_all_detected = conversionRates_all_injection[conversionRates_all_injection$name %in% inTP1$gene ,]
              conversions_all_detected_R2 = conversions_all_detected %>% dplyr::select(dplyr::contains("Inj_R2"))
              conversions_all_detected_R3 = conversions_all_detected %>% dplyr::select(dplyr::contains("Inj_R3"))
              conversions_all_detected_R2 = t(apply(conversions_all_detected_R2,1,function(x) x>0))
              conversions_all_detected_R3 = t(apply(conversions_all_detected_R3,1,function(x) x>0))
              
              conversions_all_detected_R3[which(conversions_all_detected_R3==T)] <-1
              conversions_all_detected_R3[which(conversions_all_detected_R3==F)] <-0
              conversions_all_detected_R2[which(conversions_all_detected_R2==T)] <-1
              conversions_all_detected_R2[which(conversions_all_detected_R2==F)] <-0
              
              
              tp1Expressed = conversions_all_detected_R2 + conversions_all_detected_R3
              tp1Expressed = as.data.frame(tp1Expressed)
              tpContaining = tp1Expressed[
                with(tp1Expressed, order(Inj_R2_TP9, Inj_R2_TP8,Inj_R2_TP7,Inj_R2_TP6,Inj_R2_TP5,Inj_R2_TP4,Inj_R2_TP3,Inj_R2_TP2)),
                ]
              rownames(tpContaining) = conversions_all_detected$name
              return(tpContaining)
            }

    #########################            

            ### 0.75hpf

          t1_r1 = percentageTCscore$R2 %>% dplyr::filter(Inj_R2_TP1_pass == T) %>% dplyr::mutate(replicate = 'R1') %>% dplyr::mutate(ptc_t1=Inj_R2_TP1,gene=metadata_add)  %>% dplyr::select(ptc_t1,gene,replicate)
          t1_r2 = percentageTCscore$R3 %>% dplyr::filter(Inj_R3_TP1_pass == T) %>% dplyr::mutate(replicate = 'R2')  %>% dplyr::mutate(ptc_t1=Inj_R3_TP1,gene=metadata_add)  %>% dplyr::select(ptc_t1,gene,replicate) 
          t1_total = rbind.data.frame(t1_r1,t1_r2)
          
          x = rownames(plotHeatmap(inTP1 = t1_total ,conversionRates_all_injection = conversionRates_all_injection ) %>% dplyr::filter(Inj_R2_TP9 !=0 | Inj_R2_TP8 !=0 | Inj_R2_TP7 !=0 | Inj_R2_TP6==0 ))
       
                          
          t1_total_decreasing = t1_total %>% dplyr::group_by(gene) %>% dplyr::summarise(ptc = mean(ptc_t1)) %>% dplyr::arrange(desc(ptc))
          t1_total_decreasing = t1_total_decreasing[as.character(t1_total_decreasing$gene) %in% x,]
          
          t1_total_increasing = t1_total %>% dplyr::group_by(gene) %>% dplyr::summarise(ptc = mean(ptc_t1)) %>% dplyr::arrange((ptc))
          t1_total_increasing = t1_total_increasing[as.character(t1_total_increasing$gene) %in% x,]
          
          seqs_samples = seq(50,250,by = 25)
          samples_gsea_t1ATACseq = vector('list',length(seqs_samples))
          names(samples_gsea_t1ATACseq) = seqs_samples
          
          samples_gsea_t1ATACseq_control = samples_gsea_t1ATACseq
          samples_gsea_t1h3k27acseq_control = samples_gsea_t1ATACseq
          samples_gsea_t1h3k27acseq = samples_gsea_t1ATACseq
          
          ranks <- read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Desktop/GSEA_samples/ATACseq/ATACseq.rnk",
                              header=F, colClasses = c("character", "numeric"))
          ranks <- setNames(ranks$V2, ranks$V1)
          
          
          ranks_h3k27ac <- read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Desktop/GSEA_samples/H3K27Ac//h3k27ac_ordered.rnk",
                              header=F, colClasses = c("character", "numeric"))
          ranks_h3k27ac = ranks_h3k27ac[!duplicated(ranks_h3k27ac$V1),]
          ranks_h3k27ac <- setNames(ranks_h3k27ac$V2, ranks_h3k27ac$V1)
          
          for(i in 1:length(samples_gsea_t1ATACseq)){

            pathways_ = list(as.character(t1_total_decreasing$gene[1:seqs_samples[i]]))
            names(pathways_) =seqs_samples[i]
            samples_gsea_t1ATACseq[[i]] = fgsea(pathways_, ranks, minSize=1, maxSize=5000,nperm = 1000,gseaParam=1)
            
            pathways_ = list(as.character(t1_total_increasing$gene[1:seqs_samples[i]]))
            names(pathways_) =seqs_samples[i]
             samples_gsea_t1ATACseq_control[[i]] = fgsea(pathways_, ranks, minSize=1, maxSize=5000,nperm = 1000,gseaParam=1)
            
             pathways_ = list(as.character(t1_total_decreasing$gene[1:seqs_samples[i]]))
             names(pathways_) =seqs_samples[i]
             samples_gsea_t1h3k27acseq[[i]] = fgsea(pathways_, ranks_h3k27ac, minSize=1, maxSize=5000,nperm = 1000,gseaParam=1)
             
             pathways_ = list(as.character(t1_total_increasing$gene[1:seqs_samples[i]]))
             names(pathways_) =seqs_samples[i]
             samples_gsea_t1h3k27acseq_control[[i]] = fgsea(pathways_, ranks_h3k27ac, minSize=1, maxSize=5000,nperm = 1000,gseaParam=1)
             
             
            
          }
          
          samples_gsea_t1ATACseq = do.call(rbind.data.frame,samples_gsea_t1ATACseq)
          samples_gsea_t1ATACseq$padj = -log10(samples_gsea_t1ATACseq$padj)
          samples_gsea_t1ATACseq$pathway = as.numeric(samples_gsea_t1ATACseq$pathway) 
          samples_gsea_t1ATACseq = samples_gsea_t1ATACseq %>% dplyr::mutate(tyep ='highest')
          
          samples_gsea_t1ATACseq_control = do.call(rbind.data.frame,samples_gsea_t1ATACseq_control)
          samples_gsea_t1ATACseq_control$padj = -log10(samples_gsea_t1ATACseq_control$padj)
          samples_gsea_t1ATACseq_control$pathway = as.numeric(samples_gsea_t1ATACseq_control$pathway)
          samples_gsea_t1ATACseq_control = samples_gsea_t1ATACseq_control %>% dplyr::mutate(tyep ='control')
          
          atacseq_t1 = rbind.data.frame(samples_gsea_t1ATACseq,samples_gsea_t1ATACseq_control)
          
          atacseq_t1 = ggplot(atacseq_t1 ,aes(x=pathway,y=NES,fill=tyep,alpha=padj>1.30103)) + geom_bar(position='dodge', stat = "identity")+
            scale_alpha_discrete(range = c(0.2, 0.9)) + xlab ("Top/bottom n genes") + ggtitle("T1 - ATACseq") + theme_cowplot()
          atacseq_t1 = atacseq_t1 + theme_ameres(type = 'barplot') + scale_x_continuous(breaks=c(50,75,100,125,150,175,200,225,250))
           
         ### same for h3k27ac
          
          samples_gsea_t1h3k27acseq = do.call(rbind.data.frame,samples_gsea_t1h3k27acseq)
          samples_gsea_t1h3k27acseq$padj = -log10(samples_gsea_t1h3k27acseq$padj)
          samples_gsea_t1h3k27acseq$pathway = as.numeric(samples_gsea_t1h3k27acseq$pathway) 
          samples_gsea_t1h3k27acseq = samples_gsea_t1h3k27acseq %>% dplyr::mutate(tyep ='highest')
          
          samples_gsea_t1h3k27acseq_control = do.call(rbind.data.frame,samples_gsea_t1h3k27acseq_control)
          samples_gsea_t1h3k27acseq_control$padj = -log10(samples_gsea_t1h3k27acseq_control$padj)
          samples_gsea_t1h3k27acseq_control$pathway = as.numeric(samples_gsea_t1h3k27acseq_control$pathway)
          samples_gsea_t1h3k27acseq_control = samples_gsea_t1h3k27acseq_control %>% dplyr::mutate(tyep ='control')
          
          h3k27ac_t1 = rbind.data.frame(samples_gsea_t1h3k27acseq,samples_gsea_t1h3k27acseq_control)
          
          h3k27ac_t1 = ggplot(h3k27ac_t1 ,aes(x=pathway,y=NES,fill=tyep,alpha=padj>1.30103)) + geom_bar(position='dodge', stat = "identity")+
            scale_alpha_discrete(range = c(0.2, 0.9)) + xlab ("Top/bottom n genes") + ggtitle("T1 - H3K27Ac")+ theme_cowplot()
          h3k27ac_t1 = h3k27ac_t1 + theme_ameres(type = 'barplot') + scale_x_continuous(breaks=c(50,75,100,125,150,175,200,225,250))
          
          
          
 ############# 
        
         t2_r1 = percentageTCscore$R2 %>% dplyr::filter(Inj_R2_TP2_pass == T) %>% dplyr::mutate(replicate = 'R1') %>% dplyr::mutate(ptc_t2=Inj_R2_TP2,gene=metadata_add)  %>% dplyr::select(ptc_t2,gene,replicate)
         t2_r2 = percentageTCscore$R3 %>% dplyr::filter(Inj_R3_TP2_pass == T) %>% dplyr::mutate(replicate = 'R2')  %>% dplyr::mutate(ptc_t2=Inj_R3_TP2,gene=metadata_add)  %>% dplyr::select(ptc_t2,gene,replicate) 
         t2_total = rbind.data.frame(t2_r1,t2_r2)
         
         x = rownames(plotHeatmap(inTP1 = t2_total ,conversionRates_all_injection = conversionRates_all_injection ) %>% dplyr::filter(Inj_R2_TP9 !=0 | Inj_R2_TP8 !=0 | Inj_R2_TP7 !=0 | Inj_R2_TP6 !=0))
         
         
         t2_total_decreasing = t2_total %>% dplyr::group_by(gene) %>% dplyr::summarise(ptc = mean(ptc_t2)) %>% dplyr::arrange(desc(ptc))
         t2_total_decreasing = t2_total_decreasing[as.character(t2_total_decreasing$gene) %in% x,]
         
         t2_total_increasing = t2_total %>% dplyr::group_by(gene) %>% dplyr::summarise(ptc = mean(ptc_t2)) %>% dplyr::arrange((ptc))
         t2_total_increasing = t2_total_increasing[as.character(t2_total_increasing$gene) %in% x,]
         
         seqs_samples = seq(50,250,by = 25)
         samples_gsea_t2ATACseq = vector('list',length(seqs_samples))
         names(samples_gsea_t2ATACseq) = seqs_samples
         
         samples_gsea_t2ATACseq_control = samples_gsea_t2ATACseq
         samples_gsea_t2h3k27acseq_control = samples_gsea_t2ATACseq
         samples_gsea_t2h3k27acseq = samples_gsea_t2ATACseq
         
          
         for(i in 1:length(samples_gsea_t2ATACseq)){
           
           pathways_ = list(as.character(t2_total_decreasing$gene[1:seqs_samples[i]]))
           names(pathways_) =seqs_samples[i]
           samples_gsea_t2ATACseq[[i]] = fgsea(pathways_, ranks, minSize=1, maxSize=5000,nperm = 1000,gseaParam=1)
           
           pathways_ = list(as.character(t2_total_increasing$gene[1:seqs_samples[i]]))
           names(pathways_) =seqs_samples[i]
           samples_gsea_t2ATACseq_control[[i]] = fgsea(pathways_, ranks, minSize=1, maxSize=5000,nperm = 1000,gseaParam=1)
           
           pathways_ = list(as.character(t2_total_decreasing$gene[1:seqs_samples[i]]))
           names(pathways_) =seqs_samples[i]
           samples_gsea_t2h3k27acseq[[i]] = fgsea(pathways_, ranks_h3k27ac, minSize=1, maxSize=5000,nperm = 1000,gseaParam=1)
           
           pathways_ = list(as.character(t2_total_increasing$gene[1:seqs_samples[i]]))
           names(pathways_) =seqs_samples[i]
           samples_gsea_t2h3k27acseq_control[[i]] = fgsea(pathways_, ranks_h3k27ac, minSize=1, maxSize=5000,nperm = 1000,gseaParam=1)
           
         }
        
         samples_gsea_t2ATACseq = do.call(rbind.data.frame,samples_gsea_t2ATACseq)
         samples_gsea_t2ATACseq$padj = -log10(samples_gsea_t2ATACseq$padj)
         samples_gsea_t2ATACseq$pathway = as.numeric(samples_gsea_t2ATACseq$pathway) 
         samples_gsea_t2ATACseq = samples_gsea_t2ATACseq %>% dplyr::mutate(tyep ='highest')
         
         samples_gsea_t2ATACseq_control = do.call(rbind.data.frame,samples_gsea_t2ATACseq_control)
         samples_gsea_t2ATACseq_control$padj = -log10(samples_gsea_t2ATACseq_control$padj)
         samples_gsea_t2ATACseq_control$pathway = as.numeric(samples_gsea_t2ATACseq_control$pathway)
         samples_gsea_t2ATACseq_control = samples_gsea_t2ATACseq_control %>% dplyr::mutate(tyep ='control')
         
         atacseq_t2 = rbind.data.frame(samples_gsea_t2ATACseq,samples_gsea_t2ATACseq_control)
         
         atacseq_t2 =   ggplot(atacseq_t2 ,aes(x=pathway,y=NES,fill=tyep,alpha=padj>1.30103)) + geom_bar(position='dodge', stat = "identity")+
           scale_alpha_discrete(range = c(0.2, 0.9)) + xlab ("Top/bottom n genes") + ggtitle("T2 - ATACseq") + theme_cowplot()
         atacseq_t2 = atacseq_t2 + theme_ameres(type = 'barplot')  + scale_x_continuous(breaks=c(50,75,100,125,150,175,200,225,250))
         
         ### same for h3k27ac
         
         samples_gsea_t2h3k27acseq = do.call(rbind.data.frame,samples_gsea_t2h3k27acseq)
         samples_gsea_t2h3k27acseq$padj = -log10(samples_gsea_t2h3k27acseq$padj)
         samples_gsea_t2h3k27acseq$pathway = as.numeric(samples_gsea_t2h3k27acseq$pathway) 
         samples_gsea_t2h3k27acseq = samples_gsea_t2h3k27acseq %>% dplyr::mutate(tyep ='highest')
         
         samples_gsea_t2h3k27acseq_control = do.call(rbind.data.frame,samples_gsea_t2h3k27acseq_control)
         samples_gsea_t2h3k27acseq_control$padj = -log10(samples_gsea_t2h3k27acseq_control$padj)
         samples_gsea_t2h3k27acseq_control$pathway = as.numeric(samples_gsea_t2h3k27acseq_control$pathway)
         samples_gsea_t2h3k27acseq_control = samples_gsea_t2h3k27acseq_control %>% dplyr::mutate(tyep ='control')
         
         h3k27ac_t2 = rbind.data.frame(samples_gsea_t2h3k27acseq,samples_gsea_t2h3k27acseq_control)
         
         h3k27ac_t2 = ggplot(h3k27ac_t2 ,aes(x=pathway,y=NES,fill=tyep,alpha=padj>1.30103)) + geom_bar(position='dodge', stat = "identity")+
           scale_alpha_discrete(range = c(0.2, 0.9)) + xlab ("Top/bottom n genes") + ggtitle("T2 - H3K27Ac")+ theme_cowplot()
         h3k27ac_t2 = h3k27ac_t2 + theme_ameres(type = 'barplot')  + scale_x_continuous(breaks=c(50,75,100,125,150,175,200,225,250))
         
        pdf("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure2/plots/atacSeq_h3k27AC.pdf") 
              atacseq_t1 = atacseq_t1 + ylim(c(-2,3) )
              atacseq_t2 = atacseq_t2 + ylim(c(-2,3) )
              h3k27ac_t1 = h3k27ac_t1 +  ylim(c(-2,3) )
              h3k27ac_t2 = h3k27ac_t2 +  ylim(c(-2,3) )
              (atacseq_t1 + atacseq_t2)/(h3k27ac_t1 + h3k27ac_t2)
        dev.off()
         ######### plotting properties of introns for top 100, top 200 etc.

         
         inrtonData = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure2/data//intronInfo_proteinCodingGenes.bed",
                                 sep="\t",stringsAsFactors = F)
         inrtonData$transcript = unlist(lapply(strsplit(inrtonData$V4,"_",T),function(x) x[1]))
         inrtonData = inrtonData %>% dplyr::group_by(transcript) %>% dplyr::mutate(numberOfIntrons=n()) %>% dplyr::ungroup() %>%
           dplyr::group_by(V7) %>% dplyr::mutate(maxIntrons = max(numberOfIntrons))
         inrtonData_distinct = inrtonData[!duplicated(inrtonData$V7),] %>% dplyr::select(c('maxIntrons','V7')) %>% dplyr::mutate(gene = V7)
         inrtonData_distinct = as.data.frame(inrtonData_distinct)

         total_t2 = vector('list',length(seqs_samples))
         names(total_t2) = seqs_samples
         
         total_t1 = vector('list',length(seqs_samples))
         names(total_t1) = seqs_samples
         
         
         total_t2_increasing = vector('list',length(seqs_samples))
         names(total_t2) = seqs_samples
         
         total_t1_increasing = vector('list',length(seqs_samples))
         names(total_t1) = seqs_samples
         
         
         `%!in%` = Negate(`%in%`)

         
         for(i in 1:length(seqs_samples)){
           
           topGene_t2 = t2_total_decreasing$gene[1:seqs_samples[i]]
           intronDistinct_top = inrtonData_distinct[inrtonData_distinct$gene %in% as.character(topGene_t2),] %>% dplyr::mutate(type = 'highest')
           intronDistinct_other = inrtonData_distinct[inrtonData_distinct$gene %in% classifiedGenes$gene,]  %>% dplyr::mutate(type = 'other')
           bottomGene_t2 = t2_total_increasing$gene[1:seqs_samples[i]]
           intronDistinct_bottom = inrtonData_distinct[inrtonData_distinct$gene %in% as.character(bottomGene_t2),] %>% dplyr::mutate(type = 'bottom')
           
           
           
           total_t2[[i]]  =  rbind.data.frame(intronDistinct_top,intronDistinct_bottom,intronDistinct_other) %>% dplyr::mutate(threshold= seqs_samples[i])
           
              
                ## t1
           
           topGene_t1 = t1_total_decreasing$gene[1:seqs_samples[i]]
           intronDistinct_top = inrtonData_distinct[inrtonData_distinct$gene %in% as.character(topGene_t1),] %>% dplyr::mutate(type = 'highest')
           intronDistinct_other = inrtonData_distinct[inrtonData_distinct$gene %in% classifiedGenes$gene,]  %>% dplyr::mutate(type = 'other')
           
           bottomGene_t1 = t1_total_increasing$gene[1:seqs_samples[i]]
           intronDistinct_bottom = inrtonData_distinct[inrtonData_distinct$gene %in% as.character(bottomGene_t1),] %>% dplyr::mutate(type = 'bottom')
           
           total_t1[[i]]  =  rbind.data.frame(intronDistinct_top,intronDistinct_bottom,intronDistinct_other) %>% dplyr::mutate(threshold= seqs_samples[i])
           
           
              
         }
         mycomparisons = list(c('top','other'),c('bottom','other'))


         total_t2 =  do.call(rbind.data.frame,total_t2) %>% dplyr::mutate(maxIntrons = log10(maxIntrons)) 
         total_t1 =  do.call(rbind.data.frame,total_t1) %>% dplyr::mutate(maxIntrons = log10(maxIntrons)) 
         
         total_t2_highest_other  =total_t2 %>%  dplyr::filter(type != 'bottom')
         total_t2_bottom_other  =total_t2 %>%  dplyr::filter(type != 'highest')
        
         total_t1_highest_other  =total_t1 %>%  dplyr::filter(type != 'bottom')
         total_t1_bottom_other  =total_t1 %>%  dplyr::filter(type != 'highest')
         
                   
         t2_introns =   total_t2 %>% dplyr::group_by(type,threshold)%>% dplyr::summarise(median = median(maxIntrons),iqr_75 = quantile(maxIntrons,0.75), iqr_25 = quantile(maxIntrons,0.25))
         t1_introns = total_t1 %>% dplyr::group_by(type,threshold)%>% dplyr::summarise(median = median(maxIntrons),iqr_75 = quantile(maxIntrons,0.75), iqr_25 = quantile(maxIntrons,0.25))
           
         t2_introns = reshape2::melt(t2_introns ,id.vars= c('type','threshold')) %>% dplyr::mutate(time = 'T2')
         t1_introns = reshape2::melt(t1_introns ,id.vars= c('type','threshold'))%>% dplyr::mutate(time = 'T1')
         
              
        
      intronsGenes = rbind.data.frame(t1_introns,t2_introns)
         
     intronsPlot =  ggplot() + geom_line(data=intronsGenes %>% dplyr::filter(type == 'other'),
                                            aes(x=threshold,y=value,group=variable),linetype = "dashed") + 
       geom_line(data=intronsGenes %>% 
                   dplyr::filter(type == 'highest' & variable=='median'),aes(x=threshold,y=value,color=time),size = 2) +
       theme_cowplot() + ylim(c(0,1.5)) 
     intronsPlot =  intronsPlot + theme_ameres(type = 'barplot') + 
       scale_x_continuous(limits = c(50,250), breaks = seq(50,250, by = 25)) + 
       geom_point(data=intronsGenes %>% 
                   dplyr::filter(type == 'highest' & variable=='median'),aes(x=threshold,y=value,color=time),size = 2)
     
     intronsPlot_control =  ggplot() + geom_line(data=intronsGenes %>% dplyr::filter(type == 'other'),
                                         aes(x=threshold,y=value,group=variable),linetype = "dashed") + 
       geom_line(data=intronsGenes %>% 
                   dplyr::filter(type == 'bottom' & variable=='median'),aes(x=threshold,y=value,color=time),size = 2) +
       theme_cowplot() + ylim(c(0,1.5))
     intronsPlot_control =  intronsPlot_control + theme_ameres(type = 'barplot') + 
       scale_x_continuous(limits = c(50,250), breaks = seq(50,250, by = 25))+ 
       geom_point(data=intronsGenes %>% 
                    dplyr::filter(type == 'bottom' & variable=='median'),aes(x=threshold,y=value,color=time),size = 2)
     

     pdf("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure2/plots/intonsAndCDS_basedOnPTC.pdf", 
         height=4, width = 10)
     ggpubr::ggboxplot(total_t2_highest_other,x='threshold',y='maxIntrons',fill='type') + stat_compare_means(aes(group=type),label = "p.signif", method.args = list(alternative='greater'))
     ggpubr::ggboxplot(total_t2_bottom_other,x='threshold',y='maxIntrons',fill='type') + stat_compare_means(aes(group=type),label = "p.signif",title='T2',method.args = list(alternative='greater'))
     
     ggpubr::ggboxplot(total_t1_highest_other,x='threshold',y='maxIntrons',fill='type') + stat_compare_means(aes(group=type),label = "p.signif",title='T1')
     ggpubr::ggboxplot(total_t1_bottom_other,x='threshold',y='maxIntrons',fill='type') + stat_compare_means(aes(group=type),label = "p.signif",title='T1')
     
     
        intronsPlot + intronsPlot_control
     
############## same for CDS. ################################
         
         library("biomaRt")
         listMarts()
         ensembl=useMart("ensembl")
         listDatasets(ensembl)
         ensembl = useDataset("drerio_gene_ensembl",mart=ensembl)
         
         allGenes = getBM(attributes = c("transcript_start","transcript_end","external_gene_name"),filters = "external_gene_name",mart = ensembl,values = classifiedGenes$gene)
         allGenes = allGenes %>% dplyr::mutate(transcript_length= transcript_end-transcript_start) %>% dplyr::mutate(gene = external_gene_name)
         allGenes = allGenes[complete.cases(allGenes),]
         allGenes = allGenes %>% dplyr::group_by(gene) %>% dplyr::mutate(avgLength = max(transcript_length)/1000) 
         allGenes = allGenes[!duplicated(allGenes[,c('avgLength','gene')]),]

         total_t2 = vector('list',length(seqs_samples))
         names(total_t2) = seqs_samples
         
         total_t1 = vector('list',length(seqs_samples))
         names(total_t1) = seqs_samples
         
         total_t2_increasing = vector('list',length(seqs_samples))
         names(total_t2_increasing) = seqs_samples
         
         total_t1_increasing = vector('list',length(seqs_samples))
         names(total_t1_increasing) = seqs_samples
         
         
    for(i in 1:length(seqs_samples)){
      
        
         topGene_t2 = t2_total_decreasing$gene[1:seqs_samples[i]]
         #notTop_t2 = classifiedGenes[classifiedGenes$gene %!in% topGene_t2,'gene']
         bottomGene_t2 = t2_total_increasing$gene[1:seqs_samples[i]]
         
         cds_top = as.data.frame(allGenes[allGenes$gene %in% as.character(topGene_t2),] %>% dplyr::mutate(type = 'highest'))
         cds_other = as.data.frame(allGenes[allGenes$gene %in% classifiedGenes$gene,]  %>% dplyr::mutate(type = 'other'))
         cds_bottom = as.data.frame(allGenes[allGenes$gene %in% as.character(bottomGene_t2),] %>% dplyr::mutate(type = 'bottom'))
         
           total_t2[[i]]=  rbind.data.frame(cds_bottom,cds_top,cds_other) %>% dplyr::mutate(threshold= seqs_samples[i]) %>% dplyr::mutate(avgLength = log10(avgLength))
     
           topGene_t1 = t1_total_decreasing$gene[1:seqs_samples[i]]
           #notTop_t2 = classifiedGenes[classifiedGenes$gene %!in% topGene_t2,'gene']
           bottomGene_t1 = t1_total_increasing$gene[1:seqs_samples[i]]
           
           cds_top = as.data.frame(allGenes[allGenes$gene %in% as.character(topGene_t1),] %>% dplyr::mutate(type = 'highest'))
           cds_other = as.data.frame(allGenes[allGenes$gene %in% classifiedGenes$gene,]  %>% dplyr::mutate(type = 'other'))
           cds_bottom = as.data.frame(allGenes[allGenes$gene %in% as.character(bottomGene_t1),] %>% dplyr::mutate(type = 'bottom'))
           
           total_t1[[i]]=  rbind.data.frame(cds_bottom,cds_top,cds_other) %>% dplyr::mutate(threshold= seqs_samples[i]) %>% dplyr::mutate(avgLength = log10(avgLength))
           
           
    }
         
          
      
         
         total_t2 =  rbind.data.frame(do.call(rbind.data.frame,total_t2) )
         total_t1 =  rbind.data.frame(do.call(rbind.data.frame,total_t1))
        
         total_t2_highest_other  =total_t2 %>%  dplyr::filter(type != 'bottom')
         total_t2_bottom_other  =total_t2 %>%  dplyr::filter(type != 'highest')
         
         total_t1_highest_other  =total_t1 %>%  dplyr::filter(type != 'bottom')
         total_t1_bottom_other  =total_t1 %>%  dplyr::filter(type != 'highest')
         
         ggpubr::ggboxplot(total_t2_highest_other,x='threshold',y='avgLength',fill='type') + stat_compare_means(aes(group=type),label = "p.signif",main='T2')
         ggpubr::ggboxplot(total_t2_bottom_other,x='threshold',y='avgLength',fill='type') + stat_compare_means(aes(group=type),label = "p.signif",main='T2')
         
         ggpubr::ggboxplot(total_t1_highest_other,x='threshold',y='avgLength',fill='type') + stat_compare_means(aes(group=type),label = "p.signif",main='T1')
         ggpubr::ggboxplot(total_t1_bottom_other,x='threshold',y='avgLength',fill='type') + stat_compare_means(aes(group=type),label = "p.signif",main='T1')
         
         
         t2_cds =   total_t2 %>% dplyr::group_by(type,threshold)%>% dplyr::summarise(median = median(avgLength),iqr_75 = quantile(avgLength,0.75), iqr_25 = quantile(avgLength,0.25))
         t1_cds = total_t1 %>% dplyr::group_by(type,threshold)%>% dplyr::summarise(median = median(avgLength),iqr_75 = quantile(avgLength,0.75), iqr_25 = quantile(avgLength,0.25))
         
         t2_cds = reshape2::melt(t2_cds ,id.vars= c('type','threshold')) %>% dplyr::mutate(time = 'T2')
         t1_cds = reshape2::melt(t1_cds ,id.vars= c('type','threshold'))%>% dplyr::mutate(time = 'T1')
         
         
         
        CDSGenes = rbind.data.frame(t1_cds,t2_cds)
         
         cdsPlot =  ggplot() + geom_line(data=CDSGenes %>% dplyr::filter(type == 'other'),
                                             aes(x=threshold,y=value,group=variable),linetype = "dashed") + 
           geom_line(data=CDSGenes %>% 
                       dplyr::filter(type == 'highest' & variable=='median'),aes(x=threshold,y=value,color=time),size = 2) +
           theme_cowplot() + ylim(c(0,2))
         cdsPlot =  cdsPlot + theme_ameres(type = 'barplot') + 
           scale_x_continuous(limits = c(50,250), breaks = seq(50,250, by = 25)) +  geom_point(data=CDSGenes %>% 
           dplyr::filter(type == 'highest' & variable=='median'),aes(x=threshold,y=value,color=time),size = 2) 
         
       cds_control =  ggplot() + geom_line(data=CDSGenes %>% dplyr::filter(type == 'other'),
                                                     aes(x=threshold,y=value,group=variable),linetype = "dashed") + 
           geom_line(data=CDSGenes %>% 
                       dplyr::filter(type == 'bottom' & variable=='median'),aes(x=threshold,y=value,color=time),size = 2) +
           theme_cowplot() + ylim(c(0,2))
       cds_control =  cds_control + theme_ameres(type = 'barplot') + 
           scale_x_continuous(limits = c(50,250), breaks = seq(50,250, by = 25))  + geom_point(data=CDSGenes %>% 
        dplyr::filter(type == 'bottom' & variable=='median'),aes(x=threshold,y=value,color=time),size = 2) 
         
       cdsPlot + cds_control
       
dev.off()
  
  
###########################################  
  
  total_allGenes_split = split(total_allGenes,total_allGenes$time,T)
  total_allGenes_split = lapply(total_allGenes_split,function(x) x %>% dplyr::group_by(metadata_add) %>% dplyr::mutate(ptc= mean(ptc)) %>% dplyr::distinct(metadata_add,.keep_all=T))
  
  total_allGenes_split = lapply(total_allGenes_split,function(x)  x%>% dplyr::arrange(desc(ptc) ))
  total_allGenes_split = lapply(total_allGenes_split,function(x)  plyr::join(as.data.frame(x),inrtonData_distinct))
  total_allGenes_split = lapply(total_allGenes_split,function(x) x[complete.cases(x),])
  intronLengths  = lapply(total_allGenes_split,function(x) (x[1:50,'maxIntrons']))
  intronLengths = reshape2::melt(intronLengths)
  intronLengths$value = log10(intronLengths$value)
  ggpubr::ggviolin(intronLengths,'L1','value',add = 'boxplot')
  
  total_allGenes_split = split(total_allGenes,total_allGenes$time,T)
  total_allGenes_split = lapply(total_allGenes_split,function(x) x %>% dplyr::group_by(metadata_add) %>% dplyr::mutate(ptc= mean(ptc)) %>% dplyr::distinct(metadata_add,.keep_all=T))
  
  total_allGenes_split = lapply(total_allGenes_split,function(x)  x%>% dplyr::arrange(desc(ptc) ))
  total_allGenes_split = lapply(total_allGenes_split,function(x)  plyr::join(as.data.frame(x),allGenes))
  total_allGenes_split = lapply(total_allGenes_split,function(x) x[complete.cases(x),])
  intronLengths  = lapply(total_allGenes_split,function(x) (x[1:75,'avgLength']))
  intronLengths = reshape2::melt(intronLengths)
  intronLengths$value = log10(intronLengths$value)
  ggpubr::ggviolin(intronLengths,'L1','value',add = 'boxplot')
  
  
    
  