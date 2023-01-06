
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
###

classification_slamSeq =  read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/packages/data/countingWindows_classified_conversions.txt",sep="\t",stringsAsFactors = F,header = T)

errorRates = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure2/data//errorRates_predicted_observed.txt",header = T)
conversionRates_all = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure2//data/perGenes_conversion.txt",sep="\t",stringsAsFactors = F,header=T)

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


TC_samples = read.table('/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure2/data//numberOfreadsWithMultipleTC.txt',
                        sep = "\t",stringsAsFactors = F, header = T)

allTC  = read.table('/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure2/data/numberOfreadsWithTC.txt',
                    sep="\t",stringsAsFactors = F, header = T)

TA_samples = read.table('/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures//figure2/data/numberOfreadsWithTA.txt',
                        sep="\t",stringsAsFactors = F, header = T)

TAmultiple_samples = read.table('/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures//figure2/data/numberOfreadsWithMultipleTA.txt',
                                sep="\t",stringsAsFactors = F, header = T)

totalReads = read.table('//Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure2/data//totalCounts_allCws.txt',
                        sep="\t",stringsAsFactors = F, header = T)

TC_samples = TC_samples %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Untreated_TP9), sum, na.rm = TRUE) %>% tibble::column_to_rownames('name')
TA_samples = TA_samples %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Untreated_TP9), sum, na.rm = TRUE) %>% tibble::column_to_rownames('name')
TC_all = allTC %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Untreated_TP9), sum, na.rm = TRUE) %>% tibble::column_to_rownames('name')
totalReads = totalReads %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Untreated_TP9), sum, na.rm = TRUE) %>% tibble::column_to_rownames('name')
TAmultiple_samples = TAmultiple_samples %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Untreated_TP9), sum, na.rm = TRUE) %>% tibble::column_to_rownames('name')
names_genes = rownames(TAmultiple_samples)


###### combine the replicates...

TC_all = TC_all %>% dplyr::mutate(TP1 = Inj_R2_TP1+Inj_R3_TP1,TP2 = Inj_R2_TP2 + Inj_R3_TP2,
                                  TP3 = Inj_R2_TP3 + Inj_R3_TP3 , TP4 = Inj_R2_TP4 + Inj_R3_TP4, TP5 = Inj_R2_TP5 + Inj_R3_TP5,
                                  TP6 = Inj_R2_TP6 + Inj_R3_TP6, TP7 = Inj_R2_TP7 + Inj_R3_TP7, TP8 = Inj_R2_TP8 + Inj_R3_TP8, TP9 = Inj_R2_TP9 + Inj_R3_TP9)  %>%
                                dplyr::select(c(TP1:TP9, Untreated_TP1, Untreated_TP6, Untreated_TP9))


TC_samples = TC_samples %>% dplyr::mutate(TP1 = Inj_R2_TP1+Inj_R3_TP1,TP2 = Inj_R2_TP2 + Inj_R3_TP2,
                         TP3 = Inj_R2_TP3 + Inj_R3_TP3 , TP4 = Inj_R2_TP4 + Inj_R3_TP4, TP5 = Inj_R2_TP5 + Inj_R3_TP5,
                         TP6 = Inj_R2_TP6 + Inj_R3_TP6, TP7 = Inj_R2_TP7 + Inj_R3_TP7, TP8 = Inj_R2_TP8 + Inj_R3_TP8, TP9 = Inj_R2_TP9 + Inj_R3_TP9)  %>%
                        dplyr::select(c(TP1:TP9, Untreated_TP1, Untreated_TP6, Untreated_TP9))

TA_samples = TA_samples %>% dplyr::mutate(TP1 = Inj_R2_TP1+Inj_R3_TP1,TP2 = Inj_R2_TP2 + Inj_R3_TP2,
                             TP3 = Inj_R2_TP3 + Inj_R3_TP3 , TP4 = Inj_R2_TP4 + Inj_R3_TP4, TP5 = Inj_R2_TP5 + Inj_R3_TP5,
                             TP6 = Inj_R2_TP6 + Inj_R3_TP6, TP7 = Inj_R2_TP7 + Inj_R3_TP7, TP8 = Inj_R2_TP8 + Inj_R3_TP8, TP9 = Inj_R2_TP9 + Inj_R3_TP9)  %>%
                             dplyr::select(c(TP1:TP9, Untreated_TP1, Untreated_TP6, Untreated_TP9))


TAmultiple_samples = TAmultiple_samples %>% dplyr::mutate(TP1 = Inj_R2_TP1+Inj_R3_TP1,TP2 = Inj_R2_TP2 + Inj_R3_TP2,
                             TP3 = Inj_R2_TP3 + Inj_R3_TP3 , TP4 = Inj_R2_TP4 + Inj_R3_TP4, TP5 = Inj_R2_TP5 + Inj_R3_TP5,
                             TP6 = Inj_R2_TP6 + Inj_R3_TP6, TP7 = Inj_R2_TP7 + Inj_R3_TP7, TP8 = Inj_R2_TP8 + Inj_R3_TP8, TP9 = Inj_R2_TP9 + Inj_R3_TP9)  %>%
                           dplyr::select(c(TP1:TP9, Untreated_TP1, Untreated_TP6, Untreated_TP9))


totalReads = totalReads %>% dplyr::mutate(TP1 = Inj_R2_TP1+Inj_R3_TP1,TP2 = Inj_R2_TP2 + Inj_R3_TP2,
                             TP3 = Inj_R2_TP3 + Inj_R3_TP3 , TP4 = Inj_R2_TP4 + Inj_R3_TP4, TP5 = Inj_R2_TP5 + Inj_R3_TP5,
                             TP6 = Inj_R2_TP6 + Inj_R3_TP6, TP7 = Inj_R2_TP7 + Inj_R3_TP7, TP8 = Inj_R2_TP8 + Inj_R3_TP8, TP9 = Inj_R2_TP9 + Inj_R3_TP9) %>%
                             dplyr::select(c(TP1:TP9, Untreated_TP1, Untreated_TP6, Untreated_TP9))

conversionRates_all_injection = conversionRates_all_injection %>% dplyr::mutate(TP1 = Inj_R2_TP1+Inj_R3_TP1,TP2 = Inj_R2_TP2 + Inj_R3_TP2,
                                                                                TP3 = Inj_R2_TP3 + Inj_R3_TP3 , TP4 = Inj_R2_TP4 + Inj_R3_TP4, TP5 = Inj_R2_TP5 + Inj_R3_TP5,
                                                                                TP6 = Inj_R2_TP6 + Inj_R3_TP6, TP7 = Inj_R2_TP7 + Inj_R3_TP7, TP8 = Inj_R2_TP8 + Inj_R3_TP8, TP9 = Inj_R2_TP9 + Inj_R3_TP9) %>%
  dplyr::select(c(TP1:TP9, Untreated_TP1, Untreated_TP6, Untreated_TP9,name))


### normalizing to library depth 

totalReads_scale = 1000000/colSums(totalReads[,1:ncol(totalReads)] ) 
initialFiles = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/analysis/UTRannotation/data/prePocessingStats.txt",sep="\t", header = T)
colSums(TC_all)/colSums(totalReads)

for(i in 1:length(totalReads_scale)){
  TC_samples[,i] = TC_samples[,i]* totalReads_scale[i]
  TA_samples[,i] = TA_samples[,i] * totalReads_scale[i]
  TC_all[,i] = TC_all[,i] * totalReads_scale[i]
  totalReads[,i] = totalReads[,i] * totalReads_scale[i]
  TAmultiple_samples[,i]= TAmultiple_samples[,i] * totalReads_scale[i]
}




TC_samples_untreated = TC_samples %>% dplyr::select(c(Untreated_TP1,Untreated_TP6,Untreated_TP9))
TA_samples_untreated = TA_samples %>% dplyr::select(c(Untreated_TP1,Untreated_TP6,Untreated_TP9))
TC_all_untreated = TC_all  %>% dplyr::select(c(Untreated_TP1,Untreated_TP6,Untreated_TP9))
totalReads_untreated = totalReads %>% dplyr::select(c(Untreated_TP1,Untreated_TP6,Untreated_TP9))
TA_multiple_untreated = TAmultiple_samples %>% dplyr::select(c(Untreated_TP1,Untreated_TP6,Untreated_TP9))
conversionRates_untreated = conversionRates_all_injection %>% dplyr::select(c(Untreated_TP1,Untreated_TP6,Untreated_TP9))


TC_all_untreated = TC_all_untreated[,c(1,3)]
TC_samples_untreated= TC_samples_untreated[,c(1,3)]
totalReads_untreated = totalReads_untreated[,c(1,3)]
conversionRates_untreated = conversionRates_untreated[,c(1:3)]
perTimePointReads = vector('list',2)
RPMs_allCws = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure2/data/RPM_allCws.txt",
                         sep="\t",stringsAsFactors = F, header = T)
RPMs_allCws = RPMs_allCws %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Untreated_TP9), sum, na.rm = TRUE)
RPMs_Inj = splitReplicates(RPMs_allCws,condition = "Inj",metadata_add = RPMs_allCws$name)[[3]]
RPMS_untreated  = RPMs_allCws[,c('Untreated_TP1','Untreated_TP9')]
  
pdf('/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure2/plots//correlateTA_total.pdf', height = 5, width = 3.8)


timepoints=c(1,9)

for(i in 1:length(timepoints)){
  
  perTimePointReads[[i]] = data.frame(TA=TA_samples[,timepoints[i]],
                                      allReads= RPMs_Inj[,timepoints[i]],TC=TC_all[,timepoints[i]]/0.82,TC_multiple= TC_samples[,timepoints[i]]/0.57, 
                                      TA_multiple = TAmultiple_samples[,timepoints[i]], TC_untreated = TC_all_untreated[,i], TC_samples_untreated = TC_samples_untreated[,i],
                                      totalReads_untreated = RPMS_untreated[,i],conversionRate = conversionRates_all_injection[,timepoints[i]],conversionRate_untreated = conversionRates_untreated[,i])
  colnames(perTimePointReads[[i]]) = c('TA','totalReads','TC','TC_multiple','TA_multiple','TC_untreated','TC_multiple_untreated','totalReads_Untreated','conversionRate','conversionRates_untreated')
  rownames(perTimePointReads[[i]]) = names_genes

  perTimePointReads[[i]] = (perTimePointReads[[i]][perTimePointReads[[i]]$totalReads>1 & perTimePointReads[[i]]$totalReads_Untreated>1,])
  
  
    perTimePointReads[[i]] = perTimePointReads[[i]] %>% dplyr::mutate(TC= log10(TC+1),TA= log10(TA+1), totalReads = log10(totalReads+1),
                                                                    TC_multiple = log10(TC_multiple+1), TA_multiple = log10(TA_multiple+1),
                                                                    TC_untreated = log10(TC_untreated+1),TC_multiple_untreated = log10(TC_multiple_untreated+1),
                                                                    totalReads_Untreated = log10(totalReads_Untreated+1))

  write.table(perTimePointReads[[i]],paste0("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure2/data/fig2a",timepoints[i],".txt"),
              sep="\t",quote = F,row.names = T)
  TA_multiple = perTimePointReads[[i]]
  TC_multiple = perTimePointReads[[i]] 
  TC_multiple_untreated = perTimePointReads[[i]]

  TA_single =  perTimePointReads[[i]] 
  TC_single =  perTimePointReads[[i]]
  TC_single_untreated = perTimePointReads[[i]] 
  
  nTA_multiple  = nrow(TA_multiple %>% dplyr::filter(TA_multiple>0))
  nTC_multiple = nrow(TC_multiple %>% dplyr::filter(TC_multiple>0) )
  nTC_untreated_multiple = nrow(TC_multiple_untreated %>% dplyr::filter(TC_multiple_untreated>0))
  
  nTA_single  = nrow(TA_single %>% dplyr::filter(TA>0))
  nTC_single = nrow(TC_single  %>% dplyr::filter(TC>0) )
  nTC_untreated_single = nrow(TC_single_untreated %>% dplyr::filter(TC_untreated>0) ) 
  
  cor_TC_multiple = round(cor(TC_multiple$TC_multiple,TC_multiple$totalReads),2)
  cor_TA_multiple = round(cor(TA_multiple$TA_multiple,TA_multiple$totalReads),2)
  cor_TC_untreated_multiple = round(cor(TC_multiple_untreated$TC_multiple_untreated,TC_multiple_untreated$totalReads_Untreated),2)
  
  cor_TC_single = round(cor(TC_single$TC,TC_single$totalReads),2)
  cor_TA_single = round(cor(TA_single$TA,TA_single$totalReads),2)
  cor_TC_untreated_single = round(cor(TC_single_untreated$TC_untreated,TC_single_untreated$totalReads_Untreated),2)
  
 
  TA_TC_density = ggplot(TA_multiple , aes(x =totalReads , y = TA_multiple)) +xlab('log10(total reads)') + 
    ylab('log10(Conversion containing reads)') +geom_density_2d(color='black') + 
    geom_density_2d(data = TC_multiple,aes(x=totalReads,y=TC_multiple),color='red') + 
    xlim(c(-0.5,4))  + ylim(c(-0.5,4))  + theme_cowplot() + annotate('text',label = paste0(timepoints[i] ,"\n", "TA=",nTA_multiple,"\n","TC=",nTC_multiple, "\n", "multiple conversions"),x=1,y=2) 
  TA_TC_density = TA_TC_density + theme_ameres(type = 'barplot') + annotate("text",label=paste0( "corTA=",cor_TA_multiple,"\n","corTC=",cor_TC_multiple, "\n"),x=2,y=3) + theme_ameres(type = 'barplot')
  #print(TA_TC_density)
  
  
  
  TCbg_TC_density = ggplot(TC_multiple_untreated , aes(x =totalReads_Untreated , y = TC_multiple_untreated)) +xlab('log10(total reads)') + 
    ylab('log10(Conversion containing reads)') +geom_density_2d(color='black') + 
    geom_density_2d(data = TC_multiple,aes(x=totalReads,y=TC_multiple),color='red') + 
    xlim(c(-0.5,4))  + ylim(c(-0.5,4))  + theme_cowplot() + annotate('text',label = paste0(timepoints[i]-1 ,"\n","TC=",nTC_multiple, "\n","TC untreated=",nTC_untreated_multiple, "\n", "multiple conversions"),x=1,y=3) 
  TCbg_TC_density = TCbg_TC_density + theme_ameres(type = 'barplot') + annotate("text",label=paste0( "corTC=",cor_TC_multiple, "\n","corTC untreated",cor_TC_untreated_multiple),x=3,y=3) + theme_ameres(type = 'barplot')
  
  #print(TCbg_TC_density)
  
  
  

    RPMs_multiple =  rbind.data.frame(TA_multiple %>% dplyr::select(totalReads) %>% dplyr::mutate(type = 'TA'), TC_multiple_untreated %>% 
                                              dplyr::select(totalReads)  %>% dplyr::mutate(type = 'TC background'),
                   TC_multiple %>% dplyr::select(totalReads) %>% dplyr::mutate(type = 'TC conversions') ) 
  
    totalReads_TC_bg = ggplot(RPMs_multiple %>% dplyr::filter(type != 'TA'),aes(type,totalReads)) + geom_boxplot( )+ theme_cowplot() 
    totalReads_TC_bg =  totalReads_TC_bg + theme_ameres(type = 'barplot') + 
      annotate(geom='text',label=paste0('Multiple',timepoints[i]),x=1,y=3) + ylim(c(-0.5,4))  +  ggpubr::stat_compare_means(label.y = c(3))  + coord_flip()
    #print(totalReads_TC_bg)
    
    totalReads_TC_TA = ggplot(RPMs_multiple %>% dplyr::filter(type != 'TC background'),aes(type,totalReads)) + geom_boxplot( )+ theme_cowplot() 
    totalReads_TC_TA =  totalReads_TC_TA + theme_ameres(type = 'barplot') + annotate(geom='text',label=paste0('Multiple',timepoints[i]),x=1,y=3) + ylim(c(-0.5,4))  + 
      ggpubr::stat_compare_means(label.y = c(3))  + coord_flip()
    #print(totalReads_TC_TA)
 
    
     
    multi_TA = data.frame(conversion = TA_multiple  %>% dplyr::select(TA_multiple) %>% dplyr::pull(1),type = 'TA')
    multi_TC = data.frame(conversion = TC_multiple %>% dplyr::select(TC_multiple) %>% dplyr::pull(1),type = 'TC conversions')
    multi_TC_untreated = data.frame(conversion = TC_multiple_untreated %>% dplyr::select(TC_multiple_untreated) %>% dplyr::pull(1), type = 'TC background')
    RPMs_multiple_conversion = rbind.data.frame(multi_TA,multi_TC_untreated,multi_TC)  

    
    TCbg_TC_multiple = ggplot(RPMs_multiple_conversion %>% dplyr::filter(type !="TA") ,aes(type, (conversion))) +
      geom_boxplot( )+  ggpubr::stat_compare_means(label.y = c(3,3.1,3.2)) +
      theme_cowplot() + ylab ('Multiple conversions') + ylim(c(-0.5,4))+ coord_flip()
    TCbg_TC_multiple = TCbg_TC_multiple + theme_ameres(type = 'barplot')+ annotate(geom = 'text',label=timepoints[i],x=1,y=4)
    #print(TCbg_TC_multiple)
    
 
    TA_TC_multiple = ggplot(RPMs_multiple_conversion %>% dplyr::filter(type !="TC background") ,aes(type, (conversion))) +
      geom_boxplot( )+  ggpubr::stat_compare_means(label.y = c(3,3.1,3.2)) +
      theme_cowplot() + ylab ('Multiple conversions') + ylim(c(-0.5,4))+ coord_flip()
    TA_TC_multiple = TA_TC_multiple + theme_ameres(type = 'barplot')+ annotate(geom = 'text',label=timepoints[i],x=1,y=4)
    #print(TA_TC_multiple)
      library(patchwork)
    # combined_TA_TC  =   TA_TC_multiple / totalReads_TC_TA
    # combined_TA_TC = TA_TC_density / combined_TA_TC 
    # print(combined_TA_TC)
    # 
    # combined_bg_TC  = TCbg_TC_multiple / totalReads_TC_bg
    # combined_bg_TC = TCbg_TC_density / combined_bg_TC 
    # print(combined_bg_TC)

################################### single conversions... 
        
    
    
    TA_TC_density = ggplot(TA_single , aes(x =totalReads , y = TA)) +xlab('log10(total reads)') + 
      ylab('log10(Conversion containing reads)') +geom_density_2d(color='black') + 
      geom_density_2d(data = TC_single,aes(x=totalReads,y=TC),color='red') + 
      xlim(c(-0.5,4))  + ylim(c(-0.5,4))  + theme_cowplot() + annotate('text',label = paste0(timepoints[i]-1 ,"\n","all=",nrow(perTimePointReads[[i]]), "\n","TA=",nTA_single,"\n","TC=",nTC_single, "\n", "single conversions"),x=1,y=3) 
    TA_TC_density = TA_TC_density + theme_ameres(type = 'barplot') + annotate("text",label=paste0( "corTA=",cor_TA_single,"\n","corTC=",cor_TC_single, "\n"),x=3,y=3) + theme_ameres(type = 'barplot')
    #print(TA_TC_density)
    

    
    TCbg_TC_density = ggplot(TC_single_untreated , aes(x =totalReads , y = TC_untreated)) +xlab('log10(total reads)') + 
      ylab('log10(Conversion containing reads)') +geom_density_2d(color='black') + 
      geom_density_2d(data = TC_single,aes(x=totalReads,y=TC),color='red') + 
      xlim(c(-0.5,4))  + ylim(c(-0.5,4))  + theme_cowplot() + annotate('text',label = paste0(timepoints[i]-1 ,"\n","all=",nrow(perTimePointReads[[i]]),"\n","TC=",nTC_single, "\n","TC untreated=",nTC_untreated_single, "\n", "single conversions"),x=1,y=3) 
    TCbg_TC_density = TCbg_TC_density + theme_ameres(type = 'barplot') + annotate("text",label=paste0( "corTC=",cor_TC_single, "\n","corTC untreated",cor_TC_untreated_single),x=3,y=3) + theme_ameres(type = 'barplot')
    #print(TCbg_TC_density)
    a = TA_single %>% dplyr::select(TC_untreated,TA)
    a= reshape2::melt(a)
    ggplot(a,aes(x=variable,y=value)) + geom_boxplot()
    # ggplot(TA_single , aes(x =totalReads , y = TC_untreated)) +xlab('log10(total reads)') + 
    #   ylab('log10(Conversion containing reads)') +geom_density_2d(color='black') + 
    #   geom_density_2d(data = TA_single,aes(x=totalReads,y=TA),color='red') + 
    #   xlim(c(-0.5,4))  + ylim(c(-0.5,4))  + theme_cowplot() + annotate('text',label = paste0(timepoints[i]-1 ,"\n","all=",nrow(perTimePointReads[[i]]),"\n","TC=",nTC_single, "\n","TC untreated=",nTC_untreated_single, "\n", "single conversions"),x=1,y=3) 
    
    TCbg_TC_density = TCbg_TC_density + theme_ameres(type = 'barplot') + annotate("text",label=paste0( "corTC=",cor_TC_single, "\n","corTC untreated",cor_TC_untreated_single),x=3,y=3) + theme_ameres(type = 'barplot')
    
    
    
    RPMs_single=  rbind.data.frame(TA_single %>% dplyr::select(totalReads) %>% dplyr::mutate(type = 'TA'),
                                   TC_single_untreated = TC_single_untreated %>%
                                     dplyr::select(totalReads_Untreated) %>% dplyr::mutate(totalReads=totalReads_Untreated,type = 'TC background') %>% 
                                     dplyr::select(totalReads,type), TC_single %>% dplyr::select(totalReads) %>% dplyr::mutate(type = 'TC conversions'))
    
    
    totalReads_TC_bg = ggplot(RPMs_single %>% dplyr::filter(type != 'TA'),aes(type,totalReads)) + geom_boxplot(outlier.shape = NA )+ theme_cowplot() 
    totalReads_TC_bg =  totalReads_TC_bg + theme_ameres(type = 'barplot') + 
      annotate(geom='text',label=paste0('Single',timepoints[i]),x=1,y=3) + ylim(c(-0.5,4))  +  ggpubr::stat_compare_means(label.y = c(3.5))  + coord_flip()
    #print(totalReads_TC_bg)
    
    
    totalReads_TC_TA = ggplot(RPMs_single %>% dplyr::filter(type != 'TC background'),aes(type,totalReads)) + geom_boxplot(outlier.shape = NA )+ theme_cowplot() 
    totalReads_TC_TA =  totalReads_TC_TA + theme_ameres(type = 'barplot') + annotate(geom='text',label=paste0('Single',timepoints[i]),x=1,y=3) + ylim(c(-0.5,4))  + 
      ggpubr::stat_compare_means(label.y = c(3.5))  + coord_flip()
    #print(totalReads_TC_TA)
    
    
    single_TA = data.frame(conversion = TA_single %>% dplyr::select(TA) %>% dplyr::pull(1),type = 'TA')
    single_TC = data.frame(conversion = TC_single %>% dplyr::select(TC) %>% dplyr::pull(1),type = 'TC conversions')
    single_TC_untreated = data.frame(conversion = TC_single_untreated %>% dplyr::select(TC_untreated) %>% dplyr::pull(1),type = 'TC background')
    
    RPMs_single_conversion = rbind.data.frame(single_TA,single_TC_untreated,single_TC)  
    
    
    TCbg_TC_single = ggplot(RPMs_single_conversion %>% dplyr::filter(type !="TA") ,aes(type, (conversion))) +
      geom_boxplot(outlier.shape = NA )+  ggpubr::stat_compare_means(label.y = c(3,3.1,3.2)) +
      theme_cowplot() + ylab ('Single conversions') + ylim(c(-0.5,4))+ coord_flip()
    TCbg_TC_single = TCbg_TC_single + theme_ameres(type = 'barplot')+ annotate(geom = 'text',label=timepoints[i],x=1,y=3)
    #print(TCbg_TC_single)
    
    
    
    TA_TC_single = ggplot(RPMs_single_conversion %>% dplyr::filter(type !="TC background") ,aes(type, (conversion))) +
      geom_boxplot(outlier.shape = NA )+  ggpubr::stat_compare_means(label.y = c(3,3.1,3.2)) +
      theme_cowplot() + ylab ('Single conversions') + ylim(c(-0.5,4))+ coord_flip()
    TA_TC_single = TA_TC_single + theme_ameres(type = 'barplot')+ annotate(geom = 'text',label=timepoints[i],x=1,y=3)
    #print(TA_TC_single)
    
    
    combined_TA_TC  =   TA_TC_single / totalReads_TC_TA
    combined_TA_TC = TA_TC_density / combined_TA_TC 
    print(combined_TA_TC) 
    
    combined_bg_TC  = TCbg_TC_single / totalReads_TC_bg
    combined_bg_TC = TCbg_TC_density / combined_bg_TC 
    print(combined_bg_TC)
    
      TC_pie = data.frame( val =c(nTC_single,nrow(perTimePointReads[[i]])-nTC_single),sample = c('TC_single','Other'))
      TC_untreated_pie = data.frame( val =c(nTC_untreated_single,nrow(perTimePointReads[[i]])-nTC_untreated_single),sample = c('TC_untreated','Other'))
      TA_pie = data.frame( val =c(nTA_single,nrow(perTimePointReads[[i]])-nTA_single),sample = c('TA','Other'))
    
      
    pie_TC =  ggplot(TC_pie,aes(x="",y=val,fill=sample)) + geom_bar(width = 1, stat = "identity")+ coord_polar("y", start=0)+
       geom_text(aes(label = round(val*100/sum(val,4))), size=5)+ theme(axis.text.x=element_blank()) + theme_cowplot()
     pie_TC = pie_TC + theme_ameres(type = 'barplot')
     
     pie_untreated_TC =  ggplot(TC_untreated_pie,aes(x="",y=val,fill=sample)) + geom_bar(width = 1, stat = "identity")+ coord_polar("y", start=0)+
       geom_text(aes(label = round(val*100/sum(val,4))), size=5)+ theme(axis.text.x=element_blank()) + theme_cowplot()
     pie_untreated_TC = pie_untreated_TC + theme_ameres(type = 'barplot')
    
       pie_TA =  ggplot(TA_pie,aes(x="",y=val,fill=sample)) + geom_bar(width = 1, stat = "identity")+ coord_polar("y", start=0)+
       geom_text(aes(label = round(val*100/sum(val,4))), size=5) + theme(axis.text.x=element_blank()) + theme_cowplot()
      pie_TA = pie_TA + theme_ameres(type = 'barplot')
   
     MultiplePies = pie_TC / pie_untreated_TC / pie_TA   
      print(MultiplePies)
      
      }


dev.off()

plot(density (perTimePointReads[[2]]  %>% dplyr::filter(TA>0) %>% dplyr::pull(2)),col='red')
lines(density(perTimePointReads[[1]]  %>% dplyr::filter(TA>0) %>% dplyr::pull(2)))

nrow(perTimePointReads[[1]])
