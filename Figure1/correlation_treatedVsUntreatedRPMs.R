###### fold change of genes... 

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
###################   
# RPM = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_may2019/dataTables_expandedCounting/RPM_allCws.txt",
#                  sep="\t", stringsAsFactors = F, header = T)

RPM = read.table("/Volumes/Macintosh HD/users/pooja.bhat/Downloads//RPM_allCws.txt",
                 sep="\t", stringsAsFactors = F, header = T)
RPM = RPM %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Untreated_TP9), sum, na.rm = TRUE)
RPM_zygoticGenes_inj = splitReplicates(dataFrameToSplit = RPM,condition = "Inj",metadata_add = RPM$name)[[3]]
RPM_zygoticGenes_Untreated = RPM[,grep('Untre',colnames(RPM))]
RPM_total = cbind.data.frame(RPM_zygoticGenes_inj,RPM_zygoticGenes_Untreated)
# RPM_zygoticGenes = RPM_total[RPM_total$metadata_add %in% zygoticGenes$external_gene_name,] %>% dplyr::mutate(Injected_RPM = log10(Inj_R2_TP9+1), untreatedTPM = log10(Untreated_TP9+1))
# RPM_maternalGenes = RPM_total[RPM_total$metadata_add %in% maternalGenes$external_gene_name,] %>% dplyr::mutate(Injected_RPM = log10(Inj_R2_TP1+1), untreatedTPM = log10(Untreated_TP1+1))
# 
# ggpubr::ggscatter(RPM_zygoticGenes,x='Injected_RPM',y='untreatedTPM')

##### just plotting the correlation between treated and untreated samples...

 
tmp_df = RPM_total %>% dplyr::filter(Inj_R2_TP9 > 1 & Untreated_TP9 > 1)

tmp_df = tmp_df %>% dplyr::mutate(log10TP9_inj = log10(Inj_R2_TP9),log10TP9_unt = log10(Untreated_TP9)) %>%
  mutate(densityCols_TP = densCols(log10TP9_inj, log10TP9_unt,  colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))) 
pdf("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/plots//correlation_RPMsTCreads_treatedVsUntreated.pdf")

          p = ggplot2::ggplot(tmp_df) + ggrastr::geom_point_rast(aes(x=log10TP9_inj,y=log10TP9_unt),col=tmp_df$densityCols_TP,alpha=0.3) + theme_cowplot() +xlim(c(0,4))+ylim(c(0,4))
          p = p  + theme_ameres(type = 'barplot') + ggpubr::stat_cor(aes(x=log10TP9_inj,y=log10TP9_unt)) + xlab("log10(+4SU)")  + ylab("log10(-4SU)")
          p = p + ggtitle(paste("TP9","n=",nrow(tmp_df)))
          print(p)
          
          
          tmp_df = RPM_total %>% dplyr::filter(Inj_R2_TP1 > 1 & Untreated_TP1 > 1)
          
          tmp_df = tmp_df %>% dplyr::mutate(log10TP1_inj = log10(Inj_R2_TP1),log10TP1_unt = log10(Untreated_TP1)) %>%
            mutate(densityCols_TP = densCols(log10TP1_inj, log10TP1_unt,  colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))) 
          
          
          p = ggplot2::ggplot(tmp_df) + ggrastr::geom_point_rast(aes(x=log10TP1_inj,y=log10TP1_unt),col=tmp_df$densityCols_TP,alpha=0.3) + theme_cowplot() +xlim(c(0,4))+ylim(c(0,4))
          p = p  + theme_ameres(type = 'barplot') + ggpubr::stat_cor(aes(x=log10TP1_inj,y=log10TP1_unt)) + xlab("log10(+4SU)")  + ylab("log10(-4SU)")
          p = p + ggtitle(paste("TP1","n=",nrow(tmp_df)))
          print(p)
          
          
          ### TP6 
          
          tmp_df = RPM_total %>% dplyr::filter(Inj_R2_TP6 >1 & Untreated_TP6 >1)
          
          tmp_df = tmp_df %>% dplyr::mutate(log10TP6_inj = log10(Inj_R2_TP6),log10TP6_unt = log10(Untreated_TP6)) %>%
            mutate(densityCols_TP = densCols(log10TP6_inj, log10TP6_unt,  colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))) 
          
          
          p = ggplot2::ggplot(tmp_df) + ggrastr::geom_point_rast(aes(x=log10TP6_inj,y=log10TP6_unt),col=tmp_df$densityCols_TP, alpha=0.3) + theme_cowplot() +xlim(c(0,4))+ylim(c(0,4))
          p = p  + theme_ameres(type = 'barplot') + ggpubr::stat_cor(aes(x=log10TP6_inj,y=log10TP6_unt),) + xlab("log10(+4SU)")  + ylab("log10(-4SU)")
          p = p + ggtitle(paste("TP6","n=",nrow(tmp_df)))
          print(p)
          
          
         

dev.off()
