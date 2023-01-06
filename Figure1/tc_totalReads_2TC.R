###### script plotting total reads vs tc reads to show that for a set of pre-defined purely zygotic genes there most reads have tc
### requirements ####

### input 
# 1) number of TC reads
# 2) Number of total reads
# 3) error rates 

library(dplyr)
library(ggpubr)
library(ggplot2)
library(reshape)
library(patchwork)
library(ggExtra)
library(ggjoy)
library(gridExtra)
library(ggbeeswarm)
library(ggrepel)
'%!in%' <- function(x,y)!('%in%'(x,y))

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
fun_mean <- function(x){
  return(data.frame(y=median(x),label=round(median(x,na.rm=T),2)) )}
plotFractionTC = function(dataF){
  dataF = dataF[order(dataF$type),]
  p =   ggpubr::ggboxplot(dataF ,x='type',y='fractionTC',error.plot = 'errorbar', outlier.shape=NA,palette = c('grey','grey30','red'),
                          fill = 'type',bxp.errorbar=T) + theme_ameres(type = 'barplot') + ylab ('Fraction of TC reads') + 
    stat_summary(fun.data = fun_mean, geom="text", vjust=-0.7)
  
  return(p)
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
plotdensity = function(dataF){
  dataF_split = split(dataF,dataF$type,T)
  dataF_split = lapply(dataF_split,function(x) x[x$totalReads>10,] )
  
  n_dataSplit = reshape::melt(lapply(dataF_split,function(x) nrow(x))) %>% dplyr::mutate(value,L1) %>% dplyr::mutate(n= paste(value,L1))
  n_dataDetaile = paste(n_dataSplit$n[1],"\n",n_dataSplit$n[2],"\n",n_dataSplit$n[3])
  
  cor_dataSplit = reshape::melt(lapply(dataF_split,function(x) cor(x$tcReads,x$totalReads))) %>% dplyr::mutate(cor= paste(round(value,2),L1))
  cor_details = paste(cor_dataSplit$cor[1],"\n",cor_dataSplit$cor[2],"\n",cor_dataSplit$cor[3])
  
  p = ggplot()  + 
    geom_density_2d(data = dataF %>% dplyr::filter(type == "InBetween"),aes(x=log10(tcReads+1),y=log10(totalReads+1)),col='grey',alpha=0.5)+ 
    geom_density_2d(data = dataF %>% dplyr::filter(type == "Maternal"),aes(x=log10(tcReads+1),y=log10(totalReads+1)),col='black',alpha=0.5) +
    geom_density_2d(data = dataF %>% dplyr::filter(type == "Zygotic"),aes(x=log10(tcReads+1),y=log10(totalReads+1)),col='red',alpha=0.5)
  p = p + theme_ameres(type = 'barplot') + xlim(c(0,4))+ ylim(c(0,4)) + xlab('TC reads (log10)') + ylab('Total reads (log10)')
  p= p+ggplot2::annotate(geom = 'text', label =cor_details , x = 3, y = 3) +ggplot2::annotate(geom = 'text', label =n_dataDetaile , x = 3, y = 1) 
  return(p)
  
}
plotTC_total = function(dataF,timeP){
  pattern_use = paste0("_",timeP)
  ### also plotting TC reads v.s total reads... 
  dataF  = dataF %>% dplyr::select(c(dplyr::matches(pattern_use),dplyr::contains('Unt'),'external_gene_name')) %>% 
    dplyr::mutate(logTC = log10(.[[2]]+1) ,logTotal = log10(.[[3]]+1) ) %>%
    mutate(densityCols_TP = densCols(logTotal, logTC,  colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))) 
  contour =   ggplot(dataF,aes(x=logTotal,y=logTC)) + geom_density_2d_filled(contour_var="count") + xlim(c(0,3)) + ylim(c(0,3)) + ggtitle(timeP) + theme_cowplot()
  contour = contour + theme_ameres(type = 'barplot') 
  
  dataF = dataF[order(dataF[,2],decreasing = T),]
  cor_data = round(cor(dataF$totalReads_9,dataF$TC_9),2)
  n_data= nrow(dataF)
  pVal_cor = cor.test(dataF$totalReads_9,dataF$TC_9)$p.value
  
  filledDensity    = ggplot(dataF,aes(x=logTotal,y=logTC,label=external_gene_name)) + stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
    scale_fill_distiller(palette='Greys', direction=1)  +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
      legend.position='none'
    )+ ylim(c(0,3))  + xlim(c(0,3)) + theme_cowplot()
  filledDensity = filledDensity + ggtitle(timeP) + theme_ameres(type = 'barplot') + 
    geom_density_2d(color='grey') + 
    geom_point(data = dataF,aes(x=logTotal,y=logTC),alpha=0.2,color='red',size=0.5) + ggrepel::geom_text_repel(data = dataF[1:10,])
  filledDensity =  filledDensity + ggplot2::annotate(geom = 'text', label =paste(cor_data,pVal_cor,sep = ",") , x = 2, y = 2.7) +
    ggplot2::annotate(geom = 'text', label =n_data , x = 2, y = 2.5)  + ylab('Number of TC reads (log10)') + xlab('Number of total reads (log10)')
  
  pointDensity = ggplot(dataF,aes(x=logTotal,y=logTC)) + geom_point(col=dataF$densityCols_TP)+ xlim(c(0,3)) + ylim(c(0,3))+ ggtitle(timeP) + theme_cowplot()
  pointDensity = pointDensity + theme_ameres(type = 'barplot')
  
  returnpPlots = list(contour,filledDensity,pointDensity)
  
  return(returnpPlots)
}
plotTC_dotPlot = function(dataF,timeP){
  pattern_use = paste0("_",timeP)
  ### also plotting TC reads v.s total reads... 
  dataF  = dataF %>% dplyr::select(c(dplyr::matches(pattern_use),dplyr::contains('Unt'))) %>% 
    dplyr::mutate(logTC = log10(.[[2]]+1) ,logTotal = log10(.[[3]]+1) ) %>%
    mutate(densityCols_TP = densCols(logTotal, logTC,  colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))) 
  
  pointDensity = ggplot(dataF,aes(x=logTotal,y=logTC)) + geom_point(col=dataF$densityCols_TP)+ xlim(c(0,3)) + ylim(c(0,3))+ ggtitle(timeP) + theme_cowplot()
  pointDensity = pointDensity + theme_ameres(type = 'barplot')
  
  returnpPlots = list(pointDensity)
  
  return(returnpPlots)
}
plotConverionRates_TC = function(dataF,addName){
  getConversionRates= vector('list',9)
  TCreads= vector('list',9)
  fractionTC = TCreads
  fractionALLTC = TCreads
  RPM_genes = vector('list',9)
  geneName = vector('list',9)
  tc_with_names = vector('list',9)
  names(getConversionRates) = c(0.75,2,2.5,3,3.5,4,4.5,5,5.5)
  names(TCreads) = c(0.75,2,2.5,3,3.5,4,4.5,5,5.5)
  names(fractionTC) = c(0.75,2,2.5,3,3.5,4,4.5,5,5.5)
  names(RPM_genes) =  c(0.75,2,2.5,3,3.5,4,4.5,5,5.5)
  names(fractionALLTC) = c(0.75,2,2.5,3,3.5,4,4.5,5,5.5)
  names(geneName) = c(0.75,2,2.5,3,3.5,4,4.5,5,5.5)
  names(tc_with_names) = c(0.75,2,2.5,3,3.5,4,4.5,5,5.5)
  for(i in 1:9){
    pattern_use = paste0("_",i)
    RPM_get = paste0('Inj_mean_TP',i)
    rates_rep1 = paste0('Inj_R2_TP',i)
    rates_rep2 = paste0('Inj_R3_TP',i)
    
    getConversionRates[[i]]  = dataF %>% 
      dplyr::select(c(dplyr::matches(pattern_use),dplyr::contains('Unt'),dplyr::matches(rates_rep1),dplyr::matches(rates_rep2)),
                    dplyr::matches(RPM_get)) %>% dplyr::filter(.[[8]]>0 & .[[9]]>0 ) %>% dplyr::pull(1)
    geneName[[i]] = dataF %>% 
      dplyr::select(c(dplyr::matches(pattern_use),dplyr::contains('Unt'),dplyr::matches(rates_rep1),dplyr::matches(rates_rep2)),
                    dplyr::matches(RPM_get),external_gene_name) %>% dplyr::filter(.[[8]]>0 & .[[9]]>0 ) %>% dplyr::select('external_gene_name')
    TCreads[[i]]  = dataF %>% dplyr::select(c(dplyr::matches(pattern_use),dplyr::contains('Unt')),dplyr::matches(RPM_get))  %>% dplyr::pull(2)
    fractionTC[[i]]= dataF %>% dplyr::select(c(dplyr::matches(pattern_use),dplyr::contains('Unt')),dplyr::matches(RPM_get))  %>% dplyr::mutate(fractionTC = .[[2]]/.[[3]]) %>% dplyr::pull('fractionTC')
    fractionALLTC[[i]]= dataF %>% dplyr::select(c(dplyr::matches(pattern_use),dplyr::contains('Unt')),dplyr::matches(RPM_get))  %>% dplyr::mutate(fractionTC = .[[4]]/.[[3]]) %>% dplyr::pull('fractionTC')
    RPM_genes[[i]] = dataF %>% dplyr::select(c(dplyr::matches(pattern_use),dplyr::contains('Unt')),dplyr::matches(RPM_get))  %>% dplyr::pull(8)
    tc_with_names[[i]] = cbind.data.frame(geneName[[i]],getConversionRates[[i]])
  }
  
  allRPMs =  do.call(cbind.data.frame,RPM_genes)
  allRPMs = as.data.frame(t(apply(allRPMs,1,function(x) x/max(x))))
  allRPMs$external_gene_name = dataF$external_gene_name
  allRPMs = melt(allRPMs)
  getConversionRates = reshape2::melt(getConversionRates)
  TCreads = reshape2::melt(TCreads)
  fractionTC = reshape2::melt(fractionTC)
  fractionAllTC = reshape2::melt(fractionALLTC)
  RPM_genes = reshape2::melt(RPM_genes)
  
  addThis_complete = data.frame(L1=c(0.75,2,2.5,3),value = c(0,0,0,0))
  getConversionRates = rbind.data.frame(getConversionRates,addThis_complete)
  TCreads = rbind.data.frame(TCreads,addThis_complete)
  fractionTC = rbind.data.frame(fractionTC,addThis_complete)
  fractionAllTC = rbind.data.frame(fractionAllTC,addThis_complete)
  RPM_genes = rbind.data.frame(RPM_genes,addThis_complete)
  
  conversionRate_plot = ggplot(getConversionRates,aes(x=L1,y=(value),group=L1)) + geom_violin(fill = 'red',scale = 'width',trim=T)  + theme_cowplot()+
    ylab("Conversion rate") + ggtitle(paste0(addName))+ xlab("Timepoint") 
  conversionRate_plot = conversionRate_plot  + theme_ameres(type = "barplot") + geom_boxplot(width=0.1, outlier.shape=NA)+ ylim(c(0,0.15))
  
  conversionRate_plot
  
  ### plotting a quasirandome plot
  df.summary <- getConversionRates %>%
    group_by(L1) %>%
    summarise(
      sd = sd(value, na.rm = TRUE),
      len = mean(value,na.rm=T)
    )
  df.summary[is.na(df.summary)]<-0
  df.summary = df.summary %>% dplyr::mutate(top = len+sd, bottom = len-sd)
  getConversionRates %>%  dplyr::group_by(L1) %>% summarise(n=n(), avg = mean(value,na.rm=T)) ->Summary.data
  
  quasirandom_plot = ggplot(getConversionRates,aes(L1, value)) +
    geom_quasirandom(alpha=.2,varwidth = T) + geom_errorbar(data = df.summary,aes(x = L1,y=len,ymin=len-sd,ymax=len+sd,width=0.1),col='red')
  quasirandom_plot = quasirandom_plot + theme_cowplot() + xlab('Time (hpf)') + ylab ('Conversion rate')
  quasirandom_plot = quasirandom_plot + theme_ameres(type = 'barplot')+ ylim(c(-0.01,0.15))  +
    geom_line(data = df.summary,aes(x=L1,y=len-sd,group=1),linetype = "dashed")  +
    geom_line(data = df.summary,aes(x=L1,y=len+sd,group=1),linetype = "dashed") +
    geom_text(data=Summary.data ,aes(x = L1, y = avg, label=n),color="red", fontface =1, size = 3) + ggtitle(addName)
 
  # quasirandom_plot = ggplot(getConversionRates,aes(L1, value)) +
  #   geom_quasirandom(alpha=.2,varwidth = T) + geom_errorbar(data = df.summary,aes(x = L1,y=len,ymin=len-sd,ymax=len+sd,width=0.1),col='red')
  # quasirandom_plot = quasirandom_plot + theme_cowplot() + xlab('Time (hpf)') + ylab ('Conversion rate')
  # quasirandom_plot = quasirandom_plot + theme_ameres(type = 'barplot')+ ylim(c(-0.01,0.15))  +   
  #   geom_text(data=Summary.data ,aes(x = L1, y = avg, label=n),color="red", fontface =1, size = 3) + ggtitle(addName)
  
  
   linePlot_conversions = ggplot(getConversionRates,aes(L1, value)) +
    geom_errorbar(data = df.summary,aes(x = L1,y=len,ymin=len-sd,ymax=len+sd,width=0.1),col='red')
   linePlot_conversions = linePlot_conversions + theme_cowplot() + xlab('Time (hpf)') + ylab ('Conversion rate')
   linePlot_conversions = linePlot_conversions + theme_ameres(type = 'barplot')+ ylim(c(-0.02,0.15))  + 
    geom_line(data = df.summary,aes(x=L1,y=len-sd,group=1),linetype = "dashed")  + 
    geom_line(data = df.summary,aes(x=L1,y=len+sd,group=1),linetype = "dashed") +   
    geom_text(data=Summary.data ,aes(x = L1, y = avg, label=n),color="red", fontface =1, size = 3) + ggtitle(addName)
  ############## plotting the RPMS.. 
  df.summary_RPMs <- allRPMs %>%
    group_by(variable) %>%
    summarise(
      sd = sd(value, na.rm = TRUE),
      len = mean(value,na.rm=T)
    )
  df.summary_RPMs[is.na(df.summary_RPMs)]<-0
  allRPMs %>%  dplyr::group_by(variable) %>% summarise(n=n(), avg = mean(value,na.rm=T)) ->Summary.data_RPMs
  quasirandom_plot_RPM = ggplot(allRPMs,aes(variable, value)) +
    geom_quasirandom(alpha=.2,varwidth = T) + geom_errorbar(data = df.summary_RPMs,aes(x = variable,y=len,ymin=len-sd,ymax=len+sd,width=0.1),col='red')
  quasirandom_plot_RPM = quasirandom_plot_RPM + theme_cowplot() + xlab('Time (hpf)') + ylab ('Normalized RPM')
  quasirandom_plot_RPM = quasirandom_plot_RPM + theme_ameres(type = 'barplot')  + 
    geom_line(data = df.summary_RPMs,aes(x=variable,y=len-sd,group=1),linetype = "dashed")  + 
    geom_line(data = df.summary_RPMs,aes(x=variable,y=len+sd,group=1),linetype = "dashed") + ggtitle(paste0(addName,"_",nrow(allRPMs)/9)) + ylim(c(-0.1,1.1))
  
  
line_plot_RPM = ggplot(allRPMs,aes(variable, value)) +
    geom_errorbar(data = df.summary_RPMs,aes(x = variable,y=len,ymin=len-sd,ymax=len+sd,width=0.1),col='red')+ ylim(c(-0.1,1.1))
  line_plot_RPM = line_plot_RPM + theme_cowplot() + xlab('Time (hpf)') + ylab ('Normalized RPM') 
  line_plot_RPM = line_plot_RPM + theme_ameres(type = 'barplot')  + 
    geom_line(data = df.summary_RPMs,aes(x=variable,y=len-sd,group=1),linetype = "dashed")  + 
    geom_line(data = df.summary_RPMs,aes(x=variable,y=len+sd,group=1),linetype = "dashed") + ggtitle(paste0(addName,"_",nrow(allRPMs)/9))

  FractionTC_plot = ggplot(fractionTC,aes(x=L1,y=value,group=L1)) + geom_violin(fill = 'red',scale = 'width',trim=T)  + theme_cowplot()+
    ylab("Fraction T>C containing reads") + ggtitle(paste0(addName))+ xlab("Timepoint")
  FractionTC_plot = FractionTC_plot  + theme_ameres(type = "barplot") + geom_boxplot(width=0.1, outlier.shape=NA) + ylim(c(0,1))
  FractionTC_plot
  
  
  tc_with_names = do.call(rbind.data.frame,tc_with_names)
  return_plotList = list(quasirandom_plot,quasirandom_plot_RPM,getConversionRates,line_plot_RPM,linePlot_conversions,fractionTC,fractionAllTC,tc_with_names)
  names(return_plotList) = c('quasiRandom_conversion','quasiRandom_RPM','data_conversionRates','linePlot_rpm','linePlot_conversions','fractionTC','fractionAllTC','tcPertimePoint')
  
  
  return(return_plotList)
}

totalReads = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure2//data/totalCounts_allCws.txt",
                        sep="\t", header = T, stringsAsFactors = F)

# totalReads = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures//figure2/data/totalCounts_allCws.txt",
#                         sep="\t", header = T, stringsAsFactors = F)

tcReads = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure2//data/numberOfreadsWithMultipleTC.txt",
                     sep="\t", header = T, stringsAsFactors = F)

# tcReads = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure2/data/numberOfreadsWithMultipleTC.txt",
#                      sep="\t", header = T, stringsAsFactors = F)

tcReads_3OrMore = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure2//data/numberOfreadsWith_3orMoreTC.txt",
                             sep="\t",stringsAsFactors = F, header = T)

# tcReads_3OrMore = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures//figure2/data/numberOfreadsWith_3orMoreTC.txt",
#                              sep="\t",stringsAsFactors = F, header = T)

ALLtcReads = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure2//data/numberOfreadsWithTC.txt",
                        sep="\t", header = T, stringsAsFactors = F)

# ALLtcReads = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures//figure2/data/numberOfreadsWithTC.txt",
#                      sep="\t", header = T, stringsAsFactors = F)

ALLTAReads = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure2//data/numberOfreadsWithTA.txt",
                        sep="\t", header = T, stringsAsFactors = F)

# ALLTAReads = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures//figure2/data/numberOfreadsWithTA.txt",
#                         sep="\t", header = T, stringsAsFactors = F)


totalReads =     totalReads[,-grep("Inc",colnames(totalReads))]
ALLtcReads = ALLtcReads[,-grep("Inc",colnames(ALLtcReads))]
ALLTAReads = ALLTAReads[,-grep("Inc",colnames(ALLTAReads))]

totalReads_Inj = totalReads %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Inj_R3_TP9), sum, na.rm = TRUE)
tcReads_Inj = tcReads %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Inj_R3_TP9), sum, na.rm = TRUE)
ALLtcReads_Inj = ALLtcReads %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Inj_R3_TP9), sum, na.rm = TRUE)
ALLtaReads_Inj =   ALLTAReads %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Inj_R3_TP9), sum, na.rm = TRUE)
tcReads3TC_Inj = tcReads_3OrMore %>% dplyr::group_by(Name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Inj_R3_TP9), sum, na.rm = TRUE)


totalReads_Inj_mean = splitReplicates(dataFrameToSplit = totalReads_Inj, condition = "Inj",metadata_add = totalReads_Inj$name)[[3]]
tcReads_Inj_mean = splitReplicates(dataFrameToSplit = tcReads_Inj, condition = "Inj",metadata_add = tcReads_Inj$name)[[3]]
ALLtcReads_Inj_mean = splitReplicates(dataFrameToSplit = ALLtcReads_Inj, condition = "Inj",metadata_add = ALLtcReads_Inj$name)[[3]]
ALLtaReads_Inj_mean = splitReplicates(dataFrameToSplit = ALLtaReads_Inj, condition = "Inj",metadata_add = ALLtaReads_Inj$name)[[3]]
tcReads3TC_Inj_mean = splitReplicates(dataFrameToSplit = tcReads3TC_Inj, condition = "Inj",metadata_add = tcReads3TC_Inj$Name)[[3]]


colnames(tcReads_Inj_mean) = c(paste0('TC_',c(1:9)), 'external_gene_name')
colnames(totalReads_Inj_mean) = c(paste0('totalReads_',c(1:9)), 'external_gene_name')
colnames(ALLtcReads_Inj_mean) = c(paste0('allTCReads_',c(1:9)), 'external_gene_name')
colnames(tcReads3TC_Inj_mean) = c(paste0('TCReads3_',c(1:9)), 'external_gene_name')

tcReads_Inj_mean$external_gene_name = NULL
ALLtcReads_Inj_mean$external_gene_name = NULL
tcReads_totalReads = cbind.data.frame(tcReads_Inj_mean, totalReads_Inj_mean,ALLtcReads_Inj_mean)
allGenes = read.table("/Volumes//Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure1/data//allGenes.txt",sep="\t",stringsAsFactors = F)

#allGenes = read.table("/Volumes//Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure1/data//allGenes.txt",sep="\t",stringsAsFactors = F)
tcReads_totalReads = plyr::join(tcReads_totalReads, allGenes)
preDefinedGenes = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawFigures/figure1/data/preDefGenes.txt",
                             sep="\t",stringsAsFactors = F, header = T)

# preDefinedGenes = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawFigures/figure1/data/preDefGenes.txt",
#                              sep="\t",stringsAsFactors = F, header = T)

preDefinedGenes = preDefinedGenes %>% dplyr::mutate(type_otherStudies= type ) %>% dplyr::select(external_gene_name,type_otherStudies)
tcReads_totalReads = plyr::join(tcReads_totalReads,preDefinedGenes)

classifiedGenes = read.table('/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/packages/data/countingWindows_classified_conversions.txt',
                             sep="\t", header = T,stringsAsFactors = F)
tcReads_totalReads = tcReads_totalReads[tcReads_totalReads$external_gene_name %in% classifiedGenes$gene,]
zygoticGenes = tcReads_totalReads[which(tcReads_totalReads$type_otherStudies == "Z"),] %>% dplyr::filter(Inj_mean_TP1==0) %>% dplyr::filter(totalReads_9>0)
MaternalGenes = tcReads_totalReads[which(tcReads_totalReads$type_otherStudies == "M"),] %>% dplyr::filter(Inj_mean_TP1>1)
OtherGenes = tcReads_totalReads[is.na(tcReads_totalReads$type_otherStudies),] %>% dplyr::filter(Inj_mean_TP1>1)
####### it seems that there are a lot of maternal genes whcih have T>C conversions at tp9

  ##### plotting conversion rates for these genes.. 
TCconversionRate = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other//rawfigures/figure2/data//perGenes_conversion.txt",
                              sep="\t",stringsAsFactors = F, header = T)

# TCconversionRate = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure2/data//perGenes_conversion.txt",
#                               sep="\t",stringsAsFactors = F, header = T)
errorRates = read.table("/Volumes//Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawFigures/figure2/data//errorRates_predicted_observed.txt",header = T,stringsAsFactors = F)
TCconversionRate = TCconversionRate[TCconversionRate$name  %in% classifiedGenes$gene,]


Injection_R2 =  TCconversionRate[grep("Inj_R2",names(TCconversionRate))]
Injection_R3 =  TCconversionRate[grep("Inj_R3",names(TCconversionRate))]

order_samples = c("Inj_R2_TP1","Inj_R2_TP2","Inj_R2_TP3","Inj_R2_TP4","Inj_R2_TP5","Inj_R2_TP6","Inj_R2_TP7","Inj_R2_TP8","Inj_R2_TP9")
Injection_R2 = Injection_R2[order_samples]  
order_samples = c("Inj_R3_TP1","Inj_R3_TP2","Inj_R3_TP3","Inj_R3_TP4","Inj_R3_TP5","Inj_R3_TP6","Inj_R3_TP7","Inj_R3_TP8","Inj_R3_TP9")
Injection_R3 = Injection_R3[order_samples]  



errorRates = errorRates %>% dplyr::mutate(condition = paste("Inj",replicate,time,sep="_"))
errorRates_R2 = errorRates %>% dplyr::filter(replicate == "R2")
errorRates_R3 = errorRates %>% dplyr::filter(replicate == "R3")

for(i in 1:length(Injection_R2)){
  Injection_R2[[i]] = Injection_R2[[i]] - errorRates_R2$predicted[i]
  Injection_R2[[i]][which(Injection_R2[[i]] < 0)] <-0
}


for(i in 1:length(Injection_R3)){
  Injection_R3[[i]] = Injection_R3[[i]] - errorRates_R3$predicted[i]
  Injection_R3[[i]][which(Injection_R3[[i]] < 0)] <-0
}



Injection_R2.mat = do.call(cbind.data.frame,Injection_R2)
Injection_R3.mat = do.call(cbind.data.frame,Injection_R3)

bothConversionRates = cbind.data.frame(Injection_R2.mat,Injection_R3.mat)
colnames(bothConversionRates) = paste0('ConversionRates',colnames(bothConversionRates))
colnames(Injection_R3.mat) = c(0.75,2,2.5,3,3.5,4,4.5,5,5.5)
colnames(Injection_R2.mat) = c(0.75,2,2.5,3,3.5,4,4.5,5,5.5)
mean_Inj.mat = (Injection_R2.mat + Injection_R3.mat)/2

colnames(mean_Inj.mat) = c(paste0("ConversionRates_",c(1:9)))
untreatedSamples = TCconversionRate %>% dplyr::select(dplyr::contains('Untre'))
untreatedSamples[is.na(untreatedSamples)]<-0
mean_Inj.mat = cbind.data.frame(mean_Inj.mat,bothConversionRates)
mean_Inj.mat$external_gene_name = TCconversionRate$name
mean_Inj.mat =  plyr::join(mean_Inj.mat,tcReads_totalReads)
mean_Inj.mat = cbind.data.frame(mean_Inj.mat,untreatedSamples)
mean_Inj.mat = mean_Inj.mat[mean_Inj.mat$external_gene_name %in% classifiedGenes$gene,]

M = mean_Inj.mat %>% dplyr::filter(type_otherStudies == "M")   %>% dplyr::filter(Inj_mean_TP1>1)
Z = mean_Inj.mat %>% dplyr::filter(type_otherStudies == "Z") %>% dplyr::filter(Inj_mean_TP1==0)
write.table(Z,"/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/data/finalZgenes.txt",
            sep="\t",quote = F, row.names = F,col.names = T)
M_Z = c(M$external_gene_name,Z$external_gene_name)
`%!in%` = Negate(`%in%`)
O = mean_Inj.mat[mean_Inj.mat$external_gene_name %!in% M_Z,] 




##### now alsio add thing the predefined genes based on oogenesis meiosis
preDefinedTotalGenes = read.table('/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures//figure1/data/oocyte_meiosis.txt',sep="\t",stringsAsFactors = F)
preDefinedTotalGenes = preDefinedTotalGenes %>% dplyr::filter(!is.na(GO))
Maternal_selectedGenes = mean_Inj.mat[mean_Inj.mat$external_gene_name %in% preDefinedTotalGenes$external_gene_name,]


### manually removing some genes which are also invilved in mitosiss.  defined in ⁨Users⁩ ▸ ⁨pooja.bhat⁩ ▸ ⁨Dropbox (VBC)⁩ ▸ ⁨Paperoutline⁩ ▸ ⁨figureOutline⁩ ▸ ⁨paper_ZebrafishSLAMseq⁩ ▸ ⁨figures⁩ ▸ ⁨figure1⁩ ▸ ⁨data⁩ ▸ ⁨exampleGenes⁩

Maternal_selectedGenes = Maternal_selectedGenes[Maternal_selectedGenes$external_gene_name %!in% c('ccnb1','wee2','pttg1','birc5b','rad21l1','ercc4','smad2'),]

Maternal_write = Maternal_selectedGenes %>% dplyr::select(c(external_gene_name,ConversionRates_1 ,ConversionRates_2, 
                                    ConversionRates_3, ConversionRates_4, ConversionRates_5,
                                    ConversionRates_6, ConversionRates_7,ConversionRates_8,ConversionRates_9,
                                    Untreated_TP1,Untreated_TP6,Untreated_TP9,TC_1 ,
                                    TC_2, TC_3, TC_4, TC_5, TC_6, TC_7, TC_8, TC_9,allTCReads_1 ,
                                    allTCReads_2, allTCReads_3, allTCReads_4, allTCReads_5, allTCReads_6,
                                    allTCReads_7, allTCReads_8, allTCReads_9,
                                    totalReads_1, totalReads_2, totalReads_3, totalReads_4,
                                    totalReads_5,totalReads_6, totalReads_7, totalReads_8, totalReads_9,Inj_mean_TP1,
                                    Inj_mean_TP2, Inj_mean_TP3, Inj_mean_TP4, Inj_mean_TP5, Inj_mean_TP6, Inj_mean_TP7,
                                    Inj_mean_TP8, Inj_mean_TP9 ) ) %>% dplyr::mutate(type = 'Maternal')

colnames(Maternal_write ) = c('external_gene_name',paste0('meanConversion_',c(1:9)) ,'conversion_untreated_1','conversion_untreated_6',
                              'conversion_untreated_9',paste0('multipleTC_',c(1:9)), paste0('allTCReads_',c(1:9)), paste0('totalReads_',c(1:9)),
                               paste0('RPM_',c(1:9)),'type') 

zygotic_write = Z %>% dplyr::select(c(external_gene_name,ConversionRates_1 ,ConversionRates_2, 
                                                    ConversionRates_3, ConversionRates_4, ConversionRates_5,
                                                    ConversionRates_6, ConversionRates_7,ConversionRates_8,ConversionRates_9,
                                                    Untreated_TP1,Untreated_TP6,Untreated_TP9,TC_1 ,
                                                    TC_2, TC_3, TC_4, TC_5, TC_6, TC_7, TC_8, TC_9,allTCReads_1 ,
                                                    allTCReads_2, allTCReads_3, allTCReads_4, allTCReads_5, allTCReads_6,
                                                    allTCReads_7, allTCReads_8, allTCReads_9,
                                                    totalReads_1, totalReads_2, totalReads_3, totalReads_4,
                                                    totalReads_5,totalReads_6, totalReads_7, totalReads_8, totalReads_9,Inj_mean_TP1,
                                                    Inj_mean_TP2, Inj_mean_TP3, Inj_mean_TP4, Inj_mean_TP5, Inj_mean_TP6, Inj_mean_TP7,
                                                    Inj_mean_TP8, Inj_mean_TP9 ))%>% dplyr::mutate(type = 'Zygotic')

colnames(zygotic_write ) = c('external_gene_name',paste0('meanConversion_',c(1:9)) ,'conversion_untreated_1','conversion_untreated_6',
                              'conversion_untreated_9',paste0('multipleTC_',c(1:9)), paste0('allTCReads_',c(1:9)), paste0('totalReads_',c(1:9)),
                              paste0('RPM_',c(1:9)),'type') 



TC3_zgenes = cbind.data.frame(tcReads3TC_Inj_mean[tcReads3TC_Inj_mean$external_gene_name %in% zygotic_write$external_gene_name,] %>%
                                         dplyr::select(dplyr::contains("TCReads3")),zygotic_write)
TC3_zgenes = TC3_zgenes %>% dplyr::mutate(fraction3TC = TCReads3_9/totalReads_9) %>% dplyr::select(fraction3TC) %>% dplyr::pull(1)

#### I also want to plot the fractions of T>C containing reads.. 
TCreads = totalWrite %>% dplyr::filter(type == 'Zygotic') %>% dplyr::select(c(allTCReads_9,multipleTC_9,totalReads_9))  %>% 
  dplyr::mutate(fractionAllTC=allTCReads_9/totalReads_9, fractionMultiplTC =multipleTC_9/ totalReads_9) %>% dplyr::select(c(fractionAllTC,fractionMultiplTC))
TCreads = cbind.data.frame(TCreads,TC3_zgenes)
TCreads = reshape2::melt(TCreads)

### samme to high confidence maternal genes
TCreads_hcM = totalWrite %>% dplyr::filter(type=='Maternal') %>% dplyr::select(c(allTCReads_9,multipleTC_9,totalReads_9))  %>% 
  dplyr::mutate(fractionAllTC=allTCReads_9/totalReads_9, fractionMultiplTC =multipleTC_9/ totalReads_9) %>% dplyr::select(c(fractionAllTC,fractionMultiplTC))
TC3_hcMgenes = cbind.data.frame(tcReads3TC_Inj_mean[tcReads3TC_Inj_mean$external_gene_name %in% Maternal_write$external_gene_name,] %>%
                                dplyr::select(dplyr::contains("TCReads3")),Maternal_write)
TC3_hcMgenes = TC3_hcMgenes %>% dplyr::mutate(fraction3TC = TCReads3_9/totalReads_9) %>% dplyr::select(fraction3TC) %>% dplyr::pull(1)
TCreads_hcM = cbind.data.frame(TCreads_hcM,TC3_hcMgenes)
TCreads_hcM = reshape2::melt(TCreads_hcM)


pdf("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/plots/fractionTC_singleAndMultiple.pdf", 
    height=4, width=3)
    fractionAllAndMultipleTC =   ggpubr::ggviolin(TCreads,x='variable',y='value',add='boxplot',fill='red',ylab = 'Fraction T>C reads', ylim=c(0,1),trim = T) + theme_ameres(type='barplot')
    print(fractionAllAndMultipleTC)
    fractionAllAndMultipleTC_M =   ggpubr::ggviolin(TCreads_hcM,x='variable',y='value',add='boxplot',fill='red',ylab = 'Fraction T>C reads', ylim=c(0,1),trim = T) + theme_ameres(type='barplot')
    print(fractionAllAndMultipleTC_M)
dev.off()
removeThese = which(Maternal_selectedGenes$ConversionRates_9 < Maternal_selectedGenes$Untreated_TP9 | 
                                  Maternal_selectedGenes$ConversionRates_9 < Maternal_selectedGenes$Untreated_TP6 | Maternal_selectedGenes$ConversionRates_9 < Maternal_selectedGenes$Untreated_TP1)




### combining all M and the selected M genes
M = rbind.data.frame(M,Maternal_selectedGenes)
M = M[!duplicated(M$external_gene_name),]

Mall_write = M %>% dplyr::select(c(external_gene_name,ConversionRates_1 ,ConversionRates_2, 
                                      ConversionRates_3, ConversionRates_4, ConversionRates_5,
                                      ConversionRates_6, ConversionRates_7,ConversionRates_8,ConversionRates_9,
                                      Untreated_TP1,Untreated_TP6,Untreated_TP9,TC_1 ,
                                      TC_2, TC_3, TC_4, TC_5, TC_6, TC_7, TC_8, TC_9,allTCReads_1 ,
                                      allTCReads_2, allTCReads_3, allTCReads_4, allTCReads_5, allTCReads_6,
                                      allTCReads_7, allTCReads_8, allTCReads_9,
                                      totalReads_1, totalReads_2, totalReads_3, totalReads_4,
                                      totalReads_5,totalReads_6, totalReads_7, totalReads_8, totalReads_9,Inj_mean_TP1,
                                      Inj_mean_TP2, Inj_mean_TP3, Inj_mean_TP4, Inj_mean_TP5, Inj_mean_TP6, Inj_mean_TP7,
                                      Inj_mean_TP8, Inj_mean_TP9 ))%>% dplyr::mutate(type = 'Maternal based on gene expression')

colnames(Mall_write ) = c('external_gene_name',paste0('meanConversion_',c(1:9)) ,'conversion_untreated_1','conversion_untreated_6',
                             'conversion_untreated_9',paste0('multipleTC_',c(1:9)), paste0('allTCReads_',c(1:9)), paste0('totalReads_',c(1:9)),
                             paste0('RPM_',c(1:9)),'type') 

totalWrite = rbind.data.frame(zygotic_write, Maternal_write,Mall_write)


write.table(totalWrite, '/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures//figure1/data/info_exampleMandZ.txt',
            sep = "\t",quote = F, row.names = F)


M_noConversions = M %>% dplyr::filter(TC_9 == 0)
M_Conversions = M %>% dplyr::filter(TC_9 > 0)

write.table(M_noConversions,"/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure1/data/Mgenes_noConversions.txt",
            sep="\t",quote = F, row.names = F,col.names = F)

write.table(M_Conversions,"/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure1/data/Mgenes_Conversions.txt",
            sep="\t",quote = F, row.names = F,col.names = F)


countingWindows = read.table('/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures//figure1/data/ConversionRates_allCws.txt',
                             sep="\t",stringsAsFactors = F, header = T)

countingWindows = splitReplicates( countingWindows,condition = "Inj",metadata_add =countingWindows[,c(1:6)] )[[3]]
countingWindows$Inj_R2_TP9[is.na(countingWindows$Inj_R2_TP9)]<-0
countingWindows_split = split(countingWindows,countingWindows$name,T)
countingWindows_split = lapply(countingWindows_split,function(x) x[which.max(x$Inj_R2_TP9),])
countingWindows_split = do.call(rbind.data.frame,countingWindows_split)


Z_countingWindows =  countingWindows_split[countingWindows_split$name %in% Z$external_gene_name,]  %>% 
  dplyr::select(c('chr','start','end','name','score','strand'))
M_countingWindows =  countingWindows_split[countingWindows_split$name %in% M$external_gene_name,] %>%
  dplyr::select(c('chr','start','end','name','score','strand'))

M_noConversions_countingWindows =  countingWindows_split[countingWindows_split$name %in% M_noConversions$external_gene_name,]  %>% 
  dplyr::select(c('chr','start','end','name','score','strand'))

M_Conversions_countingWindows =  countingWindows_split[countingWindows_split$name %in% M_Conversions$external_gene_name,]  %>% 
  dplyr::select(c('chr','start','end','name','score','strand'))

O_Conversions_countingWindows =  countingWindows_split[countingWindows_split$name %in% O$external_gene_name,]  %>% 
  dplyr::select(c('chr','start','end','name','score','strand'))

M_selectedGenes_countingWindows = countingWindows_split[countingWindows_split$name %in% Maternal_selectedGenes$external_gene_name,]  %>% 
  dplyr::select(c('chr','start','end','name','score','strand'))


write.table(Z_countingWindows,"/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/data/preDefinedZgenes.bed",
            sep="\t",quote = F, col.names = F,row.names = F)
write.table(M_countingWindows,"/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/data/preDefinedMgenes.bed",
            sep="\t",quote = F, col.names = F,row.names = F)

write.table(M_noConversions_countingWindows,"/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/data/preDefinedM_noConversionsgenes.bed",
            sep="\t",quote = F, col.names = F,row.names = F)
write.table(M_Conversions_countingWindows,"/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/data/preDefinedM_Conversionsgenes.bed",
            sep="\t",quote = F, col.names = F,row.names = F)

write.table(O_Conversions_countingWindows,"/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/data/preDefinedOther_Conversionsgenes.bed",
            sep="\t",quote = F, col.names = F,row.names = F)

write.table(M_selectedGenes_countingWindows,"/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/data/preDefinedMselectedGenes_Conversionsgenes.bed",
            sep="\t",quote = F, col.names = F,row.names = F)


pdf('/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/plots/fractionTC_conversionRates_NumberTC.pdf', height = 3, width =4 )
      
      Z_parameterplots = plotConverionRates_TC(dataF = Z,addName = 'Zygotic')
      MnoConv_parameterplots = plotConverionRates_TC(dataF = M %>%dplyr::filter(TC_9==0),addName = 'Maternal_noConversions')
      MConvContaining_parameterplots = plotConverionRates_TC(dataF = M %>%dplyr::filter(TC_9>0),addName = 'Maternal with conversions')
      Mall_parameterplots = plotConverionRates_TC(dataF = M ,addName = 'Maternal all')
      Other_parameterPlots = plotConverionRates_TC(dataF = O ,addName = 'other all')
      M_oogenesis_meiosis = plotConverionRates_TC(dataF = Maternal_selectedGenes ,addName = 'Maternal_oogeneesis_meiosis')
      
      Z_parameterplots[[1]]
      MnoConv_parameterplots[[1]]
      MConvContaining_parameterplots[[1]]
      Mall_parameterplots[[1]]
      Other_parameterPlots[[1]]
      M_oogenesis_meiosis[[1]]
      
      Z_parameterplots[[2]]
      Z_parameterplots[[4]]
      Z_parameterplots[[5]]

      MnoConv_parameterplots[[2]]
      MConvContaining_parameterplots[[2]]
      Mall_parameterplots[[2]]
      Mall_parameterplots[[4]]
      Mall_parameterplots[[5]]

      Other_parameterPlots[[2]]
      Other_parameterPlots[[4]]
      Other_parameterPlots[[5]]
      M_oogenesis_meiosis[[2]]

      ########### plotting violin plots ... #####################

      Z_convRates = Z_parameterplots$data_conversionRates %>% dplyr::mutate(type = 'Z')
      MnoConv_parameterplots$data_conversionRates = MnoConv_parameterplots$data_conversionRates %>% dplyr::mutate(type = 'M-noConversion')
      MConvContaining_parameterplots$data_conversionRates = MConvContaining_parameterplots$data_conversionRates %>% dplyr::mutate(type = 'M-Conversion')
      Mall_parameterplots$data_conversionRates = Mall_parameterplots$data_conversionRates %>% dplyr::mutate(type = 'M-all')
      
      allM = list(MnoConv_parameterplots$data_conversionRates,MConvContaining_parameterplots$data_conversionRates,Mall_parameterplots$data_conversionRates)
      names(allM) = c('noConversions','ConversionContaining','All')
      
      
      for(i in 1:length(allM)){
        allConvRates = rbind.data.frame(Z_convRates,allM[[i]]) 
        
        df.summary <- allConvRates %>%
          group_by(L1,type) %>%
          summarise(
            sd = sd(value, na.rm = TRUE),
            len = mean(value,na.rm=T)
          )
        df.summary[is.na(df.summary)]<-0
        
        my_comparisons <- list( c("M-noConversion", "Z") )
        allConvRates %>%  dplyr::group_by(L1,type) %>% summarise(n=n(), avg = mean(value,na.rm=T)) ->Summary.data
        
        quasirandom_plot = ggplot(allConvRates,aes(L1, value,color=type)) + geom_quasirandom(varwidth =T, alpha=.2) + scale_color_brewer(palette = 'Set1')
        quasirandom_plot = quasirandom_plot + theme_cowplot() + xlab('Time (hpf)') + ylab ('T>C conversion rate')
        quasirandom_plot = quasirandom_plot + theme_ameres(type = 'barplot')+ ylim(c(0,0.15))
        quasirandom_plot = quasirandom_plot + geom_line(data = df.summary,aes(x=L1,y=len-sd,group=type),linetype = "dashed")  + geom_line(data = df.summary,aes(x=L1,y=len+sd,group=type),linetype = "dashed") + 
          stat_compare_means(aes(group = type), label = "p.signif",method = 'wilcox.test') +   geom_text(data=Summary.data ,aes(x = L1, y = avg, label=n,group=type),color="black", fontface =1, size = 2)
        print(quasirandom_plot)
        
      }
  
dev.off()

pdf('/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/plots/fractionTC_violin.pdf', height = 4, width =2 )

    ggviolin(Z_parameterplots$fractionTC %>% dplyr::filter(L1==5.5),x='L1',y='value',add='boxplot',ylab='Fraction >1 TC containing reads',trim = F) + 
      theme_ameres(type = 'barplot') + theme_ameres(type = 'barplot')
    ggviolin(M_oogenesis_meiosis$fractionTC %>% dplyr::filter(L1==5.5),x='L1',y='value',add='boxplot',ylab='Fraction >1 TC containing reads',ylim=c(-0.1,1),trim = F) +
      theme_ameres(type = 'ameres') + theme_ameres(type = 'barplot')
    ggviolin(Mall_parameterplots$fractionTC %>% dplyr::filter(L1==5.5),x='L1',y='value',add='boxplot',ylab='Fraction >1 TC containing reads',ylim=c(-0.1,1),trim = T) +
      theme_ameres(type = 'barplot') + theme_ameres(type = 'barplot')

    ggviolin(Z_parameterplots$fractionAllTC %>% dplyr::filter(L1==5.5),x='L1',y='value',add='boxplot',ylab='Fraction all TC containing reads',ylim=c(-0.1,1),trim = F)
    allTC = Z_parameterplots$fractionAllTC %>% dplyr::mutate(category = 'All TC')
    MultiTC = Z_parameterplots$fractionTC %>% dplyr::mutate(category = 'Multi TC')
    allZfractions = rbind.data.frame(allTC,MultiTC)
    ggpubr::ggviolin(allZfractions %>% dplyr::filter(L1 == 5.5) ,x='category',y='value', add = 'boxplot',label = median(allZfractions$value),ylab = 'Fraction TC reads')   + 
  stat_summary(fun.data = function(x) data.frame(y=0.8, label = paste("Median=",round(median(x),2))), geom="text")

dev.off()



###### fractionTC and allTC 
library(cowplot)



pdf('/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure1//plots/tcreads_totalReads_2TC.pdf', height = 4,width=4)
      plotTC_total(dataF = Z,timeP = "9")[[2]]
      plotTC_dotPlot(dataF = Z,timeP = "1")
      plotTC_total(dataF =M ,timeP = "9")
      plotTC_dotPlot(dataF = M,timeP = "1")
      plotTC_total(dataF =M %>% dplyr::filter(TC_9>0) ,timeP = "9")
      plotTC_dotPlot(dataF = O,timeP = "1")
      plotTC_total(dataF = O,timeP = "9")
      plotTC_dotPlot(dataF = Maternal_selectedGenes,timeP = "9")
      plotTC_total(dataF = Maternal_selectedGenes,timeP = "9")[[2]]
      
dev.off()


##################### Lookss like there are some maternal genes which have increase T>C conversions at TP9
M_withConversions = M %>% dplyr::filter(TC_9>0)
M_withoutConversions = M %>% dplyr::filter(TC_9==0)

###################################################################################
############### the change in gene expression over time ###########################

M_highestConversion = M_withConversions[order(M_withConversions$TC_9,decreasing = T ),] 
row.names(M_highestConversion) = NULL
M_highestConversion= M_highestConversion %>% tibble::column_to_rownames('external_gene_name')
M_Conversions  = M_highestConversion[,grep('ConversionRates',colnames(M_highestConversion))]

M_highestConversion  = M_highestConversion[,grep('Inj',colnames(M_highestConversion))]
colnames(M_highestConversion) = c(0.75,2,2.5,3,3.5,4,4.5,5,5.5)
colnames(M_Conversions) = c(0.75,2,2.5,3,3.5,4,4.5,5,5.5)

kpna2 = data.frame(melt(M_highestConversion[1,]),gene = 'kpna2')
cox4i1 = data.frame(melt(M_highestConversion[2,]),gene = 'cox4i1')
cdk1 = data.frame(melt(M_highestConversion[3,]),gene = 'cdk1')

allGenes = rbind.data.frame(kpna2,cox4i1,cdk1)
allGenes = split(allGenes,allGenes$gene,T)

pdf('/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/plots/maternalGenes_withConversions.pdf', height=4, width=4)
  
    for(i in 1:nrow(M_highestConversion)){
      gene_current = melt(M_highestConversion[i,])
      p = ggpubr::ggline(gene_current,x='variable',y='value',color='black', xlab = 'Time (hpf)',ylab = 'Gene expression (RPM)',title = rownames(M_highestConversion)[i] ) + theme_ameres(type = 'barplot')       
      print(p)
    }


      ########
        
               for(i in 1:nrow(M_Conversions)){
                 gene_current = melt(M_Conversions[i,])
                 
                p = ggpubr::ggline(gene_current,x='variable',y='value',color='black', xlab = 'Time (hpf)',ylab = 'TC conversion rate',title = rownames(M_Conversions)[i]) + theme_ameres(type = 'barplot')       
                print(p)
              }

dev.off()

