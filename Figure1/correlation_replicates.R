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



# coutData_files = list.files("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_may2019//counts_expandedImplementation/",pattern = "*_perCW.tsv")
# coutData_files_path = paste0("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_may2019/counts_expandedImplementation//",coutData_files)

coutData_files = list.files("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/data//counts_expandedImplementation/",pattern = "*_perCW.tsv")
coutData_files_path = paste0("//Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/data//counts_expandedImplementation//",coutData_files)
countFiles = lapply(coutData_files_path,function(x) read.table(x,header = T,stringsAsFactors = F,sep="\t"))

timepoints = paste0("TP",c(1:9))
orderAll = (c(c(paste0("Inj_","R1_",timepoints),paste0("Inj_","R2_",timepoints),paste0("Inj_","R3_",timepoints)),c(paste0("Inc_","R1_",timepoints),paste0("Inc_","R2_",timepoints),paste0("Inc_","R3_",timepoints)),c(paste0("Untreated_",timepoints))))
orderAll = data.frame(L1=orderAll)

#sampleInfo = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/sampleInfo/barcodes_description.txt",sep="\t",stringsAsFactors = F)
sampleInfo = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/data/barcodes_description.txt",
                        sep="\t",stringsAsFactors = F)
sampleInfo = sampleInfo[order(sampleInfo$V2),]
names(countFiles) = sampleInfo$V3

sumCounts = reshape::melt(lapply(countFiles,function(x) sum(x$sumCounts)))
orderAll = plyr::join(orderAll,sumCounts)
orderAll = orderAll[complete.cases(orderAll),]
orderAll = orderAll[-grep("Inc",orderAll$L1),]

library(cowplot)
library(ggrastr)
library(ggplot2)

# RPMdata = read.table('/Volumes//groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_may2019/dataTables_expandedCounting/RPM_allCws.txt',
#                      sep="\t",stringsAsFactors = F, header = T)

RPMdata = read.table('/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/data//RPM_allCws.txt',
                                           sep="\t",stringsAsFactors = F, header = T)
                     
RPMdata =  RPMdata %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Inj_R3_TP9), sum, na.rm = TRUE)

dir.create("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/correlationAnalysis/")

#pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/correlationAnalysis/correlation_RPMs.pdf")
pdf("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/plots//correlation_RPMs.pdf")

      for(i in 1:9){
      R2_samples = paste0("Inj_R2_TP",i)
      R3_samples = paste0("Inj_R3_TP",i)
      tmp_df = RPMdata[,c(R2_samples,R3_samples,'name')]
      tmp_df = as.data.frame(tmp_df)
      tmp_df = tmp_df[which(tmp_df[,1]>1 & tmp_df[,2]>1),]
      colnames(tmp_df) = c("R2","R3",'name')
      tmp_df = tmp_df %>% dplyr::mutate(log10R2 = log10(R2),log10R3 = log10(R3))  %>%
        mutate(densityCols_TP = densCols(log10R2, log10R3,  colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))) 
      #mt  = tmp_df[grep('mt-',tmp_df$name),]
      p = ggplot2::ggplot(tmp_df) + ggrastr::geom_point_rast(aes(x=log10R2,y=log10R3),col=tmp_df$densityCols_TP) + theme_cowplot()  +xlim(c(0,5))+ylim(c(0,5))
      p = p  + theme_ameres(type = 'barplot') + ggpubr::stat_cor(aes(x=log10R2,y=log10R3)) + xlab("RPM-R1(log10)")  + ylab("RPM-R2(log10)")
      p = p + ggtitle(paste("TP",i,"n=",nrow(tmp_df)))
      #p =  p + geom_point(data = mt ,aes(x=log10R2,log10R3),col='black')
      print(p) 
        
      }
dev.off()



##################################### just read the TC containgin reads. 
library(dplyr)
# totalReads = read.table('/Volumes//groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_may2019/dataTables_expandedCounting/totalCounts_allCws.txt',
#                         sep="\t",stringsAsFactors = F, header = T)
totalReads = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure2/data/totalCounts_allCws.txt",
                        sep="\t",stringsAsFactors = F, header = T)
totalReads =  totalReads %>% dplyr::select(-dplyr::contains('Inc')  )
totalReads =  totalReads %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Untreated_TP9), sum, na.rm = TRUE)

#conversionRates  = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017//analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_may2019//dataTables_expandedCounting//perGenes_conversion.txt",sep="\t",stringsAsFactors = F,header=T)
conversionRates = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure2/data/perGenes_conversion.txt",
                             sep="\t",stringsAsFactors = F, header = T)
conversionRates =  conversionRates %>% dplyr::select(-dplyr::contains('Inc')  )
conversionRates =  conversionRates %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Untreated_TP9), sum, na.rm = TRUE)

# TCreads = read.table('/Volumes//groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_may2019/dataTables_expandedCounting/numberOfreadsWithTC.txt',
#                      sep="\t",stringsAsFactors = F, header = T)
TCreads = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure2/data/numberOfreadsWithMultipleTC.txt",
                     stringsAsFactors = F, header = T, sep = "\t")
TCreads = TCreads %>% dplyr::select(-dplyr::contains('Inc')  )
TCreads =  TCreads %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Untreated_TP9), sum, na.rm = TRUE)

reads_tc = TCreads

errorRates = read.table("/Volumes//Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure2/data//errorRates_predicted_observed.txt",header = T,stringsAsFactors = F)

#errorRates = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_may2019//errorRates_predicted_observed.txt",header = T)

###### I am not removing erors because i am considering read with 2 tc
for(i in 1:nrow(errorRates)){
 # TCreads[,i+1] = TCreads[,i+1] - (totalReads[,i+1] * errorRates$predicted[i] )
  TCreads[,i+1] = (TCreads[,i+1] * 1000000)/orderAll$value[i] ### still need to convert to rpm
}

####################################
#pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/correlationAnalysis/correlation_RPMsTCreads.pdf")
pdf("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/plots//correlation_RPMsTCreads.pdf")

for(i in 1:9){
        
        R2_samples = paste0("Inj_R2_TP",i)
        R3_samples = paste0("Inj_R3_TP",i)
        ###### putting a condition that conversion at the time point should be higher than conversion rates in untreated samples
       # consider_tp_R1 = conversionRates[which(conversionRates[,R2_samples] >conversionRates$Untreated_TP1 & conversionRates[,R2_samples]>conversionRates$Untreated_TP6 & conversionRates[,R2_samples] > conversionRates$Untreated_TP9),]
       # consider_tp_R2 = conversionRates[which(conversionRates[,R3_samples] >conversionRates$Untreated_TP1 & conversionRates[,R3_samples]>conversionRates$Untreated_TP6 & conversionRates[,R3_samples] > conversionRates$Untreated_TP9),]
       # consider_tp = intersect(consider_tp_R2$name,consider_tp_R1$name)

       #TP1_reads = TCreads[TCreads$name %in% consider_tp,]
       TP1_reads = TCreads
        tmp_df = TP1_reads[,c(R2_samples,R3_samples,'name')]
        
        tmp_df = tmp_df[which(tmp_df[,1]>1 & tmp_df[,2]>1),]
        colnames(tmp_df) = c("R2","R3","name")
        tmp_df = tmp_df %>% dplyr::mutate(log10R2 = log10(R2+1),log10R3 = log10(R3+1)) %>%
          mutate(densityCols_TP = densCols(log10R2, log10R3,  colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))) 
        mt  = tmp_df[grep('mt-',tmp_df$name),]
        
        p = ggplot2::ggplot(tmp_df) + ggrastr::geom_point_rast(aes(x=log10R2,y=log10R3),col=tmp_df$densityCols_TP) + theme_cowplot() +xlim(c(0,5))+ylim(c(0,5))
        p = p  + theme_ameres(type = 'barplot') + ggpubr::stat_cor(aes(x=log10R2,y=log10R3)) + xlab("RPM-R1(log10)")  + ylab("RPM-R2(log10)")
        p = p + ggtitle(paste("TP",i,"n=",nrow(tmp_df))) 
       
        print(p)
}
dev.off()

