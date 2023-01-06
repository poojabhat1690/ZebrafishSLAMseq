
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

conversionRates  = read.table("//Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline//paper_ZebrafishSLAMseq/other/rawfigures/figure2/data//perGenes_conversion.txt",sep="\t",stringsAsFactors = F,header=T)
conversionRates =  conversionRates %>% dplyr::select(-dplyr::contains('Inc')  )
conversionRates =  conversionRates %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Untreated_TP9), sum, na.rm = TRUE)
#conversionRates = conversionRates[conversionRates$name %in% classifiedGenes$name,]

errorRates = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures//figure2/data//errorRates_predicted_observed.txt",header = T)
conversionRates = conversionRates %>% tibble::column_to_rownames('name')
for(i in 1:ncol(conversionRates)){
  conversionRates[,i] = conversionRates[,i] - errorRates$predicted[i]
}


conversionRates = as.matrix(conversionRates)
conversionRates[which(conversionRates <0 )] <- 0
conversionRates[grep("Inf",conversionRates)]<-0
conversionRates = as.data.frame(conversionRates)
conversionRates= conversionRates %>% tibble::rownames_to_column(var = 'name')
conversionRates_untreated = conversionRates[,c('Untreated_TP1','Untreated_TP6','Untreated_TP9')]
conversionRates_untreated$maxUntreated = apply(conversionRates_untreated ,1,function(x) max(x))
conversionRates = splitReplicates(dataFrameToSplit = conversionRates,condition = "Inj",metadata_add = conversionRates$name)[[3]]
conversionRates = cbind.data.frame(conversionRates,conversionRates_untreated)

colnames(conversionRates) = c(paste0("TP",c(1:9)),'name',"Untreated_TP1","Untreated_TP6","Untreated_TP9",'maxUntreated')


TC_samples = read.table('/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures//figure2/data/numberOfreadsWithMultipleTC.txt',
                        sep = "\t",stringsAsFactors = F, header = T)

allTC  = read.table('/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures//figure2/data/numberOfreadsWithTC.txt',
                    sep="\t",stringsAsFactors = F, header = T)

TA_samples = read.table('/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure2/data/numberOfreadsWithTA.txt',
                        sep="\t",stringsAsFactors = F, header = T)

TAmultiple_samples = read.table('/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure2/data/numberOfreadsWithMultipleTA.txt',
                                sep="\t",stringsAsFactors = F, header = T)

totalReads = read.table('//Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure2/data//totalCounts_allCws.txt',
                        sep="\t",stringsAsFactors = F, header = T)


TC_samples = TC_samples %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Untreated_TP9), sum, na.rm = TRUE)
TA_samples = TA_samples %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Untreated_TP9), sum, na.rm = TRUE)
TC_all = allTC %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Untreated_TP9), sum, na.rm = TRUE)
totalReads = totalReads %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Untreated_TP9), sum, na.rm = TRUE)
TAmultiple_samples = TAmultiple_samples %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Untreated_TP9), sum, na.rm = TRUE)


normalizers = 10^6/colSums(totalReads[,c(2:ncol(totalReads))])
for(i in 1:39){
  TC_samples[,i+1] = TC_samples[,i+1] * normalizers[i]
}

TC_samples_untreated = TC_samples [,c('Untreated_TP1','Untreated_TP6','Untreated_TP9')] 


RPM_allGenes = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures//figure2/data/RPM_allCws.txt",
                          sep="\t",stringsAsFactors = F, header = T)
RPM_allGenes = RPM_allGenes %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Untreated_TP9), sum, na.rm = TRUE)
RPM_allGenes_mean = splitReplicates(dataFrameToSplit = RPM_allGenes,condition = "Inj",metadata_add =RPM_allGenes$name )[[3]]
untreated_t1 = cbind.data.frame(RPM_allGenes %>% dplyr::select(Untreated_TP1),TC_samples_untreated$Untreated_TP1)
untreated_t6 = cbind.data.frame(RPM_allGenes %>% dplyr::select(Untreated_TP6),TC_samples_untreated$Untreated_TP6)
untreated_t9 = cbind.data.frame(RPM_allGenes %>% dplyr::select(Untreated_TP9),TC_samples_untreated$Untreated_TP9)
names_untreated = c('readRPM','TCrpm')
colnames(untreated_t1) = names_untreated
colnames(untreated_t6) = names_untreated
colnames(untreated_t9) = names_untreated
untreated_t1 = untreated_t1 %>% dplyr::filter(TCrpm>1)
untreated_t6 = untreated_t6 %>% dplyr::filter(TCrpm>1)
untreated_t9 = untreated_t9  %>% dplyr::filter(TCrpm>1)
total_samples = rbind.data.frame(untreated_t1,untreated_t6,untreated_t9)

threshold = t.test(total_samples$TCrpm,conf.level = 0.99)$conf.int[2]
#threshold = quantile(total_samples$TCrpm)[4]
TC_samples = splitReplicates(dataFrameToSplit =TC_samples,condition = "Inj",metadata_add = TC_samples$name )[[3]]
TC_samples = TC_samples %>% dplyr::group_by(Inj_R2_TP1,Inj_R2_TP2,Inj_R2_TP3,Inj_R2_TP4,Inj_R2_TP5,Inj_R2_TP6,Inj_R2_TP7,Inj_R2_TP8,Inj_R2_TP9) %>% dplyr::mutate(count=n()) %>%
  dplyr::filter(count==1)  %>% dplyr::ungroup()


tp1 = TC_samples[TC_samples$Inj_R2_TP1>threshold & TC_samples$Inj_R2_TP2>threshold & TC_samples$Inj_R2_TP3>threshold & TC_samples$Inj_R2_TP4>threshold
                 & TC_samples$Inj_R2_TP5>threshold & TC_samples$Inj_R2_TP6>threshold & TC_samples$Inj_R2_TP7>threshold & TC_samples$Inj_R2_TP8>threshold & TC_samples$Inj_R2_TP9>threshold ,]
tp2 = TC_samples[TC_samples$Inj_R2_TP2>threshold & TC_samples$Inj_R2_TP3>threshold & TC_samples$Inj_R2_TP4>threshold & TC_samples$Inj_R2_TP5>threshold & TC_samples$Inj_R2_TP6>threshold & TC_samples$Inj_R2_TP7>threshold & TC_samples$Inj_R2_TP8>threshold & TC_samples$Inj_R2_TP9>threshold ,]
tp3 = TC_samples[ TC_samples$Inj_R2_TP3>threshold & TC_samples$Inj_R2_TP4>threshold & TC_samples$Inj_R2_TP5>threshold & TC_samples$Inj_R2_TP6>threshold & TC_samples$Inj_R2_TP7>threshold & TC_samples$Inj_R2_TP8>threshold & TC_samples$Inj_R2_TP9>threshold ,]
tp4 = TC_samples[ TC_samples$Inj_R2_TP4>threshold & TC_samples$Inj_R2_TP5>threshold & TC_samples$Inj_R2_TP6>threshold & TC_samples$Inj_R2_TP7>threshold & TC_samples$Inj_R2_TP8>threshold & TC_samples$Inj_R2_TP9>threshold ,]
tp5 = TC_samples[  TC_samples$Inj_R2_TP5>threshold & TC_samples$Inj_R2_TP6>threshold & TC_samples$Inj_R2_TP7>threshold & TC_samples$Inj_R2_TP8>threshold & TC_samples$Inj_R2_TP9>threshold ,]
tp6 = TC_samples[ TC_samples$Inj_R2_TP6>threshold & TC_samples$Inj_R2_TP7>threshold & TC_samples$Inj_R2_TP8>threshold & TC_samples$Inj_R2_TP9>threshold ,]
tp7 = TC_samples[ TC_samples$Inj_R2_TP7>threshold & TC_samples$Inj_R2_TP8>threshold & TC_samples$Inj_R2_TP9>threshold ,]
tp8 = TC_samples[  TC_samples$Inj_R2_TP8>threshold & TC_samples$Inj_R2_TP9>threshold ,]
tp9 = TC_samples[  TC_samples$Inj_R2_TP9>threshold ,]



conversionRates$greaterThan_t1 = apply(conversionRates ,1, function(x) length(which(x[1:9]>x[14])))
conversionRates$greaterThan_t2 = apply(conversionRates ,1, function(x) length(which(x[2:9]>x[14])))
conversionRates$greaterThan_t3 = apply(conversionRates ,1, function(x) length(which(x[3:9]>x[14])))
conversionRates$greaterThan_t4 = apply(conversionRates ,1, function(x) length(which(x[4:9]>x[14])))
conversionRates$greaterThan_t5 = apply(conversionRates ,1, function(x) length(which(x[5:9]>x[14])))
conversionRates$greaterThan_t6 = apply(conversionRates ,1, function(x) length(which(x[6:9]>x[14])))
conversionRates$greaterThan_t7 = apply(conversionRates ,1, function(x) length(which(x[7:9]>x[14])))
conversionRates$greaterThan_t8 = apply(conversionRates ,1, function(x) length(which(x[8:9]>x[14])))
conversionRates$greaterThan_t9 = apply(conversionRates ,1, function(x) length(which(x[9]>x[14])))

t1_conversions = conversionRates %>% dplyr::filter(greaterThan_t1==9)
tp1 = tp1[tp1$metadata_add %in% t1_conversions$name,]

t2_conversions = conversionRates %>% dplyr::filter(greaterThan_t2==8)
tp2 = tp2[tp2$metadata_add %in% t2_conversions$name,]

t3_conversions = conversionRates %>% dplyr::filter(greaterThan_t3==7)
tp3 = tp3[tp3$metadata_add %in% t3_conversions$name,]

t4_conversions = conversionRates %>% dplyr::filter(greaterThan_t4==6)
tp4 = tp4[tp4$metadata_add %in% t4_conversions$name,]

t5_conversions = conversionRates %>% dplyr::filter(greaterThan_t5==5)
tp5 = tp5[tp5$metadata_add %in% t5_conversions$name,]

t6_conversions = conversionRates %>% dplyr::filter(greaterThan_t6==4)
tp6 = tp6[tp6$metadata_add %in% t6_conversions$name,]

t7_conversions = conversionRates %>% dplyr::filter(greaterThan_t7==3)
tp7 = tp7[tp7$metadata_add %in% t7_conversions$name,]

t8_conversions = conversionRates %>% dplyr::filter(greaterThan_t8==2)
tp8 = tp8[tp8$metadata_add %in% t8_conversions$name,]

t9_conversions = conversionRates %>% dplyr::filter(greaterThan_t9==1)
tp9 = tp9[tp9$metadata_add %in% t9_conversions$name,]

allPassing = list(tp1,tp2,tp3,tp4,tp5,tp6,tp7,tp8,tp9)
allTps = rbind.data.frame(tp1 %>% dplyr::mutate(time = '0.75hpf') %>% dplyr::select(metadata_add,time),
tp2 %>% dplyr::mutate(time = '2hpf')%>% dplyr::select(metadata_add,time),
tp3 %>% dplyr::mutate(time = '2.5hpf')%>% dplyr::select(metadata_add,time),
tp4 %>% dplyr::mutate(time = '3hpf')%>% dplyr::select(metadata_add,time),
tp5 %>% dplyr::mutate(time = '3.5hpf')%>% dplyr::select(metadata_add,time),
tp6 %>% dplyr::mutate(time = '4hpf')%>% dplyr::select(metadata_add,time),
tp7 %>% dplyr::mutate(time = '4.5hpf')%>% dplyr::select(metadata_add,time),
tp8 %>% dplyr::mutate(time = '5hpf')%>% dplyr::select(metadata_add,time),
tp9 %>% dplyr::mutate(time = '5.5hpf')%>% dplyr::select(metadata_add,time))

colnames(allTps) = c('gene','time')
write.table(allTps,"/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure2/data/genesPertimePoint.txt",
sep="\t",quote = F)
nMT = unlist(lapply(allPassing,function(x) nrow(x[grep('mt-',x$metadata_add),])))
nonMt = unlist(lapply(allPassing,function(x) nrow(x[-grep('mt-',x$metadata_add),])))
ngenesTO = data.frame(mt  = nMT, genomic = nonMt,time = c(0.75,2,2.5,3,3.5,4,4.5,5,5.5))
ngenesTO = reshape2::melt(ngenesTO, id.vars= c('time'))

mt_fraction=c()
rpm_time= c()
rpm_mt = c()

for(i in 1:length(allPassing)){
  mt_fraction=c(mt_fraction, sum(allPassing[[i]][grep('mt-',allPassing[[i]]$metadata_add),i]/0.57) /sum(allPassing[[i]][,i]/0.57))
  rpm_time = c(rpm_time,sum(allPassing[[i]][,i])/0.57)
  rpm_mt =  c(rpm_mt,sum(allPassing[[i]][grep('mt-',allPassing[[i]]$metadata_add),i]))
  
  
}

rpm_genomic = rpm_time-rpm_mt

fractionOftranscritome = rpm_time/colSums(RPM_allGenes_mean[,c(1:9)])
fractionMaternal = 1-fractionOftranscritome 
deNovoFrac = cbind.data.frame(new=fractionOftranscritome,maternal=fractionMaternal,time=c(0.75,2,2.5,3,3.5,4,4.5,5,5.5))
deNovoFrac = melt(deNovoFrac,id.vars = c('time'))
deNovoFrac$value = deNovoFrac$value * 100
pdf('/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure2/plots/fractionGenomic_Mitochondrial.pdf', height=4, width=4)

fractionTranscriptome = ggpubr::ggbarplot(deNovoFrac,x='time',y='value',fill = 'variable',ylab='Percent contribution [%]', label = round(deNovoFrac$value,2)) + theme_ameres(type = 'barplot') 
print(fractionTranscriptome)

# fc_genomic = log2(c(rpm_genomic[2]/rpm_genomic[1],rpm_genomic[3]/rpm_genomic[2],rpm_genomic[4]/rpm_genomic[3],rpm_genomic[5]/rpm_genomic[4],
#   rpm_genomic[6]/rpm_genomic[5],rpm_genomic[7]/rpm_genomic[6],rpm_genomic[8]/rpm_genomic[7],rpm_genomic[9]/rpm_genomic[8]))
# fc_mt = log2(c(rpm_mt[2]/rpm_mt[1],rpm_mt[3]/rpm_mt[2],rpm_mt[4]/rpm_mt[3],rpm_mt[5]/rpm_mt[4],
#   rpm_mt[6]/rpm_mt[5],rpm_mt[7]/rpm_mt[6],rpm_mt[8]/rpm_mt[7],rpm_mt[9]/rpm_mt[8]))
# cumulative_sum =data.frame(genomicFc= fc_genomic,mtFC=fc_mt,time = c(2,2.5,3,3.5,4,4.5,5,5.5))
# cumulative_sum = melt(cumulative_sum,id.vars = c('time'))
# ggpubr::ggline(cumulative_sum,x='time',y='value',color='variable',ylab='FC transcriptional OP (log2)') + theme_ameres(type = 'barplot')

genomicTranscription = 1-mt_fraction
fractionTranscriptome = data.frame(mt=mt_fraction,genomic =genomicTranscription,time = c(0.75,2,2.5,3,3.5,4,4.5,5,5.5) )
fractionTranscriptome = melt(fractionTranscriptome,id.vars = c('time'))
fractionTranscriptome$value = fractionTranscriptome$value * 100
fractionTranscriptome$ngenes = ngenesTO$value
    p = ggpubr::ggbarplot(fractionTranscriptome,x='time',y='value',fill='variable',ylab='Transcriptome [%]', label = ngenesTO$value) + theme_cowplot()
    p = p + theme_ameres(type = 'barplot') 
    p
dev.off()


###### now i am looking at the output of MZ and Z genes after classification of the genes.

classified_genes = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/packages/data/countingWindows_classified_conversions.txt",sep="\t",
                              stringsAsFactors = F, header = T)
classified_genes = classified_genes %>%  dplyr::mutate(metadata_add = gene) %>% dplyr::select(c(metadata_add,class))
allPassing = lapply(allPassing,function(x) plyr::join(x, classified_genes))

MZ_time = vector('list',9)
Z_time = vector('list',9)
for(i in 1:9){
  passing_tp = allPassing[[i]][-grep('mt-',allPassing[[i]]$metadata_add),]
  MZ_time[[i]] =passing_tp%>% dplyr::filter(class=="MZ") %>% dplyr::pull(i)
  Z_time[[i]] = passing_tp %>% dplyr::filter(class=="Z") %>% dplyr::pull(i)
  
  MZ = unlist(lapply(MZ_time,function(x) sum(x/0.57)))
  Z = unlist(lapply(Z_time,function(x) sum(x/0.57)))
  
}

sum_denovo = MZ+Z
MZcontribution = data.frame(fraction=1-Z/sum_denovo,type='MZ',value =MZ)
Zcontribution  = data.frame(fraction=Z/sum_denovo ,type='Z',value=Z)

totalDenovo = rbind.data.frame(MZcontribution,Zcontribution)
totalDenovo$time = c(0.75,2,2.5,3,3.5,4,4.5,5,5.5)

nGenes_classified_Z = cbind.data.frame(type='Z',ngenes= unlist(lapply(Z_time,length)))
nGenes_classified_MZ = cbind.data.frame(type='MZ',ngenes= unlist(lapply(MZ_time,length)))
totalDenovo$nGenes = c(nGenes_classified_MZ$ngenes,nGenes_classified_Z$ngenes)

pdf("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure3/plots/transcriptionalOP_Z_MZ.pdf", height=4, width=8)
  ggpubr::ggbarplot(totalDenovo,x='time',y='value',facet.by ='type',label =totalDenovo$nGenes ) + theme_ameres(type = 'barplot') 
  totalDenovo$log10 = log10(totalDenovo$value+1)
  ggpubr::ggbarplot(totalDenovo,x='time',y='log10',facet.by ='type',label =totalDenovo$nGenes ) + theme_ameres(type = 'barplot')
dev.off()

pdf("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure3/plots/transcriptionalOP_Z_MZ_fraction.pdf", height=4, width=3)
    totalDenovo_complete = totalDenovo[complete.cases(totalDenovo),]
    ggpubr::ggbarplot(totalDenovo_complete,x='time',y='fraction',fill='type',label =totalDenovo_complete$nGenes ) + theme_ameres(type = 'barplot')
  
  dev.off()
