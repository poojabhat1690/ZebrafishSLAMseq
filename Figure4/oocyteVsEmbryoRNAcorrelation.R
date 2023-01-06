##### plotting gene expression between oocyte and tp1... 
###### script plotting total reads vs tc reads to show that for a set of pre-defined purely zygotic genes there most reads have tc

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

#### oocyte data


ActivatedOocytes = list.files("///Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure3/data/oocyteData/",pattern = "WT_Ooc_FER")
ActivatedOocytes = ActivatedOocytes[grep('.tsv',ActivatedOocytes)]
ActivatedOocytes = paste0("///Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure3/data/oocyteData/",ActivatedOocytes)
ActivatedOocytes_data  = lapply(ActivatedOocytes,function(x) read.table(x,stringsAsFactors = F, header = T, sep="\t"))
names(ActivatedOocytes_data) = c('R1','R2','R3')
ActivatedOocytes_data_RPM  = do.call(cbind.data.frame,lapply(ActivatedOocytes_data,function(x) x$ReadsCPM))
ActivatedOocytes_data_RPM$name = ActivatedOocytes_data$R1$Name

ActivatedOocytes_data = ActivatedOocytes_data_RPM %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(R1:R3), sum, na.rm = TRUE) %>%
  dplyr::mutate(mean = (R1+R2+R3)/3)
ActivatedOocytes_data  = ActivatedOocytes_data %>% dplyr::select(c('name','mean'))

RPM = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure2/data//RPM_allCws.txt",
                 sep="\t", stringsAsFactors = F, header = T)
RPM = RPM %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Untreated_TP9), sum, na.rm = TRUE)

RMP_oocyte_embryos = plyr::join(RPM, ActivatedOocytes_data,by='name')

classifiedGenes = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/packages/data//countingWindows_classified_conversions.txt",
                             sep="\t",stringsAsFactors = F,header = T)
classifiedGenes=  classifiedGenes %>% dplyr::mutate(name = gene) %>% dplyr::select(c('name','class'))
classifiedGenes = plyr::join(classifiedGenes,RMP_oocyte_embryos)


classifiedGenes$T1 = (classifiedGenes$Inj_R2_TP1 + classifiedGenes$Inj_R3_TP1)/2
classifiedGenes = classifiedGenes %>% dplyr::mutate(oocyteRPM_FER = mean) %>% dplyr::select(c('oocyteRPM_FER','T1','class'))
write.table(classifiedGenes, "//groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/paper/fig2/data/oocyteVst1.txt",sep="\t",row.names = F)

Z = classifiedGenes %>% 
  dplyr::filter(class == 'Z')
MZ = classifiedGenes %>% 
  dplyr::filter(class == 'MZ')

M = classifiedGenes %>% 
  dplyr::filter(class == 'M')

M_stable = classifiedGenes %>% 
  dplyr::filter(class == 'M-stable')

classifiedGenes = classifiedGenes %>% dplyr::mutate(log10T1 =log10(T1+1),log10Oocyte = log10(oocyteRPM_FER+1) )

pdf("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure3/plots//scatterplot_oocyteVsembryo.pdf", height = 4.5, width = 14)
    
    Z = ggplot(classifiedGenes %>% dplyr::filter(class=='Z'),aes(x=log10T1,y=log10Oocyte)) +geom_density_2d() + ylab("RPM fertilized oocytes (log10)") + 
      xlab ('RPM 0.75hpf') + theme_cowplot()+ ggpubr::stat_cor(aes(color = class))+ annotate('text',x = 4,y=4,label=nrow(classifiedGenes %>% dplyr::filter(class=='Z'))) 
    Z = Z + theme_ameres(type = 'barplot') + ylim(c(-0.5,5)) + xlim(c(-0.5,5))+ coord_fixed() + geom_abline(linetype='dashed')

    M = ggplot(classifiedGenes %>% dplyr::filter(class=='M'),aes(x=log10T1,y=log10Oocyte)) + geom_density_2d() + ylab("RPM fertilized oocytes (log10)") + 
      xlab ('RPM 0.75hpf') + theme_cowplot()+ ggpubr::stat_cor(aes(color = class))+ annotate('text',x = 4,y=4,label=nrow(classifiedGenes %>% dplyr::filter(class=='M'))) 
    M = M + theme_ameres(type = 'barplot') + ylim(c(-0.5,5)) + xlim(c(-0.5,5))+ coord_fixed() + geom_abline(linetype='dashed')
    
    MZ = ggplot(classifiedGenes %>% dplyr::filter(class=='MZ'),aes(x=log10T1,y=log10Oocyte)) + geom_density_2d() + ylab("RPM fertilized oocytes (log10)") + 
      xlab ('RPM 0.75hpf') + theme_cowplot()+ ggpubr::stat_cor(aes(color = class))+ annotate('text',x = 4,y=4,label=nrow(classifiedGenes %>% dplyr::filter(class=='MZ'))) 
    MZ = MZ + theme_ameres(type = 'barplot') + ylim(c(-0.5,5)) + xlim(c(-0.5,5))+ coord_fixed() + geom_abline(linetype='dashed')
    
    Mstable = ggplot(classifiedGenes %>% dplyr::filter(class=='M-stable'),aes(x=log10T1,y=log10Oocyte)) + geom_density_2d() + ylab("RPM fertilized oocytes (log10)") + 
      xlab ('RPM 0.75hpf') + theme_cowplot()+ ggpubr::stat_cor(aes(color = class))+ annotate('text',x = 4,y=4,label=nrow(classifiedGenes %>% dplyr::filter(class=="M-stable"))) 
    Mstable = Mstable + theme_ameres(type = 'barplot') + ylim(c(-0.5,5)) + xlim(c(-0.5,5))+ coord_fixed() + geom_abline(linetype='dashed')
    
  
       M + Mstable + MZ+Z+ plot_layout(nrow =1)
 dev.off()
# 
# pdf("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/packages/plots/figure2/boxplot_oocyteVsembryo.pdf", height = 5, width = 5)
# 
#   q = ggpubr::ggboxplot(classifiedGenes[order(classifiedGenes$class),],x='class',y='log10T1',fill='class', palette = 'Set1') + theme_ameres(type = 'barplot')
#   print(q)  
# 
#   q = ggpubr::ggboxplot(classifiedGenes[order(classifiedGenes$class),],x='class',y='log10Oocyte',fill='class', palette = 'Set1') + theme_ameres(type = 'barplot')
#   print(q)  
#   
# dev.off()