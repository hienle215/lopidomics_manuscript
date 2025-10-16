library(ggplot2)
library(mGSZ) #install.packages("mGSZ", repos="http://R-Forge.R-project.org")
library(reshape2)
library(numform)
library(ggpubr)
library(gridExtra)
library(ggplotify)
library(cowplot)
library(purrr)
library(rstatix)
library(biomaRt)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

#setwd("C:/Users/mevado/OneDrive - TUNI.fi/ISE/RNA seq CD new project/our data")

#
Reactome <- read.table("./reviewers comments/ReactomePathways.gmt/ReactomePathways.txt", header = F, sep = '\t')

Pathways_to_show <- c("Intestinal absorption", "Metabolism of lipids", "Peroxisomal lipid metabolism",  "Phospholipid metabolism","Sphingolipid metabolism",
                      "ABC transporters in lipid homeostasis","Regulation of lipid metabolism by PPARalpha","Sphingolipid de novo biosynthesis",
                      "Beta-oxidation of very long chain fatty acids","Intracellular metabolism of fatty acids regulates insulin secretion","Mitochondrial Fatty Acid Beta-Oxidation",
                      "Synthesis of very long-chain fatty acyl-CoAs","The fatty acid cycling model","Transport of fatty acids",
                      "Triglyceride metabolism","Ketone body metabolism") #"Interferon gamma signaling", "TNF signaling",
Reactome1 <- Reactome[Reactome$V1  %in% Pathways_to_show,]

mylist <- list()
for(i in 1:nrow(Reactome1)){
  tmp <- t(Reactome1[i,-c(1:2)])
  annot_human <- getBM(attributes = c("ensembl_gene_id","external_gene_name"),
                       filters = "external_gene_name", values = tmp[,1], mart = human)
  mylist[[i]] <- annot_human$ensembl_gene_id
  names(mylist)[i] <- Reactome1[i,(1)]
}
saveRDS(mylist, file="./LIPIDS_pathways.rds")

mylist <- readRDS("./LIPIDS_pathways.rds")
CD <- readRDS('./Figures for publication/Celiac biopsies DEG results.rds')
ORG <- readRDS('./Figures for publication/hDuo ZED1277 DEG results.rds')

#create theme for all plots
mytheme <- theme_classic() + theme(axis.text.x = element_text(size=8, color = "black"),
                                   axis.text.y= element_text(size=8, color = "black"),
                                   axis.title.y = element_text(size=10, color = "black"),
                                   axis.title.x = element_text(size=10, color = "black"),
                                   legend.text = element_text( size = 6, color = "black"),
                                   legend.position = "top",
                                   legend.title = element_blank(),
                                   legend.key.size = unit(3, 'mm'),
                                   plot.title = element_text(hjust = 0.5, size = 10,face = "bold"),
                                   axis.line = element_line(colour = "black", size=1/.pt),
                                   panel.background = element_rect(fill = 'White'),
                                   panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#--------------------------------------------------------------------------panel B ------------------------------------------------------------------------------------------------------------------------------
CDraw <- readRDS( "./Figures for publication/Celiac biopsies DESeqDataSet object.rds")
tran_data<- log2(counts(CDraw)+1)
metadata <- colData(CDraw)

#function for calculate Sample GSZ
Sample.GSZ <- function(gene.set,ex.data){
  GSZ <- data.frame(Samples=colnames(ex.data))
  for(i in 1:length(gene.set)){
    if(length(gene.set[[i]]) > 1){
      tmp <- sqrt(length(gene.set[[i]])) * ( colMeans( na.omit(ex.data[rownames(ex.data) %in% gene.set[[i]],] )) - colMeans( ex.data ) ) / matrixStats::colSds(as.matrix(ex.data))
    }else{
      tmp <- sqrt(length(gene.set[[i]])) * ( ex.data[rownames(ex.data) %in% gene.set[[i]],] - colMeans( ex.data ) ) / matrixStats::colSds(as.matrix(ex.data))
    }
    names(tmp) <- colnames(ex.data)
    tmp <- as.data.frame(tmp)
    colnames(tmp) <- names(gene.set)[i]
    GSZ <- cbind(GSZ,tmp)
  }
  return(GSZ)
  
}

#calculate Sample GSZ

SampleGSZ <- Sample.GSZ(gene.set=mylist,ex.data=tran_data)
GSZ_melt <- melt(SampleGSZ)
colnames(GSZ_melt)[2] <- "gene.sets"
GSZ_melt <- merge(GSZ_melt,metadata, by.x="Samples", by.y = "ID")
#write_xlsx(GSZ_melt, path="./SampleGSZ 3 IFNg gene sets.xlsx")
#abbr <- "IFNg gene sets"
#calculate groupGSZ for groups
groups <- data.frame(group1=c("GFD","GFD","PGCd"),
                     group2=c("PGCd","PGCp","PGCp"))
set.seed(10)
groupGSZ <- data.frame()
for(i in 1:nrow(groups)){
  groups2 <- metadata[metadata$label3 %in% c(groups[i,1],groups[i,2]),]$label3
  tran_data2 <- tran_data[,metadata[metadata$label3 %in% c(groups[i,1],groups[i,2]),]$ID]
  mGSZ.obj2 <- mGSZ(tran_data2, mylist, groups2, p = 100)
  
  group1 <- c(rep(groups[i,1],length(mylist)))
  group2 <- c(rep(groups[i,2],length(mylist)))
  comparison <- paste(group1,"-", group2, sep="")
  tmp <- cbind(toTable(mGSZ.obj2, n = length(mylist)), group1, group2, comparison)
  groupGSZ <- rbind(groupGSZ, tmp)
}
groupGSZ <- na.exclude(groupGSZ)
groupGSZ %>% 
  mutate(p_new = ifelse(`pvalue` > 0.01, c(paste("italic('P')~`=", f_num(`pvalue`,2), "`")), `pvalue`))%>% 
  mutate(p_new = ifelse(`pvalue` < 0.01, c(paste("italic('P')~`=", f_num(`pvalue`,3), "`")), p_new)) %>%
  mutate(p_new = ifelse(`pvalue` < 0.001, c(paste("italic('P')~`", "<.001", "`")),p_new))->groupGSZ

df <-   as.data.frame(GSZ_melt) %>% group_by(gene.sets, label3) %>% 
  dplyr::summarise(max=max(value))
groupGSZ <- merge(groupGSZ,df, by.x=c("gene.sets", "group1"), by.y=c("gene.sets", "label3"))

groupGSZ$y.position <- rep(c(1.3,1.5,1.7), length(unique(groupGSZ$gene.sets)))
groupGSZ1 <- groupGSZ %>% add_xy_position(data = groupGSZ,formula = max ~ group1)



bp<-ggboxplot(as.data.frame(GSZ_melt), x = "label3", y = "value", outlier.colour = NA, order = c("GFD","PGCd","PGCp"),
              palette = c("#00b347","#6500ff","#ff6500"),
              fill="label3")+
  labs(y = "gene set GSZ")+
  geom_jitter(size=0.5, alpha=0.5)+
  geom_hline(yintercept=c(mean(GSZ_melt$value)), linetype="dashed", color = "grey60")+
  facet_wrap(.~gene.sets, scales = "fixed")+
  mytheme+
  #stat_pvalue_manual(groupGSZ1, label = "p_new", tip.length = 0.01,parse=T)
  geom_signif(data=groupGSZ,
              aes(xmin=group1, xmax=group2, annotations=p_new, y_position=y.position),
              tip_length = 0.005,
              textsize = 6/.pt,
              manual=TRUE, parse=T, size=0.3)

B <- bp + mytheme+
  theme( legend.position = "none", axis.title.x = element_blank())
B


ggexport(B,filename = "./Figures for publication/Suppl to Figure 3.tiff",
         width = 2000, # 7 inch
         height = 2700, # 9 inch
         res = 300) #max 2100 X 2700

