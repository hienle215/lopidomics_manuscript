#setwd("C:/Users/mevado/OneDrive - TUNI.fi/ISE/Proteomics + Metabolomics/Lipidomics")
library(readxl)
library(plotly)
#library(LipidSigR)
library(dplyr)
library(ggpubr)
library(factoextra)
library(DEGreport)
library(tibble)
library("circlize")
library(ComplexHeatmap)
library(ggrepel)
library(gridExtra)
library(cowplot)


#create theme for all plots
mytheme <- theme_classic(base_family='sans') + theme(axis.text.x = element_text(size=6, color = "black"),
                                                     axis.text.y= element_text(size=6, color = "black"),
                                                     axis.title.y = element_text(size=8, color = "black"),
                                                     axis.title.x = element_text(size=8, color = "black"),
                                                     legend.text = element_text( size = 4, color = "black"),
                                                     legend.key.size = unit(2, 'mm'),
                                                     legend.position = "top",
                                                     legend.title = element_blank(),
                                                     plot.title = element_text(hjust = 0.5, size = 8,face = "bold"),
                                                     axis.line = element_line(colour = "black", linewidth=1/.pt),
                                                     strip.text.x = element_text(size = 6),
                                                     panel.background = element_rect(fill = 'White'),
                                                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                     axis.title.x.top = element_text(color='#4daf4a',size=7), 
                                                     axis.text.x.top = element_text(color='#4daf4a',size=6))

palette <- c(GFDd = "#00b347", GFDp = "#005923",PGCd = "#6500ff", PGCp = "#ff6500")
hulls_p2 <- c("#ff89b2","#89edff")
stcat_lipids_palette <- c(FA = "#ae3c60", 
                         GP = "#df473c",
                         ST = "#f3c33c",
                         SP = "#255e79",
                         GL = "#267778",
                         PR = "#82b4bb")


#open and transfrom data
raw_data <- readxl::read_xlsx("./Input_data/Lipids+Semipolar_metabolomics_raw.xlsx")
meta <- readxl::read_xlsx("./Input_data/proteome metadata.xlsx")
meta$`Sample_ID` <- meta$`Sample ID` 

lipids <- readxl::read_xlsx("./Input_data/lipids_characteristics.xlsx")

raw_data <- raw_data[raw_data$feature %in% lipids$feature,]
raw_data <- raw_data %>% mutate_at(meta$`Sample ID`, as.numeric)

raw_data[raw_data == 0] <- min(unlist(raw_data[,-c(1,96)])[unlist(raw_data[,-c(1,96)])>0])
row.names(raw_data) <- raw_data$feature

data_tr <- raw_data %>% 
  mutate_if(is.numeric,  function(x, na.rm = FALSE) log10(x / sum(x)*100))

data_tr2 <- as.data.frame(t(data_tr))
colnames(data_tr2) <-  raw_data$feature
data_tr2 <- data_tr2[-c(1,96),] %>% 
  mutate_all(as.numeric)

#----------------------------------------------------------------------------Density plot----------------------------------------------------------------------------------
#data_tr <- rownames_to_column(data_tr, "feature")
data_tr_melt <- reshape2::melt(data_tr, id = c("feature","Name"))
colnames(data_tr_melt)[3] <-"Sample_ID" 
data_tr_melt <- data_tr_melt %>% left_join(meta[, c("Sample_ID","label2", "BMI")], by = "Sample_ID") %>% right_join(lipids[, c("feature","Abbr", "structural_category", "totallength","totaldb")], by = "feature")

A <- ggplot(data_tr_melt, aes(x=value, color=label2)) +
  geom_density(alpha=0.6) +
  color_palette(palette) +
  ylab("Density") +
  xlab("log10(expression)")+mytheme
A

#Number of lipids indentified
Numbers <- lipids %>% group_by(structural_category) %>%summarise(n=n())

A1 <- ggplot(Numbers, aes(x = n, y = structural_category, fill = structural_category)) + fill_palette(stcat_lipids_palette)+
  geom_col() + 
  coord_flip()+ 
  labs(x= "Number \n of indentified lipids")+
  mytheme + theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none")
A1
#---------------------------------------------------------------------------lipid composition by class-----------------------------------------------------------------------------
percent_comp <- raw_data %>% 
  mutate_if(is.numeric,  function(x, na.rm = FALSE) x / sum(x)*100)

percent_comp <- reshape2::melt(percent_comp, id = c("feature","Name"))
colnames(percent_comp)[3] <-"Sample_ID" 
percent_comp <- percent_comp %>% full_join(meta[, c("Sample_ID","label2", "BMI")], by = "Sample_ID") %>% right_join(lipids[, c("feature","Abbr", "structural_category", "totallength","totaldb", "Sub_Class")], by = "feature")
percent_comp$value <- as.numeric(percent_comp$value)

pc1 <- ggplot(percent_comp, aes(x = Sample_ID, y = value, fill = structural_category)) + fill_palette(stcat_lipids_palette)+
  geom_col() + facet_wrap(~factor(label2, c("GFDd", "PGCd", "GFDp", "PGCp")), scales = "free_x", ncol=4) +
  #coord_flip()+ 
  labs(x= "Sample", y = "Lipid Composition, %")+
  mytheme + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+ guides(fill = guide_legend(nrow = 1))
pc1
#percent_comp <- raw_data %>% 
#  mutate_if(is.numeric,  function(x, na.rm = FALSE) x / sum(x)*100)

sum_by_group <- data.frame(feature = raw_data$feature)
for(i in unique(meta$label2)){
  tmp <- raw_data[,c("feature",meta[meta$label2==i,]$Sample_ID)]
  tmp2 <- as.data.frame(rowSums(tmp[,-1]))
  colnames(tmp2) <- i
  sum_by_group <- cbind(sum_by_group, tmp2)
}

sum_by_group_percent <- sum_by_group %>% 
  mutate_if(is.numeric,  function(x, na.rm = FALSE) x / sum(x)*100)
sum_by_group_percent <- sum_by_group_percent %>% left_join(lipids[, c("feature", "structural_category")], by = "feature")
sum_by_group_percent <- reshape2::melt(sum_by_group_percent)
sum_by_group_percent <- sum_by_group_percent %>% group_by(variable, structural_category) %>% summarise(total = sum(value))

sum_by_group_percent <- sum_by_group_percent %>% 
  mutate( cs = rev(cumsum(rev(total))), 
          pos = total/2 + lead(cs, 1),
          pos = if_else(is.na(pos), total/2, pos))
sum_by_group_percent$structural_category <- factor(sum_by_group_percent$structural_category, levels= unique(sum_by_group_percent$structural_category))

pc2 <- ggplot(sum_by_group_percent, aes(x = "" , y = total, fill = structural_category)) +
  geom_col(width = 1) +
  coord_polar(theta = "y", start = 0) +
  scale_fill_manual(values = stcat_lipids_palette) +guides(fill = guide_legend(title = "Status")) + facet_wrap(~~factor(variable, c("GFDd", "PGCd", "GFDp", "PGCp")), ncol = 4)+
  geom_text_repel(aes(y = pos, label = paste0(structural_category,"\n",round(total,1), "%")), data = sum_by_group_percent, size=2, show.legend = F, nudge_x = 1) +
  mytheme +  theme_void() + theme(legend.position = "none",strip.text.x = element_blank())
pc2

B <- ggarrange(pc1, pc2, ncol = 1, align = c("v"), heights = c(2,1))

B


#--------------------------------------------------------------------------PCA analysis------------------------------------------------------------------------------------------------
res.pca <- prcomp(data_tr2, scale = TRUE)
summary(res.pca)
resultNames <- c("PC1","PC2", "PC3")
df <- as.data.frame(res.pca$x[,c("PC1","PC2", "PC3")])
df$name <- rownames(df)

meta$BMI_group2 <- "BMI 18.7-26 "
meta[meta$BMI > 26,]$BMI_group2 <- "BMI > 26 "

df <- merge(df, meta[,c("Sample ID","label2", "Haplotype_G", "BMI_group2", "BMI", "Gender", "Country", "IELs")], by.x="name", by.y="Sample ID")
hulls <- plyr::ddply(df, "BMI_group2", function(df) df[chull(df$PC1, df$PC2),])
hulls2 <- plyr::ddply(df, "Gender", function(df) df[chull(df$PC2, df$PC3),])




C <- ggplot()+
  geom_polygon(data = hulls,aes(PC1,PC2,group=BMI_group2, fill=BMI_group2), alpha = 0.25)+
  #geom_polygon(data = hulls2,aes(PC2,PC3,group=Gender, fill=Gender), alpha = 0.25)+
  geom_point(data = df,aes(PC1, PC2,color=label2,
                           text = paste(
                             "ID: ", meta$newID, "\n",
                             "Group: ", meta$Treatment, "\n",
                             "Timepoint: ", meta$Timepoint, "\n",
                             "Pair #: ", meta$Pair, "\n",
                             "VH:CrD: ",meta$VHCrD, "\n",
                             "BMI: ",meta$BMI, "\n",
                             sep = "")), size=1)+
  scale_color_manual(values=palette)+
  scale_fill_manual(values=hulls_p2) #+
#labs(x=sprintf("PC1: %s Variance", percent(pca_obja$plot_env$percentVar[1])),
#     y=sprintf("PC2: %s Variance", percent(pca_obja$plot_env$percentVar[2])) ) 
C <- C + mytheme + theme(legend.box="vertical", legend.margin=margin()#,plot.margin = margin(0, 1.5, 1, 1.5, "cm")) #top clockwise
                            )
C

#ggplotly(A, tooltip = "text")

#PC correlations
#all genes



data_tr <- data_tr[,c("feature", meta$`Sample ID`)]
data_tr <- column_to_rownames(data_tr, var = "feature")
meta <- column_to_rownames(meta, var = "Sample ID")



res <- degCovariates(data_tr,meta[, c("Age", "Gender", "Country","BMI")],plot = T)
cov <- res$corMatrix#[res$corMatrix$pvalue < 0.05,]#[1:4,1:5]
cov <- cov[order(cov$compare),]
cov_cast <- reshape2::dcast(cov[, 1:3], covar ~ compare)
writexl::write_xlsx(cov_cast, path = "./regression of principal components.xlsx" )

#---------------------------------------------------------------------------Top features-----------------------------------------------------------------------------------------
aload <- abs(res.pca$rotation)[, c("PC1","PC2","PC3")]
aload <- as.data.frame(aload) %>% mutate_if(is.numeric,function(x, na.rm = FALSE) x / sum(x)*100)
aload <-  rownames_to_column(aload, var = "feature")
aload <- merge(aload, lipids[,c( "feature","Abbr", "structural_category")], by ="feature")

aload <- reshape2::melt(aload)

top6 <- aload %>% group_by(variable) %>% top_n(6) %>% arrange(variable, desc(value))

D <- ggbarplot(top6,"Abbr","value", fill = "structural_category", palette = stcat_lipids_palette) + 
  facet_wrap(~variable, scales ="free_x") +labs(y = "Contribution, %") +
  mytheme +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_blank())

#----------------------------------------------------------------------------combine all together ---------------------------------------------------------------------------------
filler <- ggplot()
p <- arrangeGrob(A,A1, B,C,D,
                 ncol = 2, nrow = 4, layout_matrix = rbind(c(1,2),
                                                           c(3,3),
                                                           c(3,3),
                                                           c(4,5)
                                                           ))
as_ggplot(p)


# Add labels to the arranged plots
p <- as_ggplot(p) +    
  draw_plot_label(label = c("A", "B", "C","D", "E"), size = 16,
                  x = c(0, 0.5, 0,  0, 0.5), 
                  y = c(1, 1,  0.75, 0.27,0.27)) # transform to a ggplot


p
ggexport(p,filename = "./Figures for publication/Figure 1.tiff",
         width = 2100, # 7 inch max width 2100
         height = 2100, 
         res = 300) # Figure panels should be prepared at a minimum resolution of 300 dpi and saved at a maximum width of 180 mm (7,08661 inch).  




#---------------------------------------------------------------------------Supplemental material-------------------------------------------------------------
#PC1 correlation
SM1 <- ggscatter(df, x = "BMI", y = "PC1",
          color = "label2", palette  = palette,  size = 1,# Points color, shape and size
          add = "reg.line",  # Add regressin line
          #add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          #conf.int = TRUE, # Add confidence interval
          #cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.x = 25, label.sep = "\n") 
)+ mytheme
ggexport(SM1,filename = "./Figures for publication/Figure S2.tiff",
         width = 1200, # 7 inch max width 2100
         height = 1200, 
         res = 300) # Figure panels should be prepared at a minimum resolution of 300 dpi and saved at a maximum width of 180 mm (7,08661 inch).
#PC2 correlation
SM2 <- ggscatter(df, x = "GGT_U_L", y = "PC2",
          color = "label2", palette  = palette,  # Points color, shape and size
          add = "reg.line",  # Add regressin line
          #add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          #conf.int = TRUE, # Add confidence interval
          #cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.x = 25, label.sep = "\n")
)


#top5 correlation with PCs
#from PC1
expr <- data_tr[data_tr$feature %in% top6[top6$variable == "PC1",]$feature,]
#expr <- rownames_to_column(expr, var = "feature")
expr_melt <- reshape2::melt(expr)
colnames(expr_melt)[2] <- "Sample_ID"

expr_melt <- expr_melt %>% left_join(meta[, c("Sample_ID","label2", "BMI")], by = "Sample_ID") %>% left_join(lipids[, c("feature","Abbr")], by = "feature")

SM3 <- ggscatter(expr_melt, x = "BMI", y = "value",
          color = "label2", palette  = palette,  # Points color, shape and size
          add = "reg.line", ylab = "log10 normalized concentration") + 
  facet_wrap(~Abbr , scales ="free") + mytheme


#from PC2
expr <- data_tr[data_tr$feature %in% top6[top6$variable == "PC2",]$feature,]
#expr <- rownames_to_column(expr, var = "feature")
expr_melt <- reshape2::melt(expr)
colnames(expr_melt)[2] <- "Sample_ID"

expr_melt <- expr_melt %>% left_join(meta[, c("Sample_ID","label2", "GGT_U_L")], by = "Sample_ID") %>% left_join(lipids[, c("feature","Abbr")], by = "feature")

ggscatter(expr_melt, x = "GGT_U_L", y = "value",
          color = "label2", palette  = palette,  # Points color, shape and size
          add = "reg.line", ylab = "log10 normalized concentration") + 
  facet_wrap(~Abbr , scales ="free") + mytheme


#---------------------------------------------------------------------------CORRELATION HEATMAP-------------------------------------------------------------------------------------
#1 Sample to Sample

cormat <- round(cor(data_tr, method = c("spearman")),2)

type1 <- meta$label2
col_type1 = list(`label` = palette)
top = HeatmapAnnotation(Samples = anno_block(gp = gpar(fill = palette, lty="blank")),
                        annotation_name_gp= gpar(fontsize = 6),
                        gap = unit(1, "mm"),
                        annotation_height = unit(1:1, "cm"))

row=rowAnnotation(`label` = as.character(meta$label2), col=col_type1, 
                  show_annotation_name = F,
                  annotation_legend_param = list(`label` = list(title_position = "topcenter",labels_gp = gpar(fontsize = 6),
                                                                title_gp = gpar(fontsize = 5),direction = "horizontal",nrow = 1,
                                                                legend_width=unit(1,"cm"),grid_height = unit(0.2, "cm"))),
                  
                  gap = unit(1, "mm"),
                  annotation_height = unit(1:1, "cm"))
 
#create color scheme for heatmap
newcolors <- c(#"#67001F",
  "#B2182B",
  #"#D6604D", 
  #"#F4A582",
  #"#FDDBC7", 
  "#F7F7F7", #white 
  #"#D1E5F0", 
  #"#92C5DE",
  #"#4393C3", 
  "#2166AC"  #, 
  #"#053061"
)

ht_list = Heatmap(cormat, name = "Pearson\nCorrelation", top_annotation = top, 
                  show_column_names = F,
                  show_row_names = F, 
                  col=newcolors,
                  column_order= meta$Sample_ID,
                  row_order= meta$Sample_ID,
                  cluster_rows = F,  
                  cluster_columns = F, 
                  column_split = factor(type1, levels = c("GFDd","GFDp", "PGCd","PGCp")),
                  column_gap = unit(0.5,"mm"),
                  row_split = factor(type1, levels = c("GFDd","GFDp", "PGCd","PGCp")),
                  row_title = NULL,
                  row_gap = unit(2,"mm"),
                  column_title_gp = gpar(fontsize = 0.5, col="White"),
                  heatmap_legend_param = list(
                    at = c(0:1),
                    labels = c(0:1),
                    title = "Pearson\nCorrelation",
                    labels_gp = gpar(fontsize = 5),
                    title_gp = gpar(fontsize = 5),direction = "horizontal",
                    legend_width=unit(2,"cm"),grid_height = unit(2, "mm"),
                    title_position = "topcenter"))+
  row
#draw(ht_list,merge_legend = TRUE,heatmap_legend_side = "bottom", annotation_legend_side = "bottom", row_gap = unit(0.5, "mm"),auto_adjust = T)
H <- grid.grabExpr(draw(ht_list,merge_legend = TRUE,heatmap_legend_side = "bottom", annotation_legend_side = "bottom", row_gap = unit(0.5, "mm"),auto_adjust = T))
#hm1
as_ggplot(H)


#2 Lipid to Lipid

cormat2 <- round(cor(data_tr2, method = c("spearman")),2)

type2 <- lipids$structural_category
col_type2 = list(`Category` = stcat_lipids_palette)
top = HeatmapAnnotation(Category = anno_block(gp = gpar(fill = stcat_lipids_palette, lty="blank")),
                        annotation_name_gp= gpar(fontsize = 6),
                        gap = unit(1, "mm"),
                        annotation_height = unit(1:1, "cm"))

row=rowAnnotation(`Category` = as.character(lipids$structural_category), col=col_type2, 
                  show_annotation_name = F,
                  annotation_legend_param = list(`Category` = list(title_position = "topcenter",labels_gp = gpar(fontsize = 6),
                                                                title_gp = gpar(fontsize = 5),direction = "horizontal",nrow = 1,
                                                                legend_width=unit(1,"cm"),grid_height = unit(0.2, "cm"))),
                  
                  gap = unit(1, "mm"),
                  annotation_height = unit(1:1, "cm"))

#create color scheme for heatmap
newcolors <- c(#"#67001F",
  "#B2182B",
  "#D6604D", 
  "#F4A582",
  "#FDDBC7", 
  "#F7F7F7", #white 
  "#D1E5F0", 
  #"#92C5DE",
  #"#4393C3", 
  "#2166AC" #, 
  #"#053061"
)

ht_list = Heatmap(cormat2, name = "Pearson\nCorrelation", top_annotation = top, 
                  show_column_names = F,
                  show_row_names = F, 
                  col=newcolors,
                  column_order= lipids$feature,
                  row_order= lipids$feature,
                  cluster_rows = F,  
                  cluster_columns = F, 
                  column_split = factor(type2, levels = c("FA", "GP", "ST", "SP", "GL", "PR" )),
                  column_gap = unit(0.5,"mm"),
                  row_split = factor(type2, levels = c("FA", "GP", "ST", "SP", "GL", "PR")),
                  row_title = NULL,
                  row_gap = unit(2,"mm"),
                  column_title_gp = gpar(fontsize = 0.5, col="White"),
                  heatmap_legend_param = list(
                    at = c(0:1),
                    labels = c(0:1),
                    title = "Pearson\nCorrelation",
                    labels_gp = gpar(fontsize = 5),
                    title_gp = gpar(fontsize = 5),direction = "horizontal",
                    legend_width=unit(2,"cm"),grid_height = unit(2, "mm"),
                    title_position = "topcenter"))+
  row
#draw(ht_list,merge_legend = TRUE,heatmap_legend_side = "bottom", annotation_legend_side = "bottom", row_gap = unit(0.5, "mm"),auto_adjust = T)
H2 <- grid.grabExpr(draw(ht_list,merge_legend = TRUE,heatmap_legend_side = "bottom", annotation_legend_side = "bottom", row_gap = unit(0.5, "mm"),auto_adjust = T))
#hm1
as_ggplot(H2)





