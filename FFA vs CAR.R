library(dplyr)
meta <- readxl::read_xlsx("./Input_data/proteome metadata.xlsx")

lipids <- readxl::read_xlsx("./Input_data/lipids_characteristics.xlsx")
raw_data <- readxl::read_xlsx("./Input_data/Lipids_metabolomics_raw.xlsx")
raw_data <- raw_data %>% mutate_at(meta$`Sample ID`, as.numeric)
#=============================================================================CAR vs FFA============================================================================

t_raw <- t(raw_data)
colnames(t_raw) <- t_raw[96,]
t_raw <- t_raw[-c(1,96),]

t_raw <- as.data.frame(t_raw) %>% mutate_all(as.numeric)
 
ratios <- data.frame(`Sample_ID` = row.names(t_raw))
for (i in c("2:0", "3:0","4:0", "6:0", "8:0", "14:0" , "16:0", "18:0", "18:1",  "18:2")){
  tmp_FA <- na.omit(lipids[lipids$Abbr == paste("FA", i, sep = " "),]$feature)
  tmp_Car <- unique(na.omit(lipids[lipids$Abbr == paste("CAR", i, sep = " "),]$feature))
  ratio <- t_raw[, tmp_FA]/t_raw[, tmp_Car]
  name <- paste("RATIO", i, sep = "_")
  tmp <- cbind(row.names(t_raw), t_raw[, tmp_FA],t_raw[, tmp_Car], ratio)
  colnames(tmp) <- c("Sample_ID", tmp_FA, tmp_Car,name)
  ratios <- merge(ratios, tmp, by = "Sample_ID")
}
ratios <- data.frame(`Sample ID` = row.names(t_raw))
ratios$ratioCar <- t_raw$acetylcarnitine/(t_raw$Palmitoylcarnitine + t_raw$octadecenoylcarnitine__C18_1_Car_)
ratios$ratio18_2 <- t_raw$C18_2/t_raw$octadecadienylcarnitine__C18_2_Car_
ratios$ratio14 <- t_raw$myristic_acid/t_raw$Myristoyl_L_carnitine
ratios$ratio18_1 <- t_raw$oleic_acid/t_raw$octadecenoylcarnitine__C18_1_Car_
ratios$ratio16 <- t_raw$palmitic_acid/t_raw$Palmitoylcarnitine
ratios <- reshape2::melt(ratios, id = "Sample.ID")

ratios <- merge(ratios, meta, by.x  = "Sample.ID", by.y = "Sample ID")



ggboxplot(ratios, x = "label2", y = "value", color = "label2", facet.by = "variable", scales = "free_y",
          palette = palette, add = "jitter", ylab = "RATIO", size = 0.1, add.params = list(size = 0.1))

annotate("rect", xmin = -Inf, xmax = Inf, ymin = 60, ymax = 150, 
         alpha = 0.2,fill = "#ccffd4")+
  geom_signif(data=stat_comp %>% filter(variable=="GFR_mL_min"), 
              aes(xmin=group1, xmax=group2, annotations=p_new, y_position=y.position),
              textsize = 6/.pt, 
              manual=TRUE, parse=T, size=0.3)+
  mytheme + theme(axis.title.x = element_blank())