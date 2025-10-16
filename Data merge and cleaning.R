#open metabolisms raw data and combine\
meta <- readxl::read_xlsx("./Input_data/proteome metadata.xlsx")
PATH = "C:/Users/dotsd/OneDrive - TUNI.fi/"

raw_data_1 <- readxl::read_xlsx(paste0(PATH,"ISE/Proteomics + Metabolomics/Metabolomics/Copy Final_combined_results_Viiri_2022_12_13.xlsx"), sheet = 1)
raw_data_1 <- raw_data_1[-c(1:2),-c(1:3,5)]
raw_data_1 <- raw_data_1 %>% mutate_at(meta$`Sample ID`, as.numeric)
names = make.names(raw_data_1$...4, unique=TRUE)
raw_data_1_m <- as.matrix(raw_data_1[,-c(1)])
rownames(raw_data_1_m) <- names

raw_data_1_m[is.na(raw_data_1_m) == T] <- 0
raw_data_1_m[raw_data_1_m == 0] <- min(unlist(raw_data_1_m)[unlist(raw_data_1_m)>0])

z_scores_1 <- (raw_data_1_m-mean(raw_data_1_m))/sd(raw_data_1_m)

raw_data_2 <- readxl::read_xlsx(paste0(PATH,"ISE/Proteomics + Metabolomics/Metabolomics/Copy Final_combined_results_Viiri_2022_12_13.xlsx"), sheet = 2)
raw_data_2 <- raw_data_2[-c(1:2),-c(1:3,5)]

raw_data_3 <- readxl::read_xlsx(paste0(PATH,"ISE/Proteomics + Metabolomics/Metabolomics/Copy Final_combined_results_Viiri_2022_12_13.xlsx"), sheet = 3)
raw_data_3 <- raw_data_3[-c(1:2),-c(1:3,5)]

raw_data_3 <- raw_data_3 %>% mutate_at(meta$`Sample ID`, as.numeric)
names = make.names(raw_data_3$...4, unique=TRUE)
raw_data_3_m <- as.matrix(raw_data_3[,-c(1)])
rownames(raw_data_3_m) <- names

raw_data_3_m[is.na(raw_data_3_m) == T] <- 0
raw_data_3_m[raw_data_3_m == 0] <- min(unlist(raw_data_3_m)[unlist(raw_data_3_m)>0])

z_scores_3 <- (raw_data_3_m-mean(raw_data_3_m))/sd(raw_data_3_m)

raw_data_4 <- readxl::read_xlsx(paste0(PATH,"ISE/Proteomics + Metabolomics/Metabolomics/Copy Final_combined_results_Viiri_2022_12_13.xlsx"), sheet = 4)
raw_data_4 <- raw_data_4[-c(1:2),-c(1:3,5)]

raw_data_5 <- readxl::read_xlsx(paste0(PATH,"ISE/Proteomics + Metabolomics/Metabolomics/Copy Final_combined_results_Viiri_2022_12_13.xlsx"), sheet = 5)
raw_data_5 <- raw_data_5[-c(1:2),-c(1:3,5)]

rawdata <- rbind(raw_data_1,raw_data_2,raw_data_3,raw_data_4,raw_data_5)
colnames(rawdata)[1] <- "Name"

lipids <- readxl::read_xlsx("./Input_data/lipids_characteristics.xlsx")
rawdata <- rawdata[rawdata$Name %in% lipids$name,]

rawdata$feature <- gsub("[-():/ ]", "_", rawdata$Name)
rawdata$feature = make.names(rawdata$feature, unique=TRUE)
rawdata <- rawdata %>% mutate_at(meta$`Sample ID`, as.numeric)

rawdata[is.na(rawdata) == T] <- 0
rawdata[raw_data == 0] <- min(unlist(rawdata[,-c(1,96)])[unlist(rawdata[,-c(1,96)])>0])

writexl::write_xlsx(rawdata, path = "./Input_data/Lipids_metabolomics_raw.xlsx")


data <- c(8, 7, 7, 10, 13, 14, 15, 16, 18) 
z_scores <- (data-mean(data))/sd(data)
z_scores
