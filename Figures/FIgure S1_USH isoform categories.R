
#Load packages
library(dplyr) 
library(readxl) 
library(xlsx) 
library(readr)
library(data.table)
library(splitstackshape)
library(ggplot2)


#Set your working directory
setwd("/Users/erikdevrieze/Library/CloudStorage/OneDrive-Radboudumc/z918116/UMCN/Manuscripts/2024 - USH isoseq manuscript/Analysis paper 2024")




###################################################################################################################
# quantification done from isoquant read assignment file containg all 4 samples. 
# script was also done on file for samples 1-3 confirming that results are the same between files. 
###################################################################################################################

#Read full data excel file containing read assignments from all 4 runs. 
s1_4_sorted_gene_assignments <- fread("dataset 1_CCS3/isoquant/00_aln.sorted.read_assignments.tsv") 

# change col name to "run"
colnames(s1_4_sorted_gene_assignments)[1] ="run"                                                                                               
                                                                                                 
# filter for target genes (featureID)
USHs1_4_sorted_gene_assignments <- s1_4_sorted_gene_assignments %>% dplyr::filter(gene_id %in% c("ENSG00000137474.23", "ENSG00000006611.17", "ENSG00000107736.22", 
                                                                                              "ENSG00000150275.20", "ENSG00000182040.9", "ENSG00000136425.14", 
                                                                                              "ENSG00000042781.14", "ENSG00000164199.18", "ENSG00000095397.16", 
                                                                                              "ENSG00000163646.12", "ENSG00000141337.13"))  

#save the file as backup (the original read assignment file is too large).
USHs1_4_sorted_gene_assignments <- read_xlsx("dataset 1_CCS3/USHs1_4_sorted_gene_assignments.xlsx") 
write.xlsx(USHs1_4_sorted_gene_assignments, "USHs1_4_sorted_gene_assignments.xlsx")


#rename ENSG codes to gene name (number by USH subtype)

USHs1_4_sorted_gene_assignments$gene_id[USHs1_4_sorted_gene_assignments$gene_id=="ENSG00000137474.23"]<-"01_MYO7A"
USHs1_4_sorted_gene_assignments$gene_id[USHs1_4_sorted_gene_assignments$gene_id=="ENSG00000006611.17"]<-"02_USH1C"
USHs1_4_sorted_gene_assignments$gene_id[USHs1_4_sorted_gene_assignments$gene_id=="ENSG00000107736.22"]<-"03_CDH23"
USHs1_4_sorted_gene_assignments$gene_id[USHs1_4_sorted_gene_assignments$gene_id=="ENSG00000150275.20"]<-"04_PCDH15"
USHs1_4_sorted_gene_assignments$gene_id[USHs1_4_sorted_gene_assignments$gene_id=="ENSG00000182040.9"]<-"05_SANS"
USHs1_4_sorted_gene_assignments$gene_id[USHs1_4_sorted_gene_assignments$gene_id=="ENSG00000136425.14"]<-"06_CIB2"
USHs1_4_sorted_gene_assignments$gene_id[USHs1_4_sorted_gene_assignments$gene_id=="ENSG00000042781.14"]<-"07_USH2A"
USHs1_4_sorted_gene_assignments$gene_id[USHs1_4_sorted_gene_assignments$gene_id=="ENSG00000164199.18"]<-"08_ADGRV1"
USHs1_4_sorted_gene_assignments$gene_id[USHs1_4_sorted_gene_assignments$gene_id=="ENSG00000095397.16"]<-"09_WHRN"
USHs1_4_sorted_gene_assignments$gene_id[USHs1_4_sorted_gene_assignments$gene_id=="ENSG00000163646.12"]<-"10_CLRN1"
USHs1_4_sorted_gene_assignments$gene_id[USHs1_4_sorted_gene_assignments$gene_id=="ENSG00000141337.13"]<-"11_ARSG"


#splice colum with to filter for structural categories and runs
s1_4.df.split2 <- cSplit(USHs1_4_sorted_gene_assignments, "structural_cat", "=", direction = "wide")
s1_4.df.split3 <- cSplit(s1_4.df.split2, "run", "/", direction = "wide")


#rename structural categories to short name (for graph)
s1_4.df.split3$structural_cat_2[s1_4.df.split3$structural_cat_2=="full_splice_match"]<-"05_FSM"
s1_4.df.split3$structural_cat_2[s1_4.df.split3$structural_cat_2=="incomplete_splice_match"]<-"04_ISM"
s1_4.df.split3$structural_cat_2[s1_4.df.split3$structural_cat_2=="novel_in_catalog"]<-"03_NIC"
s1_4.df.split3$structural_cat_2[s1_4.df.split3$structural_cat_2=="novel_not_in_catalog"]<-"02_NNIC"
s1_4.df.split3$structural_cat_2[s1_4.df.split3$structural_cat_2=="genic"]<-"01_genic"

#filter for the different runs and create individual dfs containing structural category counts per gene
sample1 <- s1_4.df.split3 %>% dplyr::filter(run_1 %in% c("m64167e_210819_012015"))
sample1_s <- sample1 %>%
  group_by(gene_id, structural_cat_2) %>% 
  summarise(n = n())

sample2 <- s1_4.df.split3 %>% dplyr::filter(run_1 %in% c("m64167e_210902_121011"))
sample2_s <- sample2 %>%
  group_by(gene_id, structural_cat_2) %>% 
  summarise(n = n())

sample3 <- s1_4.df.split3 %>% dplyr::filter(run_1 %in% c("m64167e_210903_121024"))
sample3_s <- sample3 %>%
  group_by(gene_id, structural_cat_2) %>% 
  summarise(n = n())

sample4 <- s1_4.df.split3 %>% dplyr::filter(run_1 %in% c("m64037e_211216_160915"))
sample4_s <- sample4 %>%
  group_by(gene_id, structural_cat_2) %>% 
  summarise(n = n())



######### make a plot for the structural cathegory findings  ########

#set color scheme

group.colors.struc <- c("05_FSM" = "#581845", 
                        "04_ISM" = "#900C3F", 
                        "03_NIC" ="#C70039", 
                        "02_NNIC" = "#FF5733",
                        "01_genic" = "#dadada")


#plot events per run and per gene _ sample 1 (m64167e_210819_012015)                          
df_1.2_struct_cat <- sample1 %>% 
  group_by(run_1, gene_id, structural_cat_2)%>% 
  summarise(n = n())

df_1.2_struct_cat %>% 
  ggplot(aes(x = n, y = gene_id, fill = structural_cat_2))+
  geom_col(aes(fill = structural_cat_2))+ #barplot with category as fill color
  facet_wrap(gene_id ~ ., ncol = 1, #multiple plots as horizontal facets
             scales = "free")+ #free scale
  theme(strip.text.x = element_blank()) + #remove facet titles
  scale_fill_manual(values = group.colors.struc)

ggsave("structural_category_counts_sample1.pdf", width = 24, height = 24, units = "cm")


#plot events per run and per gene _ sample 2 (m64167e_210902_121011)
df_2.2_struct_cat <- sample2 %>% 
  group_by(run_2, gene_id, structural_cat_2)%>% 
  summarise(n = n())

df_2.2_struct_cat %>% 
  ggplot(aes(x = n, y = gene_id, fill = structural_cat_2))+
  geom_col(aes(fill = structural_cat_2))+ #barplot with category as fill color
  facet_wrap(gene_id ~ ., ncol = 1, #multiple plots as horizontal facets
             scales = "free")+ #free scale
  theme(strip.text.x = element_blank()) + #remove facet titles
  scale_fill_manual(values = group.colors.struc)

ggsave("structural_category_counts_sample2.pdf", width = 24, height = 24, units = "cm")


#plot events per run and per gene _ sample 3 (m64167e_210903_121024)
df_3.2_struct_cat <- sample3 %>% 
  group_by(run_3, gene_id, structural_cat_2)%>% 
  summarise(n = n())

df_3.2_struct_cat %>% 
  ggplot(aes(x = n, y = gene_id, fill = structural_cat_2))+
  geom_col(aes(fill = structural_cat_2))+ #barplot with category as fill color
  facet_wrap(gene_id ~ ., ncol = 1, #multiple plots as horizontal facets
             scales = "free")+ #free scale
  theme(strip.text.x = element_blank()) + #remove facet titles
  scale_fill_manual(values = group.colors.struc)

ggsave("structural_category_counts_sample3.pdf", width = 24, height = 24, units = "cm")


#plot events per run and per gene _ sample 4 (m64037e_211216_160915)
df_4.2_struct_cat <- sample4 %>% 
  group_by(run_3, gene_id, structural_cat_2)%>% 
  summarise(n = n())

df_4.2_struct_cat %>% 
  ggplot(aes(x = n, y = gene_id, fill = structural_cat_2))+
  geom_col(aes(fill = structural_cat_2))+ #barplot with category as fill color
  facet_wrap(gene_id ~ ., ncol = 1, #multiple plots as horizontal facets
             scales = "free")+ #free scale
  theme(strip.text.x = element_blank()) + #remove facet titles
  scale_fill_manual(values = group.colors.struc)

ggsave("structural_category_counts_sample4.pdf", width = 24, height = 24, units = "cm")

#final annotations in figures were done manually in Adobe Illustrator

