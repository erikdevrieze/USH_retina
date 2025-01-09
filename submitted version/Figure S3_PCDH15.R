
#Load packages
library(dplyr) 
library(ggplot2) 
library(readxl) 
library(reshape2) 
library(ggtranscript)
library(magrittr)
library(rtracklayer)
library(xlsx) 


setwd("/Users/erikdevrieze/Library/CloudStorage/OneDrive-Radboudumc/z918116/UMCN/Manuscripts/2024 - USH isoseq manuscript/Analysis paper 2024")


#####  Visualize from data GTF  #################

# Import Gencode reference file into a temporary directory
gencode_path <- file("gencode.v45.annotation.gtf")
gencode <- rtracklayer::import(gencode_path)
gencode <- gencode %>% dplyr::as_tibble()

# Import s1-3 CCS3 isoQuant file into a temporary directory
s1_3_isoQ_path <- file.path("dataset Atlas/00_aln.sorted.transcript_models.gtf")
s1_3_isoQ <- rtracklayer::import(s1_3_isoQ_path)
s1_3_isoQ <- s1_3_isoQ %>% dplyr::as_tibble()

# Import s1-4 CCS3 isoQuant file into a temporary directory (because one unique transcript is found here)
s1_4_isoQ_path <- file.path("dataset 1_CCS3/Isoquant/00_aln.sorted.transcript_models.gtf")
s1_4_isoQ <- rtracklayer::import(s1_4_isoQ_path)
s1_4_isoQ <- s1_4_isoQ %>% dplyr::as_tibble()


#Renome Gene_ID column to gene_name
colnames(s1_3_isoQ)[10] ="gene_name"
colnames(s1_4_isoQ)[10] ="gene_name"



# filter GTF files for the gene of interest 
gene_of_interest <- "ENSG00000150275.20"       #  PCDH15  ENSG00000150275.20   

refseq_gencode <- gencode %>% 
  dplyr::filter(
    !is.na(gene_name), 
    gene_name == gene_of_interest
  ) 

PCDH15_CCS3 <- s1_3_isoQ %>% 
  dplyr::filter(
    !is.na(gene_name), 
    gene_name == gene_of_interest
  ) 

PCDH15_long <- s1_4_isoQ %>% 
  dplyr::filter(
    !is.na(gene_name), 
    gene_name == gene_of_interest
  ) 



# change df_columns for figure

#Select MANE select/clinical refseqs, add origin column & rename ENST (proir to making GIO_combi, identical ENST in multiple datasets causes an error)  
refseq_gencode_clean <- refseq_gencode %>% dplyr::filter(transcript_id %in% c("ENST00000320301.11", "ENST00000644397.2"))
refseq_gencode_clean$origin <- "gencode"
refseq_gencode_clean$transcript_id[refseq_gencode_clean$transcript_id=="ENST00000644397.2"]<-"01_Gencode_ENST00000644397.2"   # MANE select
refseq_gencode_clean$transcript_id[refseq_gencode_clean$transcript_id=="ENST00000320301.11"]<-"02_Gencode_ENST00000320301.11"  # MANE clinical


#Select published isoforms (Ahmed et al 2008). Isoforms selected from ENS transcripts that contain the exons shows in the paper.  
#add origin column & rename ENST (proir to making GIO_combi, identical ENST in multiple datasets causes an error)  
PCDH15_pub <- refseq_gencode %>% dplyr::filter(transcript_id %in% c("ENST00000320301.11", "ENST00000395445.6", "ENST00000616114.4", "ENST00000463095.2", "ENST00000373955.5"))
PCDH15_pub$origin <- "Ahmed_2008"
PCDH15_pub$transcript_id[PCDH15_pub$transcript_id=="ENST00000320301.11"]<-"03_PCDH15_CD1_ENST00000320301.11"
PCDH15_pub$transcript_id[PCDH15_pub$transcript_id=="ENST00000395445.6"]<-"04_PCDH15_CD2 ENST00000395445.6"
PCDH15_pub$transcript_id[PCDH15_pub$transcript_id=="ENST00000616114.4"]<-"05_PCDH15_CD3_ENST00000616114.4"
PCDH15_pub$transcript_id[PCDH15_pub$transcript_id=="ENST00000463095.2"]<-"06_PCDH15_Cterm_ENST00000463095.2"
PCDH15_pub$transcript_id[PCDH15_pub$transcript_id=="ENST00000373955.5"]<-"07_PCDH15_SI_ENST00000373955.5"

#IsoQ data s1-3 rename ENST and add origin column
PCDH15_CCS3$origin <- "1-3 CCS3" 
PCDH15_CCS3$transcript_id[PCDH15_CCS3$transcript_id=="ENST00000621708.4"]<-"08_tAtlas_ENST00000621708.4"
PCDH15_CCS3$transcript_id[PCDH15_CCS3$transcript_id=="ENST00000373957.7"]<-"08_tAtlas_ENST00000373957.7"
PCDH15_CCS3$transcript_id[PCDH15_CCS3$transcript_id=="ENST00000638316.1"]<-"08_tAtlas_ENST00000638316.1"
PCDH15_CCS3$transcript_id[PCDH15_CCS3$transcript_id=="ENST00000639884.1"]<-"08_tAtlas_ENST00000639884.1"

#IsoQ data s1-4 rename ENST and add origin column. Filter for the one transcript isoform that is unique to this dataset but not in the isoquant output produced by from using samples 1-3 as input
PCDH15_long <- PCDH15_long %>% dplyr::filter(transcript_id %in% c("ENST00000373955.5"))
PCDH15_long$origin <- "1-4 CCS3"
PCDH15_long$transcript_id[PCDH15_long$transcript_id=="ENST00000373955.5"]<-"09_long_ENST00000373955.5"


# combine df_columns for figure
PCDH15_fig_full <- bind_rows(refseq_gencode_clean, PCDH15_pub, PCDH15_CCS3, PCDH15_long)


# extract the required annotation columns 
PCDH15_fig <- PCDH15_fig_full %>% 
  dplyr::select(
    seqnames,
    start,
    end,
    strand,
    type,
    gene_name,
    transcript_type,
    transcript_id,
    tag, 
    origin
  )



#####  Visualize transcript type  #################


# obtain exons
GOI_exons <- PCDH15_fig %>% dplyr::filter(type == "exon")


#####  Reduce intron length  #################

GOI_rescaled <- shorten_gaps(
  exons = GOI_exons, 
  introns = to_intron(GOI_exons, "transcript_id"), 
  group_var = "transcript_id", target_gap_width =50L
)


# shorten_gaps() returns exons and introns all in one data.frame()
# let's split these for plotting 
GOI_rescaled_exons <- GOI_rescaled %>% dplyr::filter(type == "exon") 
GOI_rescaled_introns <- GOI_rescaled %>% dplyr::filter(type == "intron") 

GOI_rescaled_exons %>% 
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_id,
  ))  +
  scale_y_discrete(limits = rev)  +                              #reverts y axis
  coord_cartesian(xlim = c(max(GOI_rescaled_exons$end), (min(GOI_rescaled_exons$start))))+       #invert orientation  
  
  
  geom_range(
    aes(fill = origin)
  ) + scale_fill_manual(values=c("#F8766D", "#E68613", "#00BFC4", "#7CAE00","#CC66CC", "#999999")) +
  geom_intron(
    data = GOI_rescaled_introns,
    aes(strand = strand), 
    arrow.min.intron.length = 700
  ) + 
  theme_classic(
  ) +
  
  ggsave("PCDH15_figure v1.pdf", width = 40, height = 42, units = "cm")
 #save on fixed width of 40 with height of 2 per transcript isoform



# additional annotations were made manually

