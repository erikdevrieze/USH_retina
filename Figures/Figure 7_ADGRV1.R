
#Load packages
library(dplyr) 
library(ggplot2) 
library(readxl) 
library(reshape2) 
library(ggtranscript)
library(magrittr)
library(rtracklayer)
library(xlsx) 
library(tidyverse)
library(egg)


setwd("/Users/erikdevrieze/Library/CloudStorage/OneDrive-Radboudumc/z918116/UMCN/Manuscripts/2024 - USH isoseq manuscript/Analysis paper 2024")


#####  Visualize from data GTF  #################

# Import Gencode reference file into a temporary directory
gencode_path <- file("gencode.v45.annotation.gtf")
gencode <- rtracklayer::import(gencode_path)
gencode <- gencode %>% dplyr::as_tibble()

# Import s1-3 CCS3 isoQuant (Riepe et al) file into a temporary directory
s1_3_isoQ_path <- file.path("dataset Atlas/00_aln.sorted.transcript_models.gtf")
s1_3_isoQ <- rtracklayer::import(s1_3_isoQ_path)
s1_3_isoQ <- s1_3_isoQ %>% dplyr::as_tibble()


#Renome Gene_ID column to gene_name
colnames(s1_3_isoQ)[10] ="gene_name"


# filter GTF files for the gene of interest 
gene_of_interest <- "ADGRV1"       #  "ADGRV1" in gencode GTF &   "ENSG00000164199.18" in isoquant GTF

refseq_gencode <- gencode %>% 
  dplyr::filter(
    !is.na(gene_name), 
    gene_name == gene_of_interest
  ) 

ADGRV1_CCS3 <- s1_3_isoQ %>% 
  dplyr::filter(
    !is.na(gene_name), 
    gene_name == gene_of_interest
  ) 


# change df_columns for figure

#Select MANE select refseq, add origin column & rename ENST (proir to making GIO_combi, identical ENST in multiple datasets causes an error) 
refseq_gencode_clean <- refseq_gencode   %>% dplyr::filter(transcript_id %in% c("ENST00000405460.9"))
refseq_gencode_clean$origin <- "gencode"
refseq_gencode_clean$transcript_id[refseq_gencode_clean$transcript_id=="ENST00000405460.9"]<-"01_Gencode_ENST00000405460.9"

#produce proposed Samplix FL transcript
ADGRV1_samplix <- refseq_gencode   %>% dplyr::filter(transcript_id %in% c("ENST00000405460.9"))
ADGRV1_samplix$origin <- "samplix"
ADGRV1_samplix$transcript_id[ADGRV1_samplix$transcript_id=="ENST00000405460.9"]<-"Xsamplix_ENST00000405460.9"

#produce proposed Samplix 1e transcript. As for figure 6A, we used to manually produced xls file containg the information on the Vlgr1e transcript isoform. 
#input file is uploaded to the data folder
ADGRV1_samplix_e <- read.xlsx("USH Genes/ADGRV1/ADGRV1_Samplix_Vlgr1e.xlsx", "ADGRV1_Samplix_Vlgr1e")
ADGRV1_samplix_e$origin <- "Xx_ADGRV1_samplix_e"
ADGRV1_samplix_e$transcript_id[ADGRV1_samplix_e$transcript_id=="X_ADGRV1_Samplix_Vlgr1e"]<-"XX_ADGRV1_Samplix_Vlgr1e"

#produce proposed Pacbio manual curation isoforms that are not refseq (read from xls). Similar to the Vlgr1e trnascript, we have taken the FL ADGRV1 isoform and manually annotated exon 39A in an xls file
#input file is uploaded to the data folder
Pacbio_manual1 <- read.xlsx("USH Genes/ADGRV1/ADGRV1_Full length incl 39A.xlsx", "ADGRV1_Exon39A")
Pacbio_manual1$origin <- "Pacbio_manual1"

#produce proposed Pacbio manual curation isoforms that are not refseq (read from xls) Similar to the Vlgr1e trnascript, we have taken the FL ADGRV1 isoform and manually annotated the U2_VLGR1-c found in the BAM file. 
#input file is uploaded to the data folder
Pacbio_manual2 <- read.xlsx("USH Genes/ADGRV1/ADGRV1_Pacbio_isoC.xlsx", "U2_VLGR1-c")
Pacbio_manual2$origin <- "Pacbio_manual2"


# combine df_columns for figure
ADGRV1_fig_full <- bind_rows(refseq_gencode_clean, Pacbio_manual1, Pacbio_manual2, ADGRV1_samplix, ADGRV1_samplix_e)  


# extract the required annotation columns 
ADGRV1_fig <- ADGRV1_fig_full %>% 
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
GOI_exons <- ADGRV1_fig %>% dplyr::filter(type == "exon")


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

  
  geom_range(
    aes(fill = origin)
  ) + scale_fill_manual(values=c("#7CAE00", "#CC0033", "#CC0033", "#E68613","#E68613")) + 
  
  geom_intron(
    data = GOI_rescaled_introns,
    aes(strand = strand), 
    arrow.min.intron.length = 700
  ) + 
  theme_classic(
  ) 



# additional annotations were made manually
