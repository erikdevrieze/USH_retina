
#Load packages
library(dplyr) 
library(ggplot2) 
library(readxl) 
library(reshape2) 
library(ggtranscript)
library(magrittr)
library(rtracklayer)
library(xlsx) 


setwd("///)

#####  Visualize from data GTF  #################

# Import Gencode reference file into a temporary directory
gencode_path <- file("gencode.v45.annotation.gtf")
gencode <- rtracklayer::import(gencode_path)
gencode <- gencode %>% dplyr::as_tibble()

# Import s1-3 CCS3 isoQuant GTF file into a temporary directory
s1_3_isoQ_path <- file.path("dataset Atlas/00_aln.sorted.transcript_models.gtf")
s1_3_isoQ <- rtracklayer::import(s1_3_isoQ_path)
s1_3_isoQ <- s1_3_isoQ %>% dplyr::as_tibble()

# Import s1-4 CCS3 isoQuant GTF file into a temporary directory (dataset including the 4th long enriched library)
s1_4_isoQ_path <- file.path("dataset 1_CCS3/Isoquant/00_aln.sorted.transcript_models.gtf")
s1_4_isoQ <- rtracklayer::import(s1_4_isoQ_path)
s1_4_isoQ <- s1_4_isoQ %>% dplyr::as_tibble()

#Renome Gene_ID column to gene_name
colnames(s1_3_isoQ)[10] ="gene_name"
colnames(s1_4_isoQ)[10] ="gene_name"

# filter GTF files for the gene of interest 
gene_of_interest <- "ENSG00000137474.23"       #  MYO7A  ENSG00000137474.23   

refseq_gencode <- gencode %>% 
  dplyr::filter(
    !is.na(gene_name), 
    gene_name == gene_of_interest
  ) 

MYO7A_CCS3 <- s1_3_isoQ %>% 
  dplyr::filter(
    !is.na(gene_name), 
    gene_name == gene_of_interest
  ) 

MYO7A_s14 <- s1_4_isoQ %>% 
  dplyr::filter(
    !is.na(gene_name), 
    gene_name == gene_of_interest
  ) 


# change df_columns for figure

#Select MANE select/clinical as reference isoforms from Gencode dataset 
#rename ENSTXXX transcripts in dataset to prevent redundant names. 
refseq_gencode_clean <- refseq_gencode %>% dplyr::filter(transcript_id %in% c("ENST00000409709.9"))
refseq_gencode_clean$origin <- "gencode"
refseq_gencode_clean$transcript_id[refseq_gencode_clean$transcript_id=="ENST00000409709.9"]<-"1_Gencode_ENST00000409709.9"


#rename ENSTXXX transcripts in CCS3 dataset to prevent redundant names. 
MYO7A_CCS3$origin <- "1-3 CCS3"
MYO7A_CCS3$transcript_id[MYO7A_CCS3$transcript_id=="ENST00000409709.9"]<-"tAtlas_ENST00000409709.9"
MYO7A_CCS3$transcript_id[MYO7A_CCS3$transcript_id=="ENST00000458637.6"]<-"tAtlas_ENST00000458637.6"

#rename ENSTXXX transcripts in sample 1-4 dataset to prevent redundant names. 
MYO7A_s14$origin <- "1-4 CCS3"
MYO7A_s14$transcript_id[MYO7A_s14$transcript_id=="ENST00000409709.9"]<-"tAtlas_ENST00000409709.9"
MYO7A_s14$transcript_id[MYO7A_s14$transcript_id=="ENST00000458637.6"]<-"tAtlas_ENST00000458637.6"


# combine df_columns for figure
MYO7A_fig_full <- bind_rows(refseq_gencode_clean, MYO7A_s14, MYO7A_CCS3 )  

# extract the required annotation columns 
MYO7A_fig <- MYO7A_fig_full %>% 
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


# obtain cds
MYO7A_cds <- MYO7A_fig %>% dplyr::filter(type == "CDS")

# obtain exons
MYO7A_exons <- MYO7A_fig %>% dplyr::filter(type == "exon")


#####  Reduce intron length  #################

MYO7A_rescaled <- shorten_gaps(
  exons = MYO7A_exons, 
  introns = to_intron(MYO7A_exons, "transcript_id"), 
  group_var = "transcript_id", target_gap_width =50L
)


# shorten_gaps() returns exons and introns all in one data.frame()
# let's split these for plotting 
MYO7A_rescaled_exons <- MYO7A_rescaled %>% dplyr::filter(type == "exon") 
MYO7A_rescaled_introns <- MYO7A_rescaled %>% dplyr::filter(type == "intron") 

MYO7A_rescaled_exons %>% 
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_id,
  ))  +
  geom_range(
    aes(fill = origin)
  ) + scale_fill_manual(values=c("#F8766D", "#00BFC4","#7CAE00", "#999999")) +
  geom_intron(
    data = MYO7A_rescaled_introns,
    aes(strand = strand), 
    arrow.min.intron.length = 700
  ) + 
  theme_classic(
  ) + 
  scale_y_discrete(limits = rev)


# further processing of the plot was manually in adobe illustrator
