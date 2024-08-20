
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

# Import s1-3 CCS3 isoQuant (Riepe et al) file into a temporary directory
s1_3_isoQ_path <- file.path("dataset Atlas/00_aln.sorted.transcript_models.gtf")
s1_3_isoQ <- rtracklayer::import(s1_3_isoQ_path)
s1_3_isoQ <- s1_3_isoQ %>% dplyr::as_tibble()

#Renome Gene_ID column to gene_name
colnames(s1_3_isoQ)[10] ="gene_name"


# filter GTF files for the gene of interest 
gene_of_interest <- "USH1G"       #  USH1G  ENSG00000182040.9   

refseq_gencode <- gencode %>% 
  dplyr::filter(
    !is.na(gene_name), 
    gene_name == gene_of_interest
  ) 

SANS_CCS3 <- s1_3_isoQ %>% 
  dplyr::filter(
    !is.na(gene_name), 
    gene_name == gene_of_interest
  ) 


# change df_columns for figure

#Select MANE select refseqs, add origin column & rename ENST (proir to making GIO_combi, identical ENST in multiple datasets causes an error)  
refseq_gencode_clean <- refseq_gencode #%>% dplyr::filter(transcript_id %in% c("ENST00000614341.5"))
refseq_gencode_clean$origin <- "gencode"
refseq_gencode_clean$transcript_id[refseq_gencode_clean$transcript_id=="ENST00000614341.5"]<-"01_Gencode_ENST00000614341.5"

#IsoQ data s1-3 rename ENST and add origin column
SANS_CCS3$origin <- "1-3 CCS3" 
SANS_CCS3$transcript_id[SANS_CCS3$transcript_id=="ENST00000614341.5"]<-"tAtlas_ENST00000614341.5"
SANS_CCS3$transcript_id[SANS_CCS3$transcript_id=="ENST00000475158.1"]<-"tAtlas_ENST00000475158.1"  
SANS_CCS3$transcript_id[SANS_CCS3$transcript_id=="transcript11235.chr10.nnic"]<-"tAtlas_transcript11235.chr10.nnic"  


# combine df_columns for figure
SANS_fig_full <- bind_rows(refseq_gencode_clean, SANS_CCS3)  


# extract the required annotation columns 
SANS_fig <- SANS_fig_full %>% 
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
GOI_exons <- SANS_fig %>% dplyr::filter(type == "exon")


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
  ) + scale_fill_manual(values=c("#F8766D", "#7CAE00", "#00BFC4", "#E68613","#E68613", "#999999")) +
  geom_intron(
    data = GOI_rescaled_introns,
    aes(strand = strand), 
    arrow.min.intron.length = 700
  ) + 
  theme_classic(
  ) 



# additional annotations were made manually

