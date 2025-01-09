
#Load packages
library(dplyr)
library(ggplot2) 
library(readxl) 
library(reshape2) 
library(ggtranscript)
library(magrittr)
library(rtracklayer)
library(xlsx) 



#####  Visualize from data GTF  #################

# Import Gencode reference file into a temporary directory
gencode_path <- file("gencode.v45.annotation.gtf")
gencode <- rtracklayer::import(gencode_path)
gencode <- gencode %>% dplyr::as_tibble()

# Import s1-3 CCS3 isoQuant (Riepe et al) file into a temporary directory
s1_3_isoQ_path <- file.path("dataset Atlas/00_aln.sorted.transcript_models.gtf")
s1_3_isoQ <- rtracklayer::import(s1_3_isoQ_path)
s1_3_isoQ <- s1_3_isoQ %>% dplyr::as_tibble()


# Import s1-4 CCS3 isoQuant (Riepe et al) file into a temporary directory
s1_4_isoQ_path <- file.path("dataset 1_CCS3/Isoquant/00_aln.sorted.transcript_models.gtf")
s1_4_isoQ <- rtracklayer::import(s1_4_isoQ_path)
s1_4_isoQ <- s1_4_isoQ %>% dplyr::as_tibble()

# Import s4 CCS0 isoQuant (Riepe et al) file into a temporary directory
s4_isoQ_path <- file.path("dataset 2_CCS0/isoquant/00_aln.sorted.transcript_models.gtf")
s4_isoQ <- rtracklayer::import(s4_isoQ_path)
s4_isoQ <- s4_isoQ %>% dplyr::as_tibble()

#Renome Gene_ID column to gene_name
colnames(s1_3_isoQ)[10] ="gene_name"
colnames(s1_4_isoQ)[10] ="gene_name"
colnames(s4_isoQ)[10] ="gene_name"



# filter GTF files for the gene of interest 
gene_of_interest <- "CDH23"       #  CDH23  ENSG00000107736.22   

refseq_gencode <- gencode %>% 
  dplyr::filter(
    !is.na(gene_name), 
    gene_name == gene_of_interest
  ) 

CDH23_CCS3 <- s1_3_isoQ %>% 
  dplyr::filter(
    !is.na(gene_name), 
    gene_name == gene_of_interest
  ) 


# change df_columns for figure

#Select MANE select refseqs, add origin column & rename ENST (proir to making GIO_combi, identical ENST in multiple datasets causes an error)  
refseq_gencode_clean <- refseq_gencode %>% dplyr::filter(transcript_id %in% c("ENST00000224721.12"))
refseq_gencode_clean$origin <- "gencode"
refseq_gencode_clean$transcript_id[refseq_gencode_clean$transcript_id=="ENST00000224721.12"]<-"01_Gencode_ENST00000224721.12"


#Select published isoforms (Lagziel et al 2005). Isoforms selected from ENS transcripts that encode the proposed protein isoforms. 
#add origin column & rename ENST (proir to making GIO_combi, identical ENST in multiple datasets causes an error)  
CDH23_A1 <- read.xlsx("USH Genes/CDH23/CDH23_v1a.xlsx", "1")
CDH23_A1$origin <- "Lagziel"

#Select published isoforms (Lagziel et al 2005). Isoforms selected from ENS transcripts that encode the proposed protein isoforms. 
#add origin column & rename ENST (proir to making GIO_combi, identical ENST in multiple datasets causes an error)  
CDH23_A1h <- read.xlsx("USH Genes/CDH23/CDH23_v1a.xlsx", "1")
CDH23_A1h$origin <- "Lagziel_human"
CDH23_A2h <- read.xlsx("USH Genes/CDH23/CDH23_v1b.xlsx", "1")
CDH23_A2h$origin <- "Lagziel_human"

CDH23_C1h <- read.xlsx("USH Genes/CDH23/CDH23_v3a.xlsx", "1")
CDH23_C1h$origin <- "Lagziel_human"
CDH23_C2h <- read.xlsx("USH Genes/CDH23/CDH23_v3b.xlsx", "1")
CDH23_C2h$origin <- "Lagziel_human"

#mouse
CDH23_A1m <- read.xlsx("USH Genes/CDH23/CDH23_v1a.xlsx", "1")
CDH23_A1m$origin <- "Lagziel_mouse"
CDH23_A1m$transcript_id[CDH23_A1m$transcript_id=="Cdh23_v1a"]<-"D_Cdh23_v1a"
CDH23_A2m <- read.xlsx("USH Genes/CDH23/CDH23_v1b.xlsx", "1")
CDH23_A2m$origin <- "Lagziel_mouse"
CDH23_A2m$transcript_id[CDH23_A2m$transcript_id=="Cdh23_v1b"]<-"D_Cdh23_v1b"

CDH23_B1m <- read.xlsx("USH Genes/CDH23/CDH23_v2a.xlsx", "1")
CDH23_B1m$origin <- "Lagziel_mouse"
CDH23_B1m$transcript_id[CDH23_B1m$transcript_id=="Cdh23_v2a"]<-"D_Cdh23_v2a"
CDH23_B2m <- read.xlsx("USH Genes/CDH23/CDH23_v2b.xlsx", "1")
CDH23_B2m$origin <- "Lagziel_mouse"
CDH23_B2m$transcript_id[CDH23_B2m$transcript_id=="Cdh23_v2b"]<-"D_Cdh23_v2b"

CDH23_C1m <- read.xlsx("USH Genes/CDH23/CDH23_v3a.xlsx", "1")
CDH23_C1m$origin <- "Lagziel_mouse"
CDH23_C1m$transcript_id[CDH23_C1m$transcript_id=="Cdh23_v3a"]<-"D_Cdh23_v3a"
CDH23_C2m <- read.xlsx("USH Genes/CDH23/CDH23_v3b.xlsx", "1")
CDH23_C2m$origin <- "Lagziel_mouse"
CDH23_C2m$transcript_id[CDH23_C2m$transcript_id=="Cdh23_v3b"]<-"D_Cdh23_v3b"

#IsoQ data s1-3 rename ENST and add origin column
CDH23_CCS3$origin <- "1-3 CCS3" 
CDH23_CCS3$transcript_id[CDH23_CCS3$transcript_id=="ENST00000461841.7"]<-"tAtlas_ENST00000461841.7"
CDH23_CCS3$transcript_id[CDH23_CCS3$transcript_id=="transcript11235.chr10.nnic"]<-"tAtlas_transcript11235.chr10.nnic"
CDH23_CCS3$transcript_id[CDH23_CCS3$transcript_id=="ENST00000461841.7"]<-"tAtlas_ENST00000461841.7"




# combine df_columns for figure
CDH23_fig_full <- bind_rows(refseq_gencode_clean, CDH23_CCS3, CDH23_A1h, CDH23_A2h, CDH23_C1h, CDH23_C2h, CDH23_A1m, CDH23_A2m, CDH23_B1m, CDH23_B2m, CDH23_C1m, CDH23_C2m,)  


# extract the required annotation columns 
CDH23_fig <- CDH23_fig_full %>% 
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
GOI_exons <- CDH23_fig %>% dplyr::filter(type == "exon")


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
  ) + scale_fill_manual(values=c("#F8766D", "#7CAE00", "#00BFC4", "#999999","#E68613", "#999999")) +
  geom_intron(
    data = GOI_rescaled_introns,
    aes(strand = strand), 
    arrow.min.intron.length = 700
  ) + 
  theme_classic(
  ) 

ggsave("Figure_Sx_CDH23.pdf", width = 40, height = 30, units = "cm")

# additional annotations were made manually

