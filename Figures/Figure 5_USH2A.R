
#Load packages
library(dplyr) 
library(ggplot2)
library(readxl)
library(reshape2) 
library(ggtranscript)
library(magrittr)
library(rtracklayer)
library(xlsx) 


setwd("/Users/erikdevrieze/Library/CloudStorage/OneDrive-Radboudumc/z918116/UMCN/Manuscripts/2024 - USH isoseq manuscript/Analysis paper 2024/")


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
gene_of_interest <- "USH2A"       #  USH2A  ENSG00000042781.14   

refseq_gencode <- gencode %>% 
  dplyr::filter(
    !is.na(gene_name), 
    gene_name == gene_of_interest
  ) 

USH2A_CCS3 <- s1_3_isoQ %>% 
  dplyr::filter(
    !is.na(gene_name), 
    gene_name == gene_of_interest
  ) 

human_retina <- gencode %>% 
  dplyr::filter(
    !is.na(gene_name), 
    gene_name == gene_of_interest)


# change df_columns for figure

#Select MANE select, add origin column & rename ENST (proir to making GIO_combi, identical ENST in multiple datasets causes an error) 
refseq_gencode_clean <- refseq_gencode %>% dplyr::filter(transcript_id %in% c("ENST00000307340.8"))
refseq_gencode_clean$origin <- "gencode"
refseq_gencode_clean$transcript_id[refseq_gencode_clean$transcript_id=="ENST00000307340.8"]<-"1_Gencode_ENST00000307340.8"


#Isoform predictions based on manual curation of BAM files (cryptic exon track)
USH2A_BAMp <- read_xlsx("USH genes/USH2A/USH2A_i20 events.xlsx", "USH2A_i20")
USH2A_BAMp$origin <- "BAM"

#Annotated REFSEQs based on BAM. selected transcripts taken from gencode. 
USH2A_annotated <- refseq_gencode %>% dplyr::filter(transcript_id %in% c("ENST00000307340.8","ENST00000366942.3"))
USH2A_annotated$transcript_id[USH2A_annotated$transcript_id=="ENST00000307340.8"]<-"2_ENST00000307340.8"
USH2A_annotated$transcript_id[USH2A_annotated$transcript_id=="ENST00000366942.3"]<-"3_ENST00000366942.3"
USH2A_annotated$origin <- "BAM"

#Isoform predictions based on manual curation of BAM files (5' variation)
USH2A_BAM5 <- read_xlsx("USH genes/USH2A/USH2A_5-variation.xlsx", "Sheet1")
USH2A_BAM5$origin <- "BAM"


#Visualize Samplix FL isoform, taken from gencode.
samplix <- refseq_gencode %>% dplyr::filter(transcript_id %in% c("ENST00000307340.8"))
samplix$origin <- "samplix"
samplix$transcript_id[samplix$transcript_id=="ENST00000307340.8"]<-"ZZ_samplix_ENST00000307340.8"


# combine df_columns for figure
USH2A_fig_full <- bind_rows(refseq_gencode_clean, USH2A_annotated, USH2A_BAM5, USH2A_BAMp, samplix)  


# extract the required annotation columns 
USH2A_fig <- USH2A_fig_full %>% 
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
GOI_cds <- USH2A_fig %>% dplyr::filter(type == "CDS")

# obtain exons
GOI_exons <- USH2A_fig %>% dplyr::filter(type == "exon")


#####  Reduce intron length  #################

GOI_rescaled <- shorten_gaps(
  exons = GOI_exons, 
  introns = to_intron(GOI_exons, "transcript_id"), 
  group_var = "transcript_id", target_gap_width =150L
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
  ) +scale_fill_manual(values=c("#CC0033", "#7CAE00", "#E68613")) + 
  
  geom_intron(
    data = GOI_rescaled_introns,
    aes(strand = strand), 
    arrow.min.intron.length = 700
  ) + 
  theme_classic(
  ) 

ggsave("USH2A_figure part2.pdf", width = 40, height = 12, units = "cm")
#save on fixed width of 40 with height of 2 per transcript isoform


# additional annotations were made manually in Adobe Illustrator.

