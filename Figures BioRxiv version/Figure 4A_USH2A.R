
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

# Import isoQuant track obtained with sample 1-3 CCS3 as input into a temporary directory
s1_3_isoQ_path <- file.path("dataset Atlas/00_aln.sorted.transcript_models.gtf")
s1_3_isoQ <- rtracklayer::import(s1_3_isoQ_path)
s1_3_isoQ <- s1_3_isoQ %>% dplyr::as_tibble()

#Renome Gene_ID column to gene_name
colnames(s1_3_isoQ)[10] ="gene_name"


# filter GTF files for the gene of interest 
gene_of_interest <- "USH2A"       #  "USH2A" in Gencode GTF & "ENSG00000042781.14" in Isoquant GTF

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


#IsoQ data s1-3 rename ENST and add origin column
USH2A_CCS3$origin <- "1-3 CCS3"
USH2A_CCS3$transcript_id[USH2A_CCS3$transcript_id=="ENST00000307340.8"]<-"tAtlas_ENST00000307340.8"


# Human_retina isoform based on van Wijk et al 2003
# to produce these transcript in the figure, they are taken from the Gencode GTF and add origin column & rename ENST (proir to making GIO_combi, identical ENST in multiple datasets causes an error) 
human_retina <- human_retina %>% dplyr::filter(transcript_id %in% c("ENST00000307340.8","ENST00000366942.3"))
human_retina$origin <- "retina"
human_retina$transcript_id[human_retina$transcript_id=="ENST00000366942.3"]<-"2_USH2A_isoA"
human_retina$transcript_id[human_retina$transcript_id=="ENST00000307340.8"]<-"3_USH2A_isoB"

# combine df_columns for figure
USH2A_fig_full <- bind_rows(refseq_gencode_clean, USH2A_CCS3, human_retina)  


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
  group_var = "transcript_id", target_gap_width =100L
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
  ) +scale_fill_manual(values=c("#F8766D", "#7CAE00", "#00BFC4")) + 
  
  geom_intron(
    data = GOI_rescaled_introns,
    aes(strand = strand), 
    arrow.min.intron.length = 700
  ) + 
  theme_classic(
  ) 

ggsave("USH2A_figure part1.pdf", width = 40, height = 28, units = "cm")
#save on fixed width of 40 with height of 2 per transcript isoform


# additional annotations were made manually

