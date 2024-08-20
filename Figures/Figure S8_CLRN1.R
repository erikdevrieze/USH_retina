
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
gene_of_interest <- "ENSG00000163646.12"       #  CLRN1  ENSG00000163646.12   

refseq_gencode <- gencode %>% 
  dplyr::filter(
    !is.na(gene_name), 
    gene_name == gene_of_interest
  ) 

CLRN1_CCS3 <- s1_3_isoQ %>% 
  dplyr::filter(
    !is.na(gene_name), 
    gene_name == gene_of_interest
  ) 




# change df_columns for figure

#Select MANE select and MANE clinical as reference for annotation. add origin column & rename ENST (proir to making GIO_combi, identical ENST in multiple datasets causes an error)  
refseq_gencode_clean <- refseq_gencode %>% dplyr::filter(transcript_id %in% c("ENST00000327047.6"))
refseq_gencode_clean$origin <- "gencode"
refseq_gencode_clean$transcript_id[refseq_gencode_clean$transcript_id=="ENST00000327047.6"]<-"01_Gencode_ENST00000327047.6"


#Select published isoforms (Joensuu et al 2001). Isoforms selected from ENS transcripts that encode the proposed protein isoforms. 
#add origin column & rename ENST (proir to making GIO_combi, identical ENST in multiple datasets causes an error)  
#ENST00000485607.1 is annotated to have 4 exons (the 4th is part of the 3' UTR), but Joensuu only annoated the first 2 exons (which is altered manually in the final figure to correspond to Joensuu et al)
CLRN1_joensuu <- refseq_gencode %>% dplyr::filter(transcript_id %in% c("ENST00000485607.1", "ENST00000295911.6"))
CLRN1_joensuu$origin <- "Adato"
CLRN1_joensuu$transcript_id[CLRN1_joensuu$transcript_id=="ENST00000485607.1"]<-"02_Joensuu_canonical_ENST00000485607.1"
CLRN1_joensuu$transcript_id[CLRN1_joensuu$transcript_id=="ENST00000295911.6"]<-"03_Joensuu_isoB_ENST00000295911.6"


#Select published isoforms (Adato et al 2002). Isoforms selected from ENS transcripts that encode the proposed protein isoforms. 
#add origin column & rename ENST (proir to making GIO_combi, identical ENST in multiple datasets causes an error)  
#ENST00000485607.1 is annotated to have 4 exon (the 4th is part of the 3' UTR), but Adato only annoated the first 3 exons (which is altered manually in the final figure to correspond to Adato et al)
CLRN1_adato <- refseq_gencode %>% dplyr::filter(transcript_id %in% c("ENST00000327047.6", "ENST00000485607.1", "ENST00000295911.6"))
CLRN1_adato$origin <- "Adato"
CLRN1_adato$transcript_id[CLRN1_adato$transcript_id=="ENST00000327047.6"]<-"04_ENST00000327047.6"
CLRN1_adato$transcript_id[CLRN1_adato$transcript_id=="ENST00000295911.6"]<-"05_ENST00000295911.6"
CLRN1_adato$transcript_id[CLRN1_adato$transcript_id=="ENST00000485607.1"]<-"06_ENST00000485607.1"


#IsoQ data s1-3 rename ENST and add origin column
CLRN1_CCS3$origin <- "1-3 CCS3" 
CLRN1_CCS3$transcript_id[CLRN1_CCS3$transcript_id=="ENST00000472224.1"]<-"07_tAtlas_ENST00000472224.1"
CLRN1_CCS3$transcript_id[CLRN1_CCS3$transcript_id=="ENST00000327047.6"]<-"08_tAtlas_ENST00000327047.6"
CLRN1_CCS3$transcript_id[CLRN1_CCS3$transcript_id=="transcript28193.chr3.nic"]<-"09_tAtlas_transcript28193.chr3.nic"





# combine df_columns for figure
CLRN1_fig_full <- bind_rows(refseq_gencode_clean, CLRN1_joensuu, CLRN1_adato, CLRN1_CCS3)  


# extract the required annotation columns 
CLRN1_fig <- CLRN1_fig_full %>% 
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
GOI_exons <- CLRN1_fig %>% dplyr::filter(type == "exon")


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
  ) + scale_fill_manual(values=c("#F8766D", "#00BFC4", "#7CAE00", "#E68613","#E68613", "#999999")) +
  geom_intron(
    data = GOI_rescaled_introns,
    aes(strand = strand), 
    arrow.min.intron.length = 700
  ) + 
  theme_classic(
  ) 

ggsave("Fig Sx_CLRN1.pdf", width = 25, height = 18, units = "cm")
#save on fixed width of 10 with height of 2 per transcript isoform

# additional annotations were made manually

