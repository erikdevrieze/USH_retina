
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

#Renome Gene_ID column to gene_name
colnames(s1_3_isoQ)[10] ="gene_name"


# filter GTF files for the gene of interest 
gene_of_interest <- "USH1C"       #  USH1C  ENSG00000006611.17   

refseq_gencode <- gencode %>% 
  dplyr::filter(
    !is.na(gene_name), 
    gene_name == gene_of_interest
  ) 

USH1C_CCS3 <- s1_3_isoQ %>% 
  dplyr::filter(
    !is.na(gene_name), 
    gene_name == gene_of_interest
  ) 


# change df_columns for figure

#Select MANE select/clinical refseqs, add origin column & rename ENST (proir to making GIO_combi, identical ENST in multiple datasets causes an error)  
refseq_gencode_clean <- refseq_gencode %>% dplyr::filter(transcript_id %in% c("ENST00000005226.12", "ENST00000318024.9"))
refseq_gencode_clean$origin <- "gencode"
refseq_gencode_clean$transcript_id[refseq_gencode_clean$transcript_id=="ENST00000005226.12"]<-"01_Gencode_ENST00000005226.12"
refseq_gencode_clean$transcript_id[refseq_gencode_clean$transcript_id=="ENST00000318024.9"]<-"02_Gencode_ENST00000318024.9"


#Select published isoforms (nagel-Wolfrum 2023 and refs therein). Isoforms selected from ENS transcripts that encode the proposed protein isoforms. 
#Additional alternaitve splicing events are shown in literature, but appear to be "splicing in progress". Determining events are inclusion/exlusion middle exons, and exon 15 skipping -> FS 
#add origin column & rename ENST (proir to making GIO_combi, identical ENST in multiple datasets causes an error)  
refseq_gencode_clean <- refseq_gencode %>% dplyr::filter(transcript_id %in% c("ENST00000527720.5", "ENST00000527020.5", "ENST00000005226.12", "ENST00000318024.9", "ENST00000526313.5"))
refseq_gencode_clean$origin <- "published"
refseq_gencode_clean$transcript_id[refseq_gencode_clean$transcript_id=="ENST00000527720.5"]<-"05_Harmonin_A2_ENST00000527720.5"
refseq_gencode_clean$transcript_id[refseq_gencode_clean$transcript_id=="ENST00000527020.5"]<-"04_Harmonin_A1_ENST00000527020.5"
refseq_gencode_clean$transcript_id[refseq_gencode_clean$transcript_id=="ENST00000005226.12"]<-"07_Harmonin_B_ENST00000005226.12"
refseq_gencode_clean$transcript_id[refseq_gencode_clean$transcript_id=="ENST00000318024.9"]<-"06_Harmonin_A3_ENST00000318024.9"
refseq_gencode_clean$transcript_id[refseq_gencode_clean$transcript_id=="ENST00000526313.5"]<-"08_Harmonin_C_ENST00000526313.5"

#IsoQ data s1-3 rename ENST and add origin column
USH1C_CCS3$origin <- "1-3 CCS3" 
USH1C_CCS3$transcript_id[USH1C_CCS3$transcript_id=="ENST00000640109.1"]<-"tAtlas_ENST00000640109.1"
USH1C_CCS3$transcript_id[USH1C_CCS3$transcript_id=="ENST00000640281.1"]<-"tAtlas_ENST00000640281.1"
USH1C_CCS3$transcript_id[USH1C_CCS3$transcript_id=="ENST00000638316.1"]<-"tAtlas_ENST00000638316.1"
USH1C_CCS3$transcript_id[USH1C_CCS3$transcript_id=="ENST00000639884.1"]<-"tAtlas_ENST00000639884.1"


# combine df_columns for figure
USH1C_fig_full <- bind_rows(refseq_gencode_clean, USH1C_CCS3)  


# extract the required annotation columns 
USH1C_fig <- USH1C_fig_full %>% 
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
GOI_exons <- USH1C_fig %>% dplyr::filter(type == "exon")


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
  ) + scale_fill_manual(values=c("#F8766D", "#7CAE00", "#00BFC4", "#00BFC4","#E68613", "#999999")) +
  geom_intron(
    data = GOI_rescaled_introns,
    aes(strand = strand), 
    arrow.min.intron.length = 700
  ) + 
  theme_classic(
  ) 



# additional annotations were made manually

