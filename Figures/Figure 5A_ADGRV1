# Packages needed:
#   dplyr, ggplot2, ggtranscript, rtracklayer, readxl

# Load packages(ev"dplyr"#Load packages(every time)
library(dplyr) #include
library(ggplot2) #include
library(ggtranscript) # include
library(rtracklayer) #include
library(readxl)


setwd("/Users/erikdevrieze/Library/CloudStorage/OneDrive-Radboudumc/z918116/UMCN/Manuscripts/2024 - USH isoseq manuscript/Analysis paper 2024")


#####  Visualize from data GTF  #################

# Import Gencode reference file into a temporary directory
gencode_path <- file("gencode.v45.annotation.gtf")
gencode <- rtracklayer::import(gencode_path)
gencode <- gencode %>% dplyr::as_tibble()

# Import isoquant dataset 1 (samples 1-3, CCS3 as used in Riepe et al) file into a temporary directory
s1_3_isoQ_path <- file.path("dataset Atlas/00_aln.sorted.transcript_models.gtf")
s1_3_isoQ <- rtracklayer::import(s1_3_isoQ_path)
s1_3_isoQ <- s1_3_isoQ %>% dplyr::as_tibble()


#Renome Gene_ID column to gene_name
colnames(s1_3_isoQ)[10] ="gene_name"


# filter GTF files for the gene of interest ("ADGRV1" is used for the gencode dataset, "ENSG00000164199.18" to filter in the isoquant data)
gene_of_interest <- "ADGRV1"     

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



# Select MANE select/clinical as reference isoforms from Gencode dataset 
# rename ENSTXXX transcripts in dataset to prevent redundant names. 
refseq_gencode_clean <- refseq_gencode   %>% dplyr::filter(transcript_id %in% c("ENST00000405460.9"))
refseq_gencode_clean$origin <- "gencode"
refseq_gencode_clean$transcript_id[refseq_gencode_clean$transcript_id=="ENST00000405460.9"]<-"01_Gencode_ENST00000405460.9"

# Select human transcript isoforms published in McMillan et al (2002) - VLGR1a and VLGR1b
# rename ENSTXXX transcripts in dataset to prevent redundant names. 
ADGRV1_McMillan_ab <- refseq_gencode   %>% dplyr::filter(transcript_id %in% c("ENST00000405460.9", "ENST00000638510.1"))
ADGRV1_McMillan_ab$origin <- "McMillan_a"
ADGRV1_McMillan_ab$transcript_id[ADGRV1_McMillan_ab$transcript_id=="ENST00000405460.9"]<-"03_VLGR1-b_ENST00000405460.9"
ADGRV1_McMillan_ab$transcript_id[ADGRV1_McMillan_ab$transcript_id=="ENST00000638510.1"]<-"04_VLGR1-a_ENST00000638510.1"

#Select and create orthologs of mouse transcripts isoform published in McMillan (2002) - Vlgr1a, Vlgr1b, Vlgr1c
Mouse_ab <- refseq_gencode   %>% dplyr::filter(transcript_id %in% c("ENST00000405460.9"))
Mouse_ab$origin <- "Mouse_ab"
Mouse_ab$transcript_id[Mouse_ab$transcript_id=="ENST00000405460.9"]<-"06_1VLGR1-b_ENST00000405460.9"
Mouse_ab$transcript_id[Mouse_ab$transcript_id=="ENST00000638510.1"]<-"06_2VLGR1-a_ENST00000638510.1"
# One isoform published in McMillan did not have an immediate human orthologue in the gencode dataset. And xlsx file containing the transcript isoform was created manually. 
ADGRV1_McMillan_c <- read.xlsx("USH Genes/ADGRV1/ADGRV1_McMillan_c.xlsx", "McMillan_c")
ADGRV1_McMillan_c$origin <- "McMillan_c"

#Select transcripts to produce McMillan 2002 reference isoforms, add origin column & rename ENST (proir to making GIO_combi, identical ENST in multiple datasets causes an error) 
Mouse_ab <- refseq_gencode   %>% dplyr::filter(transcript_id %in% c("ENST00000405460.9"))
Mouse_ab$origin <- "Mouse_ab"
Mouse_ab$transcript_id[Mouse_ab$transcript_id=="ENST00000405460.9"]<-"06_1VLGR1-b_ENST00000405460.9"
Mouse_ab$transcript_id[Mouse_ab$transcript_id=="ENST00000638510.1"]<-"06_2VLGR1-a_ENST00000638510.1"


# create orthologues of mouse transcript isoforms published in Yagi et al (2005) - Vlgr1d and Vlgr1e
# both are created manually and uploaded as xlsx file
ADGRV1_Yagi_d <- read.xlsx("USH Genes/ADGRV1/Adgrv1_Yagi_d.xlsx", "Yagi_d")
ADGRV1_Yagi_d$origin <- "Yagi_d"
ADGRV1_Yagi_e <- read.xlsx("USH Genes/ADGRV1/Adgrv1_Yagi_e.xlsx", "Yagi_e")
ADGRV1_Yagi_e$origin <- "Yagi_e"

# create orthologues of mouse transcript isoforms published in Skradski 2001 orthologous - Mass1.1, Mass1.2 and Mass1.3 
# End of all Mass1 isoforms annoated based on 3'extension of exon 39 as seen in mouse
# Start site of Mass1.1 annotated based on murine 5'extension of exon 5
# Start site of Mass1.2 annotated based on EST DA391827
# Start site of Mass1.3 annotated based on murine 5'extension of exon 26
# all are created manually and uploaded as a single xlsx file
Adgrv1_Skradski <- read.xlsx("USH Genes/ADGRV1/ADGRV1_Mass1.xlsx", "Mass combined")
Adgrv1_Skradski$origin <- "Skradski"

# isoquant data
# rename ENSTXXX transcripts in dataset to prevent redundant names. 
ADGRV1_CCS3$origin <- "1-3 CCS3" 
# ADGRV1_CCS3 <- ADGRV1_CCS3 #%>% dplyr::filter(transcript_id %in% c("ENST00000638316.1", "ENST00000639884.1", "ENST00000640109.1", "ENST00000640281.1", "transcript10131.chr5.nic", "transcript10139.chr5.nic", "transcript10171.chr5.nnic"))
ADGRV1_CCS3$transcript_id[ADGRV1_CCS3$transcript_id=="ENST00000640109.1"]<-"tAtlas_ENST00000640109.1"
ADGRV1_CCS3$transcript_id[ADGRV1_CCS3$transcript_id=="ENST00000640281.1"]<-"tAtlas_ENST00000640281.1"
ADGRV1_CCS3$transcript_id[ADGRV1_CCS3$transcript_id=="ENST00000638316.1"]<-"tAtlas_ENST00000638316.1"
ADGRV1_CCS3$transcript_id[ADGRV1_CCS3$transcript_id=="ENST00000639884.1"]<-"tAtlas_ENST00000639884.1"


# combine df_columns for figure
ADGRV1_fig_full <- bind_rows(refseq_gencode_clean, ADGRV1_McMillan_ab, ADGRV1_McMillan_c, mouse_C, Mouse_ab, ADGRV1_Yagi_d, ADGRV1_Yagi_e, Adgrv1_Skradski, ADGRV1_CCS3)  


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

  
  geom_range(
    aes(fill = origin)
  ) + scale_fill_manual(values=c("#F8766D", "#7CAE00", "#00BFC4", "#00BFC4","#999999", "#999999","#D3D3D3","#999999", "#999999")) + 
  
  geom_intron(
    data = GOI_rescaled_introns,
    aes(strand = strand), 
    arrow.min.intron.length = 700
  ) + 
  theme_classic(
  ) 


# further processing of the plot was manually in adobe illustrator
