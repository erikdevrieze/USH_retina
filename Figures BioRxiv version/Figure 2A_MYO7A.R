# Packages needed:
#   dplyr, ggplot2, ggtranscript, rtracklayer

# Load packages(ev"dplyr"#Load packages(every time)
library(dplyr) #include
library(ggplot2) #include
library(ggtranscript) # include
library(rtracklayer) #include


setwd("//")


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


# filter GTF files for the gene of interest ("MYO7A" is used for the gencode dataset, "ENSG00000137474.23" to filter in the isoquant data)
gene_of_interest <- "ENSG00000137474.23"        

refseq_gencode <- gencode %>%   # creating a dataset to plot gencode reference isoform
  dplyr::filter(
    !is.na(gene_name), 
    gene_name == gene_of_interest
  ) 

MYO7A_CCS3 <- s1_3_isoQ %>%    # creating a dataset to plot isoquant isoforms
  dplyr::filter(
    !is.na(gene_name), 
    gene_name == gene_of_interest
  ) 

refseq_gilmore <- gencode %>%   # creating a dataset to plot the transcript isoforms proposed by Gilmore et al (2023) - based on gencode data)
  dplyr::filter(
    !is.na(gene_name), 
    gene_name == gene_of_interest)

refseq_msHC <- s1_3_isoQ %>%   # creating a dataset to plot the orthologous transcript isoforms of mouse auditory hair cells as published by Li et al (2020) - based on isoquant data)
  dplyr::filter(
    !is.na(gene_name), 
    gene_name == gene_of_interest)

#Select MANE select/clinical as reference isoforms from Gencode dataset 
#rename ENSTXXX transcripts in dataset to prevent redundant names. 
refseq_gencode_clean <- refseq_gencode %>% dplyr::filter(transcript_id %in% c("ENST00000409709.9"))
refseq_gencode_clean$origin <- "gencode"
refseq_gencode_clean$transcript_id[refseq_gencode_clean$transcript_id=="ENST00000409709.9"]<-"1_Gencode_ENST00000409709.9"
MYO7A_CCS3$origin <- "1-3 CCS3"
MYO7A_CCS3$transcript_id[MYO7A_CCS3$transcript_id=="ENST00000409709.9"]<-"tAtlas_ENST00000409709.9"
MYO7A_CCS3$transcript_id[MYO7A_CCS3$transcript_id=="ENST00000458637.6"]<-"tAtlas_ENST00000458637.6"

# selection of isoforms discussed in Gilmore et al from the Gencode dataset. Gilmore 2023 refers to Weil 1996 for original identification of isoforms. 
# Gilmore mentions 2 variants of full lenght MYO7A -> have to assume that IF1 and IF2 are annotated on ENST00000409709.9, implying that the short IF2 = ENST00000458637.6
# rename ENSTXXX transcripts in dataset to prevent redundant names. 
refseq_gilmore <- refseq_gilmore %>% dplyr::filter(transcript_id %in% c("ENST00000409709.9","ENST00000458637.6"))
refseq_gilmore$origin <- "gilmore"
refseq_gilmore$transcript_id[refseq_gilmore$transcript_id=="ENST00000409709.9"]<-"Gilmore IF1"
refseq_gilmore$transcript_id[refseq_gilmore$transcript_id=="ENST00000458637.6"]<-"Gilmore IF2"
 
# selection of mouse hair cell (msHC) isoforms discussed in Li et al 2020. 
# orthologous transcripts are ENST00000409709.9 (refseq) and transcript20052.chr11.nic take from the isoquant dataset
# We assume that the mouse short orthologue has the same 3' UTR as the long orthologue. Hence, we reduce the UTR lenght of the transcript20052.chr11.nic to produce the transcript for the figure.
# Rename ENSTXXX transcripts in dataset to prevent redundant names. 
refseq_msHC <- refseq_msHC %>% dplyr::filter(transcript_id %in% c("ENST00000409709.9","transcript20052.chr11.nic"))
refseq_msHC$origin <- "msHC"
refseq_msHC$transcript_id[refseq_msHC$transcript_id=="ENST00000409709.9"]<-"Myo7a-C orthologue"
refseq_msHC$transcript_id[refseq_msHC$transcript_id=="transcript20052.chr11.nic"]<-"Myo7a-S orthologue"
refseq_msHC[99,3] = 77215241
refseq_msHC[99,4] = 635



# combine df_columns for figure
MYO7A_fig_full <- bind_rows(refseq_gencode_clean, refseq_gilmore, refseq_msHC, MYO7A_CCS3 ) 


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
  ) + scale_fill_manual(values=c("#F8766D", "#7CAE00", "#00BFC4", "#999999")) +
  geom_intron(
    data = MYO7A_rescaled_introns,
    aes(strand = strand), 
    arrow.min.intron.length = 700
  ) + 
  theme_classic(
  ) + 
  scale_y_discrete(limits = rev)


# further processing of the plot was manually in adobe illustrator
