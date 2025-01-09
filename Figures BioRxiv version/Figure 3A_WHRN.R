
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
gene_of_interest <- "WHRN"       #  "WHRN" in Gencode GTF & "ENSG00000095397.16"   in IsoQuant GTF

refseq_gencode <- gencode %>% 
  dplyr::filter(
    !is.na(gene_name), 
    gene_name == gene_of_interest) 

WHRN_CCS3 <- s1_3_isoQ %>% 
  dplyr::filter(
    !is.na(gene_name), 
    gene_name == gene_of_interest) 


#Select MANE select/clinical REFSEQs, add origin column & rename ENST (proir to making GIO_combi df, identical ENST in multiple datasets causes an error) 
refseq_gencode_clean <- refseq_gencode %>% dplyr::filter(transcript_id %in% c("ENST00000362057.4")) 
refseq_gencode_clean$origin <- "gencode"
refseq_gencode_clean$transcript_id[refseq_gencode_clean$transcript_id=="ENST00000362057.4"]<-"01_Gencode_ENST00000362057.4"

#Select validated retinal isoforms for Whirlin from the Gencode dataset (selection based on van Wijk et al, 2006) 
WHRN_Wijk2006 <- refseq_gencode %>% dplyr::filter(transcript_id %in% c("ENST00000362057.4", "ENST00000674048.1"))
WHRN_Wijk2006$origin <- "WHRN_human retina Wijk2006"
WHRN_Wijk2006$transcript_id[WHRN_Wijk2006$transcript_id=="ENST00000362057.4"]<-"02_WHRN-FL_ENST00000362057.4"
WHRN_Wijk2006$transcript_id[WHRN_Wijk2006$transcript_id=="ENST00000674048.1"]<-"03_WHRN-C_(S2)_ENST00000674048.1"

#Select proposed orthologu=ous retinal isoforms for Whirlin from the Gencode dataset (based on published mouse data, Mburu et al., 2003; Belyantseva et al., 2005; Ebrahim et al., 2016) 
WHRN_proposed <- refseq_gencode %>% dplyr::filter(transcript_id %in% c("ENST00000362057.4", "ENST00000674048.1"))
WHRN_proposed$origin <- "WHRN_mouse orthologs"
WHRN_proposed$transcript_id[WHRN_proposed$transcript_id=="ENST00000362057.4"]<-"10_WHRN-FL_ENST00000362057.4"
WHRN_proposed$transcript_id[WHRN_proposed$transcript_id=="ENST00000674048.1"]<-"11_WHRN-C_(S2)_ENST00000674048.1"

#IsoQ data s1-3 rename ENST and add origin column
WHRN_CCS3$origin <- "1-3 CCS3"
WHRN_CCS3$transcript_id[WHRN_CCS3$transcript_id=="ENST00000265134.10"]<-"t3Atlas_ENST00000265134.10"
WHRN_CCS3$transcript_id[WHRN_CCS3$transcript_id=="ENST00000362057.4"]<-"t1Atlas_ENST00000362057.4"
WHRN_CCS3$transcript_id[WHRN_CCS3$transcript_id=="ENST00000374057.3"]<-"t2Atlas_ENST00000374057.3"



# combine df_columns for figure
WHRN_fig_full <- bind_rows(refseq_gencode_clean, WHRN_Wijk2006, WHRN_CCS3, WHRN_proposed )  


# extract the required annotation columns 
WHRN_fig <- WHRN_fig_full %>% 
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
GOI_exons <- WHRN_fig %>% dplyr::filter(type == "exon")


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
  ) + scale_fill_manual(values=c("#00BFC4", "#7CAE00","#F8766D", "#999999")) +
  geom_intron(
    data = GOI_rescaled_introns,
    aes(strand = strand), 
    arrow.min.intron.length = 700
  ) + 
  theme_classic(
  ) 



# additional annotations were made manually in Adobe Illustrator. 
#The mouse ortholog of FL whirlin with exon 7B was copied manually from transcript.nnic13724.chr9. 
#The mouse ortholog of whirlin N was produced by copying FL transcript and removing intron 4 and beyond

