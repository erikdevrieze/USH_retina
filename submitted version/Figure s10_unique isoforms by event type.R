
#Load packages(every time)
library(dplyr) #for data manipulation (organisation, grouping, counting, etc)
library(readxl) #read excel files
library(readr)
library(data.table)
library(splitstackshape)


#See which working directory (wd) is used atm
getwd()
#Set your working directory
setwd("/Users/erikdevrieze/Library/CloudStorage/OneDrive-Radboudumc/z918116/UMCN/Manuscripts/2024 - USH isoseq manuscript/Analysis paper 2024")


###################################################################################################################

# read count script of Figure 1E using sorted transcript model grouped counts file as input

###################################################################################################################



######## with isoquant output based on samples 1-3 ################################################################

#read atlas 00_aln.sorted.transcript_model_reads 
AT_transcript_model_counts <- fread("dataset Atlas/00_aln.sorted.transcript_model_grouped_counts.tsv") 
colnames(AT_transcript_model_counts)[1] ="transcript_id"

s1_AT_transcript_model_counts <- AT_transcript_model_counts %>% dplyr::filter(AT_transcript_model_counts$s1 >= 1) 
s2_AT_transcript_model_counts <- AT_transcript_model_counts %>% dplyr::filter(AT_transcript_model_counts$s2 >= 1) 
s3_AT_transcript_model_counts <- AT_transcript_model_counts %>% dplyr::filter(AT_transcript_model_counts$s3 >= 1) 



#Grab, count and plot ENST matching reads:
ENST_s1 <- dplyr::filter(s1_AT_transcript_model_counts, grepl("ENST", transcript_id))
ENST_s2 <- dplyr::filter(s2_AT_transcript_model_counts, grepl("ENST", transcript_id))
ENST_s3 <- dplyr::filter(s3_AT_transcript_model_counts, grepl("ENST", transcript_id))


#Grab, count and plot NNIC matching reads:
NNIC_s1 <- dplyr::filter(s1_AT_transcript_model_counts, grepl(".nnic", transcript_id))
NNIC_s2 <- dplyr::filter(s2_AT_transcript_model_counts, grepl(".nnic", transcript_id))
NNIC_s3 <- dplyr::filter(s3_AT_transcript_model_counts, grepl(".nnic", transcript_id))

#Grab, count and plot NIC matching reads:
NIC_s1 <- dplyr::filter(s1_AT_transcript_model_counts, grepl(".nic", transcript_id)) # captures both NIC and NNIC (calculated NIC manually (NIC = total .nic and .nnic combined - NNIC))
NIC_s2 <- dplyr::filter(s2_AT_transcript_model_counts, grepl(".nic", transcript_id))
NIC_s3 <- dplyr::filter(s3_AT_transcript_model_counts, grepl(".nic", transcript_id))



######## and now for the isoquant dataset generated with the short reads from both the short (samples 1-3) and long (sample 4) prep as input. ######## 


#read atlas 00_aln.sorted.transcript_model_reads 
s1_4_transcript_model_counts <- fread("dataset 1_CCS3/Isoquant/00_aln.sorted.transcript_model_grouped_counts.tsv") 
colnames(s1_4_transcript_model_counts)[1] ="transcript_id"

s1_14_AT_transcript_model_counts <- s1_4_transcript_model_counts %>% dplyr::filter(s1_4_transcript_model_counts$s1 >= 1) 
s2_14_AT_transcript_model_counts <- s1_4_transcript_model_counts %>% dplyr::filter(s1_4_transcript_model_counts$s2 >= 1) 
s3_14_AT_transcript_model_counts <- s1_4_transcript_model_counts %>% dplyr::filter(s1_4_transcript_model_counts$s3 >= 1) 
s4_14_AT_transcript_model_counts <- s1_4_transcript_model_counts %>% dplyr::filter(s1_4_transcript_model_counts$s2_long >= 1) 

#Grab, count and plot ENST matching reads:
ENST_s1_14 <- dplyr::filter(s1_14_AT_transcript_model_counts, grepl("ENST", transcript_id))
ENST_s2_14 <- dplyr::filter(s2_14_AT_transcript_model_counts, grepl("ENST", transcript_id))
ENST_s3_14 <- dplyr::filter(s3_14_AT_transcript_model_counts, grepl("ENST", transcript_id))
ENST_s4_14 <- dplyr::filter(s4_14_AT_transcript_model_counts, grepl("ENST", transcript_id))

#Grab, count and plot NNIC matching reads:
NNIC_s1_14 <- dplyr::filter(s1_14_AT_transcript_model_counts, grepl(".nnic", transcript_id))
NNIC_s2_14 <- dplyr::filter(s2_14_AT_transcript_model_counts, grepl(".nnic", transcript_id))
NNIC_s3_14 <- dplyr::filter(s3_14_AT_transcript_model_counts, grepl(".nnic", transcript_id))
NNIC_s4_14 <- dplyr::filter(s4_14_AT_transcript_model_counts, grepl(".nnic", transcript_id))


#Grab, count and plot NNIC and NIC matching reads:
NIC_s1_14 <- dplyr::filter(s1_14_AT_transcript_model_counts, grepl(".nic", transcript_id)) # captures both NIC and NNIC (calculated NIC manually (NIC = total .nic and .nnic combined - NNIC))
NIC_s2_14 <- dplyr::filter(s2_14_AT_transcript_model_counts, grepl(".nic", transcript_id))
NIC_s3_14 <- dplyr::filter(s3_14_AT_transcript_model_counts, grepl(".nic", transcript_id))
NIC_s4_14 <- dplyr::filter(s4_14_AT_transcript_model_counts, grepl(".nic", transcript_id))
