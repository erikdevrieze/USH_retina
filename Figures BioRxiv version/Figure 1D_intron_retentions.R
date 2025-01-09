
#Load packages (not all are needed)
library(dplyr) 
library(ggplot2)
library(readxl) 
library(reshape2) 
library(xlsx) 
library(readr)
library(data.table)
library(splitstackshape)


###################################################################################################################
# using sorted read assignment files from IsoQuant 
###################################################################################################################




#Read isoquant sample 1-3 read assignment file. 
s1_3_sorted_gene_assignments <- fread("dataset Atlas/00_aln.sorted.read_assignments.tsv") 
colnames(s1_3_sorted_gene_assignments)[1] ="run"
s1_3_split <- cSplit(s1_3_sorted_gene_assignments, "run", "/", direction = "wide")

#Read isoquant sample 1-4 read assignment file.  
s1_4_sorted_gene_assignments <- fread("dataset 1_CCS3/isoquant/00_aln.sorted.read_assignments.tsv") 
colnames(s1_4_sorted_gene_assignments)[1] ="run"
s1_4_split <- cSplit(s1_4_sorted_gene_assignments, "run", "/", direction = "wide")


###################################################################################################################
#Produce and save individual datasets (because the completefiles take forever to load)

s1_sorted_gene_assigments <- dplyr::filter(s1_3_sorted_gene_assignments, grepl("m64167e_210819_012015", run)) 
s2_sorted_gene_assigments <- dplyr::filter(s1_3_sorted_gene_assignments, grepl("m64167e_210902_121011", run)) 
s3_sorted_gene_assigments <- dplyr::filter(s1_3_sorted_gene_assignments, grepl("m64167e_210903_121024", run)) 
s4_sorted_gene_assigments <- dplyr::filter(s1_4_sorted_gene_assignments, grepl("m64037e_211216_160915", run))

write_csv(s1_sorted_gene_assigments, "s1_sorted_gene_assigments.csv")               
write_csv(s2_sorted_gene_assigments, "s2_sorted_gene_assigments.csv")  
write_csv(s3_sorted_gene_assigments, "s3_sorted_gene_assigments.csv")  
write_csv(s4_sorted_gene_assigments, "s4_sorted_gene_assigments.csv") 
   
#read individual datasets

s1_sorted_gene_assigments <- fread("s1_sorted_gene_assigments.csv")
s2_sorted_gene_assigments <- fread("s2_sorted_gene_assigments.csv")
s3_sorted_gene_assigments <- fread("s3_sorted_gene_assigments.csv")
s4_sorted_gene_assigments <- fread("s4_sorted_gene_assigments.csv")
                                                                            

#Grab, count intron retentions:

#sample 1
intron_ret_s1 <- dplyr::filter(s1_sorted_gene_assigments, grepl("intron_retention", assignment_events))
sum_irs1 <- intron_ret_s1 %>%
  summarise(n = n())
sum_s1tot <- s1_sorted_gene_assigments %>%
  summarise(n = n())
sum_irs1 / sum_s1tot * 100  #to get percentage

#sample 2
intron_ret_s2 <- dplyr::filter(s2_sorted_gene_assigments, grepl("intron_retention", assignment_events))
sum_irs2 <- intron_ret_s2 %>%
  summarise(n = n())
sum_s2tot <- s2_sorted_gene_assigments %>%
  summarise(n = n())
sum_irs2 / sum_s2tot * 100  #to get percentage

#sample 3
intron_ret_s3 <- dplyr::filter(s3_sorted_gene_assigments, grepl("intron_retention", assignment_events))
sum_irs3 <- intron_ret_s3 %>%
  summarise(n = n())
sum_s3tot <- s3_sorted_gene_assigments %>%
  summarise(n = n())
sum_irs3 / sum_s3tot * 100  #to get percentage

#sample 4
intron_ret_s4 <- dplyr::filter(s4_sorted_gene_assigments, grepl("intron_retention", assignment_events))
sum_irs4 <- intron_ret_s4 %>%
  summarise(n = n())
sum_s4tot <- s4_sorted_gene_assigments %>%
  summarise(n = n())
sum_irs4 / sum_s4tot * 100  #to get percentage

