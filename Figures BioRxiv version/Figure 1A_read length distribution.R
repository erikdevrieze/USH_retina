#Load packages (not all are required for this analysis)
library(dplyr) 
library(readxl) 
library(xlsx) 
library(readr)


#Set your working directory
setwd("/.../Analysis paper 2024/ccs_read_length")


###################################################################################################################

#files with the read lengths (CCS3) extracted from the mapped sequencing data are uploaded to Galaxy:
#https://usegalaxy.org/api/datasets/f9cad7b01a47213594de9d1831d703a2/display?to_ext=tabular
#https://usegalaxy.org/api/datasets/f9cad7b01a47213504a8724179d7c40d/display?to_ext=tabular
#https://usegalaxy.org/api/datasets/f9cad7b01a4721352b6f98a1097595da/display?to_ext=tabular
#https://usegalaxy.org/api/datasets/f9cad7b01a47213558c73ab11d944db9/display?to_ext=tabular

###################################################################################################################

#Sample 1 (human retina regular library prep #1)
#Read distribution in 500nt bins
s1_lenghts <- fread("ccs_length_sample1.txt") 

group_bin.1 <- s1_lenghts %>%
  group_by(group = cut(V2, breaks = seq(0, max(V2), 500))) %>%
  count()

#calculate mean read length ± SD
mean(s1_lenghts$V2)
sd(s1_lenghts$V2)
 

#Sample 2 (human retina regular library prep #2)
#Read distribution in 500nt bins
s2_lenghts <- fread("ccs_length_sample2.txt") 

group_bin.2 <- s2_lenghts %>%
  group_by(group = cut(V2, breaks = seq(0, max(V2), 500))) %>%
  count()

#calculate mean read length ± SD
mean(s2_lenghts$V2)
sd(s2_lenghts$V2)


#Sample 3 (human retina regular library prep #3)
#Read distribution in 500nt bins
s3_lenghts <- fread("ccs_length_sample3.txt") 

group_bin.3 <- s3_lenghts %>%
  group_by(group = cut(V2, breaks = seq(0, max(V2), 500))) %>%
  count()

#calculate mean read length ± SD
mean(s3_lenghts$V2)
sd(s3_lenghts$V2)


#Sample 4 (human retina long library prep, retina #2)
#Read distribution in 500nt bins
s4_lenghts <- fread("ccs_length_sample2_long.txt") 


group_bin.2long <- s4_lenghts %>%
  group_by(group = cut(V2, breaks = seq(0, max(V2), 500))) %>%
  count()

#calculate mean read length ± SD
mean(s4_lenghts$V2)
sd(s4_lenghts$V2)


