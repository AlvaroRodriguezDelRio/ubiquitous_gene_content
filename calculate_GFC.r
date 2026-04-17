#install.packages("tidytree")
library(dplyr)
library(data.table)
library(stringr)
library(dplyr)
library(tidyverse)
library(ggridges)
library(ggpubr)
library(patchwork)
library(ggplot2)
library(ggnewscale)
library(purrr)



# get ko conservation per lineage 
# this table also includes number of detections per lineage 
abs = read.table("ko_conservation.tab",sep ='\t',header = F)
names(abs) = c("samples","totsamples","tip","ko","cons","gene","desc","general")

# get lineage ubiquity
abs$ubiq = abs$samples / abs$totsamples

# load n genomes per linenage
n_genomes_per_lin = read.table("n_genomes_per_lin.tab",sep = '\t')
names(n_genomes_per_lin) = c("tip","n_genomes")
head(n_genomes_per_lin)

abs = abs %>%
  left_join(n_genomes_per_lin,by = "tip")

# filter n genomes  > 10
nrow(abs)
abs = abs %>%
  filter(n_genomes >10)

# filter genes to analyze
abs = abs %>%
  filter(!general %in% c('Dormancy','Nitrate transporter','Nitrilase','Organic N metabolism',
                         'Organic N metabolism / Ammonification','Photosynthesis',
                         'Photosynthesis - antenna proteins'))

# complete missing lin - gene combinations
head(abs)
abs_comp = abs %>%
  group_by(samples,totsamples,dn,ko,desc,general,ubiq,f2,n_genomes) %>%
  complete(tip, gene, fill = list(cov = 0,
                                  sp  = 0)) %>%
  ungroup()

######
# join perc genomes refseq data
######

ref_genomes = read.table("lin2n_refseq.tab",sep = '\t')
names(ref_genomes) = c("tip","genome_orig","n_genomes_gtdb")

ref_genomes = ref_genomes %>%
  group_by(tip) %>%
  mutate(tot = sum(n_genomes_gtdb)) %>%
  mutate(proportion_GB_genomes = n_genomes_gtdb / tot) %>%
  filter(genome_orig == 'GB') %>%
  dplyr::select(tip,proportion_GB_genomes)

abs = abs %>%
  left_join(ref_genomes,by = 'tip')

#####
# join genome completeness and contamination, per lineage
#####

comcont = read.table('comp_cont_per_lin.tab',sep = '\t')
head(comcont)
names(comcont) = c('tip','mean_compl','median_compl','mean_cont','median_cont')

comcont = comcont %>%
  dplyr::left_join(ref_genomes,by = 'tip')


comcont = comcont %>%
  mutate(uncultivated = case_when(proportion_GB_genomes == 1~'Uncultivated',
                                  .default = 'cultivated representative'))

comcont = comcont %>%
  dplyr::select(tip,mean_compl,mean_cont)

abs = abs %>%
  left_join(comcont,by = 'tip')

# correct gene conservation by completeness
abs$cons = as.numeric(abs$cons)
abs$cons_corrected = abs$cons * (100 / abs$mean_compl)


# scale gene conservation to 1
abs = abs %>%
  group_by(ko) %>%
  mutate(max_cons_corr = max(cons_corrected)) %>%
  mutate(cons_corrected = cons_corrected / max_cons_corr) %>%
  dplyr::select(!max_cons_corr)

# GFC calculation 
abs$f2_corrected = 2 * (abs$cons_corrected *abs$ubiq) / (abs$cons_corrected + abs$ubiq)

####
# export data
####

data_export = abs %>%
  dplyr::select(ubiq,tip,ko,gene,desc,general,proportion_GB_genomes,cons_corrected,GFC)

names(data_export) = c("ubiquity",'taxa','ko','conservation','gene code','gene description',
                       'general','n_genomes','proportion_uncult_genomes',
                       'conservation_corrected_by_completness','GFC')


write.csv(data_export,"ko_TGR.csv")
