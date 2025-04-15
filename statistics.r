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


# load n genomes per lin 
n_genomes_per_lin = read.table("n_genomes_per_lin.tab",sep = '\t')
names(n_genomes_per_lin) = c("tip","n_genomes")
head(n_genomes_per_lin)

# get ko coverage per lin
abs = read.table("ko_list.cov_sp.tab",sep ='\t',header = F)
names(abs) = c("samples","totsamples","tip","ko","cov","sp","gene_name","desc","general")
head(abs)

abs$ubiq = abs$samples / abs$totsamples
abs$f2 = 2 * (abs$cov *abs$ubiq) / (abs$cov + abs$ubiq)

# join n genomes data to gene coverage data 
abs = abs %>% 
  left_join(n_genomes_per_lin,by = "tip")

# join perc genomes refseq data
ref_genomes = read.table("GTDB_perc_uncultivated/lin2n_refseq.tab",sep = '\t')
names(ref_genomes) = c("tip","genome_orig","n_genomes_gtdb")

ref_genomes = ref_genomes %>% 
  group_by(tip) %>% 
  mutate(tot = sum(n_genomes_gtdb)) %>% 
  mutate(proportion_GB_genomes = n_genomes_gtdb / tot) %>% 
  filter(genome_orig == 'GB') %>% 
  dplyr::select(tip,proportion_GB_genomes)

abs = abs %>% 
  left_join(ref_genomes,by = 'tip')

# export data with ubiquity and gene content information
data_export = abs %>% 
  mutate(ubiq = samples / totsamples) %>% 
  dplyr::select(ubiq,tip,ko,cov,sp,gene,desc,general,f2,n_genomes,proportion_GB_genomes)
  
names(data_export) = c("ubiquity",'taxa','ko','conservation','exclusivity','gene code','gene description',
                       'general','TGR','n_genomes','proportion_uncult_genomes')

write.csv(data_export,"ko_TGR.csv")


################
######
# general stats
#####
###############

# species in less than 5% samples
totsp = abs %>% 
  filter(grepl("s__",tip)) %>% 
  dplyr::select(tip) %>% 
  unique() %>% 
  nrow()

sp_men_1 = abs %>% 
  filter(grepl("s__",tip) & samples / totsamples < 0.05) %>% 
  dplyr::select(tip) %>% 
  unique() %>% 
  nrow()

sp_men_1 / totsp

## no species in more than 25% soil samples
abs %>% 
  filter(grepl("s__",tip) & samples / totsamples > 0.25) %>% 
  dplyr::select(tip) %>% 
  unique() %>% 
  nrow()


# orders> 50% samples
abs %>% 
  filter(grepl("o__",tip) & !grepl("f__",tip)) %>%
  filter(samples / totsamples > 0.5) %>% 
  dplyr::select(tip) %>% 
  unique() %>% 
  nrow()

# phyla > 90% samples
abs %>% 
  filter(grepl("p__",tip) & !grepl("c__",tip)) %>%
  filter(samples / totsamples > 0.9) %>% 
  dplyr::select(tip) %>% 
  unique() %>% 
  nrow()


# Patescibacteria example
abs %>% 
  filter(grepl("p__",tip) & !grepl("c__",tip)) %>%
  filter(grepl('Patescibacteria',tip)) %>% 
  filter(samples / totsamples > 0.9) %>% 
  dplyr::select(tip) %>% 
  unique() 

abs %>% 
  filter(grepl("g__",tip)) %>%
  filter(grepl('Patescibacteria',tip)) %>% 
  filter(samples / totsamples > 0.4) %>% 
  dplyr::select(tip) %>% 
  unique() 


############
# gene content relation with perc of genomes in refseq
############

head (abs)
da = abs %>% 
  group_by(gene,general) %>% 
  filter(gene %in% c("nifH","nifD","nifK","anfG",
                     "pmo-amoA","pmo-amoB","pmo-amoC",
                     "narG, narZ, nxrA","narH, narY, nxrB","napA","napB",
                     "norB","norC","NosZ","nirS","nirK",
                     "rbcL","rbcS","prkB",
                     "K18603","K18604","abfD","aclA","aclB",
                     "mcrA","mcrB","mcrC","mcrD","fwdA","fwdB","fwdC",
                     "mmoX","mmoC","mmoY","mmoZ","pmo-amoA","pmo-amoB","pmo-amoC",
                     "coxA","coxB","coxC",
                     "coxS","coxM","coxS")) %>%
  summarise(p = cor.test(proportion_GB_genomes,cov)$p.value,
            c = cor.test(proportion_GB_genomes,cov)$estimate) %>% 
  arrange(c) %>% 
  filter(p < 0.01 & (c) > 0) %>% 
  print(n = 200)

ggplot(abs %>% 
         filter(gene %in% da$gene) %>% 
         filter(!grepl('c__',tip)))+
  geom_hex(aes(x = proportion_GB_genomes,y = cov))+
  geom_smooth(aes(x = proportion_GB_genomes,y = cov),method = 'lm')+
  facet_grid(~gene)

######
# corr genes // ubiq
######

head (abs)
da = abs %>% 
  group_by(gene,general) %>% 
  filter(gene %in% c("nifH","nifD","nifK","anfG",
                     "pmo-amoA","pmo-amoB","pmo-amoC",
                     "narG, narZ, nxrA","narH, narY, nxrB","napA","napB",
                     "norB","norC","NosZ","nirS","nirK",
                     "rbcL","rbcS","prkB",
                     "K18603","K18604","abfD","aclA","aclB",
                     "mcrA","mcrB","mcrC","mcrD","fwdA","fwdB","fwdC",
                     "mmoX","mmoC","mmoY","mmoZ","pmo-amoA","pmo-amoB","pmo-amoC",
                     "coxA","coxB","coxC",
                     "coxS","coxM","coxS")) %>%
  filter(grepl("g__",tip) & !grepl("s__",tip)) %>% 
  summarise(p = cor.test(samples,cov)$p.value,
            c = cor.test(samples,cov)$estimate) %>% 
  arrange(-c)%>% 
  filter(p < 0.01 & abs(c) > 0.2) %>% 
  print(n = 200)
da

ggplot(abs %>% 
         filter(gene %in% da$gene) %>% 
         filter(!grepl('c__',tip)))+
  geom_hex(aes(x = samples,y = cov))+
  geom_smooth(aes(x = samples,y = cov),method = 'lm')+
  facet_grid(~gene)


#######
# Carbon fixation genes
######

# rbcLS, cons, phylum
abs %>% 
  filter(n_genomes > 10) %>% 
  filter(grepl("Carbon fixataion",general)) %>%
  filter(gene %in% c('rbcL','rbcS','prkB')) %>%
  filter(!grepl('c__',tip)) %>% 
  group_by(tip,gene) %>% 
  summarise(n = cov,
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 20) 

# rbcL
abs %>% 
  filter(n_genomes > 10) %>% 
  filter(grepl("Carbon fixataion",general)) %>%
  filter(gene == 'rbcL') %>%
  group_by(tip,gene) %>% 
  summarise(n = f2,
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 20) 


# rbcL, phylum
abs %>% 
  filter(n_genomes>10) %>% 
  filter(grepl("p__",tip)) %>% 
  filter(gene == 'rbcL') %>%
  filter(!grepl("c__",tip)) %>% 
  dplyr::select(gene,general,tip,cov,sp,f1,f2,proportion_GB_genomes) %>% 
    group_by(tip,gene) %>% 
    summarise(n = sum(f2),
              proportion_GB_genomes = proportion_GB_genomes) %>% 
    arrange(desc(n))

# rbcl, genus
abs %>% 
  filter(n_genomes>10) %>% 
  filter(grepl("p__",tip)) %>% 
  filter(gene == 'rbcL') %>%
  filter(grepl("g__",tip)) %>% 
  dplyr::select(gene,general,tip,cov,sp,f1,f2,proportion_GB_genomes) %>% 
  group_by(tip,gene) %>% 
  summarise(n = sum(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n))


# rbcl, all
abs %>% 
  filter(n_genomes>10) %>% 
  filter(grepl("p__",tip)) %>% 
  filter(gene == 'rbcL') %>%
  dplyr::select(gene,general,tip,cov,sp,f1,f2,proportion_GB_genomes) %>% 
  group_by(tip,gene) %>% 
  summarise(n = sum(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 100)


# rbcl + rbcs
abs %>% 
  filter(n_genomes>10) %>% 
  filter(grepl("p__",tip)) %>% 
  filter(gene %in% c('rbcL','rbcS')) %>%
  dplyr::select(gene,general,tip,cov,sp,f1,f2,proportion_GB_genomes) %>% 
  group_by(tip,general) %>% 
  summarise(n = min(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  unique() %>% 
  arrange(desc(n)) 



########
# Carbon fixation, others
#######

# acetyl-CoA/propionyl-CoA carboxylase, cov
abs %>% 
  filter(n_genomes>10) %>% 
  filter(!grepl("o__",tip)) %>%
  filter(gene %in% c('K18604','K18903')) %>% 
  group_by(tip,gene) %>% 
  summarise(n = min(cov),
            proportion_GB_genomes = proportion_GB_genomes) %>% # sp for specificity 
  arrange(desc(n)) %>% 
  print(n = 20)

# acetyl-CoA/propionyl-CoA carboxylase
abs %>% 
  filter(n_genomes>10) %>% 
  filter(!grepl("Photosynthesis",general)) %>%
  filter(gene %in% c('K18604','K18903')) %>% 
  group_by(tip) %>% 
  summarise(n = min(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% # sp for specificity 
  arrange(desc(n)) %>% 
  print(n = 20)


# abfD
abs %>% 
  filter(n_genomes>10) %>%
  filter(gene == 'abfD') %>% 
  dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2,proportion_GB_genomes)%>% 
  group_by(tip) %>% 
  summarise(n = sum(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n  = 20)


# AclA / AclB, cov
abs %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c('aclA','aclB')) %>% 
  filter(!grepl("o__",tip)) %>%
  dplyr::select(gene,general,gene,tip,cov,sp,f1,f3,f2,proportion_GB_genomes) %>% 
  group_by(tip) %>% 
  summarise(n = min(cov),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 20)


# AclA / AclB
abs %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c('aclA','aclB')) %>% 
  dplyr::select(gene,general,gene,tip,cov,sp,f1,f3,f2,proportion_GB_genomes) %>% 
  group_by(tip,gene) %>% 
  summarise(n = min(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 20)


# cooS, acsA	anaerobic carbon-monoxide dehydrogenase catalytic subunit
abs %>% 
  filter(!grepl("Photosynthesis",general)) %>%
  filter(grepl("p__",tip)) %>% 
  filter(gene == 'cooS, acsA') %>% 
  dplyr::select(gene,general,tip,cov,sp,f1,f3,f2,proportion_GB_genomes)%>% 
  group_by(tip) %>% 
  summarise(n = sum(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n))

# K01895, acs
abs %>% 
  filter(!grepl("Photosynthesis",general)) %>%
  filter(grepl("p__",tip)) %>% 
  #  filter(f2 > 0.5) %>%
  filter(gene == 'acs') %>% 
  dplyr::select(gene,general,tip,cov,sp,f1,f3,f2,proportion_GB_genomes)%>% 
  group_by(tip) %>% 
  summarise(n = sum(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n))

# tfrB Fumarate reductase
abs %>% 
  filter(!grepl("Photosynthesis",general)) %>%
  filter(grepl("p__",tip)) %>% 
  filter(gene == 'tfrB') %>% 
  dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2,proportion_GB_genomes)%>% 
  group_by(tip) %>% 
  summarise(n = sum(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n))

# korA/B
abs %>% 
  filter(!grepl("Photosynthesis",general)) %>%
  filter(grepl("p__",tip)) %>% 
  #  filter(f2 > 0.2) %>%
  filter(gene %in% c('korA','korB')) %>% 
  dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2,proportion_GB_genomes)%>% 
  group_by(tip,gene) %>% 
  summarise(n = sum(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n))


# por
abs %>% 
  filter(!grepl("Photosynthesis",general)) %>%
  filter(grepl("p__",tip)) %>% 
  #filter(f2 > 0.7) %>%
  filter(gene == 'por') %>% 
  dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2,proportion_GB_genomes)%>% 
  group_by(tip) %>% 
  summarise(n = sum(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n))


######
# Carbon degradation
######

abs %>% 
  filter(grepl("Carbon degradation",general)) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip) %>% 
  summarise(n = sum(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n))

abs %>% 
  filter(grepl("Carbon degradation",general)) %>%
  filter(grepl("p__",tip)) %>% 
  filter(f3 > 0.3) %>%
  dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2)


abs %>% 
  filter(grepl("Carbon degradation",general)) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip) %>% 
  summarise(n = sum(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n))

# cellulases
abs %>% 
  filter(n_genomes>10) %>%
  filter(grepl(c("K01179",
                 "CELB",
                 "bcsZ",
                 "CBH1",
                 "CBH2"),gene)) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip) %>% 
  summarise(n = min(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n  = 30)


#####
# Respiration
#####

# top ranks
abs %>% 
  filter(n_genomes > 10) %>% 
  filter(grepl("Respiration",general)) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip) %>% 
  summarise(n = min(f2),proportion_GB_genomes = proportion_GB_genomes) %>% 
  unique() %>% 
  arrange(desc(n)) %>% 
  print(n = 100)

# phyla with high ubiq and low cov of cox
head(abs)
abs %>% 
  filter(grepl("Respiration",general)) %>%
  filter(grepl("p__",tip)) %>%
  filter(!grepl("c__",tip)) %>%
  filter(cov < 0.1 & ubiq > 0.5) 

#####
# Methanogenesis
#####

# mrcABC, cov
abs %>% 
  filter(n_genomes > 10) %>% 
  filter(!grepl("c__",tip)) %>% 
  filter(gene %in% c('mcrB','mcrA','mcrC')) %>% 
  group_by(tip) %>% 
  summarise(n = min(cov),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  unique() %>% 
  print(n = 20)


# mrcABC,tgr
abs %>% 
  filter(n_genomes > 10) %>% 
  filter(grepl("p__",tip)) %>% 
  filter(gene %in% c('mcrB','mcrA','mcrC')) %>% 
  group_by(tip) %>% 
  summarise(n = min(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>%
  unique() %>% 
  print(n = 100)

abs %>% 
  filter(n_genomes > 10) %>% 
  filter(grepl('Bog',tip)) %>% 
  filter(gene %in% c('mcrB','mcrA','mcrC')) 

######
# Methanotrophy
#####

# other proteins in the cluster, min does not work.. low values in general
abs %>% 
  filter(n_genomes > 10) %>% 
  filter(general == 'Methanotrophy') %>%
  filter(!gene %in% c('mmoX','mmoC','mmoY',"mmoZ")) %>% 
  filter(grepl("p__",tip)) %>% 
  group_by(tip,gene) %>% 
  summarise(n = min(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 100)



# mmoX
abs %>% 
  filter(n_genomes > 10) %>% 
  filter(grepl("p__",tip)) %>% 
  filter(gene == 'mmoX') %>% 
  group_by(tip) %>% 
  summarise(n = sum(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n))


# mmoC
abs %>% 
  filter(grepl("p__",tip)) %>% 
  filter(gene == 'mmoC') %>% 
  group_by(tip) %>% 
  summarise(n = sum(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n))

# max for all subunits 
# other proteins in the cluster
abs %>% 
  filter(general == 'Methanotrophy') %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip) %>% 
  summarise(n = sum(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 100)


n_genomes_per_lin %>% 
  filter(lin == 'd__Bacteria;p__Desulfobacterota_B;c__Binatia;o__Bin18')


# other proteins in the cluster
abs %>% 
  filter(gene == 'mmoB') %>%
  filter(grepl("p__",tip)) %>% 
  dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = f2,n1 = cov) %>% 
  arrange(desc(n1)) %>% 
  print(n = 100)


# pmo-mmo genes
abs %>% 
  filter(gene %in% c('pmo-amoA','pmo-amoB','pmo-amoC')) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip,gene) %>% 
  summarise(n = f2,n1 = cov,
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 100)


# pmo-mmo genes in bacteria
abs %>% 
  filter(n_genomes > 10) %>% 
  filter(gene %in% c('pmo-amoA','pmo-amoB','pmo-amoC')) %>%
  filter(grepl("d__Bacteria",tip)) %>% 
  group_by(tip) %>% 
  summarise(n = min(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>%
  unique() %>% 
  print(n = 100) 



n_genomes_per_lin %>% 
  filter(tip == 'd__Bacteria;p__Desulfobacterota_B;c__Binatia;o__Bin18')




###
# CO oxidation
###

# all
abs %>% 
  filter(n_genomes > 10) %>% 
  filter(gene %in% c('coxM','coxS')) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip) %>% 
  summarise(n = min(f2),
              proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 100)

abs %>% 
  filter(grepl("d__Bacteria",tip)) %>% 
  filter(grepl("CO oxidation",general)) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip,gene) %>% 
  summarise(n = sum(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 100)

# phyla
abs %>% 
  filter(grepl("CO oxidation",general)) %>%
  filter(grepl("p__",tip)) %>% 
  filter(!grepl("c__",tip)) %>% 
  group_by(tip,gene) %>% 
  summarise(n = sum(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 100)

#####################
#############
# nitrogen metabolism
############
#####################

# nitrogen fixation,  nif
abs %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("nifH","nifD","nifK")) %>%
  filter(general %in% c("Nitrogen fixation")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip) %>% 
  summarise(n = min(f2),c = sum(sp),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  unique() %>% 
  arrange(desc(n)) %>% 
  print(n = 50)


# nitrogen fixation,  nif, no cyano
abs %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("nifH","nifD","nifK")) %>%
  filter(general %in% c("Nitrogen fixation")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip) %>% 
  summarise(n = min(f2),c = sum(sp),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  unique() %>% 
  arrange(desc(n)) %>% 
  head(n = 50) %>% 
  filter(grepl('Cyano',tip))

# anfG
abs %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("anfG")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip,gene) %>% 
  summarise(n = sum(f2),
            proportion_GB_genomes = proportion_GB_genomes,
            cov = cov, sample = samples) %>% 
  arrange(desc(n)) %>% 
  print(n = 100)


# nodAC + nif genes
head(abs)
abs %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("nodA","nodC","nifH")) %>%
  filter(general %in% c("Nitrogen fixation")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip,gene) %>% 
  summarise(n = sum(f2),cov = sum(cov),samples = samples,
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 100) 


# pmos
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("pmo-amoA","pmo-amoB","pmo-amoC")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip) %>% 
  summarise(n = min(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  unique() %>% 
  arrange(desc(n)) %>% 
  print(n = 100)


# hao
abs %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("hao")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip,gene) %>% 
  summarise(n = sum(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 100)


# nrxA / nrxB // narGH
abs %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("narG, narZ, nxrA","narH, narY, nxrB")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip) %>% 
  summarise(n = min(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  unique() %>% 
  print(n = 100)



# napA / napB
abs %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("napA","napB")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip) %>% 
  summarise(n = min(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  unique() %>% 
  print(n = 100)


 # nirA
abs %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("nirA")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip,gene) %>% 
  summarise(n = sum(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 100)


# nasA
abs %>% 
  filter(n_genomes>10) %>% 
  filter(ko %in% c("K00372")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip,gene) %>% 
  summarise(n = sum(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 100)


# nasC
abs %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("nasC")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip) %>% 
  summarise(n = sum(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 100)

# nrfAH
abs %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("nrfA","nrfH")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip) %>% 
  summarise(n = sum(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 100)


# nirK 
abs %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("nirK")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip) %>% 
  summarise(n = sum(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 100)

# NirS
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("nirS")) %>%
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = sum(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 100)


# both NirS and nirK
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("nirS","nirK")) %>%
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = min(f2),n1 = min(cov),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 100)

#norBC
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("norB","norC")) %>%
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = min(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  unique() %>% 
  print(n = 50)

# nosZ
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("NosZ")) %>%
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = mean(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 100)



# nod
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("nod")) %>%
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = mean(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 100)


abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("nirB")) %>%
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = mean(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 100)



################
##########
# most important taxa for figures, uncultivated, info
##########
###################

export = abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("nifH","nifD","nifK","anfG",
                     "pmo-amoA","pmo-amoB","pmo-amoC",
                     "narG, narZ, nxrA","narH, narY, nxrB","napA","napB",
                     "norB","norC","NosZ","nirS","nirK",
                     "rbcL","rbcS","prkB",
                     "K18603","K18604","abfD","aclA","aclB",
                     "mcrA","mcrB","mcrC","mcrD","fwdA","fwdB","fwdC",
                     "mmoX","mmoC","mmoY","mmoZ","pmo-amoA","pmo-amoB","pmo-amoC",
                     "coxA","coxB","coxC",
                     "coxS","coxM","coxS")) %>%
  filter(grepl("p__",tip)) %>%
  #filter(proportion_GB_genomes == 1) %>% 
  group_by(tip,gene) %>% 
  summarise(TGR = mean(f2),
            proportion_GB_genomes = proportion_GB_genomes,
            lineage_cov = cov,
            ubiquity = samples / totsamples,
            lineage_specificity = sp) %>%
  #group_by(gene,tip) %>% 
  ungroup() %>% 
  group_by(gene) %>% 
  arrange(desc(TGR),.by_group = T) %>%
  mutate(taxa_rank_order = row_number()) %>% 
  filter(proportion_GB_genomes == 1) %>% 
  slice_head(n = 5) %>%
  print(n = 200)

#################
##########
# Sulfur
##########
#################

# sat, met3
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("sat, met3")) %>%
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = min(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 100)
  

abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("cysM","cysK")) %>%
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = min(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 100)

abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("dsrA","dsrB")) %>%
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = min(f2),
            cov = mean(cov),
            nsamples = samples,
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 100)


abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("soxB")) %>%
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = min(f2),
            cov = mean(cov),
            nsamples = samples,
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 100)


#################
##########
# Phosphorous
##########
#################

# qcd
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("qcd")) %>%
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = min(f2),
            cov = mean(cov),
            nsamples = samples,
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 100)

# 
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("phytase","phoD","phoA")) %>%
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = mean(f2),
            cov = mean(cov),
            nsamples = samples,
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 100)



#################
##########
# Taxa involved in many / few routes
##########
#################

head(abs)
table(abs$general)

# number of times in the top 100 taxa for any gene (e.g. generalists)
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(grepl("p__",tip)) %>% 
  filter(gene %in% c("pmo-amoA","rbcL","mcrA","NosZ","nifH","coxA",
                     "mmoX","narG, narZ, nxrA","napA","norB",'K01179')) %>% 
  dplyr::select(f2,gene,tip,proportion_GB_genomes) %>%
  unique() %>% 
  group_by(gene) %>% 
  arrange(desc(f2)) %>%
  slice_head(n = 100) %>% 
  group_by(tip) %>% 
  summarise(n = n(),
            concatenated_text = paste(gene, collapse = ", "),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  unique() %>% 
  arrange(desc(n)) %>% 
  print(n = 100)

# number of times in the top 100 taxa for only one gene (e.g. specialists)
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(grepl("p__",tip)) %>% 
  filter(gene %in% c("pmo-amoA","rbcL","mcrA","NosZ","nifH","coxA",
                     "mmoX","narG, narZ, nxrA","napA","norB",'K01179')) %>%   dplyr::select(f2,gene,tip) %>%
  unique() %>% 
  group_by(gene) %>% 
  arrange(desc(f2)) %>%
  slice_head(n = 100) %>% 
  group_by(tip) %>% 
  summarise(n = n(),
            concatenated_text = paste(gene, collapse = ", ")) %>% 
  arrange(desc(n)) %>% 
  filter(n == 1) %>% 
  group_by(concatenated_text) %>% 
  summarise(n = n()) %>% 
  arrange(desc(n))




abs %>%   
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(grepl("p__",tip)) %>% 
  filter(gene %in% c("pmo-amoA","rbcL","mcrA","NosZ","nifH","coxA","mmoX","narG, narZ, nxrA","napA","NorB")) %>% 
  filter(f2 > 0.3) %>% 
  dplyr::select(samples,gene,tip) %>%
  unique() %>% 
  group_by(tip,samples) %>% 
  summarise(n = length(unique(gene)),
            all = unique(gene)) %>% 
  arrange(desc(n)) %>% 
  print(n = 100)




# with methane
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(grepl("p__",tip)) %>% 
  filter(gene %in% c("pmo-amoA","rbcL","mcrA","NosZ","nifH","coxA","mmoX","narG, narZ, nxrA","napA","NorB")) %>% 
  filter(f2 > 0.3) %>% 
  dplyr::select(samples,gene,tip) %>%
  unique() %>% 
  group_by(tip,samples) %>% 
  summarise(n = length(unique(gene)),
            all = toString(unique(gene))) %>% 
  arrange(desc(n)) %>% 
  filter(grepl('mcrA',all))


# with pmoA-amoA
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(grepl("p__",tip)) %>% 
  filter(gene %in% c("pmo-amoA","rbcL","mcrA","NosZ","nifH","coxA","mmoX","narG, narZ, nxrA","napA","NorB")) %>% 
  filter(f2 > 0.3) %>% 
  dplyr::select(samples,gene,tip) %>%
  unique() %>% 
  group_by(tip,samples) %>% 
  summarise(n = length(unique(gene)),
            all = toString(unique(gene))) %>% 
  arrange(desc(n)) %>% 
  filter(grepl('pmo-amoA',all))


# specialists
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(grepl("p__",tip)) %>% 
  filter(gene %in% c("pmo-amoA","rbcL","mcrA","NosZ","nifH","coxA","mmoX","narG, narZ, nxrA","napA","NorB")) %>% 
  filter(f2 > 0.3) %>% 
  dplyr::select(samples,gene,tip) %>%
  unique() %>% 
  group_by(tip,samples) %>% 
  summarise(n = length(unique(gene)),
            all = unique(gene)) %>% 
  filter(n == 1) %>% 
  group_by(all) %>% 
  summarise(n = n())


# specialist examples
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(grepl("p__",tip)) %>% 
  filter(gene %in% c("pmo-amoA","rbcL","mcrA","NosZ","nifH","coxA","mmoX","narG, narZ, nxrA","napA","NorB")) %>% 
  filter(f2 > 0.3) %>% 
  dplyr::select(samples,gene,tip) %>%
  unique() %>% 
  group_by(tip,samples) %>% 
  summarise(n = length(unique(gene)),
            all = unique(gene)) %>% 
  filter(n == 1) %>% 
  filter(all %in% c('mcrA','mmoX','napA'))


#####
# n taxa coding for some genes
####

abs %>% 
  filter(n_genomes>10) %>% 
  filter(grepl("p__",tip)) %>% 
  filter(gene %in% c("pmo-amoA","rbcL","mcrA","NosZ","nifH","coxA","mmoX","narG, narZ, nxrA","napA","NorB")) %>% 
  filter(f2 > 0.3) %>% 
  dplyr::select(samples,gene,tip) %>%
  unique() %>% 
  group_by(gene) %>% 
  summarise(n = length(unique(tip))) %>% 
  arrange(desc(n))

