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


setwd("~/analysis/Berlin/sandpiper/GTDB tree rep/")

# # load n genomes per lin
n_genomes_per_lin = read.table("n_genomes_per_lin.tab",sep = '\t')
names(n_genomes_per_lin) = c("tip","n_genomes")

# # ####
# # # get ko cov per lin
# # ####

abs = read.table("ko_list.cov.tab",sep ='\t',header = F)
names(abs) = c("samples","totsamples","tip","cov")

# get variables
abs$ubiq = abs$samples / abs$totsamples

# #####
# # join n genomes data to f2 data
# ####
 
abs = abs %>%
   left_join(n_genomes_per_lin,by = "tip")

# filter n genomes > 10
abs = abs %>%
  filter(n_genomes >10)
 
# ######
# # join perc genomes refseq data
# ######

ref_genomes = read.table("../GTDB_perc_uncultivated/lin2n_refseq.tab",sep = '\t')
names(ref_genomes) = c("tip","genome_orig","n_genomes_gtdb")
 
 ref_genomes = ref_genomes %>%
   group_by(tip) %>%
   mutate(tot = sum(n_genomes_gtdb)) %>%
   mutate(proportion_GB_genomes = n_genomes_gtdb / tot) %>%
   filter(genome_orig == 'GB') %>%
   dplyr::select(tip,proportion_GB_genomes)
 
 abs = abs %>%
   left_join(ref_genomes,by = 'tip')

# #####
# # join genome completeness and contamination, per lineage
# #####

 comcont = read.table('../soil_MAGs/comp_cont_per_lin.tab',sep = '\t')
 names(comcont) = c('tip','mean_compl','median_compl','mean_cont','median_cont')
 
 comcont = comcont %>%
   dplyr::select(tip,mean_compl,mean_cont)
 
 
# join with abs
 abs = abs %>%
   left_join(comcont,by = 'tip')
 
# correct coverage by completeness
 abs$cov = as.numeric(abs$cov)
 abs$cov_corrected = abs$cov * (100 / abs$mean_compl)
 
 # # scale to 1, per ko
abs = abs %>% 
   group_by(ko) %>% 
   mutate(max_cov_corr = max(cov_corrected)) %>% 
   mutate(cov_corrected = cov_corrected / max_cov_corr) %>% 
   dplyr::select(!max_cov_corr)

# corrected version of GFC index
abs$f2_corrected = 2 * (abs$cov_corrected *abs$ubiq) / (abs$cov_corrected + abs$ubiq)
 
# ####
# # export data
# ####
 
 data_export = abs %>%
   mutate(ubiq = samples / totsamples) %>%
   dplyr::select(ubiq,tip,ko,cov,sp,gene,desc,general,f2,mean_compl,mean_cont,
                 n_genomes,proportion_GB_genomes,cov_corrected,f2_corrected)
 
 names(data_export) = c("ubiquity",'taxa','ko','conservation','exclusivity','gene code','gene description',
                        'general','GFR','mean_completness',
                        'mean_contamination','n_genomes','proportion_uncult_genomes',
                        'conservation_corrected_by_completness','GFC_corrected_by_completness')
 
 
 write.csv(data_export,"../Figures -V2/Supp/ko_TGR.v3.csv")

###############
####
# load data for gene coverage / tgr
####
##############

abs = read.csv('../Figures -V3/Supp/abs_v3.csv')
head(abs)


####
# cov / cov corrected dist
#####

############
# corrr genes // perc refseq
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
                     "coxS","coxM","coxS"
  )) %>%
  summarise(p = cor.test(proportion_GB_genomes,cov_corrected)$p.value,
            c = cor.test(proportion_GB_genomes,cov_corrected)$estimate) %>% 
  arrange(c) %>% 
  filter(p < 0.01 & (c) > 0) %>% 
  print(n = 200)

ggplot(abs %>% 
         filter(gene %in% da$gene) %>% 
         filter(!grepl('c__',tip)))+
  geom_hex(aes(x = proportion_GB_genomes,y = cov))+
  geom_smooth(aes(x = proportion_GB_genomes,y = cov),method = 'lm')+
  facet_grid(~gene)


abs = abs %>% 
  mutate(uncultivated = case_when(proportion_GB_genomes == 1 ~'yes',
                                  .default = 'no'))


# more coverage in uncultivated 
# cov corrected
abs %>% 
  filter(grepl('g__',tip) & !grepl('s__',tip)) %>% 
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
  group_by(gene) %>% 
  summarise(p = wilcox.test(cov_corrected~uncultivated)$p.value,
            median_uncult = median(cov_corrected[uncultivated == "yes"], na.rm = TRUE),
            median_cult = median(cov_corrected[uncultivated == "no"], na.rm = TRUE),
            mean_uncult = mean(cov_corrected[uncultivated == "yes"], na.rm = TRUE),
            mean_cult = mean(cov_corrected[uncultivated == "no"], na.rm = TRUE),
            .groups = "drop") %>% 
  filter(p < 0.01) %>% 
  arrange((p)) %>% 
  print(n = 100) %>% 
  filter(mean_cult < mean_uncult)

# more coverage in cultivated 
# cov corrected
abs %>% 
  filter(grepl('g__',tip) & !grepl('s__',tip)) %>% 
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
  group_by(gene) %>% 
  summarise(p = wilcox.test(cov_corrected~uncultivated)$p.value,
            median_uncult = median(cov_corrected[uncultivated == "yes"], na.rm = TRUE),
            median_cult = median(cov_corrected[uncultivated == "no"], na.rm = TRUE),
            mean_uncult = mean(cov_corrected[uncultivated == "yes"], na.rm = TRUE),
            mean_cult = mean(cov_corrected[uncultivated == "no"], na.rm = TRUE),
            .groups = "drop") %>% 
  filter(p < 0.01) %>% 
  filter(mean_cult > mean_uncult) %>% 
  arrange((p)) %>% 
  print(n = 20)
  


ggplot(abs %>% 
         filter(grepl('g__',tip) & !grepl('s__',tip)) %>% 
         filter(gene %in% c("nifH","nifD","nifK","anfG",
                            "pmo-amoA","pmo-amoB","pmo-amoC",
                            "narG, narZ, nxrA","narH, narY, nxrB","napA","napB",
                            "norB","norC","NosZ","nirS","nirK",
                            "rbcL","rbcS","prkB",
                            "K18603","K18604","abfD","aclA","aclB",
                            "mcrA","mcrB","mcrC","mcrD","fwdA","fwdB","fwdC",
                            "mmoX","mmoC","mmoY","mmoZ","pmo-amoA","pmo-amoB","pmo-amoC",
                            "coxA","coxB","coxC",
                            "coxS","coxM","coxS")))+
  geom_boxplot(aes(x = uncultivated, y = cov_corrected))+
  facet_wrap(~gene)




######
# corr genes // ubiq
######

head (abs)
da = abs %>% 
  filter(n_genomes > 10) %>% 
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
                     "coxS","coxM","coxS",
                     'phoA',"qcd","soxB","dsrA",
                     "K01179",
                     "CELB",
                     "bcsZ",
                     "CBH1",
                     "CBH2", 'cysM', 'cysK')) %>%
  filter(grepl("g__",tip) & !grepl("s__",tip)) %>% 
  summarise(p = cor.test(samples,cov_corrected,method = 'spearman')$p.value,
            c = cor.test(samples,cov_corrected,method = 'spearman')$estimate) %>% 
  arrange(-c)%>% 
  filter(p < 0.01) %>% 
  print(n = 100)

da


abs %>% 
  filter(n_genomes > 10) %>% 
  filter(gene %in% c("nifH","nifD","nifK","anfG",
                     "pmo-amoA","pmo-amoB","pmo-amoC",
                     "narG, narZ, nxrA","narH, narY, nxrB","napA","napB",
                     "norB","norC","NosZ","nirS","nirK",
                     "rbcL","rbcS","prkB",
                     "K18603","K18604","abfD","aclA","aclB",
                     "mcrA","mcrB","mcrC","mcrD","fwdA","fwdB","fwdC",
                     "mmoX","mmoC","mmoY","mmoZ","pmo-amoA","pmo-amoB","pmo-amoC",
                     "coxA","coxB","coxC",
                     "coxS","coxM","coxS",
                     'phoA',"qcd","soxB","dsrA",
                     "K01179",
                     "CELB",
                     "bcsZ",
                     "CBH1",
                     "CBH2", 'cysM', 'cysK')) %>% 
  filter(grepl("g__",tip) & !grepl("s__",tip)) %>% 
  group_by(gene) %>% 
  summarise(p = cor.test(samples,cov_corrected,method = 'spearman')$p.value,
            c = cor.test(samples,cov_corrected,method = 'spearman')$estimate) %>% 
  arrange(-c)

ggplot(abs %>% 
         filter(gene %in% c("nifH","nifD","nifK","anfG",
                            "pmo-amoA","pmo-amoB","pmo-amoC",
                            "narG, narZ, nxrA","narH, narY, nxrB","napA","napB",
                            "norB","norC","NosZ","nirS","nirK",
                            "rbcL","rbcS","prkB",
                            "K18603","K18604","abfD","aclA","aclB",
                            "mcrA","mcrB","mcrC","mcrD","fwdA","fwdB","fwdC",
                            "mmoX","mmoC","mmoY","mmoZ","pmo-amoA","pmo-amoB","pmo-amoC",
                            "coxA","coxB","coxC",
                            "coxS","coxM","coxS",
                            'phoA',"qcd","soxB","dsrA",
                            "K01179",
                            "CELB",
                            "bcsZ",
                            "CBH1",
                            "CBH2", 'cysM', 'cysK')) %>%
         filter(grepl("g__",tip) & !grepl("s__",tip)) %>% 
         filter(n_genomes > 10))+
  geom_point(aes(x = 100*samples / totsamples,y = cov_corrected),
             alpha = 0.3)+
  #geom_smooth(aes(x = samples / totsamples,y = cov))+
  facet_wrap(~gene,scales = 'free')+
  theme_classic()+
  ylab('Gene conservation')+
  xlab('Ubiquity (% soil samples)')+
  stat_cor(aes(x = samples / totsamples,y = cov_corrected),
           method = "spearman",
           label.x.npc = "left",
           label.y.npc = "top",
           size = 3.2,
           color = 'red',
           p.accuracy = 0.01, r.accuracy = 0.01
  ) 




#######
# Carbon fixation
######

# rbcLS and prkb, cons, phylum, cov corrected
abs %>% 
  filter(n_genomes > 10) %>% 
  filter(grepl("Carbon fixataion",general)) %>%
  filter(gene %in% c('rbcL','rbcS','prkB')) %>%
  filter(!grepl('c__',tip)) %>% 
  group_by(tip,gene) %>% 
  summarise(c = cov_corrected,
            proportion_GB_genomes = proportion_GB_genomes,
            u = samples / totsamples) %>% 
  arrange(desc(c)) %>% 
  print(n = 20) 


# rbcl, all, f2 corrected
abs %>% 
  filter(n_genomes>10) %>% 
  filter(grepl("p__",tip)) %>% 
  #filter(cov > 0.8) %>%
  filter(gene == 'rbcL') %>%
  #  dplyr::select(gene,general,tip,cov,sp,f1,f2,proportion_GB_genomes) %>% 
  group_by(tip,gene) %>% 
  summarise(n = sum(f2_corrected),
            proportion_GB_genomes = proportion_GB_genomes,
            u = samples / totsamples,
            c = cov_corrected) %>% 
  arrange(desc(n)) %>% 
  print(n = 20)


# rbcl + rbcs, f2 corrected
abs %>% 
  filter(n_genomes>10) %>% 
  filter(!grepl("Photosynthesis",general)) %>%
  filter(grepl("p__",tip)) %>% 
  #filter(cov > 0.8) %>%
  filter(gene %in% c('rbcL','rbcS')) %>%
  #  dplyr::select(gene,general,tip,cov_corrected,sp,f2_corrected,proportion_GB_genomes) %>% 
  group_by(tip,general) %>% 
  summarise(n = min(f2_corrected),
            c = min(cov_corrected),
            u = samples / totsamples,
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  unique() %>% 
  arrange(desc(n)) %>% 
  print(n = 20)


abs %>% 
  filter(gene %in% c('rbcL','rbcS')) %>%
  filter(tip == 'd__Bacteria;p__Actinomycetota;c__Thermoleophilia;o__Solirubrobacterales;f__Solirubrobacteraceae;g__CADCVQ01') %>% 
  dplyr::select(cov_corrected,samples,totsamples, gene, ubiq)

########
# Carbon fixation, others
#######

# acetyl-CoA/propionyl-CoAcarboxylase
# f2 corrected
abs %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c('K18604','K18603')) %>% 
  group_by(tip) %>% 
  summarise(n = min(f2_corrected),
            c = min(cov_corrected),
            ub = min(samples/totsamples),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  unique() %>% 
  arrange(desc(n)) %>% 
  print(n = 20)

abs %>% 
  filter(gene %in% c('K18604','K18603')) %>% 
  filter(tip == 'd__Archaea;p__Thermoproteota;c__Nitrososphaeria;o__Nitrososphaerales;f__Nitrososphaeraceae;g__TA-21') %>% 
  dplyr::select(cov_corrected, gene, ubiq)



# abfD f2, cov corrected
abs %>% 
  filter(n_genomes>10) %>%
  filter(gene == 'abfD') %>% 
  group_by(tip) %>% 
  summarise(n = sum(f2_corrected),
            c = min(cov_corrected),
            ub = min(samples/totsamples),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n  = 20)



# AclA / AclB, 
# f2 corrected
abs %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c('aclA','aclB')) %>% 
  ungroup() %>% 
  # dplyr::select(gene,general,gene,tip,cov,sp,f1,f3,f2,proportion_GB_genomes) %>% 
  group_by(tip) %>% 
  summarise(n = min(f2_corrected),
            c = min(cov_corrected),
            ub = min(samples/totsamples),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  unique() %>% 
  arrange(desc(n)) %>% 
  print(n = 40)


abs %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c('aclA','aclB')) %>% 
  filter(tip == 'd__Bacteria;p__Nitrospirota;c__Nitrospiria;o__Nitrospirales;f__Nitrospiraceae;g__Palsa-1315') %>% 
  dplyr::select(cov_corrected,ubiq,gene)


# cooS, acsA	anaerobic carbon-monoxide dehydrogenase catalytic subunit
abs %>% 
  filter(!grepl("Photosynthesis",general)) %>%
  filter(grepl("p__",tip)) %>% 
  #  filter(f2 > 0.5) %>%
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
  #  filter(f2 > 0.2) %>%
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
  #  filter(f2 > 0.9) %>%
  #  filter(!grepl("c__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = sum(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n))



# cellulases
# f2 corrected
abs %>% 
  filter(n_genomes>10) %>%
  filter(grepl(c("K01179",
                 "CELB",
                 "bcsZ",
                 "CBH1",
                 "CBH2"),gene)) %>%
  filter(grepl("p__",tip)) %>% 
  #  filter(f2 > 0.9) %>%
  #  filter(!grepl("c__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = min(f2_corrected),
            proportion_GB_genomes = proportion_GB_genomes,
            c = min(cov_corrected),
            ub = min(samples/totsamples)) %>% 
  arrange(desc(n)) %>% 
  print(n  = 30)


#####
# Respiration
#####

# top ranks
# f2 corrected completeness
abs %>% 
  filter(n_genomes > 10) %>% 
  filter(grepl("Respiration",general)) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip) %>% 
  summarise(n = min(f2_corrected),proportion_GB_genomes = proportion_GB_genomes,
            c = min(cov_corrected),
            ub = min(samples/totsamples)
  ) %>% 
  unique() %>% 
  arrange(desc(n)) %>% 
  print(n = 100)


abs %>% 
  filter(n_genomes > 10) %>% 
  filter(grepl("Respiration",general)) %>%
  filter(tip == 'd__Bacteria;p__Acidobacteriota;c__Vicinamibacteria;o__Vicinamibacterales;f__UBA2999') %>% 
  dplyr::select(cov_corrected, gene, ubiq)



# coxABC correlation with ubiq
# cov corrected
abs %>% 
  filter(n_genomes > 10) %>% 
  filter(grepl("Respiration",general)) %>%
  filter(grepl("g__",tip) & !grepl("s__",tip)) %>% 
  #  group_by(gene) %>% 
  summarise(p = cor.test(ubiq,cov_corrected)$p.value,
            c = cor.test(ubiq,cov_corrected)$estimate)


#####
# Methanogenesis
#####


# mrcABC, cov corrected
abs %>% 
  filter(n_genomes > 10) %>% 
  filter(!grepl("c__",tip)) %>% 
  filter(gene %in% c('mcrB','mcrA','mcrC')) %>% 
  group_by(tip) %>% 
  summarise(c = min(cov_corrected),
            cc = min(cov_corrected),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(c)) %>% 
  unique() %>% 
  print(n = 20)



# mrcABC,tgr
# f2 corrected
abs %>% 
  filter(n_genomes > 10) %>% 
  filter(grepl("p__",tip)) %>% 
  filter(gene %in% c('mcrB','mcrA','mcrC')) %>% 
  group_by(tip) %>% 
  summarise(n = min(f2_corrected),
            c = min(cov_corrected),
            ub = min(samples/totsamples),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>%
  unique() %>% 
  print(n = 40)

abs %>% 
  filter(n_genomes > 10) %>% 
  filter(gene %in% c('mcrB','mcrA','mcrC')) %>% 
  filter(tip == 'd__Archaea;p__Halobacteriota;c__Bog-38') %>% 
  dplyr::select(gene,cov_corrected,ubiq)

abs %>% 
  filter(grepl('UBA11358',tip))

######
# Methanotrophy
#####

# other proteins in the cluster, min does not work.. low values in general
abs %>% 
  filter(n_genomes > 10) %>% 
  filter(general == 'Methanotrophy') %>%
  filter(!gene %in% c('mmoX','mmoC','mmoY',"mmoZ")) %>% 
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip,gene) %>% 
  summarise(n = min(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 100)



# mmoX
# f2 corrected
abs %>% 
  filter(n_genomes > 10) %>% 
  #  filter(grepl("Respiration",general)) %>%
  filter(grepl("p__",tip)) %>% 
  #  filter(f2 > 0.5) %>%
  filter(gene == 'mmoX') %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = sum(f2_corrected),
            c = min(cov_corrected),
            ub = min(samples/totsamples),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n))


# mmoX, y & Z 
# f2 corrected
abs %>% 
  filter(n_genomes > 10) %>% 
  #  filter(grepl("Respiration",general)) %>%
  filter(grepl("p__",tip)) %>% 
  #  filter(f2 > 0.5) %>%
  filter(gene %in% c('mmoX', 'mmoY', 'mmoZ')) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = min(f2_corrected),
            c = min(cov_corrected),
            ub = min(samples/totsamples),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  unique() %>% 
  arrange(desc(n))


abs %>% 
  filter(n_genomes > 10) %>% 
  filter(gene %in% c('mmoX', 'mmoY', 'mmoZ')) %>% 
  filter(tip == 'd__Bacteria;p__Actinomycetota;c__Acidimicrobiia;o__IMCC26256;f__VFJJ01') %>% 
  dplyr::select(gene,cov_corrected,ubiq)



# mmoC
abs %>% 
  #  filter(grepl("Respiration",general)) %>%
  filter(grepl("p__",tip)) %>% 
  #filter(f2 > 0.4) %>%
  filter(gene == 'mmoC') %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = sum(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n))

# max for all subunits 
# other proteins in the cluster
abs %>% 
  filter(general == 'Methanotrophy') %>%
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
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
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip,gene) %>% 
  summarise(n = f2,n1 = cov,
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 100)


# pmo-mmo genes in bacteria
# f2 corrected
abs %>% 
  filter(n_genomes > 10) %>% 
  filter(gene %in% c('pmo-amoA','pmo-amoB','pmo-amoC')) %>%
  filter(grepl("d__Bacteria",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = min(f2_corrected),
            proportion_GB_genomes = proportion_GB_genomes,
            c = min(cov_corrected),
            ub = min(samples/totsamples)) %>% 
  arrange(desc(n)) %>%
  unique() %>% 
  print(n = 10) 


abs %>% 
  filter(n_genomes > 10) %>% 
  filter(gene %in% c('pmo-amoA','pmo-amoB','pmo-amoC')) %>%
  filter(tip == 'd__Bacteria;p__Chloroflexota;c__UBA6077;o__UBA6077;f__CF-72') %>% 
  dplyr::select(gene,cov_corrected,ubiq)




###
# CO oxidation
###


# all
# f2 corrected
abs %>% 
  filter(n_genomes > 10) %>% 
  #filter(grepl("d__Bacteria",tip)) %>% 
  filter(gene %in% c('coxM','coxS')) %>%
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = min(f2_corrected),
            proportion_GB_genomes = proportion_GB_genomes,
            c = min(cov_corrected),
            ub = min(samples/totsamples)) %>% 
  arrange(desc(n)) %>% 
  unique() %>% 
  print(n = 10)

abs %>% 
  filter(n_genomes > 10) %>% 
  #filter(grepl("d__Bacteria",tip)) %>% 
  filter(gene %in% c('coxM','coxS')) %>%
  filter(tip == 'd__Bacteria;p__Actinomycetota;c__Acidimicrobiia;o__IMCC26256') %>% 
  dplyr::select(gene,cov_corrected,ubiq)


#####################
#############
# nitrogen metabolism
############
#####################

# nitrogen fixation,  nif
# f2 corrected
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("nifH","nifD","nifK")) %>%
  filter(general %in% c("Nitrogen fixation")) %>%
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = min(f2_corrected),
            c = min(cov_corrected),
            ub = min(samples/totsamples),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  unique() %>% 
  arrange(desc(n)) %>% 
  print(n = 30)



abs %>% 
  filter(n_genomes > 10) %>% 
  filter(gene %in% c("nifH","nifD","nifK")) %>%
  #filter(grepl("d__Bacteria",tip)) %>% 
  filter(tip == 'd__Bacteria;p__Myxococcota;c__Polyangia;o__Polyangiales;f__Polyangiaceae;g__PMG-095') %>% 
  dplyr::select(gene,cov_corrected,ubiq)

# nitrogen fixation,  nif, no cyano
# f2 corrected
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("nifH","nifD","nifK")) %>%
  filter(general %in% c("Nitrogen fixation")) %>%
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = min(f2_corrected),
            c = min(cov_corrected),
            ub = min(samples/totsamples),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  unique() %>% 
  arrange(desc(n)) %>% 
  print(n = 60) 



# anfG
# f2 corrected
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("anfG")) %>%
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip,gene) %>% 
  summarise(n = sum(f2_corrected),
            proportion_GB_genomes = proportion_GB_genomes,
            c = min(cov_corrected),
            ub = min(samples/totsamples)) %>% 
  arrange(desc(n)) %>% 
  print(n = 100)


# nodAC + nif genes
head(abs)
abs %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("nodA","nodC","nifH")) %>%
  filter(general %in% c("Nitrogen fixation")) %>%
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip,gene) %>% 
  summarise(n = sum(f2),cov = sum(cov),samples = samples,
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 100) 


# pmos
# f2 corrected 
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("pmo-amoA","pmo-amoB","pmo-amoC")) %>%
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = min(f2_corrected),
            c = min(cov_corrected),
            ub = min(samples/totsamples),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  unique() %>% 
  arrange(desc(n)) %>% 
  print(n = 20)


abs %>% 
  filter(n_genomes > 10) %>% 
  filter(gene %in% c("pmo-amoA","pmo-amoB","pmo-amoC")) %>%
  #filter(grepl("d__Bacteria",tip)) %>% 
  filter(tip == 'd__Archaea;p__Thermoproteota;c__Nitrososphaeria;o__Nitrososphaerales;f__Nitrososphaeraceae;g__TA-21') %>% 
  dplyr::select(gene,cov_corrected,ubiq)


# hao
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("hao")) %>%
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip,gene) %>% 
  summarise(n = sum(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 100)


# nrxA / nrxB // narGH
# f2 corrected
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("narG, narZ, nxrA","narH, narY, nxrB")) %>%
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = min(f2_corrected),
            c = min(cov_corrected),
            ub = min(samples/totsamples),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  unique() %>% 
  print(n = 20)


abs %>% 
  filter(n_genomes > 10) %>% 
  filter(gene %in% c("narG, narZ, nxrA","narH, narY, nxrB")) %>%
  #filter(grepl("d__Bacteria",tip)) %>% 
  filter(tip == 'd__Bacteria;p__Actinomycetota;c__Thermoleophilia;o__Gaiellales;f__Gaiellaceae;g__GMQP-bins7') %>% 
  dplyr::select(gene,cov_corrected,ubiq)


# napA / napB
# f2 corrected
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("napA","napB")) %>%
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = min(f2_corrected),
            c = min(cov_corrected),
            ub = min(samples/totsamples),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  unique() %>% 
  print(n = 20)


abs %>% 
  filter(n_genomes > 10) %>% 
  filter(gene %in% c("napA","napB")) %>%
  #filter(grepl("d__Bacteria",tip)) %>% 
  filter(tip == 'd__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__SG8-39;g__SCGC-AG-212-J23') %>% 
  dplyr::select(gene,cov_corrected,ubiq)


# nirA
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("nirA")) %>%
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip,gene) %>% 
  summarise(n = sum(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 100)


# nasA
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(ko %in% c("K00372")) %>%
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip,gene) %>% 
  summarise(n = sum(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 100)


# nasC
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("nasC")) %>%
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = sum(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 100)

# nrfAH
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("nrfA","nrfH")) %>%
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = sum(f2),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 100)


# nirK 
# f2 corrected
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("nirK")) %>%
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = sum(f2_corrected),
            c = min(cov_corrected),
            ub = min(samples/totsamples),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 20)

abs %>% 
  filter(n_genomes > 10) %>% 
  filter(gene %in% c("nirK")) %>%
  #filter(grepl("d__Bacteria",tip)) %>% 
  filter(tip == 'd__Bacteria;p__Chloroflexota;c__UBA6077') %>% 
  dplyr::select(gene,cov_corrected,ubiq)



# NirS, f2 corrected
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("nirS")) %>%
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = sum(f2_corrected),
            c = min(cov_corrected),
            ub = min(samples/totsamples),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 20)


# both NirS and nirK
# f2 corrected
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("nirS","nirK")) %>%
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = min(f2_corrected),
            c = min(cov_corrected),
            ub = min(samples/totsamples),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  
  arrange(desc(n)) %>% 
  unique() %>% 
  print(n = 20)

#norBC, f2 corrected
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("norB","norC")) %>%
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = min(f2_corrected),
            proportion_GB_genomes = proportion_GB_genomes,
            c = min(cov_corrected),
            ub = min(samples/totsamples)) %>% 
  arrange(desc(n)) %>% 
  unique() %>% 
  print(n = 20)

# nosZ, f2 corrected
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("NosZ")) %>%
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = mean(f2_corrected),
            proportion_GB_genomes = proportion_GB_genomes,
            c = min(cov_corrected),
            ub = min(samples/totsamples)) %>% 
  arrange(desc(n)) %>% 
  print(n = 20)




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

write.table(export, "../Figures/summary figures/table_fig.tab",sep = '\t',row.names = F)

# cfix, k18603/4:TA-21, (abfD: VAZQ01?) , aclAB: Palsa-1315, rbcLS: CADCVQ01
# resp: UBA2999 
# coxMS: UBA2999
# mrcABC: c__Bog-38

# nifHDK: d__Bacteria;p__Nitrospirota;c__UBA9217, only with decent cov values, see below 
# nirk: GWC2-71-9
# nirS: UBA6960
# pmos: g__TA-21 (arch)
# pmos: d__Bacteria;p__Chloroflexota;c__UBA6077;o__UBA6077;f__CF-72 (bact, see below)
# nosZ: SCN-69-37

# nif, not well represented before 
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("nifH","nifD","nifK")) %>%
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
  #slice_head(n = 20) %>%
  arrange(desc(TGR)) %>% 
  print(n = 200)



# bacterial pmos
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("pmo-amoA","pmo-amoB","pmo-amoC")) %>%
  filter(grepl("d__Bacteria",tip)) %>%
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
  #slice_head(n = 20) %>%
  arrange(desc(TGR)) %>% 
  print(n = 200)




#################
##########
# Sulfur
##########
#################

# sat, met3
# sulfate adenylyltransferase 
# f2 corrected
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>%
  filter(gene %in% c("sat, met3")) %>%
  filter(grepl("p__",tip)) %>%
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = mean(f2_corrected),
            proportion_GB_genomes = proportion_GB_genomes,
            c = min(cov_corrected),
            ub = min(samples/totsamples)) %>% 
  arrange(desc(n)) %>% 
  print(n = 30)

# cysMK, corrected
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("cysM","cysK")) %>%
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = min(f2_corrected),
            c = min(cov_corrected),
            proportion_GB_genomes = proportion_GB_genomes,
            ub = min(samples/totsamples)) %>% 
  arrange(desc(n)) %>%
  unique() %>% 
  print(n = 20) 


# dsrAB, f2 corrected
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("dsrA","dsrB")) %>%
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = min(f2_corrected),
            c = min(cov_corrected),
            proportion_GB_genomes = proportion_GB_genomes,
            ub = min(samples/totsamples)) %>% 
  unique() %>% 
  arrange(desc(n)) %>% 
  print(n = 20)

# soxB, f2 corrected
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("soxB")) %>%
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = min(f2_corrected),
            c = min(cov_corrected),
            proportion_GB_genomes = proportion_GB_genomes,
            ub = min(samples/totsamples)) %>% 
  arrange(desc(n)) %>% 
  print(n = 10)


#################
##########
# Phosphorous
##########
#################

# qcd
# f2 corrected
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("qcd")) %>%
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = min(f2_corrected),
            c = min(cov_corrected),
            proportion_GB_genomes = proportion_GB_genomes,
            ub = min(samples/totsamples)) %>% 
  arrange(desc(n)) %>% 
  print(n = 20)

# phytases, pho, f2 corrected 
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("phytase","phoD","phoA")) %>%
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = min(f2_corrected),
            c = min(cov_corrected),
            proportion_GB_genomes = proportion_GB_genomes,
            ub = min(samples/totsamples)) %>% 
  unique() %>% 
  arrange(desc(n)) %>% 
  print(n = 20)

head(abs)
8533 / 9012

#################
##########
# Taxa involved in many / few routes
##########
#################

head(abs)
table(abs$general)

# number of times in the top 100 taxa for any gene (e.g. generalists)
# f2 corrected 
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(grepl("p__",tip)) %>% 
  filter(gene %in% c("pmo-amoA","rbcL","mcrA","NosZ","nifH","coxA",
                     "mmoX","narG, narZ, nxrA","napA","norB",'K01179',
                     'phoA',"qcd","soxB","dsrA")) %>% 
  dplyr::select(f2_corrected,gene,tip,proportion_GB_genomes) %>%
  unique() %>% 
  group_by(gene) %>% 
  arrange(desc(f2_corrected)) %>%
  slice_head(n = 100) %>% 
  group_by(tip) %>% 
  summarise(n = n(),
            concatenated_text = paste(gene, collapse = ", "),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  unique() %>% 
  arrange(desc(n)) %>% 
  print(n = 100)

# number of times in the top 100 taxa for only one gene (e.g. specialists)
# f2 corrected
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(grepl("p__",tip)) %>% 
  filter(gene %in% c("pmo-amoA","rbcL","mcrA","NosZ","nifH","coxA",
                     "mmoX","narG, narZ, nxrA","napA","norB",'K01179',
                     'phoA',"qcd","soxB","dsrA")) %>%   dplyr::select(f2_corrected,gene,tip) %>%
  unique() %>% 
  group_by(gene) %>% 
  arrange(desc(f2_corrected)) %>%
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
  filter(gene %in% c("pmo-amoA","rbcL","mcrA","NosZ","nifH","coxA","mmoX","narG, narZ, nxrA","napA","NorB",
                     'K01179',
                     'phoA',"qcd","soxB","dsrA")) %>% 
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
# f2 corrected
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
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(grepl("p__",tip)) %>% 
  filter(gene %in% c("pmo-amoA","rbcL","mcrA","NosZ","nifH","coxA","mmoX","narG, narZ, nxrA","napA","NorB")) %>% 
  filter(f2 > 0.3) %>% 
  dplyr::select(samples,gene,tip) %>%
  unique() %>% 
  group_by(gene) %>% 
  summarise(n = length(unique(tip))) %>% 
  arrange(desc(n))



#########
# f2 values total for each gene
#########

abs %>% 
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
  group_by(gene) %>% 
  summarise(TGR = sum(f2),
            lineage_cov = sum(cov)) %>%
  arrange(desc(TGR)) %>%
  print(n = 200)



###########################################################
###########################################################

#####
# GC responses
####

setwd("~/analysis/Berlin/soil 2019 metagenomes/tax_profile/singlem/")

metadata = read.table("../../metadata.tab",header = T,sep = '\t')
metadata$Tube.number
head(metadata)
metadata$sample = paste(metadata$Tube.number,sep  = '')


####
# load data
####

# load collapsed per tax data for plotting profiles
data = read.table("out_rare.tsv",header = F, sep = '\t')
head(data)
names(data) = c("sample","coverage","tax")
data$sample = tstrsplit(data$sample, "_")[[2]]
data$sample = str_replace_all(data$sample,"b",'')
head(data)


# add 0s to samples where lineages are not present
all_combinations <- expand.grid(
  sample = unique(data$sample),
  tax = unique(data$tax)
)

# Merge with the original dataframe to find missing combinations
data <- all_combinations %>%
  left_join(data, by = c("sample", "tax")) %>%
  replace_na(list(coverage = 0))


# join metadata
head(metadata)
data = data %>%
  left_join(metadata,by = "sample")


data$remark = factor(data$remark,levels = c("control","antibiotics", "copper", "drought",
                                            "microplastic", "Ndep", "salinity", "temp",
                                            "fungicide", "glyphosate","insecticide","Level 8"
))



data$tax = str_replace_all(data$tax,"Root; d","d")
data$tax = str_replace_all(data$tax,"; ",";")


d = data %>% 
  filter(tax %in% c("d__Bacteria;p__Pseudomonadota",
                    "d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;o__Acetobacterales",
                    "d__Bacteria;p__Cyanobacteriota",
                    "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__Usitatibacteraceae",
                    "d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;o__Rhizobiales;f__Xanthobacteraceae;g__Bradyrhizobium",
                    "d__Archaea;p__Thermoproteota;c__Nitrososphaeria;o__Nitrososphaerales",
                    "d__Archaea;p__Thermoproteota;c__Nitrososphaeria;o__Nitrososphaerales;f__Nitrososphaeraceae;g__Nitrososphaera",
                    "d__Archaea;p__Thermoproteota;c__Nitrososphaeria;o__Nitrososphaerales;f__Nitrososphaeraceae;g__TA-21",
                    "d__Archaea;p__Thermoproteota;c__Nitrososphaeria;o__Nitrososphaerales;f__Nitrososphaeraceae;g__Nitrosocosmicus",
                    "d__Archaea;p__Thermoproteota;c__Nitrososphaeria;o__Nitrososphaerales;f__Nitrososphaeraceae;g__TH1177",
                    "d__Bacteria;p__Nitrospirota;c__Nitrospiria;o__Nitrospirales;f__Nitrospiraceae;g__Palsa-1315",
                    "d__Bacteria;p__Nitrospirota", 
                    "d__Bacteria;p__Nitrospirota;c__Nitrospiria",
                    "d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;o__Reyranellales;f__Reyranellaceae;g__Reyranella",
                    "d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;o__Rhizobiales;f__Xanthobacteraceae;g__Z2-YC6860",
                    "d__Archaea;p__Halobacteriota", 
                    "d__Archaea;p__Methanobacteriota",
                    "d__Archaea;p__Thermoproteota;c__Nitrososphaeria;o__Nitrososphaerales;f__Nitrososphaeraceae;g__TA-21",
                    "d__Archaea;p__Thermoproteota;c__Nitrososphaeria;o__Nitrososphaerales;f__Nitrososphaeraceae;g__JAFAQB01",
                    "d__Archaea;p__Thermoproteota;c__Nitrososphaeria;o__Nitrososphaerales",
                    "d__Bacteria;p__Desulfobacterota_B;c__Binatia;o__Bin18",
                    "d__Bacteria;p__Chloroflexota;c__UBA6077;o__UBA6077;f__CF-72"))


head(d)
perform_wilcox <- function(df, value_col, factor_col, ref_level) {
  # List of levels excluding the reference level
  other_levels <- setdiff(unique(df[[factor_col]]), ref_level)
  
  # Perform Wilcoxon test for each comparison
  results <- map_dfr(other_levels, function(level) {
    subset_df <- df %>% filter(df[[factor_col]] %in% c(ref_level, level))
    test_result <- wilcox.test(subset_df[[value_col]] ~ subset_df[[factor_col]])
    
    tibble(
      reference_level = ref_level,
      compared_level = level,
      p_value = test_result$p.value
    )
  })
  
  return(results)
}


results <- d %>%
  group_by(tax) %>%
  group_modify(~ perform_wilcox(.x, value_col = "coverage", factor_col = "remark", ref_level = "control")) %>% 
  filter(p_value < 0.05) 

gc_effect_plot = ggplot(d %>% 
                          filter(tax %in% results$tax),aes(x  = remark, y = coverage))+
  geom_boxplot()+
  facet_wrap(~tax,scales = "free_x")+
  stat_compare_means(aes(group = as.factor(remark)),label = "p.signif", method = "wilcox.test",
                     ref.group = "control",hide.ns = TRUE, tip.length = 0)+
  coord_flip()+
  theme_minimal()+
  #  geom_hline(data = mm,aes(yintercept= m), linetype="dashed", 
  #             color = "grey", size=1)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")+
  geom_jitter(aes(color = as.factor(Lv)),shape=16, position=position_jitter(0.2),alpha = 0.8)+
  ylab("Coverage")+
  xlab("Number of factors")+
  labs(fill = "Number of factors",color = "Number of factors")+
  theme(legend.position = "none")  # Remove legend title


gc_effect_plot


