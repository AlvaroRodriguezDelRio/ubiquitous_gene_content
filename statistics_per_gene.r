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


###############
####
# load data for gene coverage / gfc index
####
##############

abs = read.csv('ko_TGR.v3.csv')
head(abs)

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
  summarise(p = cor.test(proportion_GB_genomes,conservation_corrected)$p.value,
            c = cor.test(proportion_GB_genomes,conservation_corrected)$estimate) %>% 
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
  summarise(p = wilcox.test(conservation_corrected~uncultivated)$p.value,
            median_uncult = median(conservation_corrected[uncultivated == "yes"], na.rm = TRUE),
            median_cult = median(conservation_corrected[uncultivated == "no"], na.rm = TRUE),
            mean_uncult = mean(conservation_corrected[uncultivated == "yes"], na.rm = TRUE),
            mean_cult = mean(conservation_corrected[uncultivated == "no"], na.rm = TRUE),
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
  summarise(p = wilcox.test(conservation_corrected~uncultivated)$p.value,
            median_uncult = median(conservation_corrected[uncultivated == "yes"], na.rm = TRUE),
            median_cult = median(conservation_corrected[uncultivated == "no"], na.rm = TRUE),
            mean_uncult = mean(conservation_corrected[uncultivated == "yes"], na.rm = TRUE),
            mean_cult = mean(conservation_corrected[uncultivated == "no"], na.rm = TRUE),
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
  geom_boxplot(aes(x = uncultivated, y = conservation_corrected))+
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
  summarise(p = cor.test(samples,conservation_corrected,method = 'spearman')$p.value,
            c = cor.test(samples,conservation_corrected,method = 'spearman')$estimate) %>% 
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
  summarise(p = cor.test(samples,conservation_corrected,method = 'spearman')$p.value,
            c = cor.test(samples,conservation_corrected,method = 'spearman')$estimate) %>% 
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
  geom_point(aes(x = 100*samples / totsamples,y = conservation_corrected),
             alpha = 0.3)+
  facet_wrap(~gene,scales = 'free')+
  theme_classic()+
  ylab('Gene conservation')+
  xlab('Ubiquity (% soil samples)')+
  stat_cor(aes(x = samples / totsamples,y = conservation_corrected),
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

head(abs)

# rbcLS and prkb, cons, phylum, conservation
abs %>% 
  filter(n_genomes > 10) %>% 
  filter(grepl("Carbon fixataion",general)) %>%
  filter(gene %in% c('rbcL','rbcS','prkB')) %>%
  filter(!grepl('c__',tip)) %>% 
  group_by(tip,gene) %>% 
  summarise(c = conservation_corrected,
            proportion_GB_genomes = proportion_GB_genomes,
            u = samples / totsamples) %>% 
  arrange(desc(c)) %>% 
  print(n = 20) 


# rbcl, all, GFC
abs %>% 
  filter(n_genomes>10) %>% 
  filter(gene == 'rbcL') %>%
  group_by(tip,gene) %>% 
  summarise(n = sum(GFC),
            proportion_GB_genomes = proportion_GB_genomes,
            u = samples / totsamples,
            c = conservation_corrected) %>% 
  arrange(desc(n)) %>% 
  print(n = 20)


# rbcl + rbcs, GFC
abs %>% 
  filter(n_genomes>10) %>% 
  filter(!grepl("Photosynthesis",general)) %>%
  filter(grepl("p__",tip)) %>% 
  filter(gene %in% c('rbcL','rbcS')) %>%
  group_by(tip,general) %>% 
  summarise(n = min(GFC),
            c = min(conservation_corrected),
            u = samples / totsamples,
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  unique() %>% 
  arrange(desc(n)) %>% 
  print(n = 20)


abs %>% 
  filter(gene %in% c('rbcL','rbcS')) %>%
  filter(tip == 'd__Bacteria;p__Actinomycetota;c__Thermoleophilia;o__Solirubrobacterales;f__Solirubrobacteraceae;g__CADCVQ01') %>% 
  dplyr::select(conservation_corrected,samples,totsamples, gene, ubiq)



########
# Carbon fixation, others
#######

# acetyl-CoA/propionyl-CoAcarboxylase
# GFC
abs %>% 
  filter(n_genomes>10) %>% 
  filter(!grepl("Photosynthesis",general)) %>%
  filter(gene %in% c('K18604','K18603')) %>% 
  group_by(tip) %>% 
  summarise(n = min(GFC),
            c = min(conservation_corrected),
            ub = min(samples/totsamples),
            proportion_GB_genomes = proportion_GB_genomes) %>% # sp for specificity 
  unique() %>% 
  arrange(desc(n)) %>% 
  print(n = 20)

abs %>% 
  filter(gene %in% c('K18604','K18603')) %>% 
  filter(tip == 'd__Archaea;p__Thermoproteota;c__Nitrososphaeria;o__Nitrososphaerales;f__Nitrososphaeraceae;g__TA-21') %>% 
  dplyr::select(conservation_corrected, gene, ubiq)



# abfD
abs %>% 
  filter(n_genomes>10) %>%
  filter(gene == 'abfD') %>% 
  group_by(tip) %>% 
  summarise(n = sum(GFC),
            c = min(conservation_corrected),
            ub = min(samples/totsamples),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n  = 20)



######
# Carbon degradation
######

# cellulases
# GFC
abs %>% 
  filter(n_genomes>10) %>%
  filter(grepl(c("K01179",
                 "CELB",
                 "bcsZ",
                 "CBH1",
                 "CBH2"),gene)) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip) %>% 
  summarise(n = min(GFC),
            proportion_GB_genomes = proportion_GB_genomes,
            c = min(conservation_corrected),
            ub = min(samples/totsamples)) %>% 
  arrange(desc(n)) %>% 
  print(n  = 30)


#####
# Respiration
#####

# top ranks
# GFC completeness
abs %>% 
  filter(n_genomes > 10) %>% 
  filter(grepl("Respiration",general)) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip) %>% 
  summarise(n = min(GFC),
            proportion_GB_genomes = proportion_GB_genomes,
            c = min(conservation_corrected),
            ub = min(samples/totsamples)
  ) %>% 
  unique() %>% 
  arrange(desc(n)) %>% 
  print(n = 100)


abs %>% 
  filter(n_genomes > 10) %>% 
  filter(grepl("Respiration",general)) %>%
  filter(tip == 'd__Bacteria;p__Desulfobacterota_B') %>% 
  dplyr::select(conservation_corrected, gene, ubiq)



# coxABC correlation with ubiq
# cons corrected
abs %>% 
  filter(n_genomes > 10) %>% 
  filter(grepl("Respiration",general)) %>%
  filter(grepl("g__",tip) & !grepl("s__",tip)) %>% 
  summarise(p = cor.test(ubiq,conservation_corrected)$p.value,
            c = cor.test(ubiq,conservation_corrected)$estimate)


#####
# Methanogenesis
#####


# mrcABC,tgr
# GFC
abs %>% 
  filter(n_genomes > 10) %>% 
  filter(grepl("p__",tip)) %>% 
  filter(gene %in% c('mcrB','mcrA','mcrC')) %>% 
  group_by(tip) %>% 
  summarise(n = min(GFC),
            c = min(conservation_corrected),
            ub = min(samples/totsamples),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>%
  unique() %>% 
  print(n = 40)

abs %>% 
  filter(n_genomes > 10) %>% 
  filter(gene %in% c('mcrB','mcrA','mcrC')) %>% 
  filter(tip == 'd__Archaea;p__Halobacteriota;c__Bog-38') %>% 
  dplyr::select(gene,conservation_corrected,ubiq)

######
# Methanotrophy
#####


# mmoX, y & Z 
# GFC
abs %>% 
  filter(n_genomes > 10) %>% 
  filter(grepl("p__",tip)) %>% 
  filter(gene %in% c('mmoX', 'mmoY', 'mmoZ')) %>% 
  group_by(tip) %>% 
  summarise(n = min(GFC),
            c = min(conservation_corrected),
            ub = min(samples/totsamples),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  unique() %>% 
  arrange(desc(n))


abs %>% 
  filter(n_genomes > 10) %>% 
  filter(gene %in% c('mmoX', 'mmoY', 'mmoZ')) %>% 
  filter(tip == 'd__Bacteria;p__Actinomycetota;c__Acidimicrobiia;o__IMCC26256;f__VFJJ01') %>% 
  dplyr::select(gene,conservation_corrected,ubiq)

# pmo-mmo genes in bacteria
# GFC
abs %>% 
  filter(n_genomes > 10) %>% 
  filter(gene %in% c('pmo-amoA','pmo-amoB','pmo-amoC')) %>%
  filter(grepl("d__Bacteria",tip)) %>% 
  group_by(tip) %>% 
  summarise(n = min(GFC),
            proportion_GB_genomes = proportion_GB_genomes,
            c = min(conservation_corrected),
            ub = min(samples/totsamples)) %>% 
  arrange(desc(n)) %>%
  unique() %>% 
  print(n = 10) 


abs %>% 
  filter(n_genomes > 10) %>% 
  filter(gene %in% c('pmo-amoA','pmo-amoB','pmo-amoC')) %>%
  filter(tip == 'd__Bacteria;p__Chloroflexota;c__UBA6077;o__UBA6077;f__CF-72') %>% 
  dplyr::select(gene,conservation_corrected,ubiq)




###
# CO oxidation
###


# all
# GFC
abs %>% 
  filter(n_genomes > 10) %>% 
  filter(gene %in% c('coxM','coxS')) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip) %>% 
  summarise(n = min(GFC),
            proportion_GB_genomes = proportion_GB_genomes,
            c = min(conservation_corrected),
            ub = min(samples/totsamples)) %>% 
  arrange(desc(n)) %>% 
  unique() %>% 
  print(n = 10)

abs %>% 
  filter(n_genomes > 10) %>% 
  filter(gene %in% c('coxM','coxS')) %>%
  filter(tip == 'd__Bacteria;p__Actinomycetota;c__Acidimicrobiia;o__IMCC26256') %>% 
  dplyr::select(gene,conservation_corrected,ubiq)


#####################
#############
# nitrogen metabolism
############
#####################

# nitrogen fixation,  nif
# GFC
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("nifH","nifD","nifK")) %>%
  filter(general %in% c("Nitrogen fixation")) %>%
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = min(GFC),
            c = min(conservation_corrected),
            ub = min(samples/totsamples),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  unique() %>% 
  arrange(desc(n)) %>% 
  print(n = 30)



abs %>% 
  filter(n_genomes > 10) %>% 
  filter(gene %in% c("nifH","nifD","nifK")) %>%
  filter(tip == 'd__Bacteria;p__Myxococcota;c__Polyangia;o__Polyangiales;f__Polyangiaceae;g__PMG-095') %>% 
  dplyr::select(gene,conservation_corrected,ubiq)

# nitrogen fixation,  nif, no cyano
# GFC
abs %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("nifH","nifD","nifK")) %>%
  filter(general %in% c("Nitrogen fixation")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip) %>% 
  summarise(n = min(GFC),
            c = min(conservation_corrected),
            ub = min(samples/totsamples),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  unique() %>% 
  arrange(desc(n)) %>% 
  print(n = 60) 



# anfG
# GFC
abs %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("anfG")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip,gene) %>% 
  summarise(n = sum(GFC),
            proportion_GB_genomes = proportion_GB_genomes,
            c = min(conservation_corrected),
            ub = min(samples/totsamples)) %>% 
  arrange(desc(n)) %>% 
  print(n = 100)


# pmos
# GFC 
abs %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("pmo-amoA","pmo-amoB","pmo-amoC")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip) %>% 
  summarise(n = min(GFC),
            c = min(conservation_corrected),
            ub = min(samples/totsamples),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  unique() %>% 
  arrange(desc(n)) %>% 
  print(n = 20)


abs %>% 
  filter(n_genomes > 10) %>% 
  filter(gene %in% c("pmo-amoA","pmo-amoB","pmo-amoC")) %>%
  filter(tip == 'd__Archaea;p__Thermoproteota;c__Nitrososphaeria;o__Nitrososphaerales;f__Nitrososphaeraceae;g__TA-21') %>% 
  dplyr::select(gene,conservation_corrected,ubiq)


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
# GFC
abs %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("narG, narZ, nxrA","narH, narY, nxrB")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip) %>% 
  summarise(n = min(GFC),
            c = min(conservation_corrected),
            ub = min(samples/totsamples),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  unique() %>% 
  print(n = 20)


abs %>% 
  filter(n_genomes > 10) %>% 
  filter(gene %in% c("narG, narZ, nxrA","narH, narY, nxrB")) %>%
  filter(tip == 'd__Bacteria;p__Actinomycetota;c__Thermoleophilia;o__Gaiellales;f__Gaiellaceae;g__GMQP-bins7') %>% 
  dplyr::select(gene,conservation_corrected,ubiq)


# napA / napB
# GFC
abs %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("napA","napB")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip) %>% 
  summarise(n = min(GFC),
            c = min(conservation_corrected),
            ub = min(samples/totsamples),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  unique() %>% 
  print(n = 20)


abs %>% 
  filter(n_genomes > 10) %>% 
  filter(gene %in% c("napA","napB")) %>%
  filter(tip == 'd__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__SG8-39;g__SCGC-AG-212-J23') %>% 
  dplyr::select(gene,conservation_corrected,ubiq)


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
# GFC
abs %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("nirK")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip) %>% 
  summarise(n = sum(GFC),
            c = min(conservation_corrected),
            ub = min(samples/totsamples),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 20)

abs %>% 
  filter(n_genomes > 10) %>% 
  filter(gene %in% c("nirK")) %>%
  filter(tip == 'd__Bacteria;p__Chloroflexota;c__UBA6077') %>% 
  dplyr::select(gene,conservation_corrected,ubiq)



# NirS, GFC
abs %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("nirS")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip) %>% 
  summarise(n = sum(GFC),
            c = min(conservation_corrected),
            ub = min(samples/totsamples),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  arrange(desc(n)) %>% 
  print(n = 20)


# both NirS and nirK
# GFC
abs %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("nirS","nirK")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip) %>% 
  summarise(n = min(GFC),
            c = min(conservation_corrected),
            ub = min(samples/totsamples),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  
  arrange(desc(n)) %>% 
  unique() %>% 
  print(n = 20)

#norBC, GFC
abs %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("norB","norC")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip) %>% 
  summarise(n = min(GFC),
            proportion_GB_genomes = proportion_GB_genomes,
            c = min(conservation_corrected),
            ub = min(samples/totsamples)) %>% 
  arrange(desc(n)) %>% 
  unique() %>% 
  print(n = 20)

# nosZ, GFC
abs %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("NosZ")) %>%
  filter(grepl("p__",tip)) %>% 
  #dplyr::select(samples,gene,general,tip,cov,sp,f1,f3,f2) %>% 
  group_by(tip) %>% 
  summarise(n = mean(GFC),
            proportion_GB_genomes = proportion_GB_genomes,
            c = min(conservation_corrected),
            ub = min(samples/totsamples)) %>% 
  arrange(desc(n)) %>% 
  print(n = 20)

#################
##########
# Sulfur
##########
#################

# sat, met3
# sulfate adenylyltransferase 
# GFC
abs %>% 
  filter(n_genomes>10) %>%
  filter(gene %in% c("sat, met3")) %>%
  filter(grepl("p__",tip)) %>%
  group_by(tip) %>% 
  summarise(n = mean(GFC),
            proportion_GB_genomes = proportion_GB_genomes,
            c = min(conservation_corrected),
            ub = min(samples/totsamples)) %>% 
  arrange(desc(n)) %>% 
  print(n = 30)

# cysMK, corrected
abs %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("cysM","cysK")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip) %>% 
  summarise(n = min(GFC),
            c = min(conservation_corrected),
            proportion_GB_genomes = proportion_GB_genomes,
            ub = min(samples/totsamples)) %>% 
  arrange(desc(n)) %>%
  unique() %>% 
  print(n = 20) 


# dsrAB, GFC
abs %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("dsrA","dsrB")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip) %>% 
  summarise(n = min(GFC),
            c = min(conservation_corrected),
            proportion_GB_genomes = proportion_GB_genomes,
            ub = min(samples/totsamples)) %>% 
  unique() %>% 
  arrange(desc(n)) %>% 
  print(n = 20)

# soxB, GFC
abs %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("soxB")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip) %>% 
  summarise(n = min(GFC),
            c = min(conservation_corrected),
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
# GFC
abs %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("qcd")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip) %>% 
  summarise(n = min(GFC),
            c = min(conservation_corrected),
            proportion_GB_genomes = proportion_GB_genomes,
            ub = min(samples/totsamples)) %>% 
  arrange(desc(n)) %>% 
  print(n = 20)

# phytases, pho, GFC 
abs %>% 
  filter(n_genomes>10) %>% 
  filter(gene %in% c("phytase","phoD","phoA")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip) %>% 
  summarise(n = min(GFC),
            c = min(conservation_corrected),
            proportion_GB_genomes = proportion_GB_genomes,
            ub = min(samples/totsamples)) %>% 
  unique() %>% 
  arrange(desc(n)) %>% 
  print(n = 20)

head(abs)
8533 / 9012

#################
##########
# Taxa involved in many 
##########
#################

head(abs)
table(abs$general)

# number of times in the top 100 taxa for any gene 
# GFC 
abs %>% 
  filter(n_genomes>10) %>% 
  filter(grepl("p__",tip)) %>% 
  filter(gene %in% c("pmo-amoA","rbcL","mcrA","NosZ","nifH","coxA",
                     "mmoX","narG, narZ, nxrA","napA","norB",'K01179',
                     'phoA',"qcd","soxB","dsrA")) %>% 
  dplyr::select(GFC,gene,tip,proportion_GB_genomes) %>%
  unique() %>% 
  group_by(gene) %>% 
  arrange(desc(GFC)) %>%
  slice_head(n = 100) %>% 
  group_by(tip) %>% 
  summarise(n = n(),
            concatenated_text = paste(gene, collapse = ", "),
            proportion_GB_genomes = proportion_GB_genomes) %>% 
  unique() %>% 
  arrange(desc(n)) %>% 
  print(n = 100)
