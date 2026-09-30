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


abs = read.csv('~/analysis/Berlin/sandpiper/Figures - rev1/Supp/ko_GFC.pruned.csv')


#######
# Carbon fixation
######

# rbcl + rbcs
pruned %>% 
  filter(grepl("p__", tip)) %>% 
  filter(gene %in% c('rbcL', 'rbcS')) %>%
  group_by(tip, general) %>% 
  slice_min(GFC, n = 1, with_ties = FALSE) %>% 
  transmute(
    GFC = GFC,               
    ee_GFC = ee_GFC,       
    c = cov,               
    ee_cov = ee_cov,
    u = samples / totsamples,
    ee_ubiq = ee_ubiq,
    p_unc_gtdb = p_unc_gtdb
  ) %>% 
  ungroup() %>% 
  unique() %>% 
  arrange(desc(GFC)) %>% 
  print(n = 20)


# acetyl-CoA/propionyl-CoAcarboxylase
pruned %>% 
  filter(gene %in% c('K18604','K18603')) %>% 
  group_by(tip, general) %>% 
  slice_min(GFC, n = 1, with_ties = FALSE) %>% 
  transmute(
    GFC = GFC,            
    ee_GFC = ee_GFC,      
    c = cov,               
    ee_cov = ee_cov,
    u = samples / totsamples,
    ee_ubiq = ee_ubiq,
    p_unc_gtdb = p_unc_gtdb
  ) %>% 
  ungroup() %>% 
  unique() %>% 
  arrange(desc(GFC)) %>% 
  print(n = 20)

# abfD
pruned %>% 
  filter(gene == 'abfD') %>% 
  group_by(tip, general) %>% 
  slice_min(GFC, n = 1, with_ties = FALSE) %>% 
  transmute(
    GFC = GFC,             
    ee_GFC = ee_GFC,      
    c = cov,             
    ee_cov = ee_cov,
    u = samples / totsamples,
    ee_ubiq = ee_ubiq,
    p_unc_gtdb = p_unc_gtdb
  ) %>% 
  ungroup() %>% 
  unique() %>% 
  arrange(desc(GFC)) %>% 
  print(n = 20)



# AclA / AclB, 
pruned %>% 
  filter(gene %in% c('aclA','aclB')) %>% 
  group_by(tip, general) %>% 
  slice_min(GFC, n = 1, with_ties = FALSE) %>% 
  transmute(
    GFC = GFC,               
    ee_GFC = ee_GFC,     
    c = cov,               
    ee_cov = ee_cov,
    u = samples / totsamples,
    ee_ubiq = ee_ubiq,
    p_unc_gtdb = p_unc_gtdb
  ) %>% 
  ungroup() %>% 
  unique() %>% 
  arrange(desc(GFC)) %>% 
  print(n = 20)


#####
# Respiration
#####

# coxABC
pruned %>% 
  filter(grepl("Respiration",general)) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip, general) %>% 
  slice_min(GFC, n = 1, with_ties = FALSE) %>% 
  transmute(
    GFC = GFC,               
    ee_GFC = ee_GFC,       
    c = cov,              
    ee_cov = ee_cov,
    u = samples / totsamples,
    ee_ubiq = ee_ubiq,
    p_unc_gtdb = p_unc_gtdb
  ) %>% 
  ungroup() %>% 
  unique() %>% 
  arrange(desc(GFC)) %>% 
  print(n = 100)




#####
# Methanogenesis
#####


# mrcABC
pruned %>% 
  filter(grepl("p__",tip)) %>% 
  filter(gene %in% c('mcrB','mcrA','mcrC')) %>% 
  group_by(tip, general) %>% 
  slice_min(GFC, n = 1, with_ties = FALSE) %>% 
  transmute(
    GFC = GFC,               
    ee_GFC = ee_GFC,       
    c = cov,               
    ee_cov = ee_cov,
    u = samples / totsamples,
    ee_ubiq = ee_ubiq,
    p_unc_gtdb = p_unc_gtdb
  ) %>% 
  ungroup() %>% 
  unique() %>% 
  arrange(desc(GFC)) %>% 
  print(n = 20)


######
# Methanotrophy
#####

# mmoX, y & Z 
pruned %>% 
  filter(grepl("p__",tip)) %>% 
  filter(gene %in% c('mmoX', 'mmoY', 'mmoZ')) %>% 
  group_by(tip, general) %>% 
  slice_min(GFC, n = 1, with_ties = FALSE) %>% 
  transmute(
    GFC = GFC,              
    ee_GFC = ee_GFC,       
    c = cov,               
    ee_cov = ee_cov,
    u = samples / totsamples,
    ee_ubiq = ee_ubiq,
    p_unc_gtdb = p_unc_gtdb
  ) %>% 
  ungroup() %>% 
  unique() %>% 
  arrange(desc(GFC)) %>% 
  print(n = 20)


# pmogenes in bacteria
pruned %>% 
  filter(gene %in% c('pmo-amoA','pmo-amoB','pmo-amoC')) %>%
  filter(grepl("d__Bacteria",tip)) %>% 
  group_by(tip, general) %>% 
  slice_min(GFC, n = 1, with_ties = FALSE) %>% 
  transmute(
    GFC = GFC,               
    ee_GFC = ee_GFC,      
    c = cov,              
    ee_cov = ee_cov,
    u = samples / totsamples,
    ee_ubiq = ee_ubiq,
    p_unc_gtdb = p_unc_gtdb
  ) %>% 
  ungroup() %>% 
  unique() %>% 
  arrange(desc(GFC)) %>% 
  print(n = 20)


###
# CO oxidation
###

# coxMLS
pruned %>% 
  filter(gene %in% c('coxM','coxS','coxL')) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip, general) %>% 
  slice_min(GFC, n = 1, with_ties = FALSE) %>% 
  transmute(
    GFC = GFC,             
    ee_GFC = ee_GFC,    
    c = cov,             
    ee_cov = ee_cov,
    u = samples / totsamples,
    ee_ubiq = ee_ubiq,
    p_unc_gtdb = p_unc_gtdb
  ) %>% 
  ungroup() %>% 
  unique() %>% 
  arrange(desc(GFC)) %>% 
  print(n = 20)


#####################
#############
# nitrogen metabolism
############
#####################

# nitrogen fixation
pruned %>% 
  filter(gene %in% c("nifH","nifD","nifK")) %>%
  filter(general %in% c("Nitrogen fixation")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip, general) %>% 
  slice_min(GFC, n = 1, with_ties = FALSE) %>% 
  transmute(
    GFC = GFC,               
    ee_GFC = ee_GFC,       
    c = cov,               
    ee_cov = ee_cov,
    u = samples / totsamples,
    ee_ubiq = ee_ubiq,
    p_unc_gtdb = p_unc_gtdb
  ) %>% 
  ungroup() %>% 
  unique() %>% 
  arrange(desc(GFC)) %>% 
  print(n = 50)


# pmos
pruned %>% 
  filter(gene %in% c("pmo-amoA","pmo-amoB","pmo-amoC")) %>%
  filter(grepl("d__Archaea",tip)) %>% 
  group_by(tip, general) %>% 
  slice_min(GFC, n = 1, with_ties = FALSE) %>% 
  transmute(
    GFC = GFC,           
    ee_GFC = ee_GFC,     
    c = cov,           
    ee_cov = ee_cov,
    u = samples / totsamples,
    ee_ubiq = ee_ubiq,
    p_unc_gtdb = p_unc_gtdb
  ) %>% 
  ungroup() %>% 
  unique() %>% 
  arrange(desc(GFC)) %>% 
  print(n = 20)



# narGH
pruned %>% 
  filter(gene %in% c("narG, narZ, nxrA","narH, narY, nxrB")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip, general) %>% 
  slice_min(GFC, n = 1, with_ties = FALSE) %>% 
  transmute(
    GFC = GFC,             
    ee_GFC = ee_GFC,      
    c = cov,               
    ee_cov = ee_cov,
    u = samples / totsamples,
    ee_ubiq = ee_ubiq,
    p_unc_gtdb = p_unc_gtdb
  ) %>% 
  ungroup() %>% 
  unique() %>% 
  arrange(desc(GFC)) %>% 
  print(n = 20)


# napA / napB
pruned %>% 
  filter(gene %in% c("napA","napB")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip, general) %>% 
  slice_min(GFC, n = 1, with_ties = FALSE) %>% 
  transmute(
    GFC = GFC,             
    ee_GFC = ee_GFC,       
    c = cov,              
    ee_cov = ee_cov,
    u = samples / totsamples,
    ee_ubiq = ee_ubiq,
    p_unc_gtdb = p_unc_gtdb
  ) %>% 
  ungroup() %>% 
  unique() %>% 
  arrange(desc(GFC)) %>% 
  print(n = 20)

# nirK 
pruned %>% 
  filter(gene %in% c("nirK")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip, general) %>% 
  slice_min(GFC, n = 1, with_ties = FALSE) %>% 
  transmute(
    GFC = GFC,              
    ee_GFC = ee_GFC,      
    c = cov,            
    ee_cov = ee_cov,
    u = samples / totsamples,
    ee_ubiq = ee_ubiq,
    p_unc_gtdb = p_unc_gtdb
  ) %>% 
  ungroup() %>% 
  unique() %>% 
  arrange(desc(GFC)) %>% 
  print(n = 20)


# NirS
pruned %>% 
  filter(gene %in% c("nirS")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip, general) %>% 
  slice_min(GFC, n = 1, with_ties = FALSE) %>% 
  transmute(
    GFC = GFC,              
    ee_GFC = ee_GFC,      
    c = cov,               
    ee_cov = ee_cov,
    u = samples / totsamples,
    ee_ubiq = ee_ubiq,
    p_unc_gtdb = p_unc_gtdb
  ) %>% 
  ungroup() %>% 
  unique() %>% 
  arrange(desc(GFC)) %>% 
  print(n = 20)

#norBC
pruned %>% 
  filter(gene %in% c("norB","norC")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip, general) %>% 
  slice_min(GFC, n = 1, with_ties = FALSE) %>% 
  transmute(
    GFC = GFC,              
    ee_GFC = ee_GFC,    
    c = cov,               
    ee_cov = ee_cov,
    u = samples / totsamples,
    ee_ubiq = ee_ubiq,
    p_unc_gtdb = p_unc_gtdb
  ) %>% 
  ungroup() %>% 
  unique() %>% 
  arrange(desc(GFC)) %>% 
  print(n = 20)


# nosZ
pruned %>% 
  filter(gene %in% c("NosZ")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip, general) %>% 
  slice_min(GFC, n = 1, with_ties = FALSE) %>% 
  transmute(
    GFC = GFC,               
    ee_GFC = ee_GFC,      
    c = cov,             
    ee_cov = ee_cov,
    u = samples / totsamples,
    ee_ubiq = ee_ubiq,
    p_unc_gtdb = p_unc_gtdb
  ) %>% 
  ungroup() %>% 
  unique() %>% 
  arrange(desc(GFC)) %>% 
  head(n = 20)


#################
##########
# Sulfur
##########
#################

# sat, met3
pruned %>% 
  filter(gene %in% c("sat, met3")) %>%
  filter(grepl("p__",tip)) %>%
  group_by(tip, general) %>% 
  slice_min(GFC, n = 1, with_ties = FALSE) %>% 
  transmute(
    GFC = GFC,          
    ee_GFC = ee_GFC,   
    c = cov,               
    ee_cov = ee_cov,
    u = samples / totsamples,
    ee_ubiq = ee_ubiq,
    p_unc_gtdb = p_unc_gtdb
  ) %>% 
  ungroup() %>% 
  unique() %>% 
  arrange(desc(GFC)) %>% 
  print(n = 20)

# cysMK
pruned %>% 
  filter(gene %in% c("cysM","cysK")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip, general) %>% 
  slice_min(GFC, n = 1, with_ties = FALSE) %>% 
  transmute(
    GFC = GFC,               
    ee_GFC = ee_GFC,       
    c = cov,               
    ee_cov = ee_cov,
    u = samples / totsamples,
    ee_ubiq = ee_ubiq,
    p_unc_gtdb = p_unc_gtdb
  ) %>% 
  ungroup() %>% 
  unique() %>% 
  arrange(desc(GFC)) %>% 
  print(n = 20)


# dsrAB
pruned %>% 
  filter(gene %in% c("dsrA","dsrB")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip, general) %>% 
  slice_min(GFC, n = 1, with_ties = FALSE) %>% 
  transmute(
    GFC = GFC,            
    ee_GFC = ee_GFC,      
    c = cov,               
    ee_cov = ee_cov,
    u = samples / totsamples,
    ee_ubiq = ee_ubiq,
    p_unc_gtdb = p_unc_gtdb
  ) %>% 
  ungroup() %>% 
  unique() %>% 
  arrange(desc(GFC)) %>% 
  print(n = 20)

# soxB
pruned %>% 
  filter(gene %in% c("soxB")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip, general) %>% 
  slice_min(GFC, n = 1, with_ties = FALSE) %>% 
  transmute(
    GFC = GFC,               
    ee_GFC = ee_GFC,     
    c = cov,              
    ee_cov = ee_cov,
    u = samples / totsamples,
    ee_ubiq = ee_ubiq,
    p_unc_gtdb = p_unc_gtdb
  ) %>% 
  ungroup() %>% 
  unique() %>% 
  arrange(desc(GFC)) %>% 
  print(n = 20)


#################
##########
# Phosphorous
##########
#################

# qcd
pruned %>% 
  filter(gene %in% c("qcd")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip, general) %>% 
  slice_min(GFC, n = 1, with_ties = FALSE) %>% 
  transmute(
    GFC = GFC,             
    ee_GFC = ee_GFC,       
    c = cov,              
    ee_cov = ee_cov,
    u = samples / totsamples,
    ee_ubiq = ee_ubiq,
    p_unc_gtdb = p_unc_gtdb
  ) %>% 
  ungroup() %>% 
  unique() %>% 
  arrange(desc(GFC)) %>% 
  print(n = 20)


# pho
pruned %>% 
  filter(gene %in% c("phoD","phoA")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip, general) %>% 
  slice_min(GFC, n = 1, with_ties = FALSE) %>% 
  transmute(
    GFC = GFC,             
    ee_GFC = ee_GFC,     
    c = cov,          
    ee_cov = ee_cov,
    u = samples / totsamples,
    ee_ubiq = ee_ubiq,
    p_unc_gtdb = p_unc_gtdb
  ) %>% 
  ungroup() %>% 
  unique() %>% 
  arrange(desc(GFC)) %>% 
  print(n = 20)


# "phytase",
pruned %>% 
  filter(gene %in% c("phytase")) %>%
  filter(grepl("p__",tip)) %>% 
  group_by(tip, general) %>% 
  slice_min(GFC, n = 1, with_ties = FALSE) %>% 
  transmute(
    GFC = GFC,             
    ee_GFC = ee_GFC,    
    c = cov,         
    ee_cov = ee_cov,
    u = samples / totsamples,
    ee_ubiq = ee_ubiq,
    p_unc_gtdb = p_unc_gtdb
  ) %>% 
  ungroup() %>% 
  unique() %>% 
  arrange(desc(GFC)) %>% 
  print(n = 20)



#################
################
# Number of  uncultivated genera among the most importnat for each 
###############
#################

head(pruned)

as.data.frame(pruned) %>% 
  ungroup() %>% 
  filter(gene %in% c("nifH","nifD","nifK","anfG",
                   "pmo-amoA","pmo-amoB","pmo-amoC",
                   "narG, narZ, nxrA","narH, narY, nxrB","napA","napB",
                   "norB","norC","NosZ","nirS","nirK",
                   "rbcL","rbcS","prkB",
                   "K18603","K18604","abfD","aclA","aclB",
                   "mcrA","mcrB","mcrC","mcrD","fwdA","fwdB","fwdC",
                   "mmoX","mmoC","mmoY","mmoZ","pmo-amoA","pmo-amoB","pmo-amoC",
                   "coxA","coxB","coxC",
                   "coxS","coxM","coxS",'sat, met3','soxB','qcd','cysM','cysK','dsrA','dsrB',
                   'phoA','phoD',"phytase")) %>%
  mutate(pahtway = case_when(gene %in% c('coxA','coxB','coxC')~'coxABC',
                             gene %in% c('nifH','nifD','nifK')~'nifHDK',
                             gene %in% c("mcrA","mcrB","mcrC")~'mcrABC',
                             gene %in% c("pmo-amoA","pmo-amoB","pmo-amoC")~'pmoABC',
                             gene %in% c("aclA","aclB")~'aclAB',
                             gene %in% c("coxS","coxM",'coxL')~'coxMSL',
                             gene %in% c("mmoX","mmoY","mmoZ")~'mmoXYZ',
                             gene %in% c("norB","norC")~'norBC',
                             gene %in% c("narG, narZ, nxrA","narH, narY, nxrB")~'narGH',
                             gene %in% c("napA","napB")~'napAB',
                             gene %in% c("rbcL","rbcS")~'rbcLS',
                             gene %in% c("K18603","K18604")~"K186034",
                             gene %in% c("fwdA","fwdB","fwdC")~'fwdABC',
                             gene %in% c('phoA','phoD',"phytase")~'phoAD_phy',
                             gene %in% c("nrfA","nrfB")~'nrfAB',
                             gene %in% c('nasB','nasC, nasA')~'nasABC',
                             gene %in% c('vnfD','vnfG','vnfK')~'vnfDGK',
                             gene %in% c('korA','korB')~'korAB',
                             gene %in% c('soxB')~'soxB',
                             gene %in% c('cysK','cysM')~'cysKM',
                             gene %in% c('dsrA','dsrB')~'dsrAB',
                             gene %in% c('dsrA','dsrB')~'dsrAB',
                             .default = gene
  )) %>% 
  filter(grepl("g__",tip) & !grepl('s__',tip)) %>% 
  group_by(tip, general,pahtway) %>% 
  slice_min(GFC, n = 1, with_ties = FALSE) %>% 
  transmute(
    GFC = GFC,              
    ee_GFC = ee_GFC,       
    c = cov,              
    ee_cov = ee_cov,
    u = samples / totsamples,
    ee_ubiq = ee_ubiq,
    p_unc_gtdb = p_unc_gtdb
  ) %>% 
  filter(GFC > 0) %>% 
  ungroup() %>%
  unique() %>% 
  group_by(pahtway) %>% 
  slice_max(GFC, n = 20, with_ties = FALSE) %>% 
  ungroup() %>% 
  filter(p_unc_gtdb == 1) %>% 
  group_by(pahtway) %>% 
  summarise(n = n()) %>% 
  print(n = 30)




#################
##########
# Taxa involved in many / few routes
##########
#################

head(abs)
table(abs$general)

# number of times in the top 100 taxa for any gene (e.g. generalists)
pruned %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(grepl("p__",tip)) %>% 
  filter(gene %in% c("pmo-amoA","rbcL","mcrA","NosZ","nifH","coxA",
                     "mmoX","narG, narZ, nxrA","napA","norB",'K01179',
                     'phoA',"qcd","soxB","dsrA")) %>% 
  dplyr::select(GFC,gene,tip,p_unc_gtdb) %>%
  unique() %>% 
  group_by(gene) %>% 
  arrange(desc(GFC)) %>%
  slice_head(n = 100) %>% 
  group_by(tip) %>% 
  summarise(n = n(),
            concatenated_text = paste(gene, collapse = ", "),
            p_unc_gtdb = p_unc_gtdb) %>% 
  unique() %>% 
  arrange(desc(n)) %>% 
  print(n = 100)

# number of times in the top 100 taxa for only one gene (e.g. specialists)
pruned %>% 
  #  filter(grepl("d__Bacteria",tip)) %>% 
  filter(n_genomes>10) %>% 
  filter(grepl("p__",tip)) %>% 
  filter(gene %in% c("pmo-amoA","rbcL","mcrA","NosZ","nifH","coxA",
                     "mmoX","narG, narZ, nxrA","napA","norB",'K01179',
                     'phoA',"qcd","soxB","dsrA")) %>%   dplyr::select(GFC,gene,tip) %>%
  unique() %>% 
  group_by(gene) %>% 
  arrange(desc(GFC)) %>%
  slice_head(n = 100) %>% 
  group_by(tip) %>% 
  summarise(n = n(),
            concatenated_text = paste(gene, collapse = ", ")) %>% 
  arrange(desc(n)) %>% 
  filter(n == 1) %>% 
  group_by(concatenated_text) %>% 
  summarise(n = n()) %>% 
  arrange(desc(n))



######
# corr genes // ubiq
######

da = abs %>% 
  group_by(gene,general) %>% 
  filter(gene %in% c("nifH","nifD","nifK","anfG",
                     "pmo-amoA","pmo-amoB","pmo-amoC",
                     "narG, narZ, nxrA","narH, narY, nxrB","napA","napB",
                     "norB","norC","NosZ","nirS","nirK",
                     "rbcL","rbcS","prkB",
                     "K18603","K18604","abfD","aclA","aclB",
                     "mcrA","mcrB","mcrC","mcrD","fwdA","fwdB","fwdC",
                     "pmo-amoA","pmo-amoB","pmo-amoC", #"mmoX","mmoC","mmoY","mmoZ",
                     "coxA","coxB","coxC",
                     "coxS","coxM","coxS",
                     'phoA',"qcd","soxB","dsrA",
                     "K01179",
                     "CELB",
                     "bcsZ",
                     "CBH1",'coxL',
                     "CBH2", 'cysM', 'cysK')) %>%
  filter(grepl("g__",tip) & !grepl("s__",tip)) %>% 
  summarise(p = cor.test(samples,cov,method = 'pearson')$p.value,
            c = cor.test(samples,cov,method = 'pearson')$estimate) %>% 
  arrange(-c) %>% 
  print(n = 100)


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
                            'phoA',"qcd","soxB","dsrA",'cysM','cysK')) %>%
         filter(grepl("g__",tip) & !grepl("s__",tip)))+
  geom_point(aes(x = 100*samples / totsamples,y = cov),
             alpha = 0.3)+
  facet_wrap(~gene,scales = 'free')+
  theme_classic()+
  ylab('Gene conservation')+
  xlab('Ubiquity (% soil samples)')+
  stat_cor(aes(x = samples / totsamples,y = cov),
           method = "pearson",
           label.x.npc = "left",
           label.y.npc = "top",
           size = 3.2,
           color = 'red',
           p.accuracy = 0.01, r.accuracy = 0.01
  ) 


