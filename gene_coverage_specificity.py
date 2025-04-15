import sys
from collections import defaultdict,Counter
import numpy as np

lin2num = Counter()
genome2t = {}
for line in open(sys.argv[1]):
    genome,tax = list(map(str.strip,line.split('\t')))
    genome2t[genome] = tax
    genome2t[genome.split('.')[0]] = tax
    t_comb = []
  
    # number of genomes per lin (for calculating coverage)
    for t in tax.split(';'):
        t_comb.append(t)
        lin2num[';'.join(t_comb)] += 1

# parse emapper annotations
lin2db2annot = defaultdict(lambda:defaultdict(lambda:defaultdict(lambda:[])))
db2annot2genomes = defaultdict(lambda:defaultdict(lambda:set()))
for line in open(sys.argv[2]):
    if line.startswith('#'):
        continue

    f = list(map(str.strip,line.split('\t')))

    genome = f[0].split('@')[0]
    if genome not in genome2t:
        continue


    all_annots = defaultdict(lambda:set())
    for i,db in enumerate(['kos','kpath','kmod']):
        field = i + 11
        annots = f[field]
        for annot in annots.split(','):
            all_annots[db].add(annot)
            db2annot2genomes[db][annot].add(genome)

    tax =  genome2t[genome]
    t_comb = []
    for t in tax.split(';'):
        t_comb.append(t)
        for db in all_annots:
            for annot in all_annots[db]:
                lin2db2annot[';'.join(t_comb)][db][annot].append(genome)


# get coverage / specificity
for lin in lin2db2annot:
    for db in lin2db2annot[lin]:
        for annot,genomes in lin2db2annot[lin][db].items():

            n_genomes = len(list(set(genomes)))
            n_copies = len(genomes)
            n_genomes_annot = len(list(db2annot2genomes[db][annot]))

            mean = n_copies / n_genomes
            lin_n = lin2num[lin]
            cov = n_genomes / lin_n
            sp = n_genomes / n_genomes_annot

            print ('\t'.join(list(map(str,[lin,db,annot,cov,sp,mean,lin_n]))))
