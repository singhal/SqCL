import re
import subprocess
import pandas as pd
import os
import numpy as np
import argparse
import gzip

parser = argparse.ArgumentParser(description="Calculate het and Fst.")
parser.add_argument('--cl', help="Cluster for which to run.")
args = parser.parse_args()
cl = args.cl

mt_file = '/Volumes/heloderma4/sonal/skink_mtDNA/cytb_alignment_30July15.fixed.aln.fa'
c_file = '/Volumes/heloderma4/sonal/eco_IBD_oz/data/clustering/Ctenotus_Lerista.clustering.revised.csv'
vcf_dir = '/Volumes/heloderma4/sonal/eco_IBD_oz/snp_calling/variants/'
out_dir = '/Volumes/heloderma4/sonal/eco_IBD_oz/data/diversity/'

def get_mt(mt_file):
	mt = {}
	id = ''
	f = open(mt_file, 'r')
	for l in f:
		if re.search('>', l):
			id = re.search('>(\S+)', l).group(1)
			mt[id] = ''
		else:
			mt[id] += l.rstrip()
	f.close()
	return mt


def get_clusters(c_file):
        d = pd.read_csv(c_file)
        d = d[d.GMYC_RAxML2.notnull()]
        d = d.groupby('GMYC_RAxML2')

        clusters = dict([(name, sorted(group['sample'].tolist())) for name, group in d])

        return clusters


def get_mt_dist(mt, ind1, ind2):
	allowed = ['A', 'T', 'C', 'G', 'a', 't', 'c', 'g']

	seq1 = mt[ind1]
	seq2 = mt[ind2]

	diff = 0
	denom = 0

	for bp1, bp2 in zip(seq1, seq2):
		if bp1 in allowed and bp2 in allowed:
			denom += 1
			if bp1.upper() != bp2.upper():
				diff += 1
	diff = diff / float(denom)
	return diff


def get_diversity(cl, inds, vcf_dir, out_file, mt):
	file = '%s%s.final.vcf.gz' % (vcf_dir, cl)

	pi = {'pi_sum': 0, 'sites': 0}
	het = {'het_sum': 0, 'sites': 0}

	allowed = ['0/0', '0/1', '1/1']

	f = gzip.open(file, 'r')
	for l in f:
		if not re.search('#', l) and not re.search('INDEL', l):
			d = re.split('\s+', l.rstrip())
			# don't mess with multiallelics
			if len(re.split(',', d[4])) == 1:
				genos = [re.search('^(\S\/\S)', x).group(1) for x in d[9:]]
				genos = [x for x in genos if x in allowed]

				# only do it if there is at least one non-missing site
				if len(genos) > 0:
					# calculate proportion hets
					het_prop = genos.count('0/1') / float(len(genos))
				 	het['het_sum'] += het_prop
					het['sites'] += 1

					alleles = []
					for geno in genos:
						alleles += re.split('/', geno)
					alleles = dict([(x, alleles.count(x)) for x in set(alleles)])
					
					if len(alleles) > 1:
						# https://binhe.org/2011/12/29/calculate-nucleotide-diversity-per-base-pair-summation-method/
						# total alleles
						n = float(np.sum(alleles.values()))
						# minor count
						j = float(np.min(alleles.values()))
						pi_prop = (2 * j * (n - j)) / (n * (n - 1))
					else:
						pi_prop = 0
					pi['pi_sum'] += pi_prop
					pi['sites'] += 1

	f.close()

	mt_diffs = []
	for ix, ind1 in enumerate(inds):
		for ind2 in inds[(ix+1):]:
			if ind1 in mt and ind2 in mt:
				mt_diffs.append(get_mt_dist(mt, ind1, ind2))				
	mt_pi = np.mean(mt_diffs)
	
	
	o = open(out_file, 'a')
	pi_val = pi['pi_sum'] / float(pi['sites'])
	het_val = het['het_sum'] / float(het['sites']) 
	o.write('%s,%s,%.6f,%s,%.6f,%s,%.6f\n' % (cl, len(inds), pi_val, pi['sites'], het_val, het['sites'], mt_pi))
	o.close()


def initialize_file(out_file):
	if not os.path.isfile(out_file):
		o = open(out_file, 'w')
		o.write('cluster,nInds,pi,pi_sites,het,het_sites,mt_pi\n')
		o.close()

out_file = '%sspecies_diversity.csv' % out_dir
initialize_file(out_file)
mt = get_mt(mt_file)
clusters = get_clusters(c_file)
get_diversity(cl, clusters[cl], vcf_dir, out_file, mt)
