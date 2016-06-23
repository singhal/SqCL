import pickle
import re
import subprocess
import pandas as pd
import os
import numpy as np
import argparse
import gzip
from geopy.distance import vincenty

parser = argparse.ArgumentParser(description="Calculate het and Fst.")
parser.add_argument('--cl', help="Cluster for which to run.")
args = parser.parse_args()
cl = args.cl

c_file = '/Volumes/heloderma4/sonal/eco_IBD_oz/data/clustering/Ctenotus_Lerista.clustering.revised.csv'
vcf_dir = '/Volumes/heloderma4/sonal/eco_IBD_oz/snp_calling/variants/'
out_dir = '/Volumes/heloderma4/sonal/eco_IBD_oz/data/divergence/'
dist_file = '/Volumes/heloderma4/sonal/eco_IBD_oz/data/metadata/individual_data_nomissing22Oct15.csv'
mt_file = '/Volumes/heloderma4/sonal/skink_mtDNA/cytb_alignment_30July15.fixed.aln.fa'


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


def get_dist_hash(dist_file):
	# data file with individuals and lats / longs
	inds = pd.read_csv(dist_file, sep=',')
	# get rid of undefined rows
	inds = inds[np.isfinite(inds.lat)]
	inds = inds[np.isfinite(inds.lon)]

	# gets a dictionary where key is sample id and value is tuple of lat / long
	latlong = dict([(x, (y, z)) for x, y, z in zip(inds.sample_id, inds.lat, inds.lon)])
	return latlong


def get_distance(latlong, ind1, ind2):
        if ind1 in latlong and ind2 in latlong:
                # these distances are calculated based on a model
                #       in which earth is an oblate spheroid
                # more typically used is the great-circle distance
                #       which assumes spherical earth but this isn't true
                dist = vincenty(latlong[ind1], latlong[ind2]).meters
                return round(dist, 2)
        else:
                return np.nan


def get_clusters(c_file):
        d = pd.read_csv(c_file)
        d = d[d.GMYC_RAxML2.notnull()]
        d = d.groupby('GMYC_RAxML2')

        clusters = dict([(name, sorted(group.sample.tolist())) for name, group in d])

        return clusters


def get_mt_dist(mt, ind1, ind2):
	allowed = ['A', 'T', 'C', 'G', 'a', 't', 'c', 'g']

	if ind1 in mt and ind2 in mt:
		seq1 = mt[ind1]
		seq2 = mt[ind2]

		diff = 0
		denom = 0

		for bp1, bp2 in zip(seq1, seq2):
			if bp1 in allowed and bp2 in allowed:
				denom += 1
				if bp1 != bp2:
					diff += 1
		diff = diff / float(denom)
		return (diff, denom)
	else:
		return (np.nan, np.nan)


def fst_estimator(counts, sample_sizes):
	'''
	modified from G. Bradford's R code in bedassle
	'calculate.pairwise.Fst'

	both inputs are arrays where each row is an individual
	and each column is a SNP
	'''

	counts = np.array(counts)
	sample_sizes = np.array(sample_sizes).astype('float')

	pop_af = counts / sample_sizes
	mean_af = np.sum(counts, axis=0) / np.sum(sample_sizes, axis = 0)

	MSP = np.sum((pop_af - mean_af) ** 2 * sample_sizes, axis=0)
	MSG = np.sum((1 - pop_af) * pop_af * sample_sizes, axis=0) \
			* (1 / np.sum(sample_sizes - 1, axis=0))
	n_c = np.sum(sample_sizes, axis = 0) - np.sum(sample_sizes ** 2, axis=0) \
			/ np.sum(sample_sizes, axis=0)

	fst = np.sum(MSP - MSG) / np.sum(MSP + (n_c - 1) * MSG)

	return fst


def fst_reich(counts, sample_sizes):
	counts1 = np.array(counts)[0]
	counts2 = np.array(counts)[1]

	sample_sizes1 = np.array(sample_sizes).astype('float')[0]
	sample_sizes2 = np.array(sample_sizes).astype('float')[1]

	h1 = counts1 * (sample_sizes1 - counts1) / (sample_sizes1 * (sample_sizes1 - 1))
	h2 = counts2 * (sample_sizes2 - counts2) / (sample_sizes2 * (sample_sizes2 - 1))
	
	N = []
	D = []

	for _a1, _a2, _n1, _n2, _h1, _h2 in zip(counts1, counts2, sample_sizes1, sample_sizes2, h1, h2):
		n = ((_a1 / _n1) - (_a2 / _n2)) ** 2 - (_h1 / _n1) - (_h2 / _n2)
		N.append(n)
		d = n + _h1 + _h2
		D.append(d)

	F = np.sum(N) / np.sum(D)

	return F


def get_divergence(cl, inds, vcf_dir, out_dir, mt, latlong):
	diff = { '0/0': {'0/1': 0.5, '1/1': 1, '0/0': 0},
	        '0/1': {'0/1': 0, '1/1': 0.5, '0/0': 0.5},
	        '1/1': {'0/1': 0.5, '1/1': 0, '0/0': 1} }
	count = {'0/0': 0, '1/1': 2, '0/1': 1, './.': np.nan}

	file = '%s%s.final.vcf.gz' % (vcf_dir, cl)

	div = {}
	for ix, ind1 in enumerate(inds):
		div[ind1] = {}
		for ind2 in inds[(ix + 1):]:
			div[ind1][ind2] = {'diff': 0, 'denom': 0}
		

	# for calculating fst
	counts = dict([(ind, []) for ind in inds])

	f = gzip.open(file, 'r')
	for l in f:
		if not re.search('#', l) and not re.search('INDEL', l):
			d = re.split('\s+', l.rstrip())
			# don't mess with multiallelics
			if len(re.split(',', d[4])) == 1:
				genos = d[9:]
				genos = [re.search('^(\S\/\S)', x).group(1) for x in genos]

				# variable site to be used in fst
				if d[4] in ['A', 'T', 'C', 'G']:
					for ind, geno in zip(inds, genos):
						counts[ind].append(count[geno])

				# get divergence data
				genos = dict(zip(inds, genos))
				for ind1 in div:
					for ind2 in div[ind1]:
						if genos[ind1] != './.' and genos[ind2] != './.':
							div[ind1][ind2]['denom'] += 1
							div[ind1][ind2]['diff'] += diff[genos[ind1]][genos[ind2]]
	f.close()


	out = '%s%s.divergence.csv' % (out_dir, cl)
	o = open(out, 'w')
	o.write('cl,ind1,ind2,geo_dist,nuc_dxy,nuc_denom,fst,fst_denom,mt_dxy,mt_denom\n')
	for ind1 in div:
		for ind2 in div[ind1]:
			dxy_denom = div[ind1][ind2]['denom']
			if dxy_denom > 0:
				dxy = div[ind1][ind2]['diff'] / float(div[ind1][ind2]['denom'])
			else:
				dxy = np.nan

			alleles = np.array([counts[ind1], counts[ind2]])
			to_mask = np.any(np.isnan(alleles), axis=0)
			alleles = alleles[:, -to_mask]
			if len(alleles[0]) > 0:
				sizes = [[2] * len(alleles[0]), [2] * len(alleles[0])]
				fst = fst_reich(alleles, sizes)
			else:
				fst = np.nan

			# get geographic distance
			geo_dist = get_distance(latlong, ind1, ind2)

			# get mito distance
			(mt_dist, mt_denom) = get_mt_dist(mt, ind1, ind2)

			o.write('%s,%s,%s,%s,%.6f,%s,%.6f,%s,%.6f,%s\n' % (cl, ind1, ind2, geo_dist, dxy, dxy_denom, fst, len(alleles[0]), mt_dist, mt_denom))

	o.close()


mt = get_mt(mt_file)
latlong = get_dist_hash(dist_file)
clusters = get_clusters(c_file)
get_divergence(cl, clusters[cl], vcf_dir, out_dir, mt, latlong)
