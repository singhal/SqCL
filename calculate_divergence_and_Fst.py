import argparse
import gzip
import os
import numpy as np
import pandas as pd
import re
import subprocess

'''
Sonal Singhal
created on 29 June 2016
'''


def get_args():
	parser = argparse.ArgumentParser(
		description="Calculate Fst and Dxy for a lineage.",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
		)

	# lineage
        parser.add_argument(
                '--lineage',
                type=str,
                default=None,
                help='Lineage for which to make calculations.'
                )

        # sample file
        parser.add_argument(
                '--file',
                type=str,
                default=None,
                help='File with sample info.'
                )
        
        # base dir
        parser.add_argument(
                '--dir',
                type=str,
                default=None,
                help="Base directory as necessary"
                     " when used with pipeline"
                )

	# output dir
        parser.add_argument(
                '--outdir',
                type=str,
                default=None,
                help='Directory to output pop gen stats, '
                     'only necessary if not running '
                     'in context of pipeline'
                )

	# vcfdir
	parser.add_argument(
		'--vcfdir',
		type=str,
		default=None,
		help='Directory with VCFs, '
		     'only necessary if not running '
                     'in context of pipeline'
		)

	return parser.parse_args()


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


def get_divergence(lineage, inds, vcf, outdir):
	diff = { '0/0': {'0/1': 0.5, '1/1': 1, '0/0': 0},
	        '0/1': {'0/1': 0, '1/1': 0.5, '0/0': 0.5},
	        '1/1': {'0/1': 0.5, '1/1': 0, '0/0': 1} }
	count = {'0/0': 0, '1/1': 2, '0/1': 1, './.': np.nan}

	div = {}
	for ix, ind1 in enumerate(inds):
		div[ind1] = {}
		for ind2 in inds[(ix + 1):]:
			div[ind1][ind2] = {'diff': 0, 'denom': 0}
		

	# for calculating fst
	counts = dict([(ind, []) for ind in inds])

	f = gzip.open(vcf, 'r')
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


	out = os.path.join(outdir, '%s_divergence.csv' % lineage)
	o = open(out, 'w')
	o.write('lineage,ind1,ind2,nuc_dxy,nuc_denom,fst,fst_denom\n')
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

			o.write('%s,%s,%s,%.6f,%s,%.6f,%s\n' % 
				(lineage, ind1, ind2, dxy, dxy_denom, fst, len(alleles[0])))

	o.close()


def get_data(args):
	lineage = args.lineage

	d = pd.read_csv(args.file)
	inds = d.ix[d.lineage == lineage, 'sample'].tolist()
	inds = sorted(inds)

	if not args.outdir:
		outdir = os.path.join(args.dir, 'pop_gen')
		vcf = os.path.join(args.dir, 'variants', 
                                   '%s.qual_filtered.cov_filtered.vcf.gz' % args.lineage)
	else:
		outdir = args.outdir
		vcf = os.path.join(args.vcfdir, 
                                   '%s.qual_filtered.cov_filtered.vcf.gz' % args.lineage)

	if not os.path.isdir(outdir):	
		os.mkdir(outdir)

	return lineage, inds, vcf, outdir


def main():
	args = get_args()
	lineage, inds, vcf, outdir = get_data(args)
	get_divergence(lineage, inds, vcf, outdir)

if __name__ == "__main__":
	main()
