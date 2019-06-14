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
		description="Calculate pi for a lineage.",
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


def get_diversity(lineage, inds, vcf, outdir):
	# keep track of it all
	all = {}
	for ind in inds:
		all[ind] = {'pi': {'sum': 0, 'sites': 0}}
	all['all'] = {'pi': {'sum': 0, 'sites': 0},
					  'het': {'sum': 0, 'sites': 0}}

	allowed = ['0/0', '0/1', '1/1']

	f = gzip.open(vcf, 'r')
	for l in f:
		l = l.decode('utf-8')
		if not re.search('#', l) and not re.search('INDEL', l):
			d = re.split('\s+', l.rstrip())
			# don't mess with multiallelics
			if len(re.split(',', d[4])) == 1:
				genos = [re.search('^(\S\/\S)', x).group(1) for x in d[9:]]
				# look at inds
				for ind, gen in zip(inds, genos):
					if gen in allowed:
						all[ind]['pi']['sites'] += 1
						if gen == '0/1':
							all[ind]['pi']['sum'] += 1
						
				genos = [x for x in genos if x in allowed]

				# only do it if there is at least one non-missing site
				if len(genos) > 0:
					# calculate proportion hets
					het_prop = genos.count('0/1') / float(len(genos))
					all['all']['het']['sum'] += het_prop
					all['all']['het']['sites'] += 1

					alleles = []
					for geno in genos:
						alleles += re.split('/', geno)
					
					if len(set(alleles)) > 1:
						n_c = 0
						n_diff = 0
						for i in alleles:
							for j in alleles:
								n_c += 1
								if i != j:
									n_diff += 1
						pi_prop = n_diff / float(n_c)
					else:
						pi_prop = 0
					all['all']['pi']['sum'] += pi_prop
					all['all']['pi']['sites'] += 1

	f.close()
		
	out = os.path.join(outdir, '%s_diversity.csv' % lineage)
	o = open(out, 'w')
	o.write('type,lineage,n_inds,ind,pi,pi_denom,het,het_denom\n')
	if all['all']['pi']['sites'] > 0:	
		pi = all['all']['pi']['sum'] / float(all['all']['pi']['sites'])
		het = all['all']['het']['sum'] / float(all['all']['het']['sites'])
	else:
		pi = np.nan
		het = np.nan	
	o.write('ALL,%s,%s,NA,%.6f,%s,%.6f,%s\n' % 
				(lineage, len(inds), pi, all['all']['pi']['sites'], het, 
		all['all']['het']['sites']))
	for ind in inds:
		if all[ind]['pi']['sites'] > 0:
			pi = all[ind]['pi']['sum'] / float(all[ind]['pi']['sites'])
			het = pi
		else:
			pi = np.nan
			het = np.nan
		o.write('IND,%s,%s,%s,%.6f,%s,%.6f,%s\n' % 
						(lineage, len(inds), ind, pi, 
						 all[ind]['pi']['sites'], het, 
						 all[ind]['pi']['sites']))
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
	get_diversity(lineage, inds, vcf, outdir)


if __name__ == "__main__":
	main()
