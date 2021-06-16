import argparse
import numpy as np
import os
import pandas as pd
import re
import gzip
import subprocess

"""
Sonal Singhal
created on 28 June 2016
updated on 16 June 2021
Written assuming:
"""

def get_args():
	parser = argparse.ArgumentParser(
		description="Looks at alignment quality.",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
		)

	# file
	parser.add_argument(
		'--file',
		type=str,
		default=None,
		help='File with sample info.'
		)

	# basedir
	parser.add_argument(
		'--dir',
		type=str,
		default=None,
		help="Full path to base dir with reads & assemblies & "
			 "everything else."
		)

	# lineage
	parser.add_argument(
		'--lineage',
		type=str,
		default=None,
		help='lineage to run'
	   )

	# samtools
	parser.add_argument(
		'--samtools',
		type=str,
		default=None,
		help='samtools executable, full path.'
	   )

	# outdir
	parser.add_argument(
		'--outdir',
		type=str,
		default=None,
		help='where should I output quality stats?'
	   )

	return parser.parse_args()


def get_data(args):
	d = pd.read_csv(args.file)
	sp = args.lineage

	inds = d.ix[d['lineage'] == sp, "sample"].tolist()

	bams = [os.path.join(args.dir, 'alignments', '%s.dup.rg.mateFixed.sorted.recal.bam' % x) for x in inds]
	vcf = os.path.join(args.dir, 'variants', '%s.raw.vcf.gz' % sp)

	outdir = args.outdir	
	if not os.path.isdir(outdir):
		os.mkdir(outdir)

	return sp, inds, bams, vcf, outdir


def get_mapped_count(args, inds, bams, stats):

	for ind, bam in zip(inds, bams):
		p = subprocess.Popen("%s flagstat %s" % (args.samtools, bam), stdout=subprocess.PIPE, shell=True)
		x = [l.rstrip() for l in p.stdout]

		stats[ind] = {}

		stats[ind]['orig_reads'] = int(re.search('^(\d+)', x[0]).group(1))
		stats[ind]['map_reads'] = round(float(re.search('([\d\.]+)%', x[4]).group(1)) / 100., 3)
		
		prop_paired = int(re.search('(^\d+)', x[8]).group(1))
		all_paired = int(re.search('(^\d+)', x[9]).group(1))		
		single = int(re.search('(^\d+)', x[10]).group(1))

		tot = all_paired + single
		if tot > 0:
			prop_paired = prop_paired / float(all_paired + single)
		else:
			prop_paired = np.nan

		map_reads = int(re.search('(^\d+)', x[4]).group(1))
		if map_reads > 0:
			dup = round(int(re.search('^(\d+)', x[3]).group(1)) / float(map_reads), 3)
		else:
			dup = np.nan

		stats[ind]['paired'] = round(prop_paired, 3)
		stats[ind]['duplicates'] = dup

	return stats

def get_coverage(args, vcf, stats):
	f = gzip.open(vcf, 'r')

	d = {}

	for l in f:
		l = l.decode('utf-8').rstrip()
		if re.search('^#CHROM', l):
			inds = re.split('\t', l)[9:]
			for ind in inds:
				d[ind] = {'num_sites': 0, 'tot_cov': 0, 
							'num_sites_10x': 0, 'tot_cov_10x': 0,
							'num_loci_10x': {}}
		elif not re.search('^#', l):
			x = re.split('\t', l)
			dpix = re.split(':', x[8]).index('DP')

			for ind, g in zip(inds, x[9:]):
				dp = int(re.split(':', g)[dpix])

				d[ind]['num_sites'] += 1
				d[ind]['tot_cov'] += dp

				if dp >= 10:
					d[ind]['num_loci_10x'][x[0]] = 1
					d[ind]['num_sites_10x'] += 1
					d[ind]['tot_cov_10x'] += dp
	f.close()

	for ind in d:
		for key in d[ind]:
			if key == 'num_loci_10x':
				stats[ind][key] = len(d[ind][key])
			else:
				stats[ind][key] = d[ind][key]

	return stats

def print_stats(outdir,  sp, stats):
	out = os.path.join(outdir, '%s.alignquality.csv' % sp)
	o = open(out, 'w')
	o.write("individual,lineage,metric,value\n")
	for ind in stats:
		for key in stats[ind]:
			o.write("%s,%s,%s,%s\n" % (ind, sp, key, stats[ind][key]))
	o.close()


def main():
	args = get_args()
	sp, inds, bams, vcf, outdir = get_data(args)
	stats = {}
	stats = get_mapped_count(args, inds, bams, stats)
	stats = get_coverage(args, vcf, stats)
	print_stats(outdir, sp, stats)

if __name__ == "__main__":
	main()
