import argparse
import numpy as np
import os
import pandas as pd
import re
import subprocess

"""
Sonal Singhal
created on 28 June 2016
Written assuming:
        * GATK 3.6
	* samtools 1.3.1
	* picard tools
"""

def get_args():
        parser = argparse.ArgumentParser(
                description="Call SNPs. Assumes GATK 3.6",
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

	# samtools
        parser.add_argument(
                '--samtools',
                type=str,
                default=None,
                help='samtools executable, full path.'
               )

	# GATK
	parser.add_argument(
                '--gatk',
                type=str,
                default=None,
                help='GATK executable, full path.'
               )

	# picard tools
	parser.add_argument(
                '--picard',
                type=str,
                default=None,
                help='picard executable, full path.'
               )

        # CPUs
        parser.add_argument(
                '--CPU',
                type=int,
                default=1,
                help='# of CPUs to use in alignment.'
               )

	# assembly dir
	parser.add_argument(
                '--gdir',
                type=str,
                default=None,
                help="Full path to trinity assembly dir if "
                     "you aren't running in context of pipeline."
                )

        # match file
        parser.add_argument(
                '--mfile',
                type=str,
                default=None,
                help="Full path to match file for assembly."
                )

	# alignment dir
	parser.add_argument(
                '--adir',
                type=str,
                default=None,
                help="Full path to dir with bam files if "
                     "you aren't running in context of pipeline."
                )

	# out dir
        parser.add_argument(
                '--outdir',
                type=str,
                default=None,
                help="Full path to dir for coverage output "
                     "if you aren't running in context of pipeline."
                )

	return parser.parse_args()


def get_data(args):
	d = pd.read_csv(args.file)
	sps = {}

	for x in d['sample']:
		lineage = d.ix[d['sample'] == x, 'lineage'].tolist()[0]

		if args.gdir:
			a = os.path.join(args.gdir, '%s.fasta' % x)
			m = os.path.join(args.mfile)
		else:
			a = os.path.join(args.dir, 'trinity_assembly', '%s.fasta' % x)
			m = os.path.join(args.dir, 'matches', '%s_matches.csv' % x)
	
		if args.adir:
			align = os.path.join(args.adir, '%s.realigned.dup.rg.mateFixed.sorted.recal.bam' % x)
		else:
			align = os.path.join(args.dir, 'alignments', '%s.realigned.dup.rg.mateFixed.sorted.recal.bam' % x)

		sps[x] = {'assembly': a, 'align': align, 'lineage': lineage, 'match': m}

	if args.outdir:
		outdir = args.outdir
	else:
		outdir = os.path.join(args.dir, 'coverage')
	
	if not os.path.isdir(outdir):
		os.mkdir(outdir)

	return sps, outdir


def get_mapped_count(args, sps, stats):
	for sp in sps:
		print('  * %s' % sp)
		bam = sps[sp]['align']
		p = subprocess.Popen("%s flagstat %s" % (args.samtools, bam), stdout=subprocess.PIPE, shell=True)
		x = [l.rstrip() for l in p.stdout]

		stats[sp]['orig_reads'] = int(re.search('^(\d+)', x[0]).group(1))
		stats[sp]['map_reads'] = round(float(re.search('([\d\.]+)%', x[4]).group(1)) / 100., 3)
		stats[sp]['paired'] = round(float(re.search('([\d\.]+)%', x[8]).group(1)) / 100., 3)
		stats[sp]['duplicates'] = round(int(re.search('^(\d+)', x[3]).group(1)) / float(stats[sp]['orig_reads']), 3)

	return stats


def get_stats(sps):
	stats = {}

	types = ['gene', 'AHE', 'uce', 'all']
	vals = ['num', 'mean_length', 'median_length',
 		'mean_cov', 'sites_10x']

	for sp in sps:
		stats[sp] = {}
		stats[sp]['orig_reads'] = 0
		stats[sp]['map_reads'] = 0
		stats[sp]['paired'] = 0
		stats[sp]['duplicates'] = 0
		stats[sp]['median_insert_size'] = 0
		for type in types:
			for val in vals:
				stats[sp][type + '_' + val] = 0

	return stats


def run_coverage(args, sps, outdir):
	for sp in sps:
		out = os.path.join(outdir, sp)
		if not os.path.isfile(out):
			subprocess.call("java -jar %s -T DepthOfCoverage -R %s -o %s -I %s -ct 5 --outputFormat csv "
                        	        "--omitPerSampleStats --omitLocusTable --omitIntervalStatistics -nt %s" % (args.gatk,
                        	        sps[sp]['prg'], out, sps[sp]['align'], args.CPU), shell=True)
		sps[sp]['cov'] = out

	return sps


def run_insert(args, sps):
        for sp in sps:
		out1 = re.sub('.bam', '_insert.txt', sps[sp]['align'])
		out2 = re.sub('.bam', '_insert.pdf', sps[sp]['align'])
		if not os.path.isfile(out1):
			subprocess.call("java -jar %s CollectInsertSizeMetrics I=%s O=%s H=%s" % 
                        	        (args.picard, sps[sp]['align'], out1, out2), shell=True)
		f = open(out1, 'r')
		for l in f:
			if re.search('^## METRICS', l):
				d = f.next()
				d = f.next()
				d = re.split('\s+', d)
				sps[sp]['median_insert_size'] = int(d[0])
				break

        return sps


def get_contig_length(args, sps, stats):
	for sp in sps:
		f = open(sps[sp]['assembly'], 'r')

		oseq = {}
		id = ''
		for l in f:
			if re.search('>', l):
				id = re.search('>(\S+)', l.rstrip()).group(1)
				oseq[id] = ''
			else:
				oseq[id] += l.rstrip()
		f.close()

		seq = {}
		d = pd.read_csv(sps[sp]['match'])
		for ix, row in d.iterrows():
			if row['status'] in ['easy_recip_match', 'complicated_recip_match']:
				seq[row['match']] = len(oseq[row['contig']])	

		res = {}
		types = ['gene', 'AHE', 'uce', 'all']
		for type in types:
			res[type] = []

		for id, s in seq.items():
			type = re.search('^([^-]+)', id).group(1)
			res[type].append(s)
			res['all'].append(s)
			stats[sp]['%s_num' % type] += 1
			stats[sp]['all_num'] += 1		

		for type in res:
			a = res[type]
			if len(a) > 1:
				stats[sp]['%s_mean_length' % type] = round(np.mean(a), 2)
				stats[sp]['%s_median_length' % type] = round(np.median(a), 2)
			else:
				stats[sp]['%s_mean_length' % type] = np.nan
                                stats[sp]['%s_median_length' % type] = np.nan
	
	return stats


def get_coverage(args, sps, stats):
	for sp in sps:
		print('  * %s' % sp)
		f = open(sps[sp]['cov'], 'r')
		head = f.next()

		res = {}
		types = ['gene', 'AHE', 'uce', 'all']
		for type in types:
			res[type] = {'cov': 0, 'sites': 0}

		for l in f:
			type = re.search('^([^-]+)', l).group(1)
			depth = int(re.search('(\d+)$', l.rstrip()).group(1))

			res[type]['cov'] += depth
			res[type]['sites'] += 1
			res['all']['cov'] += depth
			res['all']['sites'] += 1

			if depth >= 10:
				stats[sp]['%s_sites_10x' % type] += 1
				stats[sp]['all_sites_10x'] += 1
		f.close()

		for type in res:
			if res[type]['sites'] > 0:
				stats[sp]['%s_mean_cov' % type] = round(res[type]['cov'] / float(res[type]['sites']), 2)
			else:
				stats[sp]['%s_mean_cov' % type] = np.nan

	return stats


def print_stats(outdir, sps, stats):
	out = os.path.join(outdir, 'summary_statistics.csv')
	o = open(out, 'w')

	keys = sorted(stats[sps.keys()[0]].keys())
	o.write('%s,%s,%s\n' % ('sample', 'lineage', ','.join(keys)))	
	for sp in sps:
		vals = ['%s' % stats[sp][key] for key in keys]
		o.write('%s,%s,%s\n' % (sp, sps[sp]['lineage'], ','.join(vals)))
	o.close()


def main():
	args = get_args()
	sps, outdir = get_data(args)
	stats = get_stats(sps)
        print("Running contig lengths ...")
        stats = get_contig_length(args, sps, stats)
	print("Running insert ...")
	sps = run_insert(args, sps)
	print("Running coverage ...")
	sps = run_coverage(args, sps, outdir)
	print("Running map counts ...")
	stats = get_mapped_count(args, sps, stats)
	print("Running coverage ...")
	stats = get_coverage(args, sps, stats)
	print("Printing ...")
	print_stats(outdir, sps, stats)

if __name__ == "__main__":
	main()
