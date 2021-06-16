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
	* kseqtest
"""

def get_args():
	parser = argparse.ArgumentParser(
		description="Looks at read counts & lengths for original & trimmed reads. Assumes kseqtest.",
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

	# individual
	parser.add_argument(
		'--ind',
		type=str,
		default=None,
		help='individual to check'
	   )

	# kseqtest
	parser.add_argument(
		'--kseq',
		type=str,
		default=None,
		help='where is your Kseq? include path & binary name'
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

	ind = args.ind

	# information for individual
	row = d.ix[d['sample'] == ind, ].to_dict('list')
        sp = d.loc[d['sample'] == ind, 'lineage'].tolist()[0]

	# where cleaned reads are
	readdir = os.path.join(args.dir, 'trim_reads/')

	orig_reads = [row['read1'][0], row['read2'][0]]

	new_reads = ['%s%s_R1.final.fq.gz' % (readdir, ind),
             '%s%s_R2.final.fq.gz' % (readdir, ind),
             '%s%s_unpaired.final.fq.gz' % (readdir, ind)]

        outdir = args.outdir
        if not os.path.isdir(outdir):
            os.mkdir(outdir)

	return ind, sp, new_reads, orig_reads, outdir


def run_reads_sub(args, reads):
    tot = 0
    bp = 0
    for r in reads:
        p = subprocess.Popen("%s %s" % (args.kseq, r), stdout = subprocess.PIPE, shell=True)
        x = [l.rstrip() for l in p.stdout]
        tot += int(re.search('(\d+)', x[0]).group(1))
        bp += int(re.search('(\d+)', x[1]).group(1))
    return tot, bp


def run_reads(args, new_reads, orig_reads):
	stats = {}

	tot1, bp1  = run_reads_sub(args, new_reads[0:2]) 
	stats['trim_pairreads_count'] =  tot1
	stats['trim_pairreads_totlength'] =  bp1

        tot0, bp0  = run_reads_sub(args, [new_reads[2]])
        stats['trim_unpairreads_count'] =  tot0
        stats['trim_unpairreads_totlength'] =  bp0

	tot2, bp2  = run_reads_sub(args, orig_reads)
	stats['orig_reads_count'] =  tot2
	stats['orig_reads_totlength'] =  bp2

	return stats


def print_stats(outdir, ind, sp, stats):
	out = os.path.join(outdir, '%s.readquality.csv' % ind)
	o = open(out, 'w')
	o.write("individual,lineage,metric,value\n")
	for key in stats:
		o.write("%s,%s,%s,%s\n" % (ind, sp, key, stats[key]))
	o.close()


def main():
	args = get_args()
	ind, sp, new_reads, orig_reads, outdir = get_data(args)
	stats = run_reads(args, new_reads, orig_reads)
	print_stats(outdir, ind, sp, stats)


if __name__ == "__main__":
	main()
