import argparse
import re
import pandas as pd
import os
import numpy as np

"""
Sonal Singhal
created on 2 August 2016
updated on 16 June 2021
Written assuming: NOTHING!
"""

def get_args():
	parser = argparse.ArgumentParser(
		description="Summarizes assemblies & PRGs.",
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

	# out dir
	parser.add_argument(
		'--outdir',
		type=str,
		default=None,
		help="Full path to dir for quality data output"
		)

	return parser.parse_args()


def get_data(args):
	d = pd.read_csv(args.file)

	gdir = os.path.join(args.dir, 'trinity_assembly')

	pdir = os.path.join(args.dir, 'PRG')

	outdir = args.outdir
	
	if not os.path.isdir(outdir):
		os.mkdir(outdir)

	ind = args.ind
	sp = d.loc[d['sample'] == ind, 'lineage'].tolist()[0]

	return d, gdir, pdir, outdir, ind, sp

def get_seq(seqfile):
	seq = {}
	f = open(seqfile, 'r')
	for l in f:
		if re.search('>', l):
			id = re.search('>(\S+)', l).group(1)
			seq[id] = ''
		else:
			seq[id] += l.rstrip()
	f.close()

	seq2 = [len(y) for x,y in seq.items()]
	seq2 = sorted(seq2)

        tot_length = np.sum(seq2)

	run_count = 0
	for ix, val in enumerate(seq2):
		run_count += val
		if run_count >= (tot_length / 2.0):
			n50 = val
			break

	return seq, seq2, n50


def summarize_contigs(d, gdir, pdir, outdir, ind, sp):
	stats = {}

	a = os.path.join(gdir, '%s.fasta' % ind)
	seq, seqa, n50a = get_seq(a)

	stats['assembled_contig_count'] = len(seqa)
	stats['assembled_contig_totlength'] = np.sum(seqa)
	stats['assembled_contig_n50'] = n50a
	
	prg = os.path.join(pdir, '%s.fasta' % sp)
	seq, seqp, n50p = get_seq(prg)

	stats['annotated_contig_count'] = len(seqp)
	stats['annotated_contig_totlength'] = np.sum(seqp)
	stats['annotated_contig_n50'] = n50p

	loci = list(seq.keys())
	stats['annotated_AHE_count'] = len([x for x in loci if re.search("AHE", x)])
	stats['annotated_uce_count'] = len([x for x in loci if re.search("uce", x)])

        return stats

def print_stats(outdir, ind, sp, stats):
	out = os.path.join(outdir, '%s.assemblyquality.csv' % ind)
	o = open(out, 'w')
	o.write("individual,lineage,metric,value\n")
	for key in stats:
		o.write("%s,%s,%s,%s\n" % (ind, sp, key, stats[key]))
	o.close()


def main():
	args = get_args()
	d, gdir, pdir, outdir, ind, sp = get_data(args)
	stats = summarize_contigs(d, gdir, pdir, outdir, ind, sp)
        print_stats(outdir, ind, sp, stats)

if __name__ == "__main__":
	main()
