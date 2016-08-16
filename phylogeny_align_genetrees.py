import argparse
import glob
import os
import multiprocessing as mp
import pandas as pd
import re
import subprocess
import random

"""
Sonal Singhal
created on 23 June 2016
Written assuming:
	* mafft 7.294
	* RAxML 8.2.4
This script borrows heavily from:
https://github.com/faircloth-lab/phyluce/blob/master/bin/align/phyluce_align_seqcap_align
"""

def get_args():
    parser = argparse.ArgumentParser(
        description="Align, possibly trim, and infer gene trees for UCE loci.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--dir",
        type=str,
        default=None,
        help="The base directory if running the pipeline."
    )

    parser.add_argument(
        "--outdir",
        type=str,
	default=None,
        help="The directory with the phylogeny, if "
             "not using in context of a pipeline."
    )

    '''
    parser.add_argument(
        "--no-trim",
        action="store_true",
        default=False,
        help="""Align, but DO NOT trim alignments."""
    )

    parser.add_argument(
        "--window",
        type=int,
        default=20,
        help="""Sliding window size for trimming."""
    )

    parser.add_argument(
        "--proportion",
        type=float,
        default=0.65,
        help="""The proportion of taxa required to have sequence at alignment ends."""
    )

    parser.add_argument(
        "--threshold",
        type=float,
        default=0.65,
        help="""The proportion of residues required across the window in """ +
        """proportion of taxa."""
    )

    parser.add_argument(
        "--max-divergence",
        type=float,
        default=0.20,
        help="""The max proportion of sequence divergence allowed between any row """ +
        """of the alignment and the alignment consensus."""
    )

    parser.add_argument(
        "--min-length",
        type=int,
        default=100,
        help="""The minimum length of alignments to keep."""
    )
   '''

    parser.add_argument(
        "--CPU",
        type=int,
        default=1,
        help="""Process alignments in parallel using --CPU for alignment. """ +
        """This is the number of PHYSICAL CPUs."""
    )

    parser.add_argument(
	"--mafft",
	type=str,
	default=None,
	help="Full path to mafft executable."
	)

    parser.add_argument(
        "--raxml",
        type=str,
        default=None,
        help="Full path to RAxML executable."
        )

    return parser.parse_args()


def get_dir(args):
	if not args.outdir:
		outdir = os.path.join(args.dir, 'phylogeny', 'alignments')
		treedir = os.path.join(args.dir, 'phylogeny', 'gene_trees')
	else:
		outdir = os.path.join(args.outdir, 'alignments')
		treedir = os.path.join(args.outdir, 'gene_trees')	
	
	return outdir, treedir


def align(params):
	file, mafft = params

	aln_out = file.replace('.fasta', '.fasta.aln')
	proc = subprocess.call("%s --maxiterate 1000 --globalpair "
                               "--adjustdirection --quiet %s > %s" %
			       (mafft, file, aln_out), shell=True)

	os.remove(file)
	return aln_out


def run_alignments(outdir, args):
	files = glob.glob(outdir + '/*fasta')
	
	params = zip(files, [args.mafft] * len(files))
	
	if args.CPU > 1:
		pool = mp.Pool(args.CPU)
		alns = pool.map(align, params)
	
	return alns


def convert_phyml(locus_file):
	f = open(locus_file, 'r')
	phy_file = locus_file + '.phy'
	o = open(phy_file, 'w')

	seq = {}
	id = ''
	for l in f:
		if re.search('>', l):
			id = re.search('>(\S+)', l.rstrip()).group(1)
			if re.search('^_R_', id):
				id = re.sub('^_R_', '', id)
			seq[id] = ''
		else:
			seq[id] += l.rstrip()
	f.close()

	o.write(' %s %s\n' % (len(seq), len(seq.values()[0])))
	for sp, s in seq.items():
		o.write('%s   %s\n' % (sp, s))
	o.close()

	return phy_file


def sub_raxml(file, outdir, raxml):

	locus = re.sub('^.*/', '', file)
	locus = re.sub('\.fa.*', '', locus)

	os.chdir(outdir)
	subprocess.call('%s -x %s -# 100 -p %s -m GTRCAT -f a -n %s -s %s' % 
                        (raxml, random.randint(0,1000), random.randint(0,1000), 
                        locus, file), shell=True)

	orig_boot = 'RAxML_bootstrap.%s' % locus
	orig_tree = 'RAxML_bipartitions.%s' % locus

	new_boot = '%s.bootstrap.trees' % locus
	new_tree = '%s.bestTree.tre' % locus

	os.rename(orig_boot, new_boot)
	os.rename(orig_tree, new_tree)

	subprocess.call("rm RAxML_*%s" % locus, shell=True) 

	return new_tree


result_list = []
def log_result(result):
    # This is called whenever foo_pool(i) returns a result.
    # result_list is modified only by the main process, not the pool workers.
    result_list.append(result)


def run_raxml(outdir, treedir, alns, args):
	
	if not os.path.isdir(treedir):
		os.mkdir(treedir)

	if args.CPU > 1:
                pool = mp.Pool(args.CPU)
        	phys = pool.map(convert_phyml, alns)
	
		dirs = [treedir] * len(phys)
		raxml = [args.raxml] * len(phys)
		for i in range(len(phys)):
			pool.apply_async(sub_raxml, args=(phys[i], treedir, args.raxml, ), callback=log_result)
		pool.close()
		pool.join()

def main():
	args = get_args()
	outdir, treedir = get_dir(args)	
	# alns = run_alignments(outdir, args)
	alns = glob.glob("/scratch/drabosky_flux/sosi/brazil/phylogeny/alignments/*aln")
	run_raxml(outdir, treedir, alns, args)	

if __name__ == "__main__":
	main()
