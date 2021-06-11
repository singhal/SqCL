import argparse
import gzip
import os
import pandas as pd
import re
import subprocess
import random

"""
Sonal Singhal
created on 4 January 2017
Written assuming:
	* samtools 1.3.1
	* GATK 3.6
"""

def get_args():
	parser = argparse.ArgumentParser(
		description="Phase reads per lineage. "
					" Assumes samtools 1.3.1 "
					" and GATK 3.6",
			formatter_class=argparse.ArgumentDefaultsHelpFormatter
		)

	# lineage
	parser.add_argument(
		'--lineage',
		type=str,
		default=None,
		help='Lineage for which to run script.'
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
		help="Full path to base dir with reads & assemblies "
			 "everything else."
		)

	# bgzip
	parser.add_argument(
		'--bgzip',
		type=str,
		default=None,
		help='bgzip executable, full path.'
		)

	# tabix
	parser.add_argument(
		'--tabix',
		type=str,
		default=None,
		help='tabix executable, full path.'
		)

	# GATK
	parser.add_argument(
		'--gatk',
		type=str,
		default=None,
		help='GATK executable, full path.'
		)
	
	# memory
	parser.add_argument(
		'--mem',
		type=int,
		default=1,
		help='Memory available, as an int, in terms of Gb.'
	   )

	# outdir
	parser.add_argument(
		'--outdir',
		type=str,
		default=None,
		help='Output directory for alignments, only needed '
			 'if not running in context of pipeline.'
		)
				
	# bamfiles
	parser.add_argument(
		'--bamfile',
		type=str,
		default=None,
		help="Full path to file with BAM files, listed one "
			 "per line if running not in context of pipeline. "
		)

	# PRG
	parser.add_argument(
		'--prg',
		type=str,
		default=None,
		help="Full path to pseudoref genome if "
			 "you aren't running in context of pipeline."
		)

	# vcf
	parser.add_argument(
		'--vcf',
		type=str,
		default=None,
		help="Full path to VCF if "
			 "you aren't running in context of pipeline."
		)
	
	parser.add_argument(
		"--haplo",
		action="store_true",
		default=False,
		help="Will print out haplotypes if flagged.."
		)

	return parser.parse_args()


def get_files(args):
	# gets the bam files
	if args.bamfile:
		f = open(args.bamfile, 'r')
		files = []
		for l in f:
			files.append(l.rstrip())
		f.close()
	else:
		d = pd.read_csv(args.file)
		samps = d.ix[d['lineage'] == args.lineage, 'sample'].tolist()
		files = []
		for samp in samps:
			file = os.path.join(args.dir, 'alignments', 
											'%s.realigned.dup.rg.mateFixed.sorted.recal.bam' % samp)
			files.append(file)
	# makes sure the order stays consistent
	files = sorted(files)

	# gets the genome
	if args.prg:
		genome = args.prg
	else:
		genome = os.path.join(args.dir, 'PRG', '%s.fasta' % args.lineage)

	# gets the outdir
	if args.outdir:
		outdir = args.outdir
	else:
		outdir = os.path.join(args.dir, 'variants')
	
	if not os.path.isdir(outdir):
		os.mkdir(outdir)

	# gets the vcf
	if args.vcf:
		vcf = args.vcf
	else:
		vcf = os.path.join(args.dir, 'variants', '%s.qual_filtered.cov_filtered.vcf.gz' % args.lineage)

	return files, genome, vcf, outdir


def phase(args, files, genome, vcf, dir):
	out = re.sub('.vcf.gz$', '.phased.vcf.gz', vcf)
	bam = '-I ' + ' -I '.join(files)
	subprocess.call("java -Xmx%sg -jar %s -T ReadBackedPhasing -R %s %s --variant %s -o %s" %
			(args.mem, args.gatk, genome, bam, vcf, out), shell=True)
	return out


def get_seq(genome):
	seq = {}
	f = open(genome, 'r')
	id = ''
	for l in f:
		if re.search('>', l):
			id = re.search('>(\S+)', l.rstrip()).group(1)
			seq[id] = ''
		else:
			seq[id] += l.rstrip()

	for id, s in seq.items():
		seq[id] = len(s) * ['N']

	return seq
	

def get_haplo(outdir, vcf, seq, lineage):
	f = gzip.open(vcf, 'r')

	# for haplotype orientations
	orient = {}

	for l in f:
		if re.search('#CHROM', l):
			samps = re.split('\t', l.rstrip())[9:]
			haplo = {}
			for samp in samps:
				haplo['{}_1'.format(samp)] = {}
				haplo['{}_2'.format(samp)] = {}
				for s in seq:
					# set up the results
					haplo['%s_1' % samp][s] = list(seq[s])
					haplo['%s_2' % samp][s] = list(seq[s])
		elif not re.search('#', l):
			d = re.split('\t', l.rstrip())
			
			c = d[0]
			# correct for indexing
			pos = int(d[1]) - 1
			ref = d[3]
			alt = d[4]
			# if missing, set to N
			allele = {'0': ref, '1': alt, '.': 'N'}

			# only consider biallelic SNPs
			# no indels or multiallelics
			if ref in ['A', 'T', 'C', 'G'] and alt in ['.', 'A', 'T', 'C', 'G']:
				genos = d[9:]
				for ix, fullgeno in enumerate(genos):
					geno = re.search('(\S\S\S)', fullgeno).group(1)
					geno = re.split('/', geno)
					# possible phase!
					if re.search('0/1', fullgeno):
						x = re.split(':', fullgeno)
			
						# this is where the phase data live
						phase = x[4]
					
						# tells the orient phase
						phase_pos = re.search('^(\d+)-', phase).group(1)
						
						# define the orientation of the block at random
						# for new blocks only
						if c not in orient:
							orient[c] = {}
						if phase_pos not in orient[c]:
							orient[c][phase_pos] = random.choice([True, False])
		
						# define the SNP orientation within the block
						# it gets block orientation unless it is in 2 - 1
						invert = orient[c][phase_pos]
						if re.search('-2.*-1', phase):
							invert = not invert
						if invert:
							geno = [geno[1], geno[0]]

					genos[ix] = geno
				
				alleles = [[allele[geno[0]], allele[geno[1]]] for geno in genos]
			
				for samp, allele in zip(samps, alleles):
					name1 = '%s_1' % samp
					name2 = '%s_2' % samp
					haplo[name1][c][pos] = allele[0]
					haplo[name2][c][pos] = allele[1]
							
	suboutdir = os.path.join(outdir, lineage)
	if not os.path.isdir(suboutdir):
		os.mkdir(suboutdir)

	inds = sorted(list(haplo.keys()))
	for loc in seq:
		out = os.path.join(suboutdir, '{}.aln.fasta'.format(loc))
		o = open(out, 'w')
		for ind in inds:
			o.write('>%s\n%s\n' % (ind, ''.join(haplo[ind][loc])))
		o.close()
	

def prepare_vcf(args, vcf):
	unvcf = re.sub('.gz', '', vcf)
	subprocess.call("gunzip %s" % vcf, shell=True)
	subprocess.call("%s %s" % (args.bgzip, unvcf), shell=True)
	subprocess.call("%s -p vcf %s" % (args.tabix, vcf), shell=True)


def main():
	args = get_args()
	files, genome, vcf, outdir = get_files(args)
	prepare_vcf(args, vcf)
	out_vcf = phase(args, files, genome, vcf, outdir)
	if args.haplo:
		seq = get_seq(genome)
		get_haplo(outdir, out_vcf, seq, args.lineage)

if __name__ == "__main__":
	main()
