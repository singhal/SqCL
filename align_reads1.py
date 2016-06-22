import argparse
import os
import pandas as pd
import re
import subprocess

"""
Sonal Singhal
created on 21 June 2016
Written assuming:
	* bcftools 1.3.1
	* samtools 1.3.1
	* GATK 3.6
	* picard 2.4.1
"""

def get_args():
	parser = argparse.ArgumentParser(
		description="Align reads to lineage, step 1",
        	formatter_class=argparse.ArgumentDefaultsHelpFormatter
		)

	# sample
	parser.add_argument(
                '--sample',
                type=str,
                default=None,
                help='Sample for which to run script.'
                )

	# file
	parser.add_argument(
                '--file',
                type=str,
                default=None,
                help='File with sample info.'
                )

	# outdir
	parser.add_argument(
                '--outdir',
                type=str,
                default=None,
                help='Output directory for alignments.'
                )

	# basedir
	parser.add_argument(
                '--dir',
                type=str,
                default=None,
                help="Full path to base dir with reads & assemblies."
                )

	# read1
	parser.add_argument(
                '--read1',
                type=str,
                default=None,
                help="Full path to read 1 file if "
                     "you aren't running in context of pipeline."
                )

	# read2
	parser.add_argument(
		'--read2',
		type=str,
		default=None,
		help="Full path to read 2 file if "
		     "you aren't running in context of pipeline."
		)

	# read u
        parser.add_argument(
                '--un',
                type=str,
                default=None,
                help="Full path to unpaired file if "
                     "you aren't running in context of pipeline."
                )

	# PRG
	parser.add_argument(
                '--prg',
                type=str,
                default=None,
                help="Full path to pseudoref genome if "
                     "you aren't running in context of pipeline."
                )

	# bwa
	parser.add_argument(
                '--bwa',
                type=str,
                default=None,
                help='bwa executable, full path.'
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

	# picard
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

	# memory
        parser.add_argument(
                '--mem',
                type=int,
                default=1,
                help='Memory available, as an int, in terms of Gb.'
               )


	return parser.parse_args()


def get_info(args):
	if not args.read1 and not args.read2:
		read1 = os.path.join(args.dir, 'trim_reads', '%s_R1.final.fq.gz' % args.sample)
		read2 = os.path.join(args.dir, 'trim_reads', '%s_R2.final.fq.gz' % args.sample)
		un = os.path.join(args.dir, 'trim_reads', '%s_unpaired.final.fq.gz' % args.sample)
		reads = [read1, read2, un]
	else:
		reads = [args.read1, args.read2, args.un]

	d = pd.read_csv(args.file)
	lineage = d.ix[d['sample'] == args.sample, 'lineage'].tolist()[0]

	if not args.prg:
		genome = os.path.join(args.dir, 'PRG', '%s.fasta' % lineage)
	else:
		genome = args.prg

	return reads, lineage, genome


def prepare_seq(args, genome):
	if not os.path.isfile(genome + '.bwt'):
		subprocess.call("%s index %s" % (args.bwa, genome), shell=True)
	if not os.path.isfile(genome + '.fai'):
		subprocess.call("%s faidx %s" % (args.samtools, genome), shell=True)
	out = re.sub('.fa.*', '.dict', genome)
	if not os.path.isfile(out):
		subprocess.call("java -jar %s CreateSequenceDictionary R=%s O=%s" % 
                                (args.picard, genome, out), shell=True)


def align_seq(args, r, seq):
	if not os.path.isdir(args.outdir):
		os.mkdir(args.outdir)

	out1a = '%s%s_1.sam' % (args.outdir, args.sample)
	out1b = '%s%s_2.sam' % (args.outdir, args.sample)
	out2a = '%s%s.mateFixed.bam' % (args.outdir, args.sample)
	out2b = '%s%s.bam' % (args.outdir, args.sample)
	out3a = '%s%s_1.mateFixed.sorted.bam' % (args.outdir, args.sample)
	out3b = '%s%s_2.mateFixed.sorted.bam' % (args.outdir, args.sample)
	out3 = '%s%s.mateFixed.sorted.bam' % (args.outdir, args.sample)
	out4 = '%s%s.rg.mateFixed.sorted.bam' % (args.outdir, args.sample)
	intervals = '%s%s.intervals' % (args.outdir, args.sample)
	out5 = '%s%s.realigned.rg.mateFixed.sorted.bam' % (args.outdir, args.sample)
	
	tmpdir = os.path.join(args.outdir, args.sample)
	if not os.path.isdir(tmpdir):
		os.mkdir(tmpdir)

	# align
	subprocess.call("%s mem -t %s %s %s %s > %s" % (args.bwa, args.CPU, seq, r[0], r[1], out1a), shell=True)
	subprocess.call("%s mem -t %s %s %s > %s" % (args.bwa, args.CPU, seq, r[2], out1b), shell=True)
	# fixmate
	subprocess.call("%s fixmate -O bam %s %s" % (args.samtools, out1a, out2a), shell=True)
	subprocess.call("%s view -b %s > %s" % (args.samtools, out1b, out2b), shell=True)
	# sorted
	subprocess.call("%s sort -O bam -o %s -T %s %s" % (args.samtools, out3a, tmpdir, out2a), shell=True)
	subprocess.call("%s sort -O bam -o %s -T %s %s" % (args.samtools, out3b, tmpdir, out1b), shell=True)
	subprocess.call("%s merge %s %s %s" % (args.samtools, out3, out3a, out3b), shell=True)
	# readgroup
	subprocess.call("java -jar %s AddOrReplaceReadGroups INPUT=%s OUTPUT=%s RGLB=%s RGPL=Illumina RGPU=%s RGSM=%s" % 
                          (args.picard, out3, out4, args.sample, args.sample, args.sample), shell=True)
	# indel target
	subprocess.call("%s index %s" % (args.samtools, out4), shell=True)
	subprocess.call("java -Xmx%sg -jar %s -T RealignerTargetCreator -R %s -I %s -o %s -nt %s" % 
                        (args.mem, args.gatk, seq, out4, intervals, args.CPU), shell=True)
	# indel realigner
	subprocess.call("java -Xmx%sg -jar %s -T IndelRealigner -R %s -I %s -targetIntervals %s -o %s" % 
			(args.mem, args.gatk, seq, out4, intervals, out5), shell=True)
	

	# remove the files
	[os.remove(x) for x in [out1a, out1b, out2a, out2b, out3a, out3b, out3, out4, intervals, out4 + '.bai']]
	
	# remove the dir
	os.rmdir(tmpdir)


def main():
	# get arguments
	args = get_args()
	reads, lineage, genome = get_info(args)
	# prep sequence
	prepare_seq(args, genome)
	# do the alignments
	align_seq(args, reads, genome)

if __name__ == "__main__":
	main()
