import pandas as pd
import subprocess

file = '/scratch/drabosky_flux/sosi/uce_test/samples.csv'
d = pd.read_csv(file)

name = "match"
nodes = 1 
cpu = 8
mem = 5
hours = 2

# for ix, lineage in enumerate(d['lineage']):
for ix, sample in enumerate(d['sample']):
	sh_out = '%s%s.sh' % (name, ix)
	o = open(sh_out, 'w')

	o.write("#PBS -N %s%s\n" % (name, ix))
	o.write("#PBS -M sosi@umich.edu\n")
	o.write("#PBS -A drabosky_flux\n")
	o.write("#PBS -l qos=flux\n")
	o.write("#PBS -q flux\n")
	o.write("#PBS -l nodes=%s:ppn=%s,mem=%sgb\n" % (nodes, cpu, mem))
	o.write("#PBS -l walltime=%s:00:00\n" % hours)
	o.write("#PBS -j oe\n")
	o.write("#PBS -V\n")

	o.write("\n")
	o.write("module load lsa\n")
	o.write("module load java/1.8.0\n")

	o.write("python ~/squamateUCE/align_reads1.py --sample %s --file /scratch/drabosky_flux/sosi/uce_test/samples.csv --outdir /scratch/drabosky_flux/sosi/uce_test/alignments/ --dir /scratch/drabosky_flux/sosi/uce_test/ --bwa ~/bin/bwa-0.7.12/bwa --samtools ~/bin/samtools-1.3.1/samtools --gatk ~/bin/GenomeAnalysisTK.jar --picard ~/bin/picard-tools-2.4.1/picard.jar --CPU %s --mem %s" % (sample, cpu, mem))
	#o.write("python ~/squamateUCE/make_PRG.py --lineage %s --file /scratch/drabosky_flux/sosi/uce_test/samples.csv --mdir /scratch/drabosky_flux/sosi/uce_test/matches/ --adir /scratch/drabosky_flux/sosi/uce_test/trinity_assembly/ --outdir /scratch/drabosky_flux/sosi/uce_test/PRG --keep easy_recip_match" % lineage)
	# o.write("python ~/squamateUCE/match_contigs_to_probes.py --blat ~/bin/blat --sample %s --dir /scratch/drabosky_flux/sosi/uce_test/ --evalue 1e-30 --outdir /scratch/drabosky_flux/sosi/uce_test/matches/ --db /scratch/drabosky_flux/sosi/uce_test/uce-5k-probes.fasta" % (sample))
	# o.write("python /home/sosi/squamateUCE/trinity_assembly.py --trinity ~/bin/trinityrnaseq-2.2.0/Trinity --sample %s --readdir /scratch/drabosky_flux/sosi/uce_test/trim_reads --outdir /scratch/drabosky_flux/sosi/uce_test/trinity_assembly/ --mem 55 --CPU 16" % (sample))
	# o.write("python ~/squamateUCE/clean_reads.py --trimjar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar --PEAR ~/bin/pear-0.9.10/pear-0.9.10-bin-64 --outdir /scratch/drabosky_flux/sosi/uce_test/trim_reads/ --sample %s --file %s\n" % (sample, file))
	o.close()

	subprocess.call("qsub %s" % sh_out, shell=True)
