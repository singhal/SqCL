# squamateUCE - target capture pipeline
Scripts to work with UCE data from squamates.

## Areas of future improvements
- use read error correction prior to assembly
	- `http://www.genome.umd.edu/quorum.html`
- use cap3 to make PRG rather than random picking
- better trusted SNP set for BQSR
	- maybe do Cortex intersection with raw GATK call
- check SNP calling filters?
	- `https://www.broadinstitute.org/gatk/guide/article?id=3225`
	- `https://www.broadinstitute.org/gatk/guide/article?id=6925`
- phase SNPs

## Notes before you start
This pipeline is best run with these scripts all in a row, but each script can be run on its own. However, this functionality hasn't been tested. 

When running any script, you can see the script arguments using `python script.py --help`.

## Steps
1. Start with demultiplexed data.
2. Create a CSV file summarizing data; note, column order does not matter.
	- sample: sample name
	- read1: complete path to read 1 file
	- read2: complete path to read 2 file
	- adaptor1: adaptor 1 to scrub; put asterisk where barcode goes
	- adaptor2: adaptor 2 to scrub; put asterisk where barcode goes
	- barcode1: barcode 1 to scrub
	- barcode2: barcode 2 to scrub (if it exists)
	- lineage: the group to which the sample belongs; used for assembly and SNP calling
3. Clean reads
	- Remove adaptors using Trimmomatic
	- Merge reads using PEAR
	- Lightly remove low quality sequence using Trimmomatic
	- Ends up creating three files: 2 paired files with ordering retained, and 1 unpaired read file
	- Assumes: Trimmomatic v0.36 and PEAR 0.9.10 are downloaded and working
```
python ~/squamateUCE/clean_reads.py --trimjar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar --PEAR ~/bin/pear-0.9.10/pear-0.9.10-bin-64 --dir /scratch/drabosky_flux/sosi/uce_test/ --sample Anolis_carolinensis --file /scratch/drabosky_flux/sosi/uce_test/samples.csv
```
4. Assemble reads using Trinity
```
python /home/sosi/squamateUCE/trinity_assembly.py --trinity ~/bin/trinityrnaseq-2.2.0/Trinity --sample Anolis_carolinensis --dir /scratch/drabosky_flux/sosi/uce_test/ --mem 55 --CPU 16
```
5. Match assemblies to original targets 
```
python ~/squamateUCE/match_contigs_to_probes.py --blat ~/bin/blat --sample Anolis_carolinensis --dir /scratch/drabosky_flux/sosi/uce_test/ --evalue 1e-30 --db /scratch/drabosky_flux/sosi/uce_test/uce-5k-probes.fasta
```
6. Generate pseudo-reference genomes
```
python ~/squamateUCE/make_PRG.py --lineage l1 --file /scratch/drabosky_flux/sosi/uce_test/samples.csv --dir /scratch/drabosky_flux/sosi/uce_test/ --keep easy_recip_match
```
7. Align reads
```
python ~/squamateUCE/align_reads1.py --sample Mus_musculus --file /scratch/drabosky_flux/sosi/uce_test/samples.csv --dir /scratch/drabosky_flux/sosi/uce_test/ --bwa ~/bin/bwa-0.7.12/bwa --samtools ~/bin/samtools-1.3.1/samtools --gatk ~/bin/GenomeAnalysisTK.jar --picard ~/bin/picard-tools-2.4.1/picard.jar --CPU 1 --mem 1
python ~/squamateUCE/align_reads2.py --lineage l1 --file /scratch/drabosky_flux/sosi/uce_test/samples.csv --dir /scratch/drabosky_flux/sosi/uce_test/ --samtools ~/bin/samtools-1.3.1/samtools --gatk ~/bin/GenomeAnalysisTK.jar --mem 3 --dp 10 --qual 20 --CPU 4
```
8. Call and filter variants
```
python ~/squamateUCE/call_variants.py --lineage l1 --file /scratch/drabosky_flux/sosi/uce_test/samples.csv --dir /scratch/drabosky_flux/sosi/uce_test/ --gatk ~/bin/GenomeAnalysisTK.jar --mem 4 --CPU 4
```
