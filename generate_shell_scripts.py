import pandas as pd

file = '/scratch/drabosky_flux/sosi/uce_test/samples.csv'
d = pd.read_csv(file)

name = "clean"
nodes = 1 
cpu = 1
mem = 2
hours = 12

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

	o.write("python ~/squamateUCE/clean_reads.py --trimjar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar --PEAR ~/bin/pear-0.9.10/pear-0.9.10-bin-64 --outdir /scratch/drabosky_flux/sosi/uce_test/trim_reads/ --sample %s --file %s\n" % (sample, file))
	o.close()


