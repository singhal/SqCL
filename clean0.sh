#PBS -N clean1
#PBS -M sosi@umich.edu
#PBS -A drabosky_flux
#PBS -l qos=flux
#PBS -q flux
#PBS -l nodes=1:ppn=1,pmem=2gb
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -V

python ~/squamate_UCE/clean_reads.py --trimjar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar --PEAR ~/bin/pear-0.9.10/pear-0.9.10-bin-64 --outdir /scratch/drabosky_flux/sosi/uce_test/trim_reads/ --sample Mus_musculus --file /scratch/drabosky_flux/sosi/uce_test/samples.csv
