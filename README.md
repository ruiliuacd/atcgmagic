# atcgtoolkit

export PYTHONPATH=/path/to/atcgtoolkit/build:$PYTHONPATH
export PYTHON=/path/to/atcgtoolkit/build/bioinfodevelop:$PYTHONPATH
chmod +x atcgtools
run directly:
atcgtools -h
or 
python atcgtools.py


atcgtools VJtransf
	convert vcf into which format depands on -f 
	pedmap (default): map ped plink format (if -r is proviede, final map ped file only contian unlinked sites that r2 less than -r assigned value)
	genosnp:
[ most functions below require a special designed mysql/mariandb table that cantains all SNPs of the species and the state of the polymorphism in each populations/breeds for each SNPs, 
  and the outgroup genotype information to determine the ancestral allele. This can be done by the following to steps:
 atcgtools inittopleveltable

 atcgtools fillcontextNoutgroup
 ]