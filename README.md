# atcgtoolkit

export PYTHONPATH=/path/to/atcgtoolkit/build:$PYTHONPATH
export PYTHON=/path/to/atcgtoolkit/build/bioinfodevelop:$PYTHONPATH
chmod +x atcgtools
run directly:
atcgtools -h
or 
python atcgtools.py [analysistype]


atcgtools [analysistype]={countgaps,inittopleveltable,fillcontextNoutgroup,VJtransf,VJgendadi,Detectsignalacrossgenome,VCataAnno,test}

[analysistype]
 VJtransf：
 convert vcf into which format depands on -f 
 pedmap (default): map ped plink format (if -r is proviede, final map ped file only contian unlinked sites that r2 less than value -r assigned )
 genosnp:
( most functions below require a special designed mysql/mariandb table that cantains all SNPs of the species and the state of the polymorphism in each populations/breeds for each SNPs, 
  and the outgroup genotype information to determine the ancestral allele. This can be done by the following two steps:
 atcgtools inittopleveltable allindividuals.vcf
 atcgtools fillcontextNoutgroup 
 )
 
 inittopleveltable:
 	create the toplevel variants table for a species and load SNP/INDEL records of allindividuals.vcf into the table.
 	after this step, toplevel talbe in mysql databases contain the following fields:chrID,snp_pos,snpID,ref_allele,context,alt_alle1,alt_alleN,presenceflag
 	ref_allele corresponding the 4th(REF) col of the input vcf file and context contain the flank bases of both side of the SNP, empty for INDEL.
 	alt_alle1 conrresponding to the first allele of the 5th(ALT) of the input vcf file and the rest alternated allele will be stored in alt_alleN col.
 	presenceflag, which is paticularly noteworthy. it's a two bits flag. from right(least significant bit) to left, every to bits represent the variant state of the corresponding population,which is assigned in config.propoties file.
 	
 As the SNPs will vary when different populations/individuals used to call, the input vcf here better use as comprehensive an individual/population set as possible. 
 Whereas joint too many individuals/populations to call SNP may supper time consuming, So our solutions consists of two aspects:
 1) temporarily, using a vcf file contains as comprehensive SNPs this species have, then in following analyses programs,which using different vcf files including additional variant records, could add those new records to the toplevel table and fill their outgroup info(as next step).
 2) 
 after the two steps, the presenseflag still empty. It will be set when the assigned vcf files to be used in analyses programs through instance of class alterPresenceflag_Callback, alterDB() function
 presenseflag will be set (through alterDB()) for all SNPs for corresponding populations when alignmultPopSnpPos() function is called.

 fillcontextNoutgroup
 
 fillcontextNoutgroup


 Detectsignalacrossgenome
 subclasses Caculate_popPI,Caculate_popDiv,Caculate_Fst,Caculate_ABB_BAB_BBAA,Caculate_df,Caculate_Hp_master_slave,Caculate_S_ObsExp_difference were used to
 calculate pi,dxy,Fst,abbabaab,df,Hp,SDS that selection occured in ancient, respectively.
	
	e.g. atcgtools Detectsignalacrossgenome -T spotbilled.VCFBAMconfig -T ~/avian/duckpop/configfiles/mallardZJU1.VCFBAMconfig -o EarlyseletedRegion_auto -R ~/avian/duckpop/configfiles/beijing.VCFBAMconfig -R ~/avian/duckpop/configfiles/shaoxing.VCFBAMconfig -R ~/avian/duckpop/configfiles/newdomesticbreeds.VCFBAMconfig -R ~/avian/duckpop/configfiles/gy.VCFBAMconfig -R ~/avian/duckpop/configfiles/sm.VCFBAMconfig -R ~/avian/duckpop/configfiles/jd.VCFBAMconfig -R ~/avian/duckpop/configfiles/cv.VCFBAMconfig -R ~/avian/duckpop/configfiles/campbell.VCFBAMconfig -n 24 -p early -t toplevelDuck_ZJU1ref -w 20000 -s 10000 -1 /home/lrui/BDRobot/src/bioinfodevelop/slave/Detectsignalacrossgenome_producecorrelation_slave.py -2 /home/lrui/BDRobot/src/bioinfodevelop/slave/Detectsignalacrossgenome_slidewin_slave.py -c auto

 VCataAnno
	
	e.g. 
	atcgtools VCataAnno -V spotbilledshaoxingmallardbeijingnewdomesticbreedslcwhitesramplouzr147.majorwild_domestic.vcf -r ZJU1/bjduckallchr.fna -g ZJU1/genomic.gtf -o duckpop/variantanno/ -c ZJU1duckchrominfo -m 80000 -5 3000 -3 3000
	
	ffffffffff
