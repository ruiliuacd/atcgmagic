# atcgtoolkit

export PYTHONPATH=/path/to/atcgtoolkit/build:$PYTHONPATH
export PYTHON=/path/to/atcgtoolkit/build/bioinfodevelop:$PYTHONPATH
chmod +x atcgtools
run directly:
atcgtools -h
or 
python atcgtools.py [analysistype]


atcgtools [analysistype] [parameters for different analysistype]

[analysistype]={countgaps,inittopleveltable,fillcontextNoutgroup,VJtransf,VJgendadi,Detectsignalacrossgenome,VCataAnno,test}
 VJtransf：
 
 Convert vcf into which format depands on -f. 
 -f=	
 pedmap (default): produce .map and .ped file of plink format (if -r is proviede, final map ped file only contian unlinked sites that r2 less than value -r assigned )
 genosnp: output .geno file also called EIGENSTRST format, .snp and .ind file. 
 	 .geno file contains 1 line per position. Each line contains 1 character per individual:
 	 	0 means zero copies of reference allele
 	 	1 means one copy of reference allele.
 	 	2 means two copies of reference allele.
 	 	9 missing data.
 vcf: output vcf format. 
 
 -v vcf file [-v ...].  When multiple -v vcf files were provided, then the output data will contain inner/outer joined records of all vcf as described in [1](AI 5.8). It provided a flexible way to joint calling SNP and could avoid the significant time consumption of calling SNPs for all individuals simultaneously using GATK.
 e.g. The file 'gatkallchromAll147indvd.vcf' was generated by GATK for seperate chromosomes, using all individuals from all populations. Joint calling for all seven populations of 147 ducks's sample simultaneously will cost four month for even only chromosome 1. Whereas for each single population, process all chromosomes will cost at most a month. 
 (Our test vcf file were produced by GATK 4.5 HaplotypeCaller on Intel(R) Xeon(R) CPU E5-2690 v4 @ 2.60GHz machine with 56 CPU(s) 256G memory)
 
 	python atcgtools.py VJtransf -v gatkallchromAll147indvd.vcf -o gatkallchromAll147indvd -f pedmap [-r 0.5 -t 8 -c /home/lrui/avian/duckpop/configfiles/chrnameMap.txt]
 	# above command will transfer vcf into gatkallchromAll147indvd.map and gatkallchromAll147indvd.ped file. With -r 0.5, the out .map and .ped will only contain sites that filtered with r2=0.5 by plink and .prune.in for each chromosome.
 	# then the subsequent cammands could use the .prune.in file. e.g. transfer those r2 filtered sites to .geno format.  
 	python atcgtools.py VJtransf -v gatkallchromAll147indvd.vcf -o gatkallchromAll147indvd.prune0.5 -t 8 -f genosnp -c chrnameMap.txt -i gatkallchromAll147indvd.0.5prune.in -m all8popfile.indpop 3

 	# command 1
 	# to avoid the huge time consuming of joint calling for all seven populations' samples, one can use the following command to align all seven vcf files sites by sites to produce the almost same results with gatkallchromAll147indvd.vcf.
 	python atcgtools.py VJtransf -v mpl.VCFBAMconfig  -v newdomesticbreeds.VCFBAMconfig -v shaoxing.VCFBAMconfig  -v mallard.VCFBAMconfig -v spotbilled.VCFBAMconfig -v beijing.VCFBAMconfig -v ZJU1.VCFBAMconfig [-C NC_051772.1] -o align7vcffileNC_051772withZJU1 -t 8 -f pedmap [-c chrnameMap.txt][-r ...]
	# VCFBAMconfig file's first line 'vcffile=xxxx' specify the vcf file used to align with others. The sam/bam files used to calling this vcf should list in the following lines. Those bam files will be used to judge whether an absence in a vcf file is due to it fixed as reference or no coverage for the aligned site. 
	Different suffix names .vcf or .VCFBAMconfig enable different behaviours correspondingly. when multiple -v is used, should use .VCFBAMconfig as input as it is nesscery to align records.
	
 -r optional. Conflict with -i. If provided, the output will filter out those linked sites by invoking plink --indep-pairwise 100(default) 20(default) 0.5
 
 -i optional. Conflict with -r.
 
 -m group information file. 
 	
  	# in above example, all8popfile.indpop file contain three essential columns and could with additional columns to group differently. the second argument of -m specifies  which column to use as the group info.
  	#otherspecies_commonteal_fh_u5_56        U       otherspecies    otherspecies
 	#otherspecies_commonteal_pylake_u2_53    U       otherspecies    otherspecies
 	#otherspecies_falcatedteal_fh_u4_55      U       otherspecies    otherspecies
 	#otherspecies_gadwall_fh_f7_7    F       otherspecies    otherspecies
 	#shaoxing_SRR6323890     U       shaoxing        domestic
 	#shaoxing_SRR6323918     U       shaoxing        domestic
 	#beijing_41224_42        U       beijing domestic
	#beijing_41288_44        U       beijing domestic
 -t optional. Specifies how many chromosomes will be processed simultaneously using multiple processes. Better assign as many as possible to get results faster.
 
 -c optional. if provided, the out put's chromosome ID will be replace by the corresponding name according this file.
 
 -o output file prefix name
 	
 VJgendadi: similar to VJtransf, but requiring a special designed mysql table with flank sequences and outgroup genotype info of each SNPs.
( Most functions below require a special designed mysql/mariandb table that cantains all SNPs of the species and the state of the polymorphism in each populations/breeds for each SNPs, 
  and the outgroup genotype information to determine the ancestral allele. This can be done by the following two steps:
 atcgtools inittopleveltable allindividuals.vcf
 atcgtools fillcontextNoutgroup 
 )
 
 inittopleveltable:
 
 Create the toplevel variants table for a species and load SNP/INDEL records of allindividuals.vcf into the table.
 After this step, toplevel talbe in mysql databases contain the following fields:chrID,snp_pos,snpID,ref_allele,context,alt_alle1,alt_alleN,presenceflag
 ref_allele corresponding the 4th(REF) col of the input vcf file and context contain the flank bases of both side of the SNP, empty for INDEL.
 alt_alle1 conrresponding to the first allele of the 5th(ALT) of the input vcf file and the rest alternated allele will be stored in alt_alleN col.
 presenceflag, which is paticularly noteworthy. it's a two bits flag. from right(least significant bit) to left, every to bits represent the variant state of the corresponding population,which is assigned in config.propoties file.
 	
 As the SNPs will vary when different populations/individuals used to call, the input vcf here better use as comprehensive an individual/population set as possible. 
 Whereas joint too many individuals/populations to call SNP may supper time consuming, So our solutions consists of two aspects:
 1) temporarily, using a vcf file contains as comprehensive SNPs this species have, then in following analyses programs,which using different vcf files including additional variant records, could add those new records to the toplevel table and fill their outgroup info(as next step).
 2) 
 after the two steps, the presenseflag still empty. It will be set when the assigned vcf files to be used in analyses programs through instance of class alterPresenceflag_Callback, alterDB() function
 presenseflag will be set (through alterDB() function) for all SNPs for corresponding populations when alignmultPopSnpPos() function is called. For example, the after above command #1 or VJgendadi, this presenseflag will be filled as instructure of config.properties file. as the commands all called alignmultPopSnpPos()

 fillcontextNoutgroup
 
 fillcontextNoutgroup


 Detectsignalacrossgenome
 subclasses Caculate_popPI,Caculate_popDiv,Caculate_Fst,Caculate_ABB_BAB_BBAA,Caculate_df,Caculate_Hp_master_slave,Caculate_S_ObsExp_difference were used to
 calculate pi,dxy,Fst,abbabaab,df,Hp,SDS that selection occured in ancient, respectively.
	
	e.g. atcgtools Detectsignalacrossgenome -T spotbilled.VCFBAMconfig -T configfiles/mallardZJU1.VCFBAMconfig -o EarlyseletedRegion_auto -R configfiles/beijing.VCFBAMconfig -R configfiles/shaoxing.VCFBAMconfig -R configfiles/newdomesticbreeds.VCFBAMconfig -R configfiles/gy.VCFBAMconfig -R configfiles/sm.VCFBAMconfig -R configfiles/jd.VCFBAMconfig -R configfiles/cv.VCFBAMconfig -R configfiles/campbell.VCFBAMconfig -n 24 -p early -t toplevelDuck_ZJU1ref -w 20000 -s 10000 -1 build/bioinfodevelop/slave/Detectsignalacrossgenome_producecorrelation_slave.py -2 build/bioinfodevelop/slave/Detectsignalacrossgenome_slidewin_slave.py -c auto

 VCataAnno
	
	e.g. 
	atcgtools VCataAnno -V chr1.vcf -r ZJU1/bjduckallchr.fna -g ZJU1/genomic.gtf -o duckpop/variantanno/ -c ZJU1duckchrominfo -m 80000 -5 3000 -3 3000
	
	extremely slow so far
