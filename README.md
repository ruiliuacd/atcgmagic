# atcgtoolkit

export PYTHONPATH=/path/to/atcgtoolkit/build:$PYTHONPATH
#export PYTHON=/path/to/atcgtoolkit/build/bioinfodevelop:$PYTHONPATH
chmod +x atcgtools

atcgtools [analysistype] -h
or 
python atcgtools.py [analysistype] -h

run different module functions:

atcgtools [analysistype] [parameters for different analysistype]

[analysistype]={countgaps,inittopleveltable,fillcontextNoutgroup,VJtransf,VJgendadi,Detectsignalacrossgenome,VCataAnno,test}
 VJtransf：
 
 Convert vcf into which format depands on -f. 
 -f=	
 pedmap (default): produce .map and .ped file of plink format (if -r is proviede, final map ped file only contian unlinked sites that r2 less than value -r assigned )
 		pedmap1 produce loose format plink format.	
 genosnp: output .geno file also called EIGENSTRST format, .snp and .ind file. 
 	 
      .geno file contains 1 line per position. Each line contains 1 character per individual:
 	 	0 means zero copies of reference allele
 	 	1 means one copy of reference allele.
 	 	2 means two copies of reference allele.
 	 	9 missing data.
 vcf: output vcf format. 
 
 -v vcf file [-v ...].  When multiple -v vcf files were provided, then the output data will contain inner/outer joined records of all vcf as described in [1](AI 5.8). It provided a flexible way to joint calling SNP and could avoid the significant time consumption of calling SNPs for all individuals simultaneously using GATK.
 e.g. The file 'gatkallchromAll147indvd.vcf' was generated by GATK for seperate chromosomes, using all individuals from all populations. Joint calling for all seven populations of 147 ducks's sample simultaneously will cost more than four month for even only chromosome 1. Whereas for each single population, process all chromosomes will cost at most a month and take hours to joint all vcfs into one, with distinguishing no coverage or fixed as reference allele for the absense recoreds of the populations, which could instead replace GATK joint calling.
 (Our test vcf file were produced by GATK 4.5 HaplotypeCaller on Intel(R) Xeon(R) CPU E5-2690 v4 @ 2.60GHz machine with 56 CPU(s) 256G memory)
  	 
  	 
  	 python atcgtools.py VJtransf -v gatkallchromAll147indvd.vcf -o gatkallchromAll147indvd -f pedmap [-r 0.5 -t 8 -c configfiles/chrnameMap.txt]
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
 
 -m group information file. necessary when -f genosnp
 	
    
 	# in above example, all8popfile.indpop file contain three essential columns and could with additional columns to group differently. the second argument of -m specifies  which column to use as the group info.
 	otherspecies_commonteal_fh_u5_56        U       otherspecies    otherspecies
 	otherspecies_commonteal_pylake_u2_53    U       otherspecies    otherspecies
 	otherspecies_falcatedteal_fh_u4_55      U       otherspecies    otherspecies
 	otherspecies_gadwall_fh_f7_7    F       otherspecies    otherspecies
 	shaoxing_SRR6323890     U       shaoxing        domestic
 	shaoxing_SRR6323918     U       shaoxing        domestic
 	beijing_41224_42        U       beijing domestic
 	beijing_41288_44        U       beijing domestic
    
 -t optional. Specifies how many chromosomes will be processed simultaneously using multiple processes. Better assign as many as possible to get results faster.
 
 -c optional. if provided, the out put's chromosome ID will be replace by the corresponding name according this file.
 
 -C specific chromosomes IDs. optional. [-C ...] multiple -C for multiple chromosomes IDs respectively. 
 
 -o output file prefix name
 	
 VJgendadi: similar to VJtransf, but requiring a special designed mysql table with flank sequences and outgroup genotype info of each SNPs.
( Most functions below require a special designed mysql/mariandb table that cantains all SNPs of the species and the state of the polymorphism in each populations/breeds for each SNPs, 
  and the outgroup genotype information to determine the ancestral allele. This can be done by the following two steps:
 atcgtools inittopleveltable allindividuals.vcf
 atcgtools fillcontextNoutgroup 
 )
 
 inittopleveltable:
 
 Create the toplevel variants table for a species and load SNP/INDEL records of allindividuals.vcf into the table.
 After this step, toplevel talbe in mysql databases contain the following fields:
 chrID,snp_pos,snpID,ref_allele,context,alt_alle1,alt_alleN,presenceflag
 'ref_allele' corresponding the 4th(REF) col of the input vcf file and 'context' contain the flank bases of both side of the SNP, empty for INDEL.
 'alt_alle1' conrresponding to the first allele of the 5th(ALT) of the input vcf file and the rest alternated allele will be stored in 'alt_alleN' col.
 'presenceflag' is paticularly noteworthy. It's a two bits flag, from right(least significant bit) to left, every two bits represent the variant state of the corresponding population,which is assigned in config.propoties file.
 	
 	
 	00 the SNP site is not covered in this population samples
 	01 fixed as reference allele
 	10 fixed as alternate allele
 	11 SNP
    
 As the SNPs will vary when different populations/individuals used to call, the input vcf here better use as comprehensive an individual/population set as possible. Whereas joint calling on too many individuals/populations can be extremely time consuming, so our solutions consists of two aspects:
 1) temporarily, using a vcf ,e.g. one population vcf file, file contains as comprehensive SNPs as possible. Then in following analyses programs,which using different vcf files including additional variant records, could add new records to the toplevel table and fill their outgroup info(as next step).
 2) run any type of anaysis that would joint multiple vcf files by invoking alignmultPopSnpPos(...,jointmode="i"/"o"/"l",...) function inside.
 	This function joint the variant records from all input VCF files in three modes according 'jointmode' parameter:
	
	1. i (inner join): only retain the variants records that all combined populations have and have the same alleles.

	2. o (outer join): variants that present in any one of the combined populations will be retained, and program will test the origin BAM file whether this site isn’t covered or how many reads support it as a reference allele. And the record the depth or ‘no cover’ into top-level table.

	3. l (left join): as long as the ‘selected population’, referred as ‘ref-pop’, has the variants, then retain this variants and if any of the other populations’ VCF don’t have this records then search BAM file for depth information as 2.
	
 alignmultPopSnpPos() functon return the joint records of all input vcfs and can be output in any format (vcf/pedmap/genosnp), and will automaticlly update the toplevel variants table filling new records and update presenceflag, context, outgroup info and all fields.
 The 'presenceflag' will be set (through alterDB() function) for all SNPs for corresponding populations when alignmultPopSnpPos() function is called. For example, the after above command #1 or VJgendadi, this presenseflag will be filled as instructure of config.properties file. And the 'presenceflag' for a population/group would be re-set when some analysis include new vcf files of this population/group: e.g.
 	
    
 	old presenceflag | presenceflag for this vcf
 	00 | 00	=	00
 	00 | 01	=	01
 	01 | 01	=	01
 	10 | 01	=	11
 	11 | 01	=	11
    
 One can also leave the 'presenseflag' unchanged, by seting 'callback' parameter to be None. Otherwise use a instance of class alterPresenceflag_Callback to set this mysql table field. Use this flag, one can flexibly count some statistics.


 fillcontextNoutgroup:
 
 With '-r Reference_genome_file.fa flanksequence_length' provided to fill 'context' column of toplevel mysql table by extracting flank seqs from Reference_genome_file.fa.
 
 With '-a xxx.VCFBAMconfig [[-a yyy.VCFBAMconfig]...]' provided to add columns correspondingly in toplevel variants table and fill the depth information for both ref and alt alleles that could serve as archic/outgroup populations to provide ancestral allele information.
 
 If both -r and -a is provided, do the both. each -a recive one VCFBAMconfig file as value and write/append this file name to config.properties file 'outgroupvcfbamconfig_zju1ref=outgroup1.VCFBAMconfig;outgroup1.VCFBAMconfig;...'. One can also manully write this file and wait for the program automaticlly update the corresponding column info of toplevel table in subsequent analyses process that would invoke dynamicInsertUpdateAncestralContext() function. e.g. VJgendadi, Detectsignalacrossgenome -p early that would use the infromation of ancestral allele.
 
 -C chromlistfilename. Required. Provide a file with chromosomes name and range/length each line.


 Detectsignalacrossgenome:
 
 -p assign which signal to calculate in sliding window, the corresponding calculation program were implemented in subclasses Calculate_popPI,Calculate_popDiv,Calculate_Fst,Calculate_ABB_BAB_BBAA,Calculate_df,Calculate_Hp_master_slave,Calculate_S_ObsExp_difference.
 Different'subclass' instances of Calculator were implemented to calculate statistics:pi, dxy/ Fst/ABBABAAB/  df(fixed different)/ Hp/ seletion in ,for example, ancient wild population after divergence from domestic lineage(SDS), respectively.

	
  	 e.g. atcgtools Detectsignalacrossgenome -T spotbilled.VCFBAMconfig -T configfiles/mallardZJU1.VCFBAMconfig -o EarlyseletedRegion_auto -R configfiles/beijing.VCFBAMconfig -R configfiles/shaoxing.VCFBAMconfig -R configfiles/newdomesticbreeds.VCFBAMconfig -R configfiles/gy.VCFBAMconfig -R configfiles/sm.VCFBAMconfig -R configfiles/jd.VCFBAMconfig -R configfiles/cv.VCFBAMconfig -R configfiles/campbell.VCFBAMconfig -n 24 -p early -t toplevelDuck_ZJU1ref -w 20000 -s 10000 -1 build/bioinfodevelop/slave/Detectsignalacrossgenome_producecorrelation_slave.py -2 build/bioinfodevelop/slave/Detectsignalacrossgenome_slidewin_slave.py -c auto

 VCataAnno:
 	
 Along with extracting protein-coding genes from the genome.fa and translating them into amino acid sequences according the GTF file, the diallele variants will be annotated based on their locations in CDS/UTR/intron/intergenic region respectively and output into different files.
 Our program specially take overlaped genes into consideration. e.g,
	
    
  	 atcgtools VCataAnno -V wildchr1.vcf -r ZJU1/bjduckallchr.fna -g ZJU1/genomic.gtf -o duckpop/variantanno/ -c ZJU1duckchrominfo -m 80000 -5 3000 -3 3000
	
  	 Above command output .cds/.intron/.utr/.intergenic/.mutcds/.mutaa/.refaa outfiles
  	 The .cds outfile contain those variants that locate in CDS regions, when overlaped transcripts exist, the output will be like this:
  	 #CHROM  POS     REF     ALT     trscptID        geneName  strand  cdsidx  refcodon        refaa   altcodon        altaa   INFO    FORMAT  mallard_SRR6323906      mallard_SRR6323939      mallard_SRR7091422   ... [all samples same as vcf's samples columns]
  	 1	26118287	A	G	XM_038179863.1;XM_038179877.1;XM_038179871.1;XM_038179885.1;XM_038179867.1;XM_038179893.1	CFTR;CFTR;CFTR;CFTR;CFTR;CFTR	-;-;-;-;-;-	5;3;3;3;3;2	ttg;ttg;ttg;ttg;ttg;ttg	L;L;L;L;L;L	ctg;ctg;ctg;ctg;ctg;ctg	L;L;L;L;L;L	AC=1;AF=0.017;AN=58;BaseQRankSum=-0.363;DP=401;ExcessHet=0.0000;FS=0.000;InbreedingCoeff=-0.0241;MLEAC=1;MLEAF=0.017;MQ=60.00;MQRankSum=0.000;QD=18.56;ReadPosRankSum=-1.880;SOR=0.507	GT:AD:DP:GQ:PL	['0/0:10,0:10:30:0,30,450', '0/0:6,0:6:18:0,18,260',  omit...]
  	 1	26118332	G	A	XM_038179863.1;XM_038179877.1;XM_038179871.1;XM_038179885.1;XM_038179867.1;XM_038179893.1	CFTR;CFTR;CFTR;CFTR;CFTR;CFTR	-;-;-;-;-;-	5;3;3;3;3;2	cta;cta;cta;cta;cta;cta	L;L;L;L;L;L	tta;tta;tta;tta;tta;tta	L;L;L;L;L;L	AC=38;AF=0.655;AN=58;BaseQRankSum=-0.369;DP=402;ExcessHet=2.5144;FS=11.193;InbreedingCoeff=-0.0625;MLEAC=38;MLEAF=0.655;MQ=60.00;MQRankSum=0.000;QD=33.73;ReadPosRankSum=-0.681;SOR=0.605	GT:AD:DP:GQ:PL	['1/1:0,13:13:39:570,39,0', '1/1:0,7:7:21:260,21,0', omit...]
  	 1	26119325	G	A	XM_038179863.1;XM_038179877.1;XM_038179871.1;XM_038179885.1;XM_038179867.1	CFTR;CFTR;CFTR;CFTR;CFTR	-;-;-;-;-	4;2;2;2;2	tgc;tgc;tgc;tgc;tgc	C;C;C;C;C	tgt;tgt;tgt;tgt;tgt	C;C;C;C;C	AC=1;AF=0.017;AN=58;BaseQRankSum=-3.100;DP=362;ExcessHet=0.0000;FS=0.000;InbreedingCoeff=-0.0242;MLEAC=1;MLEAF=0.017;MQ=60.00;MQRankSum=0.000;QD=10.66;ReadPosRankSum=-0.617;SOR=0.609	GT:AD:DP:GQ:PL	['0/0:12,0:12:36:0,36,522', '0/0:16,0:16:48:0,48,625', omit...]
  	 1	26119376	T	C	XM_038179863.1;XM_038179877.1;XM_038179871.1;XM_038179885.1;XM_038179867.1	CFTR;CFTR;CFTR;CFTR;CFTR	-;-;-;-;-	4;2;2;2;2	cca;cca;cca;cca;cca	P;P;P;P;P	ccg;ccg;ccg;ccg;ccg	P;P;P;P;P	AC=16;AF=0.276;AN=58;BaseQRankSum=-0.028;DP=370;ExcessHet=5.3165;FS=1.853;InbreedingCoeff=-0.2133;MLEAC=16;MLEAF=0.276;MQ=60.00;MQRankSum=0.000;QD=19.44;ReadPosRankSum=-0.307;SOR=0.802	GT:AD:DP:GQ:PL	['0/0:11,0:11:33:0,33,465', '1/1:0,12:12:36:497,36,0', omit...]
  	 1	26123590	C	T	XM_038179867.1	CFTR	-	1	cga	R	caa	Q	AC=5;AF=0.086;AN=58;BaseQRankSum=0.758;DP=433;ExcessHet=0.8126;FS=2.383;InbreedingCoeff=-0.0979;MLEAC=5;MLEAF=0.086;MQ=60.00;MQRankSum=0.000;QD=14.99;ReadPosRankSum=-0.877;SOR=0.939	GT:AD:DP:GQ:PL	['0/0:8,0:8:24:0,24,352', '0/0:8,0:8:24:0,24,345', omit...]
  	 1	26123614	G	C	XM_038179867.1	CFTR	-	1	gct	A	ggt	G	AC=2;AF=0.034;AN=58;BaseQRankSum=-1.393;DP=446;ExcessHet=0.0769;FS=0.000;InbreedingCoeff=-0.0372;MLEAC=2;MLEAF=0.034;MQ=60.00;MQRankSum=0.000;QD=17.10;ReadPosRankSum=-0.209;SOR=0.627	GT:AD:DP:GQ:PL	['0/0:7,0:7:21:0,21,294', '0/0:9,0:9:27:0,27,392', omit...]
  	 1	26125427	A	G	XM_038179863.1	CFTR	-	3	ctt	L	ctc	L	AC=13;AF=0.224;AN=58;BaseQRankSum=-1.121;DP=400;ExcessHet=0.0586;FS=2.187;InbreedingCoeff=0.2970;MLEAC=13;MLEAF=0.224;MQ=60.00;MQRankSum=0.000;QD=24.50;ReadPosRankSum=1.311;SOR=0.551	GT:AD:DP:GQ:PL	['0/0:12,0:12:36:0,36,509', '0/1:4,3:7:73:73,0,155', omit...]
  	 1	26125430	A	G	XM_038179863.1	CFTR	-	3	gct	A	gcc	A	AC=24;AF=0.414;AN=58;BaseQRankSum=0.313;DP=406;ExcessHet=4.0053;FS=0.392;InbreedingCoeff=-0.1088;MLEAC=24;MLEAF=0.414;MQ=60.00;MQRankSum=0.000;QD=22.46;ReadPosRankSum=0.658;SOR=0.661	GT:AD:DP:GQ:PL	['0/1:4,8:12:99:210,0,107', '0/1:3,4:7:99:137,0,114', omit...]
  	 1	26126862	T	C	XM_038179863.1	CFTR	-	2	cga	R	cgg	R	AC=1;AF=0.017;AN=58;BaseQRankSum=3.083;DP=403;ExcessHet=0.0000;FS=9.064;InbreedingCoeff=-0.0192;MLEAC=1;MLEAF=0.017;MQ=60.00;MQRankSum=0.000;QD=10.85;ReadPosRankSum=-0.709;SOR=3.158	GT:AD:DP:GQ:PL	['0/0:12,0:12:36:0,36,517', '0/0:9,0:9:27:0,27,332', omit...]
  	 1	26159830	G	A	XM_027457737.2;XM_038179909.1;XM_038179905.1	ASZ1;ASZ1;ASZ1	+;+;+	1;1;1	gcg;gcg;gcg	A;A;A	gca;gca;gca	A;A;A	AC=4;AF=0.087;AN=46;BaseQRankSum=-1.278;DP=143;ExcessHet=0.0067;FS=0.000;InbreedingCoeff=0.2720;MLEAC=4;MLEAF=0.087;MQ=60.00;MQRankSum=0.000;QD=13.54;ReadPosRankSum=-0.269;SOR=0.496	GT:AD:DP:GQ:PL	['0/0:5,0:5:15:0,15,225', '0/0:6,0:6:18:0,18,270', omit...]
  	 1	26159857	T	C	XM_027457737.2;XM_038179909.1;XM_038179905.1	ASZ1;ASZ1;ASZ1	+;+;+	1;1;1	gat;gat;gat	D;D;D	gac;gac;gac	D;D;D	AC=19;AF=0.365;AN=52;BaseQRankSum=1.580;DP=165;ExcessHet=0.7722;FS=1.019;InbreedingCoeff=-0.0781;MLEAC=19;MLEAF=0.365;MQ=60.00;MQRankSum=0.000;QD=17.78;ReadPosRankSum=1.089;SOR=0.947	GT:AD:DP:GQ:PL	['0/1:2,4:6:72:162,0,72', '0/0:8,0:8:24:0,24,303', omit...]
  	 1	26159878	C	T	XM_027457737.2;XM_038179909.1;XM_038179905.1	ASZ1;ASZ1;ASZ1	+;+;+	1;1;1	ggc;ggc;ggc	G;G;G	ggt;ggt;ggt	G;G;G	AC=9;AF=0.173;AN=52;BaseQRankSum=-0.265;DP=182;ExcessHet=0.6065;FS=1.680;InbreedingCoeff=0.0620;MLEAC=10;MLEAF=0.192;MQ=60.00;MQRankSum=0.000;QD=18.38;ReadPosRankSum=1.340;SOR=0.366	GT:AD:DP:GQ:PL	['0/1:2,4:6:72:115,0,72', '0/0:6,0:6:18:0,18,244', omit...]
	
	and similar for .intron outfile
	
	
-c assign a mysql table (ZJU1duckchrominfo in example) containning chrID and chrlength in first two coloumn.

-m requires the at least length for chromosomes to be annotated.

-5/-3 define how long the extending nearby regions output the variants to .utr outfile. When a transcript without UTR ahead/after the first/end codon (not necessary to be M/*) in GTF file, collect variants between -5/-3 distance to the first/end codon as 5UTR or 3UTR, respectively, output to .utr file. Otherwise, any side of UTR exist, the corresponding -5/-3 parameter defined UTR is not used.
	the inter-regions between interrupted UTR will also output to .utr outfile with different tag.
	
It's not optimized to speed up so far, So extremely slow for the whole genome. Use regioned vcf file for efficiency.

# Additional:
Basically, python analysisAppEntry/toolkit/calculateLDUsePlinkForScaffold.py repeatly extract vcf records for every chromosome and transfer into .map/.ped file and then revoke plink to use all SNPs in the file to produce a plink_part[i].ld file.
e.g. plink_part1.ld plink_part2.ld ...... plink_part7897.ld
then one can use 
binld.py -d plink_partmyNtosub.ld 2 5 1 7897 -i 0 500000 100 -m 7 -o outnamepre
to repeatly collect above plink_part1.ld ~ plink_part7897.ld one by one into bins that assigned by -i, according to value of column 2 - column 5 . The 7th column's value for each bin will be averaged. 

casualmutationdetil.py was used to select SNPs located in bed (candidate selected regions producing from findTrscpt.py/testmakemht_hist.py) and sort by delt allele frequencies (AF) between two groups of vcf files that represent two populations, e.g. wilds and domestics.
CatalogSNPincandidateRegion.py, classify SNPs in bed file by cds,utr,intron. And blast conserved sequence to identify it's range in reference genome and extrate contained SNPs, then sort by delt AF.

python ./build/bioinfodevelop/analysisAppEntry/toolkit/testmakemht_hist.py -n ancientselectionsignal_auto_sp_ma_be_sh_ne_gy_sm_jd_cv_ca.earlypostiveselected.zscorefile20000_1000088 -0.3_picturetitle outfilename -X winvalue -u 1000 -d 1000 -w 20000 -s 10000 -N 1 -o multipleMhtInOne_outfilenameprefix
to run testmakemht_hist.py to drawpicture, you may need to install rpy2 and some R packages to ./com./Rpackages first,  
- e.g.
R CMD INSTALL  -l ./com/Rpackages/ ./com/Rpackages/gap_1.1-3.tar.gz
textshaping, ragg
and adjuest r('.libPaths("...../com/Rpackages")') acorrdingly.

When -g ensBiomartGO.table is provided, it will output the GO enrichments results using ensemble GO information  according their hypergeometric distribution.
ensBiomartGO.table file can be download from ensemble biomart with the following columns:
Ensembl Gene ID Ensembl Transcript ID   EntrezGene ID   UniProt/TrEMBL Accession        GO Term Accession       GO Term Evidence Code   GO domain       Associated Gene Name    Associated Transcript Name	GO Term Name	...

Note: As our work place greater emphasis on enhencing the continuity of user thought, by providing engineering support that aligns with theories and methodologies, rather than focus on improving the speed of program execution, our program may unable run as fast as some other tools that finish some similar funtions in minutes currently. But the hours-level speed for those single steps are also tolerable and could be easy update into high speed as minutes level.
ruiliugenetic@nwu.edu.cn
cite: https://doi.org/10.1101/2020.02.03.933069
Our this article last time was reject due to a reviewer claims that our findings have been reported in somewhere else for one reason, which is not true and they published later than ours preprint version. The improved version can be obtained via email, and I am available to work with anyone interested in the associated technology. The source code for this project, which without big innovation and high skill, was intended to be open and once opened and cloned by many. However, due to some misconduct that damaged the collaborative integrity environment and our interests, access to the source code is now restricted to user who we know.