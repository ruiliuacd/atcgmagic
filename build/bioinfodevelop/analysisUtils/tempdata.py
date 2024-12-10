'''
Created on 2024年10月28日

@author: pc
'''
import pysam
fakepoolVcfmap={"/home/lrui/avian/duckpop/domesticpool/cvbellpool.vcf":['vcfobj',pysam.AlignmentFile('/home/lrui/avian/duckpop/domesticpool/cvpool/cvpool_I2481249I238.sort.uniqmap.dedup.bam','rb')],
            "/home/lrui/avian/duckpop/domesticpool/jdpool.vcf":['vcfobj',pysam.AlignmentFile('/home/lrui/avian/duckpop/domesticpool/jdpool/jdpool_1249I248.sort.uniqmap.dedup.bam','rb')],
            "/home/lrui/avian/duckpop/domesticpool/campbellpool.vcf":['vcfobj',pysam.AlignmentFile('/home/lrui/avian/duckpop/domesticpool/campbellpool/campbellpool_1249I248.sort.uniqmap.dedup.bam','rb')],
            "/home/lrui/avian/duckpop/domesticpool/smpool.vcf":['vcfobj',pysam.AlignmentFile('/home/lrui/avian/duckpop/domesticpool/smpool/smpool_1248I249.sort.uniqmap.dedup.bam','rb')],
            "/home/lrui/avian/duckpop/domesticpool/gypool.vcf":['vcfobj',pysam.AlignmentFile('/home/lrui/avian/duckpop/domesticpool/gypool/gypool_1249I248.sort.uniqmap.dedup.bam','rb')],
            '/srv/lruiwork/sramplouzr8.indvd.vcf':[""],
            '/srv/lruiwork/newdomesticbreedslcwhitelcwhitelcwhite23.indvd.vcf':[""],
            '/srv/lruiwork/shaoxing33.indvd.vcf':[""],
            '/srv/lruiwork/mallard29.indvd.vcf':[""],
            '/srv/lruiwork/spotbilled21.indvd.vcf':[""],
            '/srv/lruiwork/beijing33.indvd.vcf':[""]
            }