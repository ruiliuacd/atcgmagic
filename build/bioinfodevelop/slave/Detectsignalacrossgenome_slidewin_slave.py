'''
Created on 2015-8-21

@author: liurui
'''
from optparse import OptionParser
import re,  fractions, copy, os, pysam,config
from bioinfodevelop.analysisUtils import Calculators,variantUtils,Utils
import bioinfodevelop.analysisUtils.DBManager as dbm


parser = OptionParser()
parser.add_option("-c", "--chromlistfilename", dest="chromlistfilename", help="early,pairfst,pbs,lsbl,is")
parser.add_option("-p", "--typeOfcalculate", dest="typeOfcalculate")
parser.add_option("-t", "--topleveltablejudgeancestral", dest="topleveltablejudgeancestral", help="assigned only if -p early")
parser.add_option("-T", "--targetpopvcfconfig", dest="targetpopvcfconfig", action="append", help="treat as P1")
parser.add_option("-U", "--P2popvcfconfig", dest="P2popvcfconfig", action="append", help="treat as P2")
parser.add_option("-V", "--P3popvcfconfig", dest="P3popvcfconfig", action="append", help="treat as P3")
parser.add_option("-R", "--refpopvcffileconfig", dest="refpopvcffileconfig", action="append", help="treat as O")
parser.add_option("-w", "--winwidth", dest="winwidth", help="default infile1_infile2")  #
parser.add_option("-s", "--slideSize", dest="slideSize", help="default infile2_infile1")  #
parser.add_option("-C", "--correlationfile", dest="correlationfile", default=None, help="conflit with numberofindvdoftargetpop_todividintobin")
parser.add_option("-o", "--outfileprewithpath", dest="outfileprewithpath")
parser.add_option("-m", "--masterpid", dest="masterpid")
parser.add_option("-b", "--bedlikefile", dest="bedlikefile", help="conflict with -c ")
(options, args) = parser.parse_args()
mindeptojudgefix = 20  # for pool only
minSNPs=1
extendsize = 500000
windowWidth = int(options.winwidth)
slideSize = int(options.slideSize)
if __name__ == '__main__':
    print("runSlave_slidewin process ID", os.getpid(), "start")
    chromlistOrBedRegionList = []
    if options.chromlistfilename != None and options.bedlikefile == None:
        chromlistfile = open(options.chromlistfilename, "r")
        
        for rec in chromlistfile:
            reclist = re.split(r'\s+', rec.strip())
            chromlistOrBedRegionList.append((reclist[0].strip(), int(reclist[1].strip())))
        chrlistfilewithoutpath = re.search(r"[^/]*$", options.chromlistfilename).group(0)
        chromlistfile.close()
    elif options.chromlistfilename == None and options.bedlikefile != None:
        bedreclistfile = open(options.bedlikefile, "r")
        for rec in bedreclistfile:
            reclist = re.split(r"\s+", rec.strip())
            chromlistOrBedRegionList.append((reclist[0].strip(), (int(reclist[1]), int(reclist[2]))))
        chrlistfilewithoutpath = re.search(r"[^/]*$", options.bedlikefile).group(0)
    else:
        print("options.chromlistfilename and options.bedlikefile conflict")
#     genomedbtools = dbm.DBTools(Util.ip, Util.username, Util.password, Util.genomeinfodbname)
    print(chromlistOrBedRegionList)
    if options.typeOfcalculate == "early":
        if  options.chromlistfilename != None and options.bedlikefile == None:
            aaaaaaa = options.chromlistfilename
        else:
            aaaaaaa = options.bedlikefile
        flankseqfafilename = aaaaaaa + str(os.getpid()) + "snpflankseq.fa"
        dbvariantstools = dbm.DBTools(config.ip, config.username, config.password, config.variantsdbname)
        dynamicIU_toptable_obj = variantUtils.dynamicInsertUpdateAncestralContext(dbvariantstools, config.beijingreffa, options.topleveltablejudgeancestral)
        obsexpcaculator = Calculators.Calculate_S_ObsExp_difference(mindeptojudgefix, options.targetpopvcfconfig, options.refpopvcffileconfig, dbvariantstools, options.topleveltablejudgeancestral, options.outfileprewithpath)
        obsexpcaculator.dynamicIU_toptable_obj = dynamicIU_toptable_obj
        obsexpcaculator.flankseqfafile = open(flankseqfafilename, "w")
        plainname = re.search(r"[^/]*$", obsexpcaculator.outputname).group(0)
        if len(plainname) >= 228:
            obsexpcaculator.outputname = obsexpcaculator.outputname[:-(len(plainname) - 228)]
        outputname = obsexpcaculator.outputname
        outfile = open(outputname + "." + options.typeOfcalculate + str(windowWidth) + "_" + str(slideSize) + str(os.getpid()) + chrlistfilewithoutpath, 'w')
        print("chrNo\twinNo\tfirstsnppos\tlastsnppos\tnoofsnp\twinvalue\tzvalue", file=outfile)
        freq_correlation_configFileName = options.outfileprewithpath + ".freq_correlation_merged" 
        if options.correlationfile != freq_correlation_configFileName:
            print("warning !", freq_correlation_configFileName, " is not equal to ", options.correlationfile)
        freq_correlation_config = open(options.correlationfile, "r")
        final_freq_xaxisKEY_yaxisVALUERelation = {}
        for line in freq_correlation_config:
            if line.split():
                linelist = re.split(r"\t", line.strip())
                a = float(linelist[0]);b = float(linelist[1]);yaxisfreq = float(linelist[2])
                final_freq_xaxisKEY_yaxisVALUERelation[(a, b)] = yaxisfreq
        obsexpcaculator.freq_xaxisKEY_yaxisVALUERelation = final_freq_xaxisKEY_yaxisVALUERelation
        freq_correlation_config.close()
    elif options.typeOfcalculate=="D":
        obsexpcaculator=Calculators.Calculate_ABB_BAB_BBAA("no",options.targetpopvcfconfig,options.P2popvcfconfig,options.P3popvcfconfig,options.refpopvcffileconfig,options.outfileprewithpath,os.getpid(),30)    
        outputname=options.outfileprewithpath
        outfile=open(outputname+"."+options.typeOfcalculate + str(windowWidth) + "_" + str(slideSize) + str(os.getpid()) + chrlistfilewithoutpath,"w")
        print("chrNo\twinNo\tfirstsnppos\tlastsnppos\tnoofsnp\twinvalue\tzvalue", file=outfile)
    elif options.typeOfcalculate == "pairfst" or options.typeOfcalculate == "dxy":
#         obsexpcaculator = Caculators.Caculate_pairFst(mindeptojudgefix, options.targetpopvcfconfig, options.refpopvcffileconfig)
        obsexpcaculator = Calculators.Calculate_popDiv("no",options.targetpopvcfconfig, options.refpopvcffileconfig,options.outfileprewithpath)
        plainname = re.search(r"[^/]*$", obsexpcaculator.outputname).group(0)
        if len(plainname) >= 228:
            outputname = obsexpcaculator.outputname[:-(len(plainname) - 228)]
        outputname = obsexpcaculator.outputname
        outfile = open(outputname + "." + options.typeOfcalculate + str(windowWidth) + "_" + str(slideSize) + str(os.getpid()) + chrlistfilewithoutpath, 'w')
        print("chrNo\twinNo\tfirstsnppos\tlastsnppos\tnoofsnp\twinvalue\tzvalue", file=outfile)
    
    elif options.typeOfcalculate == "is":
        obsexpcaculator = Calculators.Calculate_IS(mindeptojudgefix / 2, options.targetpopvcfconfig, options.refpopvcffileconfig)
        plainname = re.search(r"[^/]*$", options.outfileprewithpath).group(0)
        if len(plainname) >= 228:
            outputname = options.outfileprewithpath[:-(len(plainname) - 228)]
        outputname = options.outfileprewithpath
        outfile = open(outputname + "." + options.typeOfcalculate + str(windowWidth) + "_" + str(slideSize) + str(os.getpid()) + chrlistfilewithoutpath, 'w')
        print("chrNo\twinNo\tfirstsnppos\tlastsnppos", *obsexpcaculator.vcfname_combination, sep="\t", file=outfile)
    elif options.typeOfcalculate == "hp" or options.typeOfcalculate == "pi":
        obsexpcaculator = Calculators.Calculate_Hp_master_slave(options.targetpopvcfconfig, options.outfileprewithpath, minsnps=0) if options.typeOfcalculate == "hp" else Calculators.Calculate_popPI(options.targetpopvcfconfig, options.outfileprewithpath, minsnps=0)
        plainname = re.search(r"[^/]*$", options.outfileprewithpath).group(0)
        if len(plainname) >= 228:
            outputname = options.outfileprewithpath[:-(len(plainname) - 228)]
        outputname = options.outfileprewithpath
        outfile = open(outputname + "." + options.typeOfcalculate + str(windowWidth) + "_" + str(slideSize) + str(os.getpid()) + chrlistfilewithoutpath, 'w')
        print("chrNo\twinNo\tfirstsnppos\tlastsnppos\tnoofsnp\twinvalue\tzvalue", file=outfile)
    
    if options.bedlikefile != None:
        obsexpcaculator.minsnps = minSNPs#7
    elif options.chromlistfilename != None:
        obsexpcaculator.minsnps = minSNPs#8
    # aaaa = open(options.outfileprewithpath + ".slidwin_filelist" + options.masterpid, 'a')
    Utils.safe_write(options.outfileprewithpath + ".slidwin_filelist" + options.masterpid,outputname + "." + options.typeOfcalculate + str(windowWidth) + "_" + str(slideSize) + str(os.getpid()) + chrlistfilewithoutpath)
    # print(outputname + "." + options.typeOfcalculate + str(windowWidth) + "_" + str(slideSize) + str(os.getpid()) + chrlistfilewithoutpath, file=aaaa)
    # aaaa.close()
    win = Utils.Window()
    win.outfilehandle=open(outputname + "." + options.typeOfcalculate + str(windowWidth) + "_" + str(slideSize) + str(os.getpid()) + chrlistfilewithoutpath+"_log",'w')
    obsexpsignalmapbychrom = {}
    if options.bedlikefile != None :genomedbtools = dbm.DBTools(config.ip, config.username, config.password, config.genomeinfodbname);mysqlchromtable = config.pekingduckchromtable
    
    for currentchrID, currentchrLenOrRegion in chromlistOrBedRegionList:
        if isinstance(currentchrLenOrRegion, tuple):
            currentchrLen = int(genomedbtools.operateDB("select", "select * from " + mysqlchromtable + " where chrID='" + currentchrID + "' ")[0][1])
        else:
            currentchrLen = currentchrLenOrRegion
        try:
            for vcfname in obsexpcaculator.vcfnamelist:
                vcfobj = obsexpcaculator.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname][0]
    #         for vcfobj in poplist:
                if currentchrID in vcfobj.VcfIndexMap:
                    break
            else:
                print("this chr doesn't exist in anypop")
                fillNA = [(0,0,*obsexpcaculator.getResult())]  # [(0,0,0,'NA')]
                if isinstance(currentchrLenOrRegion, tuple) and currentchrLenOrRegion[1] + extendsize < currentchrLen:
                    fillsize_End = currentchrLenOrRegion[1] + extendsize
                else:
                    fillsize_End = currentchrLen
                if isinstance(currentchrLenOrRegion, tuple) and currentchrLenOrRegion[0] - extendsize > 0:
                    fillsize_Start = currentchrLenOrRegion[0] - extendsize
                else:
                    fillsize_Start = 0
                for i in range(int((fillsize_End - fillsize_Start) / slideSize)):
                    noofsnp,t = obsexpcaculator.getResult()
                    s = i * slideSize + fillsize_Start
                    e = i * slideSize + windowWidth + fillsize_Start
                    fillNA.append((s,e,noofsnp,t))  # (0,0,0,'NA')
                obsexpsignalmapbychrom[(currentchrID, currentchrLenOrRegion)] = fillNA
                continue
            # this chr exist in one of the vcffile,then alinmultPopSnpPos
    #         for vcfobj_idx in range(len(poplist)):
            for vcfobj_idx in range(len(obsexpcaculator.vcfnamelist)):
                obsexpcaculator.listOfpopvcfRecsmapByAChr[vcfobj_idx] = {}
                vcfobj = obsexpcaculator.vcfnameKEY_vcfobj_pyBAMfilesVALUE[obsexpcaculator.vcfnamelist[vcfobj_idx]][0]
                print(obsexpcaculator.vcfnamelist[vcfobj_idx], "getvcf")
                if isinstance(currentchrLenOrRegion, tuple):  # bedfile
                    obsexpcaculator.listOfpopvcfRecsmapByAChr[vcfobj_idx][currentchrID] = vcfobj.getVcfListByChrom(currentchrID, currentchrLenOrRegion[0] - extendsize, currentchrLenOrRegion[1] + extendsize)
                else:  # chrom
                    obsexpcaculator.listOfpopvcfRecsmapByAChr[vcfobj_idx][currentchrID] = vcfobj.getVcfListByChrom(currentchrID)
            target_ref_SNPs = variantUtils.alignmultPopSnpPos(obsexpcaculator.listOfpopvcfRecsmapByAChr, "o",None,False)
            obsexpcaculator.currentchrID = currentchrID
            if options.typeOfcalculate == "early":
                obsexpcaculator.dynamicIU_toptable_obj.currentchrLen = currentchrLen
                obsexpcaculator.alignedSNP_absentinfo = {}
                obsexpcaculator.alignedSNP_absentinfo[currentchrID] = []
            ##########
            print(len(target_ref_SNPs[currentchrID]))
            if isinstance(currentchrLenOrRegion, tuple):  # bedfile
                print(currentchrID, currentchrLenOrRegion[0] - extendsize, currentchrLenOrRegion[1] + extendsize)
                win.slidWindowOverlap(target_ref_SNPs[currentchrID], min(currentchrLenOrRegion[1] + extendsize, currentchrLen), windowWidth, slideSize, obsexpcaculator, max(0, currentchrLenOrRegion[0] - extendsize))
            else:
                win.slidWindowOverlap(target_ref_SNPs[currentchrID], currentchrLen, windowWidth, slideSize, obsexpcaculator)
    
            obsexpsignalmapbychrom[(currentchrID, currentchrLenOrRegion)] = copy.deepcopy(win.winValueL)
        except Exception as e:
            print(e)
            variantUtils.pickle.dump(win.winValueL,open("breakrescue", 'wb'))

    
    for currentchrID, currentchrLenOrRegion in chromlistOrBedRegionList:
        if isinstance(currentchrLenOrRegion, tuple):
            currentchrLen = int(genomedbtools.operateDB("select", "select * from " + mysqlchromtable + " where chrID='" + currentchrID + "' ")[0][1])
        else:
            currentchrLen = currentchrLenOrRegion
        if (currentchrID, currentchrLenOrRegion) in obsexpsignalmapbychrom:
            for i in range(len(obsexpsignalmapbychrom[(currentchrID, currentchrLenOrRegion)])):
                if options.typeOfcalculate == "early" or options.typeOfcalculate=="D":
                    try:
                        if i==0:#print log
                            print(currentchrID,"winvalue is the correct value,zvalue is the not very appropriate value I used before",obsexpsignalmapbychrom[(currentchrID, currentchrLenOrRegion)][i],type(obsexpsignalmapbychrom[(currentchrID, currentchrLenOrRegion)][i][3][0]),type(obsexpsignalmapbychrom[(currentchrID, currentchrLenOrRegion)][i][3][1]))
                        print(currentchrID + "\t" + str(i) + "\t" + str(obsexpsignalmapbychrom[(currentchrID, currentchrLenOrRegion)][i][0]) + "\t" + str(obsexpsignalmapbychrom[(currentchrID, currentchrLenOrRegion)][i][1]) + "\t" + str(obsexpsignalmapbychrom[(currentchrID, currentchrLenOrRegion)][i][2]) + "\t" + '%.15f' % (obsexpsignalmapbychrom[(currentchrID, currentchrLenOrRegion)][i][3][0]) + "\t" + '%.12f' % (obsexpsignalmapbychrom[(currentchrID, currentchrLenOrRegion)][i][3][1]), file=outfile)
                    except :
                        print(currentchrID + "\t" + str(i) + "\t" + str(obsexpsignalmapbychrom[(currentchrID, currentchrLenOrRegion)][i][0]) + "\t" + str(obsexpsignalmapbychrom[(currentchrID, currentchrLenOrRegion)][i][1]) + "\t" + str(obsexpsignalmapbychrom[(currentchrID, currentchrLenOrRegion)][i][2]) + "\t" + obsexpsignalmapbychrom[(currentchrID, currentchrLenOrRegion)][i][3][0] + "\t" + obsexpsignalmapbychrom[(currentchrID, currentchrLenOrRegion)][i][3][1], file=outfile)
                elif options.typeOfcalculate == "pairfst" or options.typeOfcalculate == "dxy":
                    try :
                        print(currentchrID + "\t" + str(i) + "\t" + str(obsexpsignalmapbychrom[(currentchrID, currentchrLenOrRegion)][i][0]) + "\t" + str(obsexpsignalmapbychrom[(currentchrID, currentchrLenOrRegion)][i][1]) + "\t" + str(obsexpsignalmapbychrom[(currentchrID, currentchrLenOrRegion)][i][2]) + "\t" + '%.15f' % (obsexpsignalmapbychrom[(currentchrID, currentchrLenOrRegion)][i][3]) + '\t%.15f' % (obsexpsignalmapbychrom[(currentchrID, currentchrLenOrRegion)][i][3]) , file=outfile)
                    except:
                        print(currentchrID + "\t" + str(i) + "\t" + str(obsexpsignalmapbychrom[(currentchrID, currentchrLenOrRegion)][i][0]) + "\t" + str(obsexpsignalmapbychrom[(currentchrID, currentchrLenOrRegion)][i][1]) + "\t" + str(obsexpsignalmapbychrom[(currentchrID, currentchrLenOrRegion)][i][2]) + "\t" + "NA" + "\t" + "NA", file=outfile)
                elif options.typeOfcalculate == "is":
                    print(currentchrID + "\t" + str(i) + "\t" + str(obsexpsignalmapbychrom[(currentchrID, currentchrLenOrRegion)][i][0]) + "\t" + str(obsexpsignalmapbychrom[(currentchrID, currentchrLenOrRegion)][i][1]), *obsexpsignalmapbychrom[(currentchrID, currentchrLenOrRegion)][i][3], sep="\t", file=outfile)  # + "\t"+str()
                elif options.typeOfcalculate == "hp" or  options.typeOfcalculate == "pi":
                    try:
                        print(currentchrID + "\t" + str(i) + "\t" + str(obsexpsignalmapbychrom[(currentchrID, currentchrLenOrRegion)][i][0]) + "\t" + str(obsexpsignalmapbychrom[(currentchrID, currentchrLenOrRegion)][i][1]) + "\t" + str(obsexpsignalmapbychrom[(currentchrID, currentchrLenOrRegion)][i][2]) + "\t" + '%.15f' % (obsexpsignalmapbychrom[(currentchrID, currentchrLenOrRegion)][i][3]) + "\t0" , file=outfile)
                    
                    except:
                        print(currentchrID + "\t" + str(i) + "\t" + str(obsexpsignalmapbychrom[(currentchrID, currentchrLenOrRegion)][i][0]) + "\t" + str(obsexpsignalmapbychrom[(currentchrID, currentchrLenOrRegion)][i][1]) + "\t" + str(obsexpsignalmapbychrom[(currentchrID, currentchrLenOrRegion)][i][2]) + "\t" + 'NA' + "\t0" , file=outfile)
    outfile.close()
    win.outfilehandle.close()
    if options.topleveltablejudgeancestral != None:
        dbvariantstools.disconnect()
#     genomedbtools.disconnect()
    print("runSlave_slidewin process ID", os.getpid(), "done")
