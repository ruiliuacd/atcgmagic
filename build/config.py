'''
Created on 2023年2月6日

@author: RuiLiu
'''

import platform, sys, inspect, os, configparser
import string,random
import secure


if not hasattr(sys.modules[__name__], '__file__'):
    __file__ = inspect.getfile(inspect.currentframe())
    
currentpath=os.path.realpath(__file__)
prjpath=currentpath[:currentpath.find("atcgtoolkit")]#os.path.join("/opt/BDengine","build")
cfparser = configparser.ConfigParser()
if 'Windows' in platform.system():
    cfparser.read(prjpath+"atcgtoolkit\\com\\config.properties")
else:
    cfparser.read(os.path.join(prjpath,"atcgtoolkit/com/config.properties"))

ip=cfparser.get("driveanalysis","ip")

pathtoPython=cfparser.get("driveanalysis", "pathtoPython")

username="root"
password=secure.DBrPW

flag_OrderedLookup = []
flag_rLookup={}
index = 1  # 从 1 开始
section = "variantdata"
while True:
    key = f"corrs_presenceflag{index}"
    if cfparser.has_option(section, key):
        print(key)
        value = cfparser.get(section, key).strip()  
        if value:  #
            vcfnameOrcantinedsamples=value.split(';')
            if len(vcfnameOrcantinedsamples)==1:
                flag_OrderedLookup.append(value.strip())
                flag_rLookup[value.strip()]=index
            elif ',' not in vcfnameOrcantinedsamples[1]:#multple ; 
                flag_OrderedLookup.append(vcfnameOrcantinedsamples)
                for elem in vcfnameOrcantinedsamples:
                    flag_rLookup[elem.strip()]=index
            elif ',' in vcfnameOrcantinedsamples[1]:
                flag_values = [v.strip() for v in vcfnameOrcantinedsamples[1].split(',')]
                flag_OrderedLookup.append({vcfnameOrcantinedsamples[0].strip(): flag_values})
                if vcfnameOrcantinedsamples[0].strip() not in flag_rLookup:
                    flag_rLookup[vcfnameOrcantinedsamples[0].strip()]=[index]
                else:
                    flag_rLookup[vcfnameOrcantinedsamples[0].strip()].append(index)

    else:
        break  # 
    index += 1  # 
# flag_OrderedLookup=[element.strip() for element in "flag ; Ordered;Lookup".split(';')]
print(*flag_OrderedLookup,sep="\n")
vport=cfparser.get(section,"port")
variantsdbname=cfparser.get(section,"vdbname")
topleveltableofvdb=cfparser.get(section,'topleveltableofvdb')
beijingreffa=cfparser.get(section,"refgenomefa_ZJU1")
outgroupVCFBAMconfig_ZJU1ref=cfparser.get(section,"outgroupVCFBAMconfig_ZJU1ref").strip()

def random_str(randomlength=8):
    a = list(string.ascii_letters)
    random.shuffle(a)
    return ''.join(a[:randomlength])