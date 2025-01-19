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
prjpath=currentpath[:currentpath.find("atcgmagic")]#os.path.join("/opt/BDengine","build")
cfparser = configparser.ConfigParser()
if 'Windows' in platform.system():
    cfparser.read(prjpath+"atcgmagic\\com\\config.properties")
else:
    cfparser.read(os.path.join(prjpath,"atcgmagic/com/config.properties"))

ip=cfparser.get("driveanalysis","ip")
webdbname=cfparser.get("driveanalysis","webdbname")
scriptdir=cfparser.get("driveanalysis","scriptdir")
pathtoPython=cfparser.get("driveanalysis", "pathtoPython")

username="root"
password=secure.DBrPW

print("loading config",scriptdir,webdbname)#currentpath[:currentpath.find("life/src")]+"life/com/config.properties")
flag_OrderedLookup = []
flag_rLookup={}
index = 1  # 从 1 开始
section = "variantdata"
while True:
    key = f"corrs_presenceflag{index}"
    if cfparser.has_option(section, key):
        # print(key)
        value = cfparser.get(section, key).strip()  # 去除前后空白字符
        if value:  # 仅在值不为空时添加到列表中
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
            # for elem in vcfnameOrcantinedsamples:
            #     if "," not in elem:
            #         flag_rLookup[elem.strip()]=index
            #         flag_OrderedLookup.append(elem.strip())
            #     else:
            #         flag_OrderedLookup.append([sm.strip() for sm in elem.split(',')])
    else:
        break  # 如果配置文件中不存在这个键，则停止循环
    index += 1  # 索引递增，继续下一个键
# flag_OrderedLookup=[element.strip() for element in "flag ; Ordered;Lookup".split(';')]
for xyx in range(len(flag_OrderedLookup)):
    print(f"corrs_presenceflag bit {xyx+1}",flag_OrderedLookup[xyx])
# print(*flag_OrderedLookup,sep="\n")
vport=cfparser.get(section,"port")
variantsdbname=cfparser.get(section,"vdbname")
topleveltableofvdb=cfparser.get(section,'topleveltableofvdb')
beijingreffa=cfparser.get(section,"refgenomefa_ZJU1")
outgroupVCFBAMconfig_ZJU1ref=cfparser.get(section,"outgroupVCFBAMconfig_ZJU1ref").strip()
genomeinfodbname=section = "genomeinfo"
#cfparser.get(section,"genomeinfodbname")
pekingduckchromtable=cfparser.get(section,'ZJU1duckchromtable')#cfparser.get(section,"BGI1duckchromtable")

ghostdbname=cfparser.get(section,"ghostdbname")
# vcfdbname=cfparser.get("mysqldatabase","vcfdbname")
TranscriptGenetable=cfparser.get(section,"ZJU1TranscriptGenetable")
# D2Bduckchromtable=cfparser.get("mysqldatabase","D2Bduckchromtable")
# username=cfparser.get("mysqldatabase","username")
# password=cfparser.get("mysqldatabase","password")

UPLOAD_FOLDER="static/video"
def random_str(randomlength=8):
    a = list(string.ascii_letters)
    random.shuffle(a)
    return ''.join(a[:randomlength])