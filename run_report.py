# QC1
# QC2
# Mapping
# Small vars
# SVs
import argparse
import os
from Cyclone_mod import report,report_nortrg
from mergeQC_report import longqc_report,nanoplot_report,merge_report

parse=argparse.ArgumentParser(prog="CycloneQC",usage="%(prog)s -i <fq>  -n <sample_name>  -o <outputFile>", description="cyclone fq QC workflow",epilog="radiomumm")
parse.add_argument('-n',"--sample_name",type=str,help="sample name")
parse.add_argument('-o', '--outputFile', help = "output file path,Note that you need to add the \"/\" symbol at the end of the folder,eg:'/home/','/home/result/'")
parse.add_argument('-qctype', '--QCtype',type=str,default="NanoStat",help = "Select the software used in the QC step. Choose between NanoPlot , NanoStat or LongQC.(optional)")
parse.add_argument('-standard', '--standardtype',type=str,default="HG002",help = "Select the standard sample file . Choose between HG002 or HG001.(optional)")
parse.add_argument('-ref', '--REF',type=str,default="hg38",help = "Select the reference genome used in this analysis. Choose between hg19 or hg38.(optional)")
parse.add_argument('-callSNP_indel',"--callsmallVariants",type=str,default="None",help="The process can autonomously choose whether or not to perform SNP and InDel testing and evaluation.Choose between None or True.(optional)")
result=parse.parse_args()
samplename=result.sample_name
outputp=result.outputFile
qctype=result.QCtype
ref_fa=result.REF
standard=result.standardtype
callvar=result.callsmallVariants

output="%s%s_Cyclone"%(outputp,samplename)
rawsamplename="%s_rawdata"%samplename
cleansamplename="%s_cleandata"%samplename

#report module
if qctype == "LongQC":
    file1=longqc_report(rawsamplename,output)
    file2=longqc_report(cleansamplename,output)
    merge_report(samplename,output,file1,file2)
elif qctype == "NanoPlot" or "NanoStat":
    file1=nanoplot_report(rawsamplename,output)
    file2=nanoplot_report(cleansamplename,output)
    merge_report(samplename,output,file1,file2)
else:
    print("The input parameters do not meet the requirements.")

if callvar == "TRUE": 
    report(samplename,output,ref_fa,standard)
elif callvar == "None":
    report_nortrg(samplename,output,ref_fa,standard)