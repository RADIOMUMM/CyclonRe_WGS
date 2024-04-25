import argparse
import os
from Cyclone_mod import cycloneQC,writeShell,Cmdruner,getConfig,cyclonemapping,cyclonecallVar,cyclonecallSV
import time

parse=argparse.ArgumentParser(prog="run_cyclone_work",usage="%(prog)s -i <fq>  -n <sample_name>  -o <outputFile>", description="cyclone fq QC workflow",epilog="radiomumm")
parse.add_argument('-i',"--fastq",type=str,help="input fsatq")
parse.add_argument('-n',"--sample_name",type=str,help="sample name")
parse.add_argument('-o', '--outputFile', help = "output file path,Note that you need to add the \"/\" symbol at the end of the folder,eg:'/home/','/home/result/'")

parse.add_argument('-callSNP_indel',"--callsmallVariants",type=str,default="None",help="The process can autonomously choose whether or not to perform SNP and InDel testing and evaluation.Choose between None or True.(optional)")
parse.add_argument('-t',"--qcthread",type=int,default=20,help="Threads required for the QC program.(optional)")
parse.add_argument('-topcrop',"--topcrop",type=int,default=100,help="Trim n nucleotides from start of read.(optional)")
parse.add_argument('-tailcrop', '--tailcrop',type=int,default=50,help = "Trim n nucleotides from end of read.(optional)")
parse.add_argument('-q',"--trim_quality",type=int,default=10,help="Filter on a minimum average read quality score.(optional)")
parse.add_argument('-qctype', '--QCtype',type=str,default="NanoStat",help = "Select the software used in the QC step. Choose between NanoPlot or NanoStat.(optional)")
parse.add_argument('-standard', '--standardtype',type=str,default="HG002",help = "Select the standard sample file . Choose between HG002 or HG001.(optional)")
parse.add_argument('-ref', '--REF',type=str,default="hg19",help = "Select the reference genome used in this analysis. Choose between hg19 or hg38.(optional)")
parse.add_argument('-t2', '--mappingthread',type=int,default=20,help = "Threads required for the Mapping program.(optional)")
parse.add_argument('-t3', '--samtoolsthread',type=int,default=20,help = "Threads required for the samtools program.(optional)")
parse.add_argument('-t4', '--Callvariantthread',type=int,default=20,help = "Threads required for the Clair3 program.(optional)")
result=parse.parse_args()

FQ_PATH=result.fastq
samplename=result.sample_name
outputp=result.outputFile

qcthread=result.qcthread
topcrop=result.topcrop
tailcrop=result.tailcrop
trim_quality=result.trim_quality
qctype=result.QCtype
ref_fa=result.REF
callvar=result.callsmallVariants
mappingthread=result.mappingthread
samtoolsthread=result.samtoolsthread
varthread=result.Callvariantthread
standard=result.standardtype
output="%s%s_Cyclone"%(outputp,samplename)
os.system("mkdir -p %s/shell/01.QC %s/result/01.QC"%(output,output))
rawsamplename="%s_rawdata"%samplename
cleansamplename="%s_cleandata"%samplename
raw_data=cycloneQC(rawsamplename,output,trim_quality,topcrop,tailcrop,qcthread)
clean_data=cycloneQC(cleansamplename,output,trim_quality,topcrop,tailcrop,qcthread)
conda_path=getConfig("config","conda","path")
Conda_env="source %s/bin/activate Cyclon_humen"%conda_path
#QC module
qcfinshfile="%s/result/01.QC/%s.QC_info.xlsx"%(output,samplename)
if os.path.exists(qcfinshfile):
    shell1="echo \"QC step is already done, skip\""
    shell2="echo \"QC step is already done, skip\""
    shell3="echo \"QC step is already done, skip\""
    pass
else:
    if qctype == "LongQC":
        longqc="%s/result/01.QC/%s_cleandata/web_summary.html"%(output,samplename)
        if os.path.exists(longqc):
            shell1="echo \"QC step is already done, skip\""
            shell2="echo \"QC step is already done, skip\""
            shell3="echo \"QC step is already done, skip\""
        else:
            shell1=raw_data.run_LongQC(FQ_PATH)
            shell2=clean_data.run_trim(FQ_PATH)
            trim_fq="%s/result/01.QC/%s_trimmed-reads.fastq.gz"%(output,cleansamplename)
            shell3=clean_data.run_LongQC(trim_fq)
    elif qctype == "NanoPlot":
        nanoplot="%s/result/01.QC/%s_cleandata/NanoPlot-report.html"%(output,samplename)
        if os.path.exists(nanoplot):
            shell1="echo \"QC step is already done, skip\""
            shell2="echo \"QC step is already done, skip\""
            shell3="echo \"QC step is already done, skip\"" 
        else:
            shell1=raw_data.run_NanoPlot(FQ_PATH)
            shell2=clean_data.run_trim(FQ_PATH)
            fail_fq="%s/result/01.QC/%s_trimmed-reads.fastq.gz"%(output,cleansamplename)
            shell3=clean_data.run_NanoPlot(fail_fq)
    elif qctype == "NanoStat":
        nanostat="%s/result/01.QC/%s_cleandata/NanoStats.txt"%(output,samplename)
        if os.path.exists(nanostat):
            shell1="echo \"QC step is already done, skip\""
            shell2="echo \"QC step is already done, skip\""
            shell3="echo \"QC step is already done, skip\""
        else:
            shell1=raw_data.run_NanoStat(FQ_PATH)
            shell2=clean_data.run_trim(FQ_PATH)
            fail_fq="%s/result/01.QC/%s_trimmed-reads.fastq.gz"%(output,cleansamplename)
            shell3=clean_data.run_NanoStat(fail_fq)   
        pass
    else:
        print("The input parameters do not meet the requirements.")

#mapping module
mappingfinshfile="%s/result/02.Mapping/%s.minimap2.sample_bamstat.txt"%(output,samplename)
if os.path.exists(mappingfinshfile):
    shell4="echo \"Mapping step is already done, skip\""
    shell5="echo \"Mapping step is already done, skip\""
    shell6="echo \"Mapping step is already done, skip\""
    shell7="echo \"Mapping step is already done, skip\""
    pass
else:
    shelldir="%s/shell/02.Mapping/%s/"%(output,samplename)
    mappingdata=cyclonemapping(samplename,output,mappingthread,samtoolsthread,ref_fa)
    shell4,shell5,shell6,shell7=mappingdata.run_mminimap2()
#call var module
Varfinshfile="%s/result/03.Variants/merge_output.vcf.gz.tbi"%(output)
if os.path.exists(Varfinshfile):
    shell8="echo \"Call var step is already done, skip\""
    shell8_1="echo \"Evaluation step is already done, skip\""
    shell8_2="echo \"Evaluation step is already done, skip\""
    shell8_3="echo \"Evaluation step is already done, skip\""
    shell8_4="echo \"Evaluation step is already done, skip\""
    pass
else:
    Vardata=cyclonecallVar(samplename,output,varthread,ref_fa)
    if standard == "HG001":
        shell8,shell8_1,shell8_2=Vardata.run_callVar()
        shell8_3,shell8_4=Vardata.run_rtgtools(standard)
    elif standard == "HG002":
        shell8,shell8_1,shell8_2=Vardata.run_callVar()
        shell8_3,shell8_4=Vardata.run_rtgtools(standard)
    else:
       shell8,shell8_1,shell8_2=Vardata.run_callVar()
       shell8_3="echo \"no need evaluation\""
       shell8_4="echo \"no need evaluation\""
       pass
#call SVs module
VSfinshfile="%s/result/04.SVs/SV_BASE/summary.json"%(output)
if os.path.exists(VSfinshfile):
    shell9="echo \"call SVs step is already done, skip\""
    shell9_1="echo \"Evaluation step is already done, skip\""
    pass
else:
    SVdata=cyclonecallSV(samplename,output,varthread,ref_fa)
    if standard == "HG001":
        shell9=SVdata.run_callSV()
        shell9_1=SVdata.run_SVcheck(standard)
    elif standard == "HG002":
        shell9=SVdata.run_callSV()
        shell9_1=SVdata.run_SVcheck(standard)
    else:
        shell9=SVdata.run_callSV()
        shell9_1="echo \"no need evaluation\""
        pass
CyclonRe_dir=getConfig("config","CyclonRe","path")
report_list=["%s && python3 %s/run_cyclone_report.py -n %s -o %s -ref %s -standard %s -callSNP_indel %s \n"%(Conda_env,CyclonRe_dir,samplename,outputp,ref_fa,standard,callvar)]
writeShell("report_model",report_list,"%s/shell/"%(output))
shelllist=[shell1,
           shell2,shell3,shell4,
           shell5,
           shell6,shell7,shell8,shell9,
           shell8_1,shell8_2,
           shell8_3,shell8_4,
           shell9_1]
newshelllist=[]
for shell in shelllist:
    shell = "%s && %s"%(Conda_env,shell)
    newshelllist.append(shell)
    
writeShell("trim_CyLR_WGS_step1",[newshelllist[1]],"%s/shell/"%(output))
writeShell("mapping_CyLR_WGS_step2",[newshelllist[0],newshelllist[2],newshelllist[3]],"%s/shell/"%(output))
writeShell("makebai_CyLR_WGS_step3",[newshelllist[4]],"%s/shell/"%(output))
mainshell1="%s/shell/trim_CyLR_WGS_step1.sh"%(output)
mainshell2="%s/shell/mapping_CyLR_WGS_step2.sh"%(output)
mainshell3="%s/shell/makebai_CyLR_WGS_step3.sh"%(output)
if callvar == "True":
    writeShell("call_Var_CyLR_WGS_step4",[newshelllist[5],newshelllist[6],newshelllist[7],newshelllist[8]],"%s/shell/"%(output))
    writeShell("SV_evaluation_CyLR_WGS_step5",[newshelllist[9],newshelllist[10],newshelllist[13]],"%s/shell/"%(output))
    writeShell("SNP_evaluation_step6",[newshelllist[11],newshelllist[12]],"%s/shell/"%(output))
    mainshell4="%s/shell/call_Var_CyLR_WGS_step4.sh"%(output)
    mainshell5="%s/shell/SV_evaluation_CyLR_WGS_step5.sh"%(output)
    mainshell6="%s/shell/SNP_evaluation_step6.sh"%(output)
    mainshell7="sh %s/shell/report_model.sh"%(output)
    mainshelllist=[mainshell1,mainshell2,mainshell3,mainshell4,mainshell5,mainshell6]
elif callvar == "None":
    writeShell("call_SVs_CyLR_WGS_step4",[newshelllist[5],newshelllist[6],newshelllist[8]],"%s/shell/"%(output))
    writeShell("SV_evaluation_CyLR_WGS_step5",[newshelllist[13]],"%s/shell/"%(output))
    mainshell4="%s/shell/call_SVs_CyLR_WGS_step4.sh"%(output)
    mainshell5="%s/shell/SV_evaluation_CyLR_WGS_step5.sh"%(output)
    mainshell7="sh %s/shell/report_model.sh"%(output)
    mainshelllist=[mainshell1,mainshell2,mainshell3,mainshell4,mainshell5]

time_start_all = time.time()

for shell in mainshelllist:
    qsub="%s/qsub_sge_plus.pl --interval 30 --convert no --queue all.q --resource vf=100G --cpu 12 --lines 1 %s"%(CyclonRe_dir,shell)
    time_start = time.time()
    Cmdruner(qsub)
    time_end = time.time()
    time_sum = time_end - time_start
    min = time_sum/60
    print("%s run %.2f min"%(shell,min))
Cmdruner(mainshell7)
time_end_all = time.time()
time_sum_all = time_end_all - time_start_all
min_all = time_sum_all/60
print("The total time the program runs is %.2f min"%min_all)