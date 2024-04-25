#QC moudle
#software ： LongQC/NanoPlot/LRAPqc

#mapping moudle
#software ： NGLMR/minimap2

#call variants moudle
#software ： Clair3/deepvariant

#call SV moudle
#software ： Sniffles/cuteSV
import os
import configparser
import subprocess
import sys
import base64
from jinja2 import Environment, FileSystemLoader
from mergeQC_report import report_tmp,get_rtg,get_SVcheck,report_tmp_var

def pic2base64(pic_path):
    with open(pic_path,'rb') as f:
        image_base64 = base64.b64encode(f.read())
        image_base64="<p align=\"center\"><img src=\"data:image/png;base64,"+str(image_base64).split("\'")[1]+"\" width=\"80%\" height=\"400px\"></p>"
    return image_base64

def Cmdruner(cmdname):
    cmd=cmdname
    project=subprocess.call(cmd,shell=True)
    print(cmd)
    if project==0:
        print("step successful")
    else:
        print("step failed")
        sys.exit(1)
def getConfig(filename,section,option):
    if getattr(sys, 'frozen', False):
        application_path = os.path.dirname(sys.executable)
    elif __file__:
        application_path = os.path.dirname(__file__)
        pass
    configPath = os.path.join(application_path, filename)
    conf = configparser.ConfigParser()
    conf.read(configPath,encoding='utf-8')
    config = conf.get(section, option)
    return config
def writeShell(stepname,cmd,shelldir):
    with open("%s%s.sh"%(shelldir,stepname),"w") as shell:
        for i in cmd:
            cmdtxt=str(i)+"\n"
            shell.write(cmdtxt)

class cycloneQC:
    def __init__(self,samplename,outputpath,quality,headcrop,tailcrop,qc_thread):
        self.samplename=samplename
        self.outputpath=outputpath
        self.quality=quality
        self.headcrop=headcrop
        self.tailcrop=tailcrop
        self.qc_thread=qc_thread
    def run_LongQC(self,fqpath):
        '''使用LongQC进行QC'''
        shelldir="%s/shell/01.QC/%s/"%(self.outputpath,self.samplename)
        Cmdruner("mkdir -p %s"%shelldir)
        # print("mkdir -p %s"%samplenamedir)
        cmdlist=[]
        LongQC_path=getConfig("config","LongQC","path")
        cmd="python3 %s sampleqc  -x ont-ligation -o %s/result/01.QC/%s %s"%(LongQC_path,self.outputpath,self.samplename,fqpath)
        cmdlist.append(cmd)
        writeShell("cycloneQC_LongQC_%s"%self.samplename,cmdlist,shelldir)
        sh="sh %scycloneQC_LongQC_%s.sh"%(shelldir,self.samplename)
        return sh
    def run_NanoPlot(self,fqpath):
        '''使用NanoPlot进行QC'''
        shelldir="%s/shell/01.QC/%s/"%(self.outputpath,self.samplename)
        Cmdruner("mkdir -p %s "%shelldir)
        cmdlist=[]
        NanoPlot_path=getConfig("config","NanoPlot","path")
        # file_size = os.stat(fqpath).st_size
        # if file_size > 2*1024*1024*1024:
        #     cmd=" %s -t %s --huge --N50 --dpi 300 --fastq %s --title %s -o %s/result/01.QC/%s"%(NanoPlot_path,self.qc_thread,fqpath,self.samplename,self.outputpath,self.samplename)
        # else:
        #     cmd=" %s -t %s --N50 --dpi 300 --fastq %s --title %s -o %s/result/01.QC/%s"%(NanoPlot_path,self.qc_thread,fqpath,self.samplename,self.outputpath,self.samplename)
        cmd=" %s -t %s --huge --N50 --dpi 300 --fastq %s --title %s -o %s/result/01.QC/%s"%(NanoPlot_path,self.qc_thread,fqpath,self.samplename,self.outputpath,self.samplename)
        cmdlist.append(cmd)
        writeShell("cycloneQC_NanoPlot_%s"%self.samplename,cmdlist,shelldir)
        sh="sh %scycloneQC_NanoPlot_%s.sh"%(shelldir,self.samplename)
        return sh
    def run_trim(self,fqpath):
        '''使用NanoFilt进行数据过滤'''
        shelldir="%s/shell/01.QC/%s/"%(self.outputpath,self.samplename)
        Cmdruner("mkdir -p %s "%(shelldir))
        cmdlist=[]
        NanoFilt_path=getConfig("config","NanoFilt","path")
        cmd="gunzip -c %s | %s -q %s --headcrop %s --tailcrop %s |pigz >%s/result/01.QC/%s_trimmed-reads.fastq.gz"%(fqpath,NanoFilt_path,self.quality,self.headcrop,self.tailcrop,self.outputpath,self.samplename)
        cmdlist.append(cmd)
        writeShell("cycloneQC_trim_%s"%self.samplename,cmdlist,shelldir)
        sh="sh %scycloneQC_trim_%s.sh"%(shelldir,self.samplename)
        return sh
    def run_NanoStat(self,fqpath):
        '''使用NanoStat数据质控'''
        shelldir="%s/shell/01.QC/%s/"%(self.outputpath,self.samplename)
        resultdir="%s/result/01.QC/%s/"%(self.outputpath,self.samplename)
        Cmdruner("mkdir -p %s %s"%(shelldir,resultdir))
        cmdlist=[]
        NanoStat_path=getConfig("config","NanoStat","path")
        cmd1="%s --fastq %s -o %s -n NanoStats.txt"%(NanoStat_path,fqpath,resultdir)
        cmdlist.append(cmd1)
        writeShell("cycloneQC_NanoStat_%s"%self.samplename,cmdlist,shelldir)
        sh="sh %scycloneQC_NanoStat_%s.sh"%(shelldir,self.samplename)
        return sh        
class cyclonemapping:
    def __init__(self,samplename,outputpath,mapping_thread,samtools_thread,fa_ref) -> None:
        self.samplename=samplename
        self.outputpath=outputpath
        self.mapping_thread=mapping_thread
        self.samtools_thread=samtools_thread
        self.fa_ref=fa_ref   
        pass  
    def run_mminimap2(self):
        '''使用minimap2进行比对'''
        if self.fa_ref == "hg19":
            fa=getConfig("config","ref","hg19")
        elif self.fa_ref == "hg38":
            fa=getConfig("config","ref","hg38")
        elif self.fa_ref == "ecoli":
            fa=getConfig("config","ref","ecoli")
        else:
            print("The input ref do not meet the requirements.")

        shelldir="%s/shell/02.Mapping/%s/"%(self.outputpath,self.samplename)
        outputdir="%s/result/02.Mapping/"%(self.outputpath)
        Cmdruner("mkdir -p %s %s"%(shelldir,outputdir))
        cmdlist1=[]
        cmdlist2=[]
        cmdlist3=[]
        cmdlist4=[]
        minimap2_path=getConfig("config","minimap2","path")
        samtools_path=getConfig("config","samtools","path")
        trimfq="%s/result/01.QC/%s_cleandata_trimmed-reads.fastq.gz"%(self.outputpath,self.samplename)
        # cmd1="%s -t %s -ax map-ont  %s %s |%s sort -@ %s -m 5G -o %s/%s.minimap2.sort.bam"%(minimap2_path,self.mapping_thread,fa,trimfq,samtools_path,self.samtools_thread,outputdir,self.samplename)
        cmd1="%s -t %s -ax map-ont  --secondary=no %s %s |%s sort -@ %s -m 5G -o %s/%s.minimap2.sort.bam"%(minimap2_path,self.mapping_thread,fa,trimfq,samtools_path,self.samtools_thread,outputdir,self.samplename)
        cmd2="%s index %s/%s.minimap2.sort.bam"%(samtools_path,outputdir,self.samplename)
        cmd3="%s flagstat %s/%s.minimap2.sort.bam >%s/%s.minimap2.sort.bam.falgstat"%(samtools_path,outputdir,self.samplename,outputdir,self.samplename)   
        cmd4="%s stats --threads %s %s/%s.minimap2.sort.bam >%s/%s.minimap2.sample_bamstat.txt"%(samtools_path,self.samtools_thread,outputdir,self.samplename,outputdir,self.samplename)
        cmdlist1.append(cmd1)
        cmdlist2.append(cmd2)
        cmdlist3.append(cmd3)
        cmdlist4.append(cmd4)
        writeShell("cycloneMapping_step1_%s"%self.samplename,cmdlist1,shelldir)
        sh1="sh %scycloneMapping_step1_%s.sh"%(shelldir,self.samplename)
        
        writeShell("cycloneMapping_step2_%s"%self.samplename,cmdlist2,shelldir)
        sh2="sh %scycloneMapping_step2_%s.sh"%(shelldir,self.samplename)
        
        writeShell("cycloneMapping_step3_%s"%self.samplename,cmdlist3,shelldir)
        sh3="sh %scycloneMapping_step3_%s.sh"%(shelldir,self.samplename)
        
        writeShell("cycloneMapping_step4_%s"%self.samplename,cmdlist4,shelldir)
        sh4="sh %scycloneMapping_step4_%s.sh"%(shelldir,self.samplename)
        return sh1,sh2,sh3,sh4
        # writeShell("cycloneMapping_%s"%self.samplename,cmdlist,shelldir)
        # sh="sh %scycloneMapping_%s.sh"%(shelldir,self.samplename)
        # return sh
        
class cyclonecallVar:
    def __init__(self,samplename,outputpath,align_thread,fa_ref) -> None:
        self.samplename=samplename
        self.outputpath=outputpath
        self.align_thread=align_thread
        self.fa_ref=fa_ref
        pass
    def run_callVar(self):
        if self.fa_ref == "hg19":
            fa=getConfig("config","ref","hg19")
        elif self.fa_ref == "hg38":
            fa=getConfig("config","ref","hg38")
        elif self.fa_ref == "ecoli":
            fa=getConfig("config","ref","ecoli")        
        else:
            print("The input ref do not meet the requirements.")
        shelldir="%s/shell/03.Variants/%s/"%(self.outputpath,self.samplename)
        outputdir="%s/result/03.Variants/"%(self.outputpath)
        Cmdruner("mkdir -p %s %s"%(shelldir,outputdir))
        cmdlist1=[]
        cmdlist2=[]
        cmdlist3=[]
        
        Clair3_path=getConfig("config","Clair3","path")
        Clair3_module_path=getConfig("config","Clair3","moudle")
        bcftools_path=getConfig("config","bcftools","path")
        BAMfile="%s/result/02.Mapping/%s.minimap2.sort.bam"%(self.outputpath,self.samplename)
        # cmd1="source /software/miniconda3/bin/activate clair3"
        # cmdlist.append(cmd1)
        cmd2="%s --bam_fn=%s --ref_fn=%s --threads=%s --platform=ont --model_path=%s --output=%s"%(Clair3_path,BAMfile,fa,self.align_thread,Clair3_module_path,outputdir)
        if self.fa_ref == "ecoli":
            cmd3="%s --include_all_ctgs"%cmd2
        else:
            cmd3=cmd2
        cmdlist1.append(cmd3)
        cmd5="%s stats %s/merge_output.vcf.gz >%s/%s.variant.stat.txt"%(bcftools_path,outputdir,outputdir,self.samplename)
        cmdlist1.append(cmd5)
        # filterSNP
        gatk_path=getConfig("config","gatk","path")
        gatk_thread=getConfig("config","gatk","thread")
        tabix_path=getConfig("config","tabix","path")
        cmd6="%s --java-options  %s SelectVariants -R %s -V %s/merge_output.vcf.gz  --select-type SNP -O %s/%s.raw.snp.vcf "%(gatk_path,gatk_thread,fa,outputdir,outputdir,self.samplename)
        cmd7="%s  --java-options %s  VariantFiltration -R %s -V %s/%s.raw.snp.vcf --filter-expression \"QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" --filter-name 'SNP_filter' -O %s/%s.filter.snp.vcf"%(gatk_path,gatk_thread,fa,outputdir,self.samplename,outputdir,self.samplename)
        cmd8="%s --java-options %s  SelectVariants  -R %s -V %s/%s.filter.snp.vcf --exclude-filtered  -O %s/%s.filtered.snp.vcf"%(gatk_path,gatk_thread,fa,outputdir,self.samplename,outputdir,self.samplename)
        SNPvcf_file="%s/%s.filtered.snp.vcf"%(outputdir,self.samplename)
        cmd14="%s view %s -Oz -o %s.gz"%(bcftools_path,SNPvcf_file,SNPvcf_file)
        cmd15="%s -p vcf %s.gz"%(tabix_path,SNPvcf_file)
        cmdlist2.append(cmd6)
        cmdlist2.append(cmd7)
        cmdlist2.append(cmd8)
        cmdlist2.append(cmd14)
        cmdlist2.append(cmd15)
        # filterindel
        cmd9="%s  --java-options  %s SelectVariants -R %s -V %s/merge_output.vcf.gz  --select-type INDEL -O %s/%s.raw.indel.vcf "%(gatk_path,gatk_thread,fa,outputdir,outputdir,self.samplename)
        cmd10="%s  --java-options %s  VariantFiltration -R %s -V %s/%s.raw.indel.vcf --filter-expression \"QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" --filter-name 'INDEL_filter' -O %s/%s.filter.indel.vcf"%(gatk_path,gatk_thread,fa,outputdir,self.samplename,outputdir,self.samplename)
        cmd11="%s --java-options %s  SelectVariants  -R %s -V %s/%s.filter.indel.vcf --exclude-filtered  -O %s/%s.filtered.indel.vcf"%(gatk_path,gatk_thread,fa,outputdir,self.samplename,outputdir,self.samplename)
        InDelvcf_file="%s/%s.filtered.indel.vcf"%(outputdir,self.samplename)
        cmd16="%s view %s -Oz -o %s.gz"%(bcftools_path,InDelvcf_file,InDelvcf_file)
        cmd17="%s -p vcf %s.gz"%(tabix_path,InDelvcf_file)
        cmdlist3.append(cmd9)
        cmdlist3.append(cmd10)
        cmdlist3.append(cmd11)
        cmdlist3.append(cmd16)
        cmdlist3.append(cmd17)
        writeShell("cycloneCallVar_step1_%s"%self.samplename,cmdlist1,shelldir)
        writeShell("cycloneCallVar_step2_%s"%self.samplename,cmdlist2,shelldir)
        writeShell("cycloneCallVar_step3_%s"%self.samplename,cmdlist3,shelldir)
        sh1="sh %scycloneCallVar_step1_%s.sh"%(shelldir,self.samplename)
        sh2="sh %scycloneCallVar_step2_%s.sh"%(shelldir,self.samplename)
        sh3="sh %scycloneCallVar_step3_%s.sh"%(shelldir,self.samplename)
        return sh1,sh2,sh3
    def run_rtgtools(self,standard):
        """https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/"""
        shelldir="%s/shell/03.Variants/%s/"%(self.outputpath,self.samplename)
        outputdir="%s/result/03.Variants/"%(self.outputpath)
        cmdlist1=[]
        cmdlist2=[]
        rtgtool_path=getConfig("config","rtgtool","path")
        if self.fa_ref == "hg38":
            if standard == "HG001":
                snp_baseline=getConfig("config","HG001","hg38_snp_baseline")
                indel_baseline=getConfig("config","HG001","hg38_indel_baseline")
                SDF=getConfig("config","HG001","hg38_SDF")
                evaluationbed=getConfig("config","HG001","hg38_evaluationbed")
            elif standard == "HG002":
                snp_baseline=getConfig("config","HG002","hg38_snp_baseline")
                indel_baseline=getConfig("config","HG002","hg38_indel_baseline")
                SDF=getConfig("config","HG002","hg38_SDF")
                evaluationbed=getConfig("config","HG002","hg38_evaluationbed")                
        elif self.fa_ref == "hg19":
            if standard == "HG001":
                snp_baseline=getConfig("config","HG001","hg19_snp_baseline")
                indel_baseline=getConfig("config","HG001","hg19_indel_baseline")
                SDF=getConfig("config","HG001","hg19_SDF")
                evaluationbed=getConfig("config","HG001","hg19_evaluationbed")
            elif standard == "HG002":
                snp_baseline=getConfig("config","HG002","hg19_snp_baseline")
                indel_baseline=getConfig("config","HG002","hg19_indel_baseline")
                SDF=getConfig("config","HG002","hg19_SDF")
                evaluationbed=getConfig("config","HG002","hg19_evaluationbed")           
            
        SNPvcf_file="%s/%s.filtered.snp.vcf"%(outputdir,self.samplename)
        InDelvcf_file="%s/%s.filtered.indel.vcf"%(outputdir,self.samplename)
        cmd12="%s vcfeval -c %s.gz -b %s --evaluation-regions %s -t %s -T %s -o %s/rtgtools_snp"%(rtgtool_path,SNPvcf_file,snp_baseline,evaluationbed,SDF,self.align_thread,outputdir)
        cmd13="%s vcfeval -c %s.gz -b %s --evaluation-regions %s -t %s -T %s -o %s/rtgtools_indel"%(rtgtool_path,InDelvcf_file,indel_baseline,evaluationbed,SDF,self.align_thread,outputdir)
        cmdlist1.append(cmd12)
        cmdlist2.append(cmd13)
        writeShell("cycloneRTG_step1%s"%self.samplename,cmdlist1,shelldir)
        sh1="sh %scycloneRTG_step1%s.sh"%(shelldir,self.samplename)
        writeShell("cycloneRTG_step2%s"%self.samplename,cmdlist2,shelldir)
        sh2="sh %scycloneRTG_step2%s.sh"%(shelldir,self.samplename)
        return sh1,sh2
        # writeShell("cycloneRTG_%s"%self.samplename,cmdlist,shelldir)
        # sh="sh %scycloneRTG_%s.sh"%(shelldir,self.samplename)
        # return sh       
class cyclonecallSV(cyclonecallVar):
    def __init__(self,samplename,outputpath,align_thread,fa_ref) -> None:
        super().__init__(samplename,outputpath,align_thread,fa_ref)
        pass
    def run_callSV(self):
        shelldir="%s/shell/04.SVs/%s/"%(self.outputpath,self.samplename)
        outputdir="%s/result/04.SVs/"%(self.outputpath)
        Cmdruner("mkdir -p %s %s"%(shelldir,outputdir))
        cmdlist=[]
        Sniffles_path=getConfig("config","Sniffles","path")
        BAMfile="%s/result/02.Mapping/%s.minimap2.sort.bam"%(self.outputpath,self.samplename)
        cmd1="%s --input %s  -t %s --vcf %s/%s.minimap2.sniffles.SVs.vcf.gz"%(Sniffles_path,BAMfile,self.align_thread,outputdir,self.samplename)
        cmd2="zcat  %s/%s.minimap2.sniffles.SVs.vcf.gz |sed '/^##/d' >%s/%s.minimap2.sniffles.SVs.txt"%(outputdir,self.samplename,outputdir,self.samplename)
        cmdlist.append(cmd1)
        cmdlist.append(cmd2)
        writeShell("cycloneCallSV_%s"%self.samplename,cmdlist,shelldir)
        sh="sh %scycloneCallSV_%s.sh"%(shelldir,self.samplename)
        return sh
    def run_SVcheck(self,standard):
        shelldir="%s/shell/04.SVs/%s/"%(self.outputpath,self.samplename)
        outputdir="%s/result/04.SVs/"%(self.outputpath)
        Cmdruner("mkdir -p %s %s"%(shelldir,outputdir))
        cmdlist=[]
        Sniffles_file="%s/%s.minimap2.sniffles.SVs.vcf.gz"%(outputdir,self.samplename)
        truvari_path=getConfig("config","truvari","path")
        if self.fa_ref == "hg38":
            fa=getConfig("config","ref","hg38")
            if standard == "HG001":
                SVbaseline=getConfig("config","HG001","hg38_SVbaseline")
                cmd="%s bench -b %s -c %s -o %s/SV_BASE --reference %s  -r 1000 -p 0.00"%(truvari_path,SVbaseline,Sniffles_file,outputdir,fa)
            elif standard == "HG002":
                SVbaseline=getConfig("config","HG002","hg38_SVbaseline")
                SVbed=getConfig("config","HG002","hg38_SVbed")
                cmd="%s bench -b %s -c %s -o %s/SV_BASE --reference %s --includebed %s -r 1000 -p 0.00"%(truvari_path,SVbaseline,Sniffles_file,outputdir,fa,SVbed)
        elif self.fa_ref == "hg19":
            fa=getConfig("config","ref","hg19")
            if standard == "HG001":
                cmd=""
                pass
            elif standard == "HG002":
                SVbaseline=getConfig("config","HG002","hg19_SVbaseline")
                SVbed=getConfig("config","HG002","hg19_SVbed")
                cmd="%s bench -b %s -c %s -o %s/SV_BASE --reference %s --includebed %s -r 1000 -p 0.00 --passonly"%(truvari_path,SVbaseline,Sniffles_file,outputdir,fa,SVbed)
        cmd0="rm -rf %s/SV_BASE"%outputdir
        cmdlist.append(cmd0)
        cmdlist.append(cmd)
        writeShell("cycloneCallSVcheck_%s"%self.samplename,cmdlist,shelldir)
        sh="sh %scycloneCallSVcheck_%s.sh"%(shelldir,self.samplename)
        return sh

        
def report(samplename,output,ref_fa,standard):
    outputdir="%s/result/05.report/"%(output)
    Cmdruner("mkdir -p  %s"%(outputdir))
    #qc report
    qcfinshfile="%s/result/01.QC/%s.QC_info.xlsx"%(output,samplename)
    rtgtools_snp="%s/result/03.Variants/rtgtools_snp/summary.txt"%(output)
    rtgtools_indel="%s/result/03.Variants/rtgtools_indel/summary.txt"%(output)
    SV_checkfile="%s/result/04.SVs/SV_BASE/summary.json"%(output)
    if standard == "HG002":
        SVcheck=get_SVcheck(SV_checkfile)
        rtg_snp=get_rtg(rtgtools_snp)
        rtg_indel=get_rtg(rtgtools_indel)
        # rtg_snp_table="\n\t<h3>2.3.1 SNP位点评价</h3>\n\t<div class=\"all_t\" style=\"width: 90%;\">\n\t\t<p align=\"center\" style=\"font-size: 16px;\"><b>表4 SNP评估指标统计表</b></p>\n\t\t<table>\n\t\t<tr>\n\t\t\t<th>Threshold</th>\n\t\t\t<th>True pos baseline</th>\n\t\t\t<th>True pos call</th>\n\t\t\t<th>False pos</th>\n\t\t\t<th>False neg</th>\n\t\t\t<th>Precision</th>\n\t\t\t<th>Sensitivity</th>\n\t\t\t<th>F measure</th>\n\t\t</tr>\n\t\t<tr align='center'>"+"\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t</tr>\n\t\t</table>\n\t</div>\n\t<hr>"%(rtg_snp[0]["Threshold"],rtg_snp[0]["True_pos_baseline"],rtg_snp[0]["True_pos_call"],rtg_snp[0]["False_pos"],rtg_snp[0]["False_neg"],rtg_snp[0]["Precision"],rtg_snp[0]["Sensitivity"],rtg_snp[0]["F_measure"])
        # rtg_indel_table="\n\t<h3>2.3.2 InDEL位点评价</h3>\n\t<div class=\"all_t\" style=\"width: 90%;\">\n\t\t<p align=\"center\" style=\"font-size: 16px;\"><b>表5 InDel评估指标统计表</b></p>\n\t\t<table>\n\t\t<tr>\n\t\t\t<th>Threshold</th>\n\t\t\t<th>True pos baseline</th>\n\t\t\t<th>True pos call</th>\n\t\t\t<th>False pos</th>\n\t\t\t<th>False neg</th>\n\t\t\t<th>Precision</th>\n\t\t\t<th>Sensitivity</th>\n\t\t\t<th>F measure</th>\n\t\t</tr>\n\t\t<tr align='center'>"+"\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t</tr>\n\t\t</table>\n\t</div>"%(rtg_indel[0]["Threshold"],rtg_indel[0]["True_pos_baseline"],rtg_indel[0]["True_pos_call"],rtg_indel[0]["False_pos"],rtg_indel[0]["False_neg"],rtg_indel[0]["Precision"],rtg_indel[0]["Sensitivity"],rtg_indel[0]["F_measure"])
        rtg_snp_table="\n\t<h3>2.3.1 SNP位点评价</h3>\n\t<div class=\"all_t\" style=\"width: 90%;\">\n\t\t<p align=\"center\" style=\"font-size: 16px;\"><b>表4 SNP评估指标统计表</b></p>\n\t\t<table>\n\t\t<tr>\n\t\t\t<th>Threshold</th>\n\t\t\t<th>True pos baseline</th>\n\t\t\t<th>True pos call</th>\n\t\t\t<th>False pos</th>\n\t\t\t<th>False neg</th>\n\t\t\t<th>Precision</th>\n\t\t\t<th>Sensitivity</th>\n\t\t\t<th>F measure</th>\n\t\t</tr>\n\t\t<tr align='center'>"+"\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t</tr>\n\t\t<tr align='center'>"%(rtg_snp[0]["Threshold"],rtg_snp[0]["True_pos_baseline"],rtg_snp[0]["True_pos_call"],rtg_snp[0]["False_pos"],rtg_snp[0]["False_neg"],rtg_snp[0]["Precision"],rtg_snp[0]["Sensitivity"],rtg_snp[0]["F_measure"])+"\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t</tr>\n\t\t</table>\n\t</div>\n\t<hr>"%(rtg_snp[1]["Threshold"],rtg_snp[1]["True_pos_baseline"],rtg_snp[1]["True_pos_call"],rtg_snp[1 ]["False_pos"],rtg_snp[1]["False_neg"],rtg_snp[1]["Precision"],rtg_snp[1]["Sensitivity"],rtg_snp[1]["F_measure"])
        rtg_indel_table="\n\t<h3>2.3.2 InDEL位点评价</h3>\n\t<div class=\"all_t\" style=\"width: 90%;\">\n\t\t<p align=\"center\" style=\"font-size: 16px;\"><b>表5 InDel评估指标统计表</b></p>\n\t\t<table>\n\t\t<tr>\n\t\t\t<th>Threshold</th>\n\t\t\t<th>True pos baseline</th>\n\t\t\t<th>True pos call</th>\n\t\t\t<th>False pos</th>\n\t\t\t<th>False neg</th>\n\t\t\t<th>Precision</th>\n\t\t\t<th>Sensitivity</th>\n\t\t\t<th>F measure</th>\n\t\t</tr>\n\t\t<tr align='center'>"+"\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t</tr>\n\t\t<tr align='center'>"%(rtg_indel[0]["Threshold"],rtg_indel[0]["True_pos_baseline"],rtg_indel[0]["True_pos_call"],rtg_indel[0]["False_pos"],rtg_indel[0]["False_neg"],rtg_indel[0]["Precision"],rtg_indel[0]["Sensitivity"],rtg_indel[0]["F_measure"])+"\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t</tr>\n\t\t</table>\n\t</div>"%(rtg_indel[1]["Threshold"],rtg_indel[1]["True_pos_baseline"],rtg_indel[1]["True_pos_call"],rtg_indel[1]["False_pos"],rtg_indel[1]["False_neg"],rtg_indel[1]["Precision"],rtg_indel[1]["Sensitivity"],rtg_indel[1]["F_measure"])
        rtg_result=rtg_snp_table+rtg_indel_table
        SV_check="\n\t<h3>2.4.1 SV位点评价</h3>\n\t<div class=\"all_t\" style=\"width: 90%;\">\n\t\t<p align=\"center\" style=\"font-size: 16px;\"><b>表4.1 SV位点评估指标统计表</b></p>\n\t\t<table>\n\t\t<tr>\n\t\t\t<th>TP-base</th>\n\t\t\t<th>FP</th>\n\t\t\t<th>FN</th>\n\t\t\t<th>precision</th>\n\t\t\t<th>recall</th>\n\t\t\t<th>f1</th>\n\t\t\t<th>base cnt</th>\n\t\t\t<th>comp cnt</th>\n\t\t\t<th>TP-comp_TP-gt</th>\n\t\t\t<th>TP-comp_FP-gt</th>\n\t\t\t<th>gt_concordance</th>\n\t\t</tr>\n\t\t<tr align='center'>"+"\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%.4f</td>\n\t\t\t<td>%.4f</td>\n\t\t\t<td>%.4f</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%.4f</td>\n\t\t</tr>\n\t\t</table>\n\t</div>\n\t<hr>"%(SVcheck["TP_comp"],SVcheck["FP"],SVcheck["FN"],float(SVcheck["precision"]),float(SVcheck["recall"]),float(SVcheck["f1"]),SVcheck["base cnt"],SVcheck["comp cnt"],SVcheck["TP_comp_TP_gt"],SVcheck["TP_comp_FP_gt"],float(SVcheck["gt_concordance"]))
        standard = standard
    elif standard == "HG001":
        rtg_snp=get_rtg(rtgtools_snp)
        rtg_indel=get_rtg(rtgtools_indel)
        rtg_snp_table="\n\t<h3>2.3.1 SNP位点评价</h3>\n\t<div class=\"all_t\" style=\"width: 90%;\">\n\t\t<p align=\"center\" style=\"font-size: 16px;\"><b>表4 SNP评估指标统计表</b></p>\n\t\t<table>\n\t\t<tr>\n\t\t\t<th>Threshold</th>\n\t\t\t<th>True pos baseline</th>\n\t\t\t<th>True pos call</th>\n\t\t\t<th>False pos</th>\n\t\t\t<th>False neg</th>\n\t\t\t<th>Precision</th>\n\t\t\t<th>Sensitivity</th>\n\t\t\t<th>F measure</th>\n\t\t</tr>\n\t\t<tr align='center'>"+"\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t</tr>\n\t\t<tr align='center'>"%(rtg_snp[0]["Threshold"],rtg_snp[0]["True_pos_baseline"],rtg_snp[0]["True_pos_call"],rtg_snp[0]["False_pos"],rtg_snp[0]["False_neg"],rtg_snp[0]["Precision"],rtg_snp[0]["Sensitivity"],rtg_snp[0]["F_measure"])+"\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t</tr>\n\t\t</table>\n\t</div>\n\t<hr>"%(rtg_snp[1]["Threshold"],rtg_snp[1]["True_pos_baseline"],rtg_snp[1]["True_pos_call"],rtg_snp[1 ]["False_pos"],rtg_snp[1]["False_neg"],rtg_snp[1]["Precision"],rtg_snp[1]["Sensitivity"],rtg_snp[1]["F_measure"])
        rtg_indel_table="\n\t<h3>2.3.2 InDEL位点评价</h3>\n\t<div class=\"all_t\" style=\"width: 90%;\">\n\t\t<p align=\"center\" style=\"font-size: 16px;\"><b>表5 InDel评估指标统计表</b></p>\n\t\t<table>\n\t\t<tr>\n\t\t\t<th>Threshold</th>\n\t\t\t<th>True pos baseline</th>\n\t\t\t<th>True pos call</th>\n\t\t\t<th>False pos</th>\n\t\t\t<th>False neg</th>\n\t\t\t<th>Precision</th>\n\t\t\t<th>Sensitivity</th>\n\t\t\t<th>F measure</th>\n\t\t</tr>\n\t\t<tr align='center'>"+"\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t</tr>\n\t\t<tr align='center'>"%(rtg_indel[0]["Threshold"],rtg_indel[0]["True_pos_baseline"],rtg_indel[0]["True_pos_call"],rtg_indel[0]["False_pos"],rtg_indel[0]["False_neg"],rtg_indel[0]["Precision"],rtg_indel[0]["Sensitivity"],rtg_indel[0]["F_measure"])+"\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t</tr>\n\t\t</table>\n\t</div>"%(rtg_indel[1]["Threshold"],rtg_indel[1]["True_pos_baseline"],rtg_indel[1]["True_pos_call"],rtg_indel[1]["False_pos"],rtg_indel[1]["False_neg"],rtg_indel[1]["Precision"],rtg_indel[1]["Sensitivity"],rtg_indel[1]["F_measure"])
        rtg_result=rtg_snp_table+rtg_indel_table
        standard = standard
        if ref_fa == "hg19":
            SV_check = ""
        elif  ref_fa == "hg38":
            SVcheck=get_SVcheck(SV_checkfile)
            SV_check="\n\t<h3>2.4.1 SV位点评价</h3>\n\t<div class=\"all_t\" style=\"width: 90%;\">\n\t\t<p align=\"center\" style=\"font-size: 16px;\"><b>表4.1 SV位点评估指标统计表</b></p>\n\t\t<table>\n\t\t<tr>\n\t\t\t<th>TP-base</th>\n\t\t\t<th>FP</th>\n\t\t\t<th>FN</th>\n\t\t\t<th>precision</th>\n\t\t\t<th>recall</th>\n\t\t\t<th>f1</th>\n\t\t\t<th>base cnt</th>\n\t\t\t<th>comp cnt</th>\n\t\t\t<th>TP-comp_TP-gt</th>\n\t\t\t<th>TP-comp_FP-gt</th>\n\t\t\t<th>gt_concordance</th>\n\t\t</tr>\n\t\t<tr align='center'>"+"\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%.4f</td>\n\t\t\t<td>%.4f</td>\n\t\t\t<td>%.4f</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%.4f</td>\n\t\t</tr>\n\t\t</table>\n\t</div>\n\t<hr>"%(SVcheck["TP_comp"],SVcheck["FP"],SVcheck["FN"],float(SVcheck["precision"]),float(SVcheck["recall"]),float(SVcheck["f1"]),SVcheck["base cnt"],SVcheck["comp cnt"],SVcheck["TP_comp_TP_gt"],SVcheck["TP_comp_FP_gt"],float(SVcheck["gt_concordance"]))
            SV_check = SV_check
        else:
            SV_check = ""
    else:
        rtg_result=""
        standard = "非标准品"
        pass
    
    Cmdruner("cp %s %s"%(qcfinshfile,outputdir))
    #mapping report
    mappingfinshfile="%s/result/02.Mapping/%s.minimap2.sample_bamstat.txt"%(output,samplename)
    Cmdruner("cp %s %s"%(mappingfinshfile,outputdir))
    #variant report
    Varfinshfile="%s/result/03.Variants/%s.variant.stat.txt"%(output,samplename)
    Cmdruner("cp %s %s"%(Varfinshfile,outputdir))
    SVfinshfile="%s/result/04.SVs/%s.minimap2.sniffles.SVs.txt"%(output,samplename)
    Cmdruner("cp %s %s"%(SVfinshfile,outputdir))
     
    qc,mapping,var,sv=report_tmp(samplename,output)
    CyclonRe_dir=getConfig("config","CyclonRe","path")
    env = Environment(loader=FileSystemLoader(CyclonRe_dir))
    template = env.get_template('report1.html')
    os.chdir(outputdir)
    outreport="%s/%s_Cyclone_WGS_report.html"%(outputdir,samplename)
    results = {}
    results.update({'sample_name': samplename,
                    'sample_sp':ref_fa,
                    'sample_standard':standard,
                    'sample_rtgresult':rtg_result,
                    'sample_SVcheck':SV_check,
                        
                    'itemQCs': qc,
                    'itemmappings': mapping,
                    'itemvars':var,
                    'itemSVs':sv
                    })
    with open(outreport, 'w+', encoding='utf-8') as f:
        out = template.render(
                            sample_name=results['sample_name'],
                            sample_sp=results['sample_sp'],
                            sample_rtgresult=results['sample_rtgresult'],
                            sample_standard=results['sample_standard'],
                            sample_SVcheck=results['sample_SVcheck'],

                            itemQCs = results['itemQCs'],
                            itemmappings = results['itemmappings'],
                            itemvars = results['itemvars'],
                            itemSVs = results['itemSVs']
                            )
        f.write(out)
        f.close()
def report_nortrg(samplename,output,ref_fa,standard):
    outputdir="%s/result/05.report/"%(output)
    Cmdruner("mkdir -p  %s"%(outputdir))
    #qc report
    qcfinshfile="%s/result/01.QC/%s.QC_info.xlsx"%(output,samplename)
    SV_checkfile="%s/result/04.SVs/SV_BASE/summary.json"%(output)
    if standard == "HG002":
        SVcheck=get_SVcheck(SV_checkfile)
        SV_check="\n\t<h3>2.3.1 SV位点评价</h3>\n\t<div class=\"all_t\" style=\"width: 90%;\">\n\t\t<p align=\"center\" style=\"font-size: 16px;\"><b>表3.1 SV位点评估指标统计表</b></p>\n\t\t<table>\n\t\t<tr>\n\t\t\t<th>TP-base</th>\n\t\t\t<th>FP</th>\n\t\t\t<th>FN</th>\n\t\t\t<th>precision</th>\n\t\t\t<th>recall</th>\n\t\t\t<th>f1</th>\n\t\t\t<th>base cnt</th>\n\t\t\t<th>comp cnt</th>\n\t\t\t<th>TP-comp_TP-gt</th>\n\t\t\t<th>TP-comp_FP-gt</th>\n\t\t\t<th>gt_concordance</th>\n\t\t</tr>\n\t\t<tr align='center'>"+"\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%.4f</td>\n\t\t\t<td>%.4f</td>\n\t\t\t<td>%.4f</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%.4f</td>\n\t\t</tr>\n\t\t</table>\n\t</div>\n\t<hr>"%(SVcheck["TP_comp"],SVcheck["FP"],SVcheck["FN"],float(SVcheck["precision"]),float(SVcheck["recall"]),float(SVcheck["f1"]),SVcheck["base cnt"],SVcheck["comp cnt"],SVcheck["TP_comp_TP_gt"],SVcheck["TP_comp_FP_gt"],float(SVcheck["gt_concordance"]))
        standard = standard
    elif standard == "HG001":
        standard = standard
        if ref_fa == "hg19":
            SV_check = ""
        elif  ref_fa == "hg38":
            SVcheck=get_SVcheck(SV_checkfile)
            SV_check="\n\t<h3>2.3.1 SV位点评价</h3>\n\t<div class=\"all_t\" style=\"width: 90%;\">\n\t\t<p align=\"center\" style=\"font-size: 16px;\"><b>表3.1 SV位点评估指标统计表</b></p>\n\t\t<table>\n\t\t<tr>\n\t\t\t<th>TP-base</th>\n\t\t\t<th>FP</th>\n\t\t\t<th>FN</th>\n\t\t\t<th>precision</th>\n\t\t\t<th>recall</th>\n\t\t\t<th>f1</th>\n\t\t\t<th>base cnt</th>\n\t\t\t<th>comp cnt</th>\n\t\t\t<th>TP-comp_TP-gt</th>\n\t\t\t<th>TP-comp_FP-gt</th>\n\t\t\t<th>gt_concordance</th>\n\t\t</tr>\n\t\t<tr align='center'>"+"\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%.4f</td>\n\t\t\t<td>%.4f</td>\n\t\t\t<td>%.4f</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%s</td>\n\t\t\t<td>%.4f</td>\n\t\t</tr>\n\t\t</table>\n\t</div>\n\t<hr>"%(SVcheck["TP_comp"],SVcheck["FP"],SVcheck["FN"],float(SVcheck["precision"]),float(SVcheck["recall"]),float(SVcheck["f1"]),SVcheck["base cnt"],SVcheck["comp cnt"],SVcheck["TP_comp_TP_gt"],SVcheck["TP_comp_FP_gt"],float(SVcheck["gt_concordance"]))
            SV_check = SV_check
        else:
            SV_check = ""
    else:
        SV_check = ""
        standard = "非标准品"
        pass
    
    
    Cmdruner("cp %s %s"%(qcfinshfile,outputdir))
    #mapping report
    mappingfinshfile="%s/result/02.Mapping/%s.minimap2.sample_bamstat.txt"%(output,samplename)
    Cmdruner("cp %s %s"%(mappingfinshfile,outputdir))
    SVfinshfile="%s/result/04.SVs/%s.minimap2.sniffles.SVs.txt"%(output,samplename)
    Cmdruner("cp %s %s"%(SVfinshfile,outputdir))
     
    qc,mapping,sv=report_tmp_var(samplename,output)
    CyclonRe_dir=getConfig("config","CyclonRe","path")
    env = Environment(loader=FileSystemLoader(CyclonRe_dir))
    template = env.get_template('report2.html')
    os.chdir(outputdir)
    outreport="%s/%s_Cyclone_WGS_report.html"%(outputdir,samplename)
    results = {}
    results.update({'sample_name': samplename,
                    'sample_sp':ref_fa,
                    'sample_standard':standard,
                    'sample_SVcheck':SV_check,
                        
                    'itemQCs': qc,
                    'itemmappings': mapping,
                    'itemSVs':sv
                    })
    with open(outreport, 'w+', encoding='utf-8') as f:
        out = template.render(
                            sample_name=results['sample_name'],
                            sample_sp=results['sample_sp'],
                            sample_standard=results['sample_standard'],
                            sample_SVcheck=results['sample_SVcheck'],

                            itemQCs = results['itemQCs'],
                            itemmappings = results['itemmappings'],
                            itemSVs = results['itemSVs']
                            )
        f.write(out)
        f.close()
      