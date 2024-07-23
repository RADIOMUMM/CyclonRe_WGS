# CyclonRe_WGS

## 1 Introduction

CyclonRe-WGS software is dedicated to MGI-cyclone/ONT triple sequencing technology, the software analyzes human sequencing data. Functions include data quality control, mapping, small variant detection, structural variant detection, and variant evaluation.


## 1.1 Software Workflow

<div align=center><img src="https://pic.imgdb.cn/item/669f5055d9c307b7e9bac487.png" alt="图片alt" title="图2" width="50%"></div>

## 1.2 Hardware/Software requirements

- x86-64 compatible processors.
- require at least 50GB of RAM and 4 CPU.
- centos 7.x 64-bit operating system (Linux kernel 3.10.0, compatible with higher software and hardware configuration).

# 2 Getting Started

## 2.1 Installation

The CyclonRe-WGS software is recommended to be installed using conda.

```
unzip CyclonRe_WGS_V1.0.1.zip
cd CyclonRe_WGS_V1.0.1
conda env create -f CyclonRe_WGS_environment.yml
```
**requirement**

<table>
    <tr>
        <td>software</td> 
        <td>version</td> 
        <td>link</td> 
   </tr>
    <tr>
  		<td>NanoStat</td> 
        <td>1.6.0</td> 
        <td>https://github.com/wdecoster/nanostat</td> 
    </tr>
    <tr>
        <td>NanoFilt</td> 
        <td>2.8.0</td> 
        <td>https://github.com/wdecoster/nanofilt</td> 
    </tr>
    <tr>
        <td>NanoPlot</td> 
        <td>1.42.0</td> 
        <td>https://github.com/wdecoster/NanoPlot</td> 
    </tr>
    <tr>
        <td>Sniffles2</td> 
        <td>2.2</td> 
        <td>https://github.com/fritzsedlazeck/Sniffles</td> 
    </tr>
    <tr>
        <td>minimap2</td> 
        <td>2.17</td> 
        <td>https://github.com/lh3/minimap2</td> 
    </tr>
    <tr>
        <td>samtools</td> 
        <td>1.7</td> 
        <td>https://github.com/samtools/samtools</td> 
    </tr>
    <tr>
        <td>clair3</td> 
        <td>1.5</td> 
        <td>https://github.com/HKU-BAL/Clair3</td> 
    </tr>
    <tr>
        <td>pigz</td> 
        <td>2.7</td> 
        <td>https://github.com/madler/pigz</td> 
    </tr>
    <tr>
        <td>bcftools</td> 
        <td>1.5</td> 
        <td>https://github.com/samtools/bcftools</td> 
    </tr>
    <tr>
        <td>tabix</td> 
        <td>0.2.6</td> 
        <td>https://github.com/samtools/tabix</td> 
    </tr>
    <tr>
        <td>truvari</td> 
        <td>1.5</td> 
        <td>https://github.com/ACEnglish/truvari</td> 
    </tr>
    <tr>
        <td>gatk</td> 
        <td>4.2</td> 
        <td>https://github.com/broadinstitute/gatk</td> 
    </tr>
    <tr>
        <td>rtg-tools</td> 
        <td>3.11</td> 
        <td>https://github.com/RealTimeGenomics/rtg-tools</td> 
    </tr>
</table>

## 2.2 Preparation

### 2.2.1 Download reference genome

The software supports WGS analysis of 3 reference genomes (hg19/hg38/ecoli), please download the reference genomes in advance.

```
##download hg38 
wget https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
##download hg19
wget https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.toplevel.fa.gz
```

**Build index for reference genome**

For example:

```
samtools faidx  Homo_sapiens.GRCh38.dna.toplevel.fa.gz
gatk CreateSequenceDictionary -R Homo_sapiens.GRCh38.dna.toplevel.fa.gz
```

### 2.2.2 Download standard baseline

While downloading the vcf.gz file, you need to download the corresponding bed and tbi files. 

**Download HG001 baseline**

```
##download HG001 small-variants baseline
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh37/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz

```

**Download HG002 baseline**

```
##download HG002 small-variants baseline
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh37/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
##download HG002 SV baseline
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.012-20231107/GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz
```

### 2.2.3 Modify the config file

Replace the software and file paths in the config file with local paths.

```
vim config
```
# 3 RUN

**Usage**

```
usage: run_cyclone_work -i <fq>  -n <sample_name>  -o <outputFile>

cyclone fq QC workflow

optional arguments:
  -h, --help            show this help message and exit
  -i FASTQ, --fastq FASTQ
                        input fsatq
  -n SAMPLE_NAME, --sample_name SAMPLE_NAME
                        sample name
  -o OUTPUTFILE, --outputFile OUTPUTFILE
                        output file path,Note that you need to add the "/" symbol at the end of the
                        folder,eg:'/home/','/home/result/'
  -callSNP_indel {None,True}, --callsmallVariants {None,True}
                        The process can autonomously choose whether or not to perform SNP and InDel testing
                        and evaluation.Choose between None or True.(optional)
  -t QCTHREAD, --qcthread QCTHREAD
                        Threads required for the QC program.(optional)
  -topcrop TOPCROP, --topcrop TOPCROP
                        Trim n nucleotides from start of read.(optional)
  -tailcrop TAILCROP, --tailcrop TAILCROP
                        Trim n nucleotides from end of read.(optional)
  -q TRIM_QUALITY, --trim_quality TRIM_QUALITY
                        Filter on a minimum average read quality score.(optional)
  -qctype {NanoPlot,NanoStat}, --QCtype {NanoPlot,NanoStat}
                        Select the software used in the QC step. Choose between NanoPlot or
                        NanoStat.(optional)
  -standard STANDARDTYPE, --standardtype STANDARDTYPE
                        Select the standard sample file . Choose between HG002 or HG001.(optional)
  -ref REF, --REF REF   Select the reference genome used in this analysis. Choose between hg19 or
                        hg38.(optional)
  -t2 MAPPINGTHREAD, --mappingthread MAPPINGTHREAD
                        Threads required for the Mapping program.(optional)
  -t3 SAMTOOLSTHREAD, --samtoolsthread SAMTOOLSTHREAD
                        Threads required for the samtools program.(optional)
  -t4 CALLVARIANTTHREAD, --Callvariantthread CALLVARIANTTHREAD
                        Threads required for the Clair3 program.(optional)

```

**Example: Use the example data to initiate the analysis.**

- Performs detection of small-variants
- Data quality requirements Q15
```
bin/run_cyclone_workflow 
-i /test/TB2000B73B-202403261655301_read.cut.fq.gz 
-n demo_test 
-o /result/ 
-callSNP_indel True 
-standard HG002 
-ref hg19
-q 15
```

# 4 Performance

<h3>Small variant evaluation</h3>

The software uses the HG002 Cyclone data  for variant detection, comparing small variants (SNP/InDel) to baseline in the giab database (https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/) to obtain SNP /InDel evaluation values.

<h4>SNPs evaluation</h4>
<table>
    <tr>
        <td>coverage</td> 
        <td>SNP number</td> 
        <td>Precision</td> 
        <td>Sensitivity</td>
        <td>F-measure</td>
   </tr>
    <tr>
        <td>10X</td> 
        <td>3850830</td> 
        <td>0.9719</td> 
        <td>0.9258</td>
        <td>0.9483</td>
   </tr>
    <tr>
        <td>15X</td> 
        <td>3923097</td> 
        <td>0.9761</td> 
        <td>0.9485</td>
        <td>0.9621</td>
   </tr>
    <tr>
        <td>30X</td> 
        <td>4011627</td> 
        <td>0.9888</td> 
        <td>0.9807</td>
        <td>0.9847</td>
   </tr>
</table>

<h4>InDel evaluation</h4>
<table>
    <tr>
        <td>coverage</td> 
        <td>InDel number</td> 
        <td>Precision</td> 
        <td>Sensitivity</td>
        <td>F-measure</td>
   </tr>
    <tr>
        <td>10X</td> 
        <td>710718</td> 
        <td>0.7175</td> 
        <td>0.3746</td>
        <td>0.4922</td>
   </tr>
    <tr>
        <td>15X</td> 
        <td>841515</td> 
        <td>0.6596</td> 
        <td>0.4147</td>
        <td>0.5092</td>
   </tr>
    <tr>
        <td>30X</td> 
        <td>909464</td> 
        <td>0.6844</td> 
        <td>0.5005</td>
        <td>0.5782</td>
   </tr>
</table>

<h3>SVs evaluation</h3>
<table>
    <tr>
        <td>coverage</td> 
        <td>SVs number</td> 
        <td>Precision</td> 
        <td>recall</td>
        <td>f1</td>
   </tr>
    <tr>
        <td>10X</td> 
        <td>35866</td> 
        <td>0.7114</td> 
        <td>0.905</td>
        <td>0.7966</td>
   </tr>
    <tr>
        <td>15X</td> 
        <td>37366</td> 
        <td>0.8022</td> 
        <td>0.9265</td>
        <td>0.8599</td>
   </tr>
    <tr>
        <td>30X</td> 
        <td>32681</td> 
        <td>0.9018</td> 
        <td>0.9467</td>
        <td>0.9237</td>
   </tr>
</table>
