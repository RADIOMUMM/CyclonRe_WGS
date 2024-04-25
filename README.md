# CyclonRe_WGS

---

<h2> Introduction </h2>
CyclonRe-WGS software is dedicated to MGI-cyclone/ONT triple sequencing technology, the software analyzes human sequencing data. Functions include data quality control, mapping, small variant detection, structural variant detection, and variant evaluation.

<hr>
<div align=center><img src="https://pic.imgdb.cn/item/662a2a6b0ea9cb14037ebad5.jpg" alt="图片alt" title="图2"></div>

---
! Getting Started

---

! Installation

CyclonRe-WGS软件推荐使用conda安装。

```Shel
conda env create -f CyclonRe-WGS_environment.yml
```

!! requirement
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
        <td>https://github.com/RealTimeGenomics/rtg-tools</td> 
    </tr>
    <tr>
        <td>rtg-tools</td> 
        <td>3.11</td> 
        <td>https://github.com/RealTimeGenomics/rtg-tools</td> 
    </tr>
</table>

---
! Performance

!! Small variant evaluation

软件对Cyclone下机的HG002数据进行变异检测，将小变异（SNP/InDel）与giab数据库 (https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/) 中的baseline进行比对，得到SNP/InDel评估值。

| ''SNP'' |
| ''coverage'' | ''SNP number'' | ''Precision'' | ''Sensitivity'' | ''F-measure'' |
| 10X | 3850830 | 0.9719 | 0.9258 | 0.9483 |
| 15X | 3923097 | 0.9761 | 0.9485 | 0.9621 |
| 30X | 4011627 | ''<font color="Salmon">0.9888</font>'' | ''<font color="Salmon">0.9807</font>'' | ''<font color="Salmon">0.9847</font>'' |


| ''InDel'' |
| ''coverage'' | ''SNP number'' | ''Precision'' | ''Sensitivity'' | ''F-measure'' |
| 10X | 710718 | ''<font color="Salmon">0.7175</font>'' | 0.3746 | 0.4922 |
| 15X | 841515 | 0.6596 | 0.4147 | 0.5092 |
| 30X | 909464 | 0.6844 | ''<font color="Salmon">0.5005</font>'' | ''<font color="Salmon">0.5782</font>'' |
