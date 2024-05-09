# CyclonRe_WGS
<hr>
<h2> Introduction </h2>
CyclonRe-WGS software is dedicated to MGI-cyclone/ONT triple sequencing technology, the software analyzes human sequencing data. Functions include data quality control, mapping, small variant detection, structural variant detection, and variant evaluation.


<div align=center><img src="https://pic.imgdb.cn/item/662a2a6b0ea9cb14037ebad5.jpg" alt="图片alt" title="图2"></div>

---
<h2> Getting Started </h2>

---

<h2> Installation </h2>

CyclonRe-WGS软件推荐使用conda安装。

```Shel
conda env create -f CyclonRe-WGS_environment.yml
```

<h3> requirement </h3>
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
<h2>Performance</h2> 

<h3>Small variant evaluation</h3>

软件对Cyclone下机的HG002数据进行变异检测，将小变异（SNP/InDel）与giab数据库 (https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/) 中的baseline进行比对，得到SNP/InDel评估值。
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

