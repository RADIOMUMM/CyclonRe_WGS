# CyclonRe_WGS
CyclonRe-WGS software is dedicated to MGI-cyclone/ONT triple sequencing technology, the software analyzes human-sequencing data. Functions include data quality control, mapping, small variant calling, structural variant calling, and variant evaluation.
! Introduction

CyclonRe-WGS软件是专用于MGI-cyclone/ONT三代测序技术的分析软件，软件对人的重测序数据进行分析。功能包括数据质控、比对、小变异检测、结构变异检测、变异评估等。

CyclonRe-WGS software is dedicated to MGI-cyclone/ONT triple sequencing technology, the software analyzes human sequencing data. Functions include data quality control, mapping, small variant detection, structural variant detection, and variant evaluation.



软件流程图

<div align=center><img src="https://pic.imgdb.cn/item/660bdcb69f345e8d039464a8.png" alt="图片alt" title="图2"></div>


---
! Getting Started



---

! Installation

CyclonRe-WGS软件推荐使用conda安装。

```Shel
conda env create -f CyclonRe-WGS_environment.yml
```

!! requirement

| ''software'' | ''version'' | ''link'' |
|NanoStat |1.6.0 |https://github.com/wdecoster/nanostat |
|NanoFilt |2.8.0 |https://github.com/wdecoster/nanofilt |
|NanoPlot |1.42.0 |https://github.com/wdecoster/NanoPlot |
|Sniffles2 |2.2 |https://github.com/fritzsedlazeck/Sniffles |
|minimap2 |2.17 |https://github.com/lh3/minimap2 |
|samtools |1.7 |https://github.com/samtools/samtools |
|clair3 |1.5 |https://github.com/HKU-BAL/Clair3 |
|pigz |2.7 |https://github.com/madler/pigz |
|bcftools |1.5 |https://github.com/samtools/bcftools |
|tabix |0.2.6 |https://github.com/samtools/tabix |
|truvari |1.5 |https://github.com/ACEnglish/truvari |
|gatk |4.2 |https://github.com/RealTimeGenomics/rtg-tools |
|rtg-tools |3.11 |https://github.com/RealTimeGenomics/rtg-tools |


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
