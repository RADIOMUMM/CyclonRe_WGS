import json
import re
import pandas as pd

def longqc_report(samplename,output):
    with open("%s/result/01.QC/%s/QC_vals_longQC_sampleqc.json"%(output,samplename),"r")as longqc:
        json_data=json.load(longqc)
        merge=[]
        merge.append("Sample\t%s\n"%samplename)
        merge.append("Total Bases\t%s\n"%json_data["Yield"])
        merge.append("Num of reads\t%s\n"%json_data["Num_of_reads"])
        merge.append("Q7 bases\t%s\n"%json_data["Q7 bases"])
        merge.append("Longest read\t%s\n"%json_data["Longest_read"])
        merge.append("Length_stats\ngamma_params\t%s\n"%json_data["Length_stats"]["gamma_params"])
        merge.append("Mean read length\t%s\n"%json_data["Length_stats"]["Mean_read_length"])
        merge.append("N50 read length\t%s\n"%json_data["Length_stats"]["N50_read_length"])
        merge.append("GC_stats\nMean_GC_content\t%s\n"%json_data["GC_stats"]["Mean_GC_content"])
        merge.append("SD_GC_content\t%s\n"%json_data["GC_stats"]["SD_GC_content"])
        # merge.append("Stats_for_adapter5\nNum_of_trimmed_reads_5\t%s\n"%json_data["Stats_for_adapter5"]["Num_of_trimmed_reads_5"])
        # merge.append("Max_identity_adp5\t%s\n"%json_data["Stats_for_adapter5"]["Max_identity_adp5"])
        # merge.append("Average_position_from_5_end\t%s\n"%json_data["Stats_for_adapter5"]["Average_position_from_5_end"])
        merge.append("Stats_for_adapter3\nNum_of_trimmed_reads_3\t%s\n"%json_data["Stats_for_adapter3"]["Num_of_trimmed_reads_3"])
        merge.append("Max_identity_adp3\t%s\n"%json_data["Stats_for_adapter3"]["Max_identity_adp3"])
        merge.append("Average_position_from_3_end\t%s\n"%json_data["Stats_for_adapter3"]["Average_position_from_3_end"])
        merge.append("Coverage_stats\nEstimated non-sense read fraction\t%s\n"%json_data["Coverage_stats"]["Estimated non-sense read fraction"])
        merge.append("Mean_coverage\t%s\n"%json_data["Coverage_stats"]["Mean_coverage"])
        merge.append("SD_coverage\t%s\n"%json_data["Coverage_stats"]["SD_coverage"])
        merge.append("Estimated crude Xome size\t%s\n"%json_data["Coverage_stats"]["Estimated crude Xome size"])
    qcfile="%s/result/01.QC/%s/cyclone_qc.txt"%(output,samplename)
    with open(qcfile,"w")as www:
        for i in merge:
            www.write(i)
    return qcfile

def nanoplot_report(samplename,output):
    with open("%s/result/01.QC/%s/NanoStats.txt"%(output,samplename),"r")as NanoPlot:
        content=NanoPlot.readlines()
        content.remove(content[0])
        newlist=[" \t%s\n"%samplename]
        for i in content:
            i=re.sub(": +","\t",i)
            i=re.sub(":\t","\t",i)
            newlist.append(i)
    qcfile="%s/result/01.QC/%s/cyclone_qc.txt"%(output,samplename)
    with open(qcfile,"w")as www:
        for i in newlist:
            www.write(i)
    return qcfile

def merge_report(samplename,output,input1,input2):
    list1=[input1,input2]
    df_list=[]
    for excel_name in list1:
        #获取每个excel的路径名称
        # print(excel_name)
        df_split=pd.read_csv(excel_name,sep="\t")
        df_list.append(df_split)
        pass
    df_concat=pd.concat(df_list,axis=1)
    df_concat= df_concat.T.drop_duplicates().T
    df_concat.to_excel("%s/result/01.QC/%s.QC_info.xlsx"%(output,samplename),index=False)

def get_rtg(iputfile):
    with open(iputfile,"r") as Var_rgt_stat:
        content=Var_rgt_stat.readlines()
        title=re.sub("\n","",content[0])
        title=re.sub("-","_",title)
        title=title.split("  ")
        con1=re.sub("^\s+|\n","",content[2])
        con1=re.sub(" +",",",con1).split(",")
        con2=re.sub("^\s+|\n","",content[3])
        con2=re.sub(" +",",",con2).split(",")
        dict1={}
        dict2={}
        for i in range(0,len(title)):
            dict1[title[i]]=con1[i]
            dict2[title[i]]=con2[i]
    return [dict1,dict2]

def report_tmp(samplename,output):
    qcfinshfile="%s/result/01.QC/%s.QC_info.xlsx"%(output,samplename)
    mappingfinshfile="%s/result/02.Mapping/%s.minimap2.sample_bamstat.txt"%(output,samplename)
    Varfinshfile="%s/result/03.Variants/%s.variant.stat.txt"%(output,samplename)
    VSfinshfile="%s/result/04.SVs/%s.minimap2.sniffles.SVs.txt"%(output,samplename)
    qc_Data=pd.read_excel(qcfinshfile)
    qc_Data.columns=["Sample","Rawdata","Cleandata"]
    data1 = qc_Data.to_dict('records')
        
    with open(mappingfinshfile,"r") as bamstat:
        content=bamstat.read().splitlines()
        datalit=[]
        for i in content:
            if re.match("SN",i):
                strinfo1 = re.compile(r'\t# .+')
                strinfo2 = re.compile(r'SN\t')
                i=strinfo1.sub('', i)
                i=strinfo2.sub('', i)
                i=str(i).replace(" ","_")
                i=str(i).replace("_(%)","")
                i=str(i).split(":\t")
                datalit.append(i[0])
                datalit.append(i[1])
    coveragedict={}
    for num in range(0,len(datalit),2):
        coveragedict[datalit[num]] = datalit[num+1]
    errorrate = coveragedict['error_rate']
    mapped_rate = float(coveragedict['reads_mapped'])/float(coveragedict['raw_total_sequences'])
    coveragedict['mapped_rate'] = round(float(mapped_rate)*100,1)
    coveragedict['error_rate'] = round(float(errorrate)*100,1)
    mappingqcdict = [coveragedict]
    
    with open(Varfinshfile,"r") as Varstat:
        content=Varstat.read().splitlines()
        datalit=[]
        for i in content:
            if re.match("SN",i):
                strinfo1 = re.compile(r'\t# .+')
                strinfo2 = re.compile(r'SN\t0\t')
                i=strinfo1.sub('', i)
                i=strinfo2.sub('', i)
                i=str(i).replace(" ","_")
                i=str(i).replace("_(%)","")
                i=str(i).split(":\t")
                datalit.append(i[0])
                datalit.append(i[1])
            elif re.match("TSTV",i):
                strinfo2 = re.compile(r'TSTV\t0\t')
                i=strinfo2.sub('', i)
                i=str(i).split(":\t")
                i=str(i[0]).split("\t")
                datalit.append("Ti")
                datalit.append(i[0])
                datalit.append("Tv")
                datalit.append(i[1])
                datalit.append("TiTv")
                datalit.append(i[2])       
    coveragedict={}
    for num in range(0,len(datalit),2):
        coveragedict[datalit[num]]=datalit[num+1]
    coveragedict["number_of_no_ALTs"]=coveragedict["number_of_no-ALTs"]
    Varqcdict=[coveragedict]
    
    dataSV=pd.read_csv(VSfinshfile,sep="\t")
    # dataSV.loc[:,"LEN"]=dataSV["INFO"].str.split(";",expand=True)[2].str.replace("SVLEN=", "").astype("int32").abs()
    dataSV['LEN'] = dataSV['INFO'].str.extract(r'(SVLEN=.{1}\d+)')
    dataSV = dataSV[dataSV["LEN"].str.contains("NaN") == False]
    dataSV.loc[:,"LEN"]=dataSV["LEN"].str.replace("SVLEN=","").astype("int32").abs()
    dataSV['SV_type'] = dataSV['INFO'].str.extract(r'(SVTYPE=[A-Z]+)')
    dataSV.loc[:,"SV_type"]=dataSV["SV_type"].str.replace("SVTYPE=","")
    # dataSV.loc[:,"SV_type"]=dataSV["ID"].str.split(".",expand=True)[1]
    dataSV.loc[:,"GT"]=dataSV["SAMPLE"].str.split(":",expand=True)[0]
    def get_GT(df):
        if df["GT"] == "1/1":
            return "ALT"
        elif df["GT"] == "0/0":
            return "REF"
        elif df["GT"] == "0/1":
            return "Hybrid"
        else:
            return "Other_GT"
    dataSV.loc[:,"GT_type"]=dataSV.apply(get_GT,axis=1)
    x=dataSV["GT_type"].value_counts()
    def get_LEN(df):
        if df["LEN"] >= 50:
            return "more50bp"
        else:
            return "less50bp"
    dataSV.loc[:,"LEN_type"]=dataSV.apply(get_LEN,axis=1)
    y=dataSV["LEN_type"].value_counts()
    def get_SVtype(df):
        if df["SV_type"] == "DUP":
            return "DUP"
        elif df["SV_type"] == "INS":
            return "INS"
        elif df["SV_type"] == "DEL":
            return "DEL"
        elif df["SV_type"] == "INV":
            return "INV"
        else:
            return df["SV_type"]
    dataSV.loc[:,"SV"]=dataSV.apply(get_SVtype,axis=1)
    z=dataSV["SV"].value_counts()
    SVqcdict=pd.DataFrame(x).T.to_dict("records")
    if "ALT" in SVqcdict[0].keys():
        pass
    else:
        SVqcdict[0]["ALT"]=0
    if "REF" in SVqcdict[0].keys():
        pass
    else:
        SVqcdict[0]["REF"]=0
    if "Hybrid" in SVqcdict[0].keys():
        pass
    else:
        SVqcdict[0]["Hybrid"]=0             
    y=pd.DataFrame(y).T.to_dict("records")
    if 'more50bp' in y[0].keys():
        pass
    else:
        y[0]['more50bp']=0
    if 'less50bp' in y[0].keys():
        pass
    else:
        y[0]['less50bp']=0
    z=pd.DataFrame(z).T.to_dict("records")
    if 'DEL' in z[0].keys():
        pass
    else:
        z[0]['DEL']=0
    if 'INS' in z[0].keys():
        pass
    else:
        z[0]['INS']=0
    if 'INV' in z[0].keys():
        pass
    else:
        z[0]['INV']=0
    if 'DUP' in z[0].keys():
        pass
    else:
        z[0]['DUP']=0
    SVqcdict[0].update(y[0])
    SVqcdict[0].update(z[0])

    return data1,mappingqcdict,Varqcdict,SVqcdict

def report_tmp_var(samplename,output):
    qcfinshfile="%s/result/01.QC/%s.QC_info.xlsx"%(output,samplename)
    mappingfinshfile="%s/result/02.Mapping/%s.minimap2.sample_bamstat.txt"%(output,samplename)
    VSfinshfile="%s/result/04.SVs/%s.minimap2.sniffles.SVs.txt"%(output,samplename)
    qc_Data=pd.read_excel(qcfinshfile)
    qc_Data.columns=["Sample","Rawdata","Cleandata"]
    data1 = qc_Data.to_dict('records')
        
    with open(mappingfinshfile,"r") as bamstat:
        content=bamstat.read().splitlines()
        datalit=[]
        for i in content:
            if re.match("SN",i):
                strinfo1 = re.compile(r'\t# .+')
                strinfo2 = re.compile(r'SN\t')
                i=strinfo1.sub('', i)
                i=strinfo2.sub('', i)
                i=str(i).replace(" ","_")
                i=str(i).replace("_(%)","")
                i=str(i).split(":\t")
                datalit.append(i[0])
                datalit.append(i[1])
    coveragedict={}
    for num in range(0,len(datalit),2):
        coveragedict[datalit[num]] = datalit[num+1]
    errorrate = coveragedict['error_rate']
    mapped_rate = float(coveragedict['reads_mapped'])/float(coveragedict['raw_total_sequences'])
    coveragedict['mapped_rate'] = round(float(mapped_rate)*100,1)
    coveragedict['error_rate'] = round(float(errorrate)*100,1)
    mappingqcdict = [coveragedict]
    
    dataSV=pd.read_csv(VSfinshfile,sep="\t")
    # dataSV.loc[:,"LEN"]=dataSV["INFO"].str.split(";",expand=True)[2].str.replace("SVLEN=", "").astype("int32").abs()
    dataSV['LEN'] = dataSV['INFO'].str.extract(r'(SVLEN=.{1}\d+)')
    dataSV = dataSV[dataSV["LEN"].str.contains("NaN") == False]
    dataSV.loc[:,"LEN"]=dataSV["LEN"].str.replace("SVLEN=","").astype("int32").abs()
    dataSV['SV_type'] = dataSV['INFO'].str.extract(r'(SVTYPE=[A-Z]+)')
    dataSV.loc[:,"SV_type"]=dataSV["SV_type"].str.replace("SVTYPE=","")
    # dataSV.loc[:,"SV_type"]=dataSV["ID"].str.split(".",expand=True)[1]
    dataSV.loc[:,"GT"]=dataSV["SAMPLE"].str.split(":",expand=True)[0]
    def get_GT(df):
        if df["GT"] == "1/1":
            return "ALT"
        elif df["GT"] == "0/0":
            return "REF"
        elif df["GT"] == "0/1":
            return "Hybrid"
        else:
            return "Other_GT"
    dataSV.loc[:,"GT_type"]=dataSV.apply(get_GT,axis=1)
    x=dataSV["GT_type"].value_counts()
    def get_LEN(df):
        if df["LEN"] >= 50:
            return "more50bp"
        else:
            return "less50bp"
    dataSV.loc[:,"LEN_type"]=dataSV.apply(get_LEN,axis=1)
    y=dataSV["LEN_type"].value_counts()
    def get_SVtype(df):
        if df["SV_type"] == "DUP":
            return "DUP"
        elif df["SV_type"] == "INS":
            return "INS"
        elif df["SV_type"] == "DEL":
            return "DEL"
        elif df["SV_type"] == "INV":
            return "INV"
        else:
            return df["SV_type"]
    dataSV.loc[:,"SV"]=dataSV.apply(get_SVtype,axis=1)
    z=dataSV["SV"].value_counts()
    SVqcdict=pd.DataFrame(x).T.to_dict("records")
    if "ALT" in SVqcdict[0].keys():
        pass
    else:
        SVqcdict[0]["ALT"]=0
    if "REF" in SVqcdict[0].keys():
        pass
    else:
        SVqcdict[0]["REF"]=0
    if "Hybrid" in SVqcdict[0].keys():
        pass
    else:
        SVqcdict[0]["Hybrid"]=0             
    y=pd.DataFrame(y).T.to_dict("records")
    if 'more50bp' in y[0].keys():
        pass
    else:
        y[0]['more50bp']=0
    if 'less50bp' in y[0].keys():
        pass
    else:
        y[0]['less50bp']=0
    z=pd.DataFrame(z).T.to_dict("records")
    if 'DEL' in z[0].keys():
        pass
    else:
        z[0]['DEL']=0
    if 'INS' in z[0].keys():
        pass
    else:
        z[0]['INS']=0
    if 'INV' in z[0].keys():
        pass
    else:
        z[0]['INV']=0
    if 'DUP' in z[0].keys():
        pass
    else:
        z[0]['DUP']=0
    SVqcdict[0].update(y[0])
    SVqcdict[0].update(z[0])

    return data1,mappingqcdict,SVqcdict

def get_SVcheck(iputfile):
    with open(iputfile)as SVck:
        new_dict={}
        json_data=json.load(SVck)
        for key,value in json_data.items():
            i=re.sub("-","_",key)
            new_dict[i]=value
    return new_dict