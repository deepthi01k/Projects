#!usr/bin/env python
import os
import re
import pandas as pd
reference="/home/alphabox0006/Deepthi/Feb-Apr/TB-pipeline/Pipeline/MTB_H37Rv.fasta"
input_fasta="/home/alphabox0006/Deepthi/Feb-Apr/TB-pipeline/Samples" 
out_path="/home/alphabox0006/Deepthi/Feb-Apr/TB-pipeline/Pipeline/output2"
sample_name="ERR4796257_VAL_HSOT"
ref=pd.read_table("/home/alphabox0006/Deepthi/Deepthi_Final/Catalogue/Output/Catalog_version4.tsv")

def indexing(reference):
    print("Indexing of ref fasta")
    os.system(f"bwa index {reference} ")
    print("Indexing done\n")
indexing(reference)

def mapping (reference):
    print("map using bwa mem")
    os.system(f"bwa mem -t 12 {reference} {input_fasta}/{sample_name}_R1.fastq.gz {input_fasta}/{sample_name}_R2.fastq.gz >>{out_path}/{sample_name}.sam")
    os.system(f"samtools sort  -@ 12  {out_path}/{sample_name}.sam -o {out_path}/{sample_name}_sorted.bam")
    os.system(f"rm {out_path}/*.sam")
    os.system(f"samtools flagstat -@ 12 {out_path}/{sample_name}_sorted.bam > {out_path}/{sample_name}_flagstat.txt ")
    os.system(f"cat {out_path}/{sample_name}_flagstat.txt ") 
    print("mapping step done!") 
mapping(reference)

def variant_calling():
   print("Starting mpileup")
   os.system(f"bcftools mpileup  -O b --threads 12 -f {reference} {out_path}/{sample_name}_sorted.bam   -o {out_path}/{sample_name}_raw.bcf ")
   print("Mpileup finished")
   os.system(f"bcftools call  --threads 12 --ploidy 1 -m -v -o {out_path}/{sample_name}_raw.vcf {out_path}/{sample_name}_raw.bcf")
   print("variant calling started")
   os.system(f"rm {out_path}/*.bcf")
   print("variant calling done")
variant_calling()

def annotation():
    print("Start annotation")
    os.system(f" bcftools annotate -x FORMAT/GT,FORMAT/PL,^INFO/DP,INFO/INDEL {out_path}/{sample_name}_raw.vcf >{out_path}/{sample_name}_raw_edited.vcf")
    os.system(f"bcftools annotate -a MTB_annotation_1.tab.gz -h header.hdr -c CHROM,FROM,TO,-,INFO/Locus_tag,-,-,INFO/Gene,INFO/Protein_Name {out_path}/{sample_name}_raw_edited.vcf>  {out_path}/{sample_name}_annotated.vcf")
    os.system(f" bcftools filter -e  '%QUAL<20 | DP<10' {out_path}/{sample_name}_annotated.vcf>{out_path}/{sample_name}_edit.vcf")
    os.system(f"rm {out_path}/*_edited.vcf {out_path}/*_annotated.vcf " )
    print("Annotation done")
annotation()

print("Beginning mutation calling")

# Reading the vcf file without comments “#”
def read_vcf():
    df1 =pd.read_csv(out_path+"/"+sample_name+"_edit.vcf",sep="\t",header=None,comment="#",names=["CHR","POS","ID","WT base","Var. base","QUAL","FILTER","INFO","FORMAT","FORMAT2"]).drop(columns=["FILTER","FORMAT","FORMAT2"])
    df1["WT base"]=df1["WT base"].str.lower()
    df1["Var. base"]=df1["Var. base"].str.lower()
    return df1
df1=read_vcf()

print("starting with variant analysis step")
def cleanup_vcf(df1):
    Type=[]
    Locus_tag=[]
    Gene=[]
    Protein_Name=[]
    DP=[]
    for i in range(len(df1)):
   
    # IF INDEL and has locus_type then split the string 5 times
        if bool(re.findall(r"INDEL",df1["INFO"][i]))==True:
            if bool(re.findall(r"Locus_Type",df1["INFO"][i]))==True:	 
                Type.append(df1["INFO"][i].split(";")[0])
                DP.append(df1["INFO"][i].split(";")[1].replace("DP=",""))
                Locus_tag.append(df1["INFO"][i].split(";")[2].replace("Locus_tag=",""))
                Gene.append(df1["INFO"][i].split(";")[3].replace("Gene=",""))
                Protein_Name.append(df1["INFO"][i].split(";")[4].replace("Protein_Name=",""))
                
            # ELSE  split 2 times and fill other 3 as empty
 
            else:
                Type.append(df1["INFO"][i].split(";")[0])
                DP.append(df1["INFO"][i].split(";")[1].replace("DP=",""))
                Locus_tag.append(" ")
                Gene.append(" ")
                Protein_Name.append(" ")
        elif bool(re.findall(r"Locus_tag",df1["INFO"][i]))==True:
            Type.append("SNP")
            DP.append(df1["INFO"][i].split(";")[0].replace("DP=",""))
            Locus_tag.append(df1["INFO"][i].split(";")[1].replace("Locus_tag=",""))
            Gene.append(df1["INFO"][i].split(";")[2].replace("Gene=",""))
            Protein_Name.append(df1["INFO"][i].split(";")[3].replace("Protein_Name=",""))
            
        else :
            Type.append("SNP")
            DP.append(df1["INFO"][i].split(";")[0].replace("DP=",""))
            Locus_tag.append(" ")
            Gene.append(" ")
            Protein_Name.append(" ")

    # append into df

    df1.insert(6," Var. type",Type)
    df1.insert(8,"Locus_Tag",Locus_tag)
    df1.insert(9,"Gene Name",Gene)
    df1.insert(10,"Product",Protein_Name)
    df1.insert(7,"DP",DP)
    df1.drop(columns=["INFO","CHR","ID"],inplace=True)
    print("Variant analysis done")
    return df1
vcf =cleanup_vcf(df1)

def called(vcf,ref):
    print("Start with mutation calling step !")
    reference=ref[['Variant position genome start', 'Variant position genome stop','Gene Name','Var. type','WT base', 'Var. base','AA change','Antibiotic']]
    
    # merge with catalog based on pos , ref and alt allele and gene name
    
    result=pd.merge(vcf,reference,how="left",left_on=["POS","WT base","Var. base","Gene Name"],right_on=["Variant position genome start","WT base","Var. base","Gene Name"])
    
    #cleanup by dropping columns not required

    result.drop(columns=["Variant position genome start","Variant position genome stop","Var. type"],inplace=True)
    result.to_csv(out_path+"/"+sample_name+"_dst.tab",sep ="\t",index=False)
    print("Mutation calling done and output tab file generated !")
    return result
called(vcf,ref)
