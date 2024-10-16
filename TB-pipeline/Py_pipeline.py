#!usr/bin/env python
#source /home/alphabox0006/Deepthi/Work_Done_Nov-Feb/snpit/snpit/.venv/bin/activate
# conda activate variant_clone
#python Py_pipeline.py -f1 /home/alphabox0006/Deepthi/Feb-Apr/TB-pipeline/Pipeline/Sample/sam_1 -m /home/alphabox0006/Deepthi/Feb-Apr/TB-pipeline/Pipeline/Script/batch_001_metadata.tsv  -C /home/alphabox0006/Deepthi/Feb-Apr/TB-pipeline/Pipeline/Script/Catalog_version4.tsv
import os
import numpy as np
import glob
import re
import argparse
import pandas as pd
from Bio.Seq import Seq
import logging
import datetime
import gzip
parser = argparse.ArgumentParser(description='TB_pipeline')
parser.add_argument('-f1', '--fastq',action ='store',dest='fastq1', help='Specify file path for fastq',default = ".")
parser.add_argument('-m', '--metadata',action ='store',dest='meta', help='Specify the reference fasta ',default ='.')
parser.add_argument('-C', '--catalog',action ='store',dest='catalog', help='Specify the name of catalog ',default = ".")
args = parser.parse_args()
#input files being used 

input_fastq = args.fastq1
catalog_path=args.catalog
metadata_table=args.meta
Rfasta = "M._tuberculosis_H37Rv_2015-11-1.fasta"
dir=os.path.basename(input_fastq).split("/")[0]
def find(name, path):
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)
        
reference=find(Rfasta,os.path.expanduser('~'))
#metadata_table=find(dir+"_metadata.tsv",os.path.expanduser('~'))

def logs():
    #function for generating log files
        logging.basicConfig(filename=dir+".log",format='%(asctime)s  %(message)s',datefmt='%m/%d/%Y %I:%M:%S %p',level=logging.DEBUG,filemode='a')
        logger = logging.getLogger()
        logger.setLevel(logging.DEBUG)
        return logger


def indexing(reference):
    logger.debug("Indexing of ref fasta")
    os.system(f"bwa index {reference} ")
    logger.debug("Indexing done\n")
    


def values(sample,dir):
    # getting sample name 
    sample_name=os.path.basename(sample).split("_R1")[0] 
    #setting output path
    parent_dir=os.path.abspath(os.path.join(os.getcwd(), os.pardir))
    out_path=os.path.join(os.path.join(parent_dir,'Staging','Analysis',dir,sample_name))
    result_out=os.path.join(os.path.join(parent_dir,'Staging','Result',dir,sample_name))

    #make output directory if does not exist
    os.makedirs(out_path,exist_ok=True)
    os.makedirs(result_out,exist_ok=True)

    fastq=sample.split("_R1")[0]
    return sample_name,out_path,result_out,fastq


        

def mapping (reference):  
        # mapping fastq with ref file            
        print("Start Bwa mem")
        logger.debug("map using bwa mem")
        os.system(f"bwa mem -t 12 {reference} {fastq}_R1.fastq.gz {fastq}_R2.fastq.gz  >{out_path}/{sample_name}.sam")
        os.system(f"samtools sort  -@ 12  {out_path}/{sample_name}.sam -o {out_path}/{sample_name}_sorted.bam")
        os.system(f"rm {out_path}/*.sam")
        os.system(f"samtools flagstat -@ 12 {out_path}/{sample_name}_sorted.bam > {out_path}/{sample_name}_flagstat.txt ")
        os.system(f"cat {out_path}/{sample_name}_flagstat.txt ") 
        logger.debug("mapping step done!") 
        print("Bwa mem done")


def variant_calling():
        #variant calling
        logger.debug("Starting mpileup")
        os.system(f"bcftools mpileup  -O b --threads 12 -f {reference} {out_path}/{sample_name}_sorted.bam   -o {out_path}/{sample_name}_raw.bcf ")
        logger.debug("Mpileup finished")
        os.system(f"bcftools call  --threads 12 --ploidy 1 -m -A -o {out_path}/{sample_name}_raw.vcf {out_path}/{sample_name}_raw.bcf")
        logger.debug("variant calling started")
        os.remove(f"{out_path}/{sample_name}_raw.bcf")
        os.system(f''' bcftools view {out_path}/{sample_name}_raw.vcf --exclude 'ALT="."'>{out_path}/{sample_name}.vcf''') # remove those entries where alt column is empty
        os.remove(f"{out_path}/{sample_name}_raw.vcf") #remove the big vcf file
        logger.debug("variant calling done")


def lineage_calling():
        # call lineage using snpit
        os.system(f"snpit --input {out_path}/{sample_name}.vcf --output {out_path}/{sample_name}_lineage.tsv ")
        classification=pd.read_csv(f"{out_path}/{sample_name}_lineage.tsv",sep="\t")
        classification.columns=["FullID","Species","Lineage_number",'Sublineage',"Lineage","Percentage"]
        # creating a new df with desired column names and filtering out if sample is MTB or not
        classification.loc[0,"FullID"]=os.path.basename(sample_name)
        if classification.loc[0,"Percentage"]<60:
            classification.loc[0,"Species"]=classification.loc[0,"Species"].replace("M. tuberculosis","Not detected")
            classification.to_csv(f"{result_out}/{sample_name}_classification.csv",index=False)
        else :
            classification.to_csv(f"{result_out}/{sample_name}_classification.csv",index=False)
        logger.debug(f"Lineage determined and file stored in {out_path}")

        os.remove(f"{out_path}/{sample_name}_lineage.tsv")
        return classification


def annotation():
        # annotating vcf file
        logger.debug("Start annotation")
        os.system(f" bcftools annotate -x FORMAT/GT,FORMAT/PL,^INFO/DP,INFO/INDEL,INFO/DP4 {out_path}/{sample_name}.vcf >{out_path}/{sample_name}_raw_edited.vcf")
        os.system(f"bcftools annotate -a MTB_annotation.tab.gz -h header.hdr -c CHROM,FROM,TO,-,INFO/Locus_tag,-,-,INFO/Gene,INFO/Protein_Name {out_path}/{sample_name}_raw_edited.vcf>  {out_path}/{sample_name}_annotate.vcf")
        os.system(f" bcftools filter -e  '%QUAL<20' {out_path}/{sample_name}_annotate.vcf>{out_path}/{sample_name}_annotated.vcf")
        os.system(f"rm {out_path}/*_edited.vcf " )
        os.remove(f"{out_path}/{sample_name}_annotate.vcf")
        logger.debug("Annotation done")

def read_files():
        logging.debug("Read files")
        # Reading the vcf file without comments “#”
        df1 =pd.read_csv(f"{out_path}/{sample_name}_annotated.vcf",sep="\t",header=None,comment="#",names=["CHR","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","FORMAT2"]).drop(columns=["FILTER","FORMAT","FORMAT2"])
        ref=pd.read_table(catalog_path)
        # lineage=pd.read_csv(f"{out_path}/{sample_name}_lineage.tsv",sep="\t")
        met=pd.read_csv(metadata_table,sep="\t")
        logging.debug("Read files done!")
        return df1,ref,met

def edit_catalog(ref):
        for i in range(len(ref)):
            ref.loc[[i],["WT base"]]=ref["WT base"][i].upper()
            ref.loc[[i],["Var. base"]]=ref["Var. base"][i].upper()
            
            #finding reverse complement of genes on negative strand
            if ref["Dir."][i]=="-":
                refb=ref["WT base"][i]
                altb=ref["Var. base"][i]
                ref.loc[[i],["WT base"]]=refb.replace(refb,str(Seq(refb).reverse_complement()))
                ref.loc[[i],["Var. base"]]=altb.replace(altb,str(Seq(altb).reverse_complement()))
            else:
                pass
        ref["AA_codon"]=ref["AA change"]+" "+ "("+ref["Codon change"]+")"
        catalog=ref[~ref["Antibiotic"].str.contains("phylo",regex=True)]
        return catalog

def cleanup_vcf(split):
        Type=[]
        Locus_tag=[]
        Gene=[]
        Protein_Name=[]
        DP=[]
        DP4=[]

        #Splitting ALT columns with 2 alt alleles by ","  and explode df to get individual rows for each allele
        df1["ALT"]=df1["ALT"].str.split(",")
        split=split.explode("ALT")
        # spliting the info column on basis of presence of INDEL and Locus_tag
        for i in split["INFO"]:
            if (i.find("INDEL") !=-1) & (i.find("Locus_tag")!=-1):
                Type.append(i.split(";")[0])
                DP.append(i.split(";")[1].replace("DP=",""))
                DP4.append(i.split(";")[2])
                Locus_tag.append(i.split(";")[3].replace("Locus_tag=",""))
                Gene.append(i.split(";")[4].replace("Gene=",""))
                Protein_Name.append(i.split(";")[5].replace("Protein_Name=",""))
            elif (i.find("INDEL") !=-1) & (i.find("Locus_tag")==-1):
                Type.append(i.split(";")[0])
                DP.append(i.split(";")[1].replace("DP=",""))
                DP4.append(i.split(";")[2])
                Locus_tag.append("")
                Gene.append("")
                Protein_Name.append(" ")
            elif (i.find("INDEL") ==-1) & (i.find("Locus_tag")!=-1):
                Type.append("SNP")
                DP.append(i.split(";")[0].replace("DP=",""))
                DP4.append(i.split(";")[1])
                Locus_tag.append(i.split(";")[2].replace("Locus_tag=",""))
                Gene.append(i.split(";")[3].replace("Gene=",""))
                Protein_Name.append(i.split(";")[4].replace("Protein_Name=",""))
            else :
                Type.append("SNP")
                DP.append(i.split(";")[0].replace("DP=",""))
                DP4.append(i.split(";")[1])
                Locus_tag.append("")
                Gene.append("")
                Protein_Name.append(" ")

        split.drop(columns=["INFO","CHR","ID"],inplace=True)
        #assigning values todf columns
        split["Var. type"]=Type
        split["DP"]=DP
        split["DP4"]=DP4
        split["Gene_ID"]=Locus_tag
        split["Gene_Name"]=Gene
        split["Product"]=Protein_Name

        # function to calculate frequency of each allele
        def calculate_freq(vcf):            
            vcf["DP4"] = vcf['DP4'].str.extract("(\d+,\d+,\d+,\d+)")
            vcf[["a","b","c","d"]] = vcf["DP4"].str.split(",", expand=True)
            vcf["Freq"]=(((vcf["c"].astype(int) + vcf["d"].astype(int))/(vcf["a"].astype(int) + vcf["b"].astype(int) + vcf["c"].astype(int) + vcf["d"].astype(int)))*100).round(decimals =3)
            vcf["95_confidence"]=vcf["Freq"].astype(int)>=95
            vcf.drop(columns=["DP4","a","b","c","d"], inplace=True)
            return vcf
        calculate_freq(split)  
        split_vcf=split[split["Freq"]>20].reset_index(drop=True)
        # os.remove(out_path+sample_name+"_edit.vcf") 
        print("editing of vcf complete!")
        logger.debug("Cleaning up vcf complete!!")
        return split_vcf

def merge_vcf_catalog():
        logger.debug("Start with mutation calling step !")
        
        catalog=ref.loc[:,['Variant position genome start','Var. type', 'Number', 'WT base', 'Var. base', 'Gene ID','Gene Name','AA_codon','Antibiotic','Reference PMID','High Confidence SNP']]
        catalog.rename({'Variant position genome start':'POS', 'WT base':'REF','Var. base':'ALT'},axis=1,inplace=True)
        
        #merge catalog with vcf based on 4 columns
        result=vcf.merge(catalog,how="left",on=['POS','ALT','REF','Var. type'])
        result.to_csv(result_out+"/"+sample_name+"_dst.tab",sep ="\t",index=False)
        
        logger.debug("Mutation calling done and output tab file generated !")
        return result

def Resistance_check(a, x):
        if any(re.search(a, x) for x in x):
            b = "Resistant"
        else:
            b = "No Resistance detected"
        return b

def dst_table(result):
        if (classi_file["Species"].str.contains("M. tuberculosis")).bool()==True:
            print("If executed")  
        #creating dst table by cheking which all mutations are present
            Antibiotic = result[result['Antibiotic'].notnull()& (result['Antibiotic'] != ' ')]['Antibiotic'].unique()
            drugs_s = {}
            drugs_s['Isoniazid (INH)'] = Resistance_check("isoniazid", Antibiotic)
            drugs_s[ 'Rifampicin (RMP)'] = Resistance_check("rifampicin", Antibiotic)
            drugs_s['Streptomycin (SM)'] = Resistance_check("streptomycin", Antibiotic)
            drugs_s['Ethambutol (EMB)'] = Resistance_check("ethambutol", Antibiotic)
            drugs_s['Pyrazinamide (PZA)'] = Resistance_check("pyrazinamide", Antibiotic)
            drugs_s['Ofloxacin (OFX)'] = Resistance_check("fluoroquinolones", Antibiotic)
            drugs_s['Moxifloxacin (MOX)'] = Resistance_check("fluoroquinolones", Antibiotic)
            drugs_s['Gatifloxacin (GAT)'] = Resistance_check("fluoroquinolones", Antibiotic)
            drugs_s['Amikacin (AMK)'] = Resistance_check("amikacin", Antibiotic)
            drugs_s['Capreomycin (CAP)'] = Resistance_check("capreomycin", Antibiotic)
            drugs_s['Ethionamide (ETO)'] = Resistance_check("ethionamide", Antibiotic)
            drugs_s['Kanamycin (KAN)'] = Resistance_check("kanamycin", Antibiotic)
            drugs_s['Linezolid (LZD)'] = Resistance_check("linezolid", Antibiotic)
            drugs_s['Para-aminosalicylic acid (PAS)'] = Resistance_check("para-aminosalicylic acid", Antibiotic)
            drugs_s['Bedaquiline (BDQ)'] = Resistance_check("bedaquiline", Antibiotic)
            drugs_s['Clofazimine (CLO)'] = Resistance_check("clofazimine", Antibiotic)
            drugs_s['Delaminid (DMD)'] = Resistance_check("delamanid", Antibiotic)
            drugs_s['Pretomanid (PTM)'] = Resistance_check("pretomanid", Antibiotic)
        else :
            print("Else executed")
            drugs_s = {}
            drugs_s['Isoniazid (INH)'] = "No Resistance detected"
            drugs_s[ 'Rifampicin (RMP)'] =  "No Resistance detected"
            drugs_s['Streptomycin (SM)'] =  "No Resistance detected"
            drugs_s['Ethambutol (EMB)'] =  "No Resistance detected"
            drugs_s['Pyrazinamide (PZA)'] = "No Resistance detected"
            drugs_s['Ofloxacin (OFX)'] =  "No Resistance detected"
            drugs_s['Moxifloxacin (MOX)'] =  "No Resistance detected"
            drugs_s['Gatifloxacin (GAT)'] =  "No Resistance detected"
            drugs_s['Amikacin (AMK)'] =  "No Resistance detected"
            drugs_s['Capreomycin (CAP)'] =  "No Resistance detected"
            drugs_s['Ethionamide (ETO)'] =  "No Resistance detected"
            drugs_s['Kanamycin (KAN)'] =  "No Resistance detected"
            drugs_s['Linezolid (LZD)'] =  "No Resistance detected"
            drugs_s['Para-aminosalicylic acid (PAS)'] =  "No Resistance detected"
            drugs_s['Bedaquiline (BDQ)'] =  "No Resistance detected"
            drugs_s['Clofazimine (CLO)'] =  "No Resistance detected"
            drugs_s['Delaminid (DMD)'] = "No Resistance detected"
            drugs_s['Pretomanid (PTM)'] =  "No Resistance detected"
              

        dst=pd.DataFrame.from_dict(data=[drugs_s],orient='columns',dtype='str')
        dst.to_csv(result_out+"/"+sample_name+"_dst_table.csv",index=False)
        logger.debug("Mutations detected !")
        return dst




def Clinical_report(d):
        # sample_name=os.path.basename(sample).split("_R1")[0]
# if classification file has MTB then get the dst table of XDR ,Pre-XDR ,MDR type
    if (classi_file["Species"].str.contains("M. tuberculosis")).bool()==True:
                print("Detected")
                tb_det=["Detected"]
                lineage=classi_file["Lineage"]
        
                

                if (d.iloc[0,0] == "No Resistance detected" and d.iloc[0,1] == "No Resistance detected"\
                    and d.iloc[0,2] == "No Resistance detected" and d.iloc[0,3] == "No Resistance detected"\
                    and d.iloc[0,4] == "No Resistance detected" and d.iloc[0,5] == "No Resistance detected"\
                    and d.iloc[0,6] == "No Resistance detected" and d.iloc[0,7] == "No Resistance detected"\
                    and d.iloc[0,8] == "No Resistance detected" and d.iloc[0,9] == "No Resistance detected"\
                    and d.iloc[0,10] == "No Resistance detected" and d.iloc[0,11] == "No Resistance detected"\
                    and d.iloc[0,12] == "No Resistance detected" and d.iloc[0,13] == "No Resistance detected"\
                    and d.iloc[0,14] == "No Resistance detected" and d.iloc[0,15] == "No Resistance detected"\
                    and d.iloc[0,16] == "No Resistance detected"and d.iloc[0,17] == "No Resistance detected"):
                        dst = "Not Drug Resistance"
                elif(d.iloc[0,1] == "Resistant"):

                        if(d.iloc[0,5] == "No Resistance detected" and  d.iloc[0,6] == "No Resistance detected" 
                    and d.iloc[0,7] == "No Resistance detected" and d.iloc[0,0] == "Resistant"):
                            dst = "Multi Drug-Resistant tuberculosis (MDR-TB)"
                        elif(d.iloc[0,5] == "Resistant" or  d.iloc[0,6] == "Resistant" or d.iloc[0,7] == "Resistant"
                    and d.iloc[0,12] == "No Resistance detected" and d.iloc[0,14] == "No Resistance detected"):
                            dst = "Pre-extensively Drug-Resistant tuberculosis (Pre-XDR-TB)"
                        elif(d.iloc[0,5] == "Resistant" or  d.iloc[0,6] == "Resistant" or d.iloc[0,8] == "Resistant"
                    and d.iloc[0,12] == "Resistant" or d.iloc[0,14] == "Resistant"):
                            dst = "Extensively Drug-Resistant tuberculosis (XDR-TB)"

                        else:
                            dst = "Drug-Resistant tuberculosis (DR-TB)"
                else:
                        dst = "Drug-Resistant tuberculosis (DR-TB)"
                
    else: #if MTB not present then following details
            tb_det=["Not detected"]
            lineage="NA"
            dst="NA"
            print("Not Detected")


    Data=pd.DataFrame()
        # Data["Sample"]=[sample_name]
    Data["Mycobacterium Tuberculosis detected"]=tb_det
    Data["Coverage"] =""
    Data["Lineage"]=lineage
    Data["Proportion of genome covered"]=" "
    Data["Genomic dst profile"]=dst
    print("Clinical Report Done")
    Data.to_csv(f"{result_out}/{sample_name}_Clinical_report.csv",index=False)
    return Data

def mutation_table():
        # get all mutations from mutation  table
        flag=not np.any(result["Antibiotic"])
        if not flag :
            mutation_table=result[['POS', 'REF', 'ALT','Var. type','Freq','AA_codon','Gene ID', 'Gene Name','Antibiotic','Reference PMID','High Confidence SNP']][result["Antibiotic"].notnull()].reset_index(drop=True)
            mutation_table.columns=['Position','Ref Allele','Alt Allele','Type','Frequency','Substitution','Gene Symbol','Gene Name','Drug','PMID','High confidence SNP']
            # mutation_table['Antibiotic']= mutation_table['Antibiotic'].str.findall(pat=r"\s\(([A-Z]*)\)").str.join(sep="|")
            mutation_table.to_csv(f"{result_out}/{sample_name}_mutation_table.csv",sep ="\t",index=False)
        else :
            Result=["No mutations detected"]
            mutation_table=pd.DataFrame(Result,columns=["Result"]).to_csv(result_out+"/"+sample_name+"_mutation_table.csv",sep ="\t",index=False)
        logger.debug("mutation_table done!")

def genome_analysis():# genome analysis file 
        df_gaas = {'Application':["OmegaTB New_v1"],
            'Run_Date':[datetime.date.today()],
            'QC_check': [""], # add this info
            'Application_run_by':["HAPL"]}
        gen_analysis=pd.DataFrame(df_gaas)
        gen_analysis.to_csv(f"{result_out}/{sample_name}_Genome_analysis_summary.csv",index=False)
        return gen_analysis

def Genome_sequence_summary(i): # genome summary file 
                logger.debug("Starting genome sequence summary")

                # getting flowcell Id from fastq file itself
                with gzip.open(f"{fastq}_R1.fastq.gz",'rt') as fh:
                    k=fh.readline()
                    flowcell=k.split(":")[0]+"-"+k.split(":")[2]

                #getting other data from metadata file
                seq_summary={"Library type":[met.loc[i,"Library_Type"]]    ,
                "Library date": [met.loc[i,"Library_Date"]],
                "Library qc": [met.loc[i,"Library_QC"]]  ,
                "Genome Sequencer": [met.loc[i,"Genome_Sequencer"]],
                "Run date": [met.loc[i,"Run_Date"]],
                "MachineID-FlowcellID":[flowcell]} 
                Genome_seq=pd.DataFrame(seq_summary)
                Genome_seq.to_csv(f"{result_out}/{sample_name}_Genome_sequence_summary.csv",index=False)
                logger.debug("Genome sequence summary done")
                return Genome_seq

# getting all data from metadata file
def sample_summary(i):
    df_samsum = {
        "Sample_ID": [met.loc[i,"Sample_ID"]],
        "Client_ID": [met.loc[i,"Client_ID"]],
        "Source": "OmegaTB",
        "Sample_Type": [met.loc[i,"Sample_Type"]],
        "Sample_Receipt_Date": [met.loc[i,"Sample_Receipt_Date"]],
        "Sample_Temperature_on_Reciept_celsius": [met.loc[i,"Sample_Temperature_on_Reciept_celsius"]],
        "Sample_QC_Date": [met.loc[i,"Sample_QC_Date"]],
        "Sample_QC": [met.loc[i,"Sample_QC"]],
        "Shipping_Date": [met.loc[i,"Shipping_Date"]],
        "Lab_Name": [met.loc[i,"Lab_Name"]],
        "Contact": [met.loc[i,"Contact"]],
        "Registered": [met.loc[i,"Registered"]],
        "Collection_date": [met.loc[i,"Collection_date"]],
        "Conc": [met.loc[i,"Conc"]],
        "QC_260": [met.loc[i,"QC_260"]],
        "Sequenced_From": [met.loc[i,"Sequenced_From"]],
    }
    sam_sum=pd.DataFrame(df_samsum)
    sam_sum.to_csv(f"{result_out}/{sample_name}_sample_summary.csv",index=False)
    return sam_sum

def run_bedtools(sample_name): # running bedtools for getting depth file
    logger.debug("Starting bedtools")
    depth = os.system(f"bedtools genomecov -ibam {out_path}/{sample_name}_sorted.bam -d > {result_out}/{sample_name}_depth.txt")
    logger.debug("Depth files created !")
    print("Depth files created !")
    return depth

files=glob.glob(os.path.join(input_fastq,"**/*_R1.fastq.gz"),recursive=True)
for sample in files :
    
    # Derive sample names from input_path
    sample_name,out_path,result_out,fastq=values(sample,dir)
    logger=logs()
    logger.debug(f'{sample_name}')
    logger.debug('The pipeline starts here')

    # perform fastqc per sample
    def fastqc ():
        os.system(f"fastqc -t 12 {fastq}_R1.fastq.gz {fastq}_R2.fastq.gz -o {out_path}")
        os.system(f"unzip {out_path}/{sample_name}_R1_fastqc.zip -d {out_path}")
        logger.debug("Fastqc completed")
    fastqc()
    print("Fastqc done!!")

    print("Indexing ...") # indexing
    indexing(reference)
    print("Indexing done !")

    print("mapping")
    mapping(reference) #mapping
    print("Mapping done !")

    print("Variant calling...")
    variant_calling() # variant calling
    print("Variant calling done !")

    print("lineage calling ...")
    classi_file=lineage_calling()
    print("Lineage calling done !")

    print("Annotation ....")
    annotation()
    print("Annotation Done !")

    logger.debug("Beginning mutation calling")

   
    #read required files
    df1,ref,met=read_files()

    logger.debug("Beginning editing catalog")
    print("Starting catalog manipulation... ")
    ref=edit_catalog(ref)
    print("Starting catalog manipulation... ")

    logger.debug("starting with variant analysis step")
    # cleaning up of vcf to get a cleaner dataframe
    vcf= cleanup_vcf(df1)
    logger.debug("Cleaning up vcf")

    #merge vcf with catalog to get dst results
    result=merge_vcf_catalog()
    print("Getting the dst results !!")

    #function to check dst profile
    logger.debug("Dst preparing")
    # using above function find dst 

    mut=dst_table(result)
    print("Preparing dst table")
    # make a table with all mutation information

    mutation_table()
    logger.debug("Final result compilation")
    print("Preparing mutation table")
    # function to generate sample summary file

    Clinical_report(mut)
    print("Preparing Clinical report")

    genome_analysis()
    print("Genome analysis file")

# finding out if our sample name is in the metadata and then only generating genome sequence summary and sample summary
    for i in range(len(met)):
         if met.loc[i,"File_name"]==sample_name:
            print(met.loc[i,"File_name"])
            gen_seq=Genome_sequence_summary(i)
            sample_sum=sample_summary(i)
    print("Genome Summary file done !")
    print("Sample Summary file done !")

    run_bedtools(sample_name)
    print("Creating depth files using bedtools")

    logger.debug("Analysis done !")
    logger.debug("\n")
print("All analysis complete")
logger.debug("The full run is over")
