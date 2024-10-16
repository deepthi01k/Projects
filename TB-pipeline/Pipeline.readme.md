TB pipeline
Input :directory containing the sample files
header.hd and MTB_annotation.tab.gz should be in same directory as pipeline
Output: in 2 folders -Analysis containing the bam, vcf file,fastwc report
Result folder having all the csv files , reports etc
Requirement : 1. SNPIT tool installed , bcftools ver>1.16 ,logtime python
command :
Activate snpit if inside venv 
 conda activate variant_clone
#python Py_pipeline.py -f1 /home/alphabox0006/Deepthi/Feb-Apr/TB-pipeline/Pipeline/Sample/sam_1 -m /home/alphabox0006/Deepthi/Feb-Apr/TB-pipeline/Pipeline/Script/batch_001_metadata.tsv  -C /home/alphabox0006/Deepthi/Feb-Apr/TB-pipeline/Pipeline/Script/Catalog_version4.tsv
