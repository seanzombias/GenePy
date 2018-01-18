#!/bin/bash



######## WELCOME TO GenePy, Gene Patogenicity model ########
#	Provide a gene name as STDIN e.g. $> sh GenePy_1.0.sh NOD2
# In the same folder of GenePy the following files must be present:
#   @ ALL_GENES_FINAL.bed : an annotated .bed file containing all the target regions to be analyzed
#   @ FINAL_GTYPED_BIALLELIC.vcf.gz : a .vcf.gz (and its index) file with all the samples to be analyzed
#############################################################

GENE=$1

echo "... EXTRACTING GENE VARIANTS ..."
# generate scripts for extracting specific gene variants from all samples
python variant_extractor_3.0.py $GENE

# create the interval file
awk -v pat=^"$GENE"$ '$4 ~ pat{print $1":"$2"-"$3}' ALL_GENES_FINAL.bed > "$GENE".intervals

#extract variants using GATK from FINAL_GTYPED_BIALLELIC.vcf.gz
sh "$GENE"_IBD.sh
sh "$GENE"_Lev.sh
sh "$GENE"_FULL.sh

## Filter vcf outputs
module load bcftools/1.2.1
module load samtools/1.2.1
vcfutils.pl varFilter -d 4 -1 0 -2 0 -3 0 -4 0 "$GENE"_IBD.vcf > "$GENE"_IBD_filt.vcf
vcfutils.pl varFilter -d 4 -1 0 -2 0 -3 0 -4 0 "$GENE"_Lev.vcf > "$GENE"_Lev_filt.vcf
vcfutils.pl varFilter -d 4 -1 0 -2 0 -3 0 -4 0 "$GENE"_FULL.vcf > "$GENE"_FULL_filt.vcf


echo "CONVERTING AND ANNOTATING WITH ANNOVAR..."
## Convert full filtered vcf to annovar for annotation
/home/em1c14/annovar/convert2annovar.pl \
  -format vcf4 "$GENE"_FULL_filt.vcf \
  -outfile "$GENE"_FULL_annovar.input \
  -allsample \
  -withfreq \
  -include 2>"$GENE"_annovar.log  

/home/em1c14/annovar/convert2annovar.pl \
  -format vcf4 "$GENE"_IBD_filt.vcf \
  -outfile "$GENE"_IBD_annovar.input \
  -allsample \
  -withfreq \
  -include 2>"$GENE"_annovar.log

/home/em1c14/annovar/convert2annovar.pl \
  -format vcf4 "$GENE"_Lev_filt.vcf \
  -outfile "$GENE"_Lev_annovar.input \
  -allsample \
  -withfreq \
  -include 2>"$GENE"_annovar.log

cut -f 18- "$GENE"_IBD_annovar.input > "$GENE"_tmp
grep '^#CHR' "$GENE"_IBD_filt.vcf | cut -f 10- > "$GENE"_tmp2
cat "$GENE"_tmp2 "$GENE"_tmp > "$GENE"_IBD_pseudo.vcf

cut -f 18- "$GENE"_Lev_annovar.input > "$GENE"_tmp
grep '^#CHR' "$GENE"_Lev_filt.vcf | cut -f 10- > "$GENE"_tmp2
cat "$GENE"_tmp2 "$GENE"_tmp > "$GENE"_Lev_pseudo.vcf



## Annotate the full_selected vcf 
/home/em1c14/annovar/table_annovar.pl \
  "$GENE"_FULL_annovar.input \
  /home/em1c14/annovar/humandb/ \
  -buildver hg38 \
  -out "$GENE" \
  -remove \
  -protocol refGene,dbnsfp33a,exac03,1000g2015aug_eur,avsnp147 \
  -operation g,f,f,f,f \
  -nastring . >>"$GENE"_annovar.log

echo "MAKING CUSTOM VCFs..."
## Merge VCFs and Annovar
python merge_vcf_pseudo_3.0.py "$GENE"_IBD_pseudo.vcf "$GENE".hg38_multianno.txt "$GENE"_IBD.meta
python merge_vcf_pseudo_3.0.py "$GENE"_Lev_pseudo.vcf "$GENE".hg38_multianno.txt "$GENE"_Lev.meta


## Run model and test
#	script takes: IBD, Levin, GENE inputs
echo "FINAL STEP, MAKE MATRICES..."
python make_scores_mat_3.0.py "$GENE"_IBD.meta "$GENE"_Lev.meta "$GENE" > "$GENE".pval

mv "$GENE".pval PVALUES/.

rm "$GENE"_* "$GENE".hg38_multianno.txt "$GENE".intervals
