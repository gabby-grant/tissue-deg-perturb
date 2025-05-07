# tissue-deg-perturb
**Download Data:**
Download uterine tissue data prepared by [Wang et al](https://figshare.com/articles/dataset/Data_record_1/5330539)

### DESeq2
For preprocessing data and DESeq2 workflow follow https://github.com/gabby-grant/deseq-workflow

# GEMDiff 
### processing gene matrix 
```bash
git clone https://github.com/feltus/gembuild
# merge matrices
python gembuild/gem_merge.py uterus-rsem-fpkm-gtex.txt ucec-rsem-fpkm-tcga-t.txt UTERN_UCECT.txt
# transpose
bash gembuild/transpose_gem.sh UTERN_UCECT.txt > UTERN_UCECT_transpose.txt
# split into test and training
bash gembuild/split_gem.sh UTERN_UCECT_transpose.txt UTERN_UCECT_train.txt UTERN_UCECT_test.txt
# convert to tab delimited
cat UTERN_UCECT_train.txt | awk '{$1=$1}1' OFS='\t'  > UTERN_UCECT.train; cat UTERN_UCECT_test.txt | awk '{$1=$1}1' OFS='\t'  > UTERN_UCECT.test
# make labels
bash gembuild/make_labels.sh UTERN_UCECT.train; bash gembuild/make_labels.sh UTERN_UCECT.test
```
check distribuition of data by running a histogram on the datasets
```bash
head -n 20 UTERN_UCECT.txt | awk '{print $1}' > gene_names_UTERN.txt
python gembuild/gem_histogram.py -e UTERN_UCECT.test -g gene_names_UTERN.txt -o test-output-histogram.png
python gembuild/gem_histogram.py -e UTERN_UCECT.train -g gene_names_UTERN.txt -o train-output-histogram.png
```
## git gemdiff
### training script
### perturb script 
### output 
## gemdif with gene set
gene set from selected studies 
1. train model on new gene set
2. perturb model on new gene set
