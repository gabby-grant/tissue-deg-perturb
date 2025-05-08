# tissue-deg-perturb
### Download Data
Download uterine tissue data prepared by [Wang et al](https://figshare.com/articles/dataset/Data_record_1/5330539)

### DESeq2
For preprocessing data and DESeq2 workflow follow https://github.com/gabby-grant/deseq-workflow

# GEMDiff 
### Preprocessing Data with `gembuild` 
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
### check distribuition of data by running and analyzing distributions on a histogram
```bash
head -n 1 UTERN_UCECT.test > test-gene-names.txt
head -n 1 UTERN_UCECT.train > train-gene-names.txt
python gembuild/gem_histogram.py -e UTERN_UCECT.test -g test-gene-names.txt -o test-output-histogram.png
python gembuild/gem_histogram.py -e UTERN_UCECT.train -g train-gene-names.txt -o train-output-histogram.png
```
## git gemdiff
 
### gene set
gene set from selected studies 
1. train model on new gene set
2. perturb model on new gene set
### training script
### perturb script 
### output
