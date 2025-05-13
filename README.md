# tissue-deg-perturb
### Download Data
Download uterine tissue data prepared by [Wang et al](https://figshare.com/articles/dataset/Data_record_1/5330539)

### DESeq2
For preprocessing data and DESeq2 workflow follow https://github.com/gabby-grant/deseq-workflow

# GEMDiff 
## Preprocessing Data with `gembuild` 
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
Check distribuition of data by running and analyzing distributions on a histogram.
```bash
# transfor data to log2 FC distribuition
python gembuild/log2_transform_gem.py -i UTERN_UCECT.test -o UTERN_UCECT.test.log2
python gembuild/log2_transform_gem.py -i UTERN_UCECT.train -o UTERN_UCECT.train.log2
head -n 1 UTERN_UCECT.test > test-gene-names.txt
head -n 1 UTERN_UCECT.train > train-gene-names.txt
python gembuild/gem_histogram.py -e UTERN_UCECT.test.log2 -g test-gene-names.txt -o test-output-histogram.png -l
python gembuild/gem_histogram.py -e UTERN_UCECT.train.log2 -g train-gene-names.txt -o train-output-histogram.png -l
```
## GEMdiff Setup

 activate the created conda environment create before running script
 
### Gene Set of Interest
This gene set is selected from: 
Bianco, B., Barbosa, C. P., Trevisan, C. M., Laganà, A. S., & Montagna, E. (2020). Endometrial cancer: A genetic point of view. Translational Cancer Research, 9(12), 7706–7715. https://doi.org/10.21037/tcr-20-2334

```
Endometrial-Cancer-Bianco	ARID1A	ELAPOR1	C1QTNF12	MTOR	GADD45A	GSTM1	LEPR	MACF1	MTHFR	NOTCH2	NRAS	ADIPOR1	AGT	AKT3	CAPN9	CD247	CFHR4	CRP	DNAH14	FASLG	IL24	MDM4	MUC1	PTGS2	RHEX	ZBTB7B	CYP1B1	DNMT3A	EIF2AK3	FABP1	MEIS1	MSH2	MSH6	BARD1	BCL2L11	CASP8	CXCR4	ERBB4	ERCC3	IL1A	IL1R2	IMP4	MCM6	UGT1A1	UGT1A8	CTDSPL	CTNNB1	OGG1	MLH1	NPRL2	GFBR2	TLR9	VGLL4	XPC	ADIPOQ	AGTR1	FOXL2	IGSF10	PIK3CA	PIK3CB	TNK2	FGFR3	JAKMIP1	MSX1	STX18	BMP3	KIT	CXCL3	FBXW7	GC	KIT	NFKB1	PDGFRA	SULT1E1	TLR2	UGT2B4	CLPTM1L	IL7R	PRLR	TERT	APC	CMYA5	CYFIP2	HSD17B4	MSH3	PCDHGB4	PDGFRB	PIK3R1	RAD50	TNFAIP8	CDKN1A	DDR1	DST	HFE	CDKN1A	VEGFA	YIPF3	ARID1B	ESR1	IGF2R	RNASET2	SYNE1	TBP	AHR	C7orf50	EGFR	IL6	NOD1	PMS2	TWIST1	BRAF	CAPZA2	DPP6	LEP	MET	POT1	SERPINE1	XRCC2	BMP1	MSR1	PPP2R2A	CCAT1	CNBD1	CPQ	LAPTM4B	MYC	NCOA2	PCMTD1	PRKDC	PSCA	PTK2	PTP4A3	SOX17	TERF1	CDKN2A	LINGO2	CD274	PDCD1LG2	TLN1	DNM1	NIPSNAP3B	TGFBR1	TLR4	TSC1	XPA	TRDMT1	NUDT5	ADAM12	C10orf88	CYP2E1	SLF2	FAZ	FGFR2	HHEX	MBL2	MGMT	PLAU	PRKG1	PTEN	RNLS	CXCL12	SIRT1	TNKS2	CYP2R1	HRAS	IGF2	P2RX3	CAVIN3	ASRGL1	ATM	BIRC2	DHCR7	GAB2	HMBS	MMP3	MMP7	P2RX3	PGR	SCGB2A1	ADIPOR2	CHD4	KRAS	CDK2	IGF1	MDM2	POLE	SH2B3	UBC	VDR	ZNF605	BRCA2	DCT	ERCC5	IRS2	PCDH17	RB1	AKT1	BMP4	CIDEB	DICER1	ESR2	HIF1A	MLH3	SNX6	XRCC3	B2M	BLM	CYP19A1	CYP1A1	IGF1R	RAD51	SYNM	TICRR	TJP1	TTC23	VPS13C	ZSCAN29	CIITA	CREBBP	ERCC4	HCFC1R1	IL32	PALB2	PDPK1	SULT1A1	CDH1	CTCF	E2F4	FTO	HSD17B2	MMP2	NOD2	TERF2	WWOX	ZFHX3	SHBG	SMYD4	SREBF1	TP53	ACE	BIRC5	BRCA1	BRIP1	CDK12	ERBB2	CYGB	HNF1B	MAPT	NF1	NME1	RAD51C	RAD51D	SPOP	BCL2	C18orf21	MC4R	MOCOS	NEDD4L	SMAD4	MUC16	DNMT1	INSR	KEAP1	NFIC	PIK3R2	RETN	SMARCA4	AKT2	APOE	BAX	CCNE1	DYRK1B	ERCC1	ERCC2	LIG1	PAK4	NOP53	PLAUR	POLD1	PPP2R1A	TGFB1	URI1	XRCC1	BMP2	ASXL1	AURKA	BMP7	DNMT3B	EYA2	MMP9	SNAI1	ZNF217	GABPA	TFF3	BIK	CBY1	CHEK2	CHEK2	COMT	GSTT1	ATP6AP2	FOXP3	ARMCX4	DACH2	HPRT1	TAF1	XIST	MT-CO2	MT-ND1
```

### Training

```bash
chmod +x training.sh
chmod +x training-gene-set.sh
sbatch training.sh
sbatch training-gene-set.sh
```

### Perturb 
```bash
chmod +x perturb.sh
chmod +x perturb-gene-set.sh
sbatch perturb.sh
sbatch perturb-gene-set.sh
```
Slurm Output:
```
visulize the perturbed data and real data
The mmd score is:0.001988342836775514
filter the perturbed gene -- 1 std
The real data mean [ 0.09740872 -0.6523077  -0.4439191  -0.3074455  -0.02253245  0.10641403
 -0.4433017  -0.42429057 -0.22459188  0.12604919 -0.2148573  -0.24754134
 -0.0072906  -0.6662008  -0.38933805 -0.16208525  0.13997407 -0.22331703
 -0.07359176  0.11816438]-- script_util
The real data std [0.05424322 0.09054648 0.18595758 0.15146753 0.12660989 0.08707687
 0.23069184 0.14454235 0.07099457 0.08162941 0.15559958 0.09192944
 0.03602436 0.09241644 0.12002992 0.10605033 0.13680498 0.07040253
 0.19098465 0.08978657]-- script_util
The perturb data mean [ 0.09319391 -0.62166816 -0.46046373 -0.46447265  0.01576618  0.02810708
 -0.57444614 -0.15925665 -0.19866227 -0.02563703 -0.1243693  -0.16826321
  0.01072228 -0.56983507 -0.53796774 -0.18446513 -0.01035629 -0.18531884
  0.1295753   0.08703887]-- script_util
The perturb data std [0.05536953 0.07519238 0.18305118 0.12617977 0.121355   0.05413768
 0.22112286 0.13434984 0.05387673 0.07199974 0.09225128 0.08927536
 0.04224531 0.11621472 0.15383069 0.10799643 0.06942936 0.0512172
 0.21273297 0.08261545]-- script_util
The differences between real and perturb data [0.02867874 0.05520887 0.05367832 0.16876146 0.0701891  0.08549368
 0.14985566 0.26503396 0.04979801 0.15173659 0.0931955  0.09258598
 0.03192927 0.10325675 0.15843865 0.06299698 0.15038326 0.0503264
 0.21243359 0.05858151] -- script_util
The standard deviation between real and perturb data data 0.06262937933206558 -- script_util
The mean between real and perturb data data 0.10462810844182968 -- script_util
The real data [-0.3074455  -0.42429057 -0.07359176]-- script_util
The perturb data [-0.46447265 -0.15925665  0.1295753 ]-- script_util
The perturbation percentages between real and perturb data data [-0.548915  -0.624652  -2.8866491]-- script_util
The indentified genes are: Index(['IL1R2', 'IRS2', 'ESR1'], dtype='object') -- 1 standard deviation of the perturbation among all 20 gene
```
#### Visualization
UMAP Plot output
![image](https://github.com/user-attachments/assets/0f4ccd8c-57b4-4e22-9a6f-a56c7f52e360)


# Gene Network 

## STRING Network Generator 

This script generates a gene interaction network file from the STRING database,
suitable for use with the Perturb visualization tool.

Usage:
```
  python generate_string_network.py --genes perturbed_genes.txt --degs deseq_results.tsv --output string_network.tsv
```
Run the script with your perturbed genes and DEG file:
```
python generate_string_network.py --genes "TMEM184A,PDE1C,GPR146,FAM110B" --degs your_degs_file.tsv --output gemdiff_network.tsv
```
Alternatively, if you have your perturbed genes in a file (one per line):
```
python generate_string_network.py --genes perturbed_genes.txt --degs your_degs_file.tsv --output gemdiff_network.tsv
```
You can adjust several parameters to fine-tune your network:
```
--score: Confidence threshold (0-1000) for STRING interactions (higher = more confident)
--additional: Number of additional interacting genes to include (0 for only direct interactions)
--species: Change to 10090 for mouse, 10116 for rat, etc.
```
## Perturbed Genes Visualization
Perturbed genes visualization with differentially expressed genes
This script visualizes perturbed genes identified by GEMDiff in a network context,
highlighting connections between perturbed genes and top differentially expressed genes.

