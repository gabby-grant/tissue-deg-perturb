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

```
chmod +x training.sh
chmod +x training-gene-set.sh
```

### Perturb 
#### Visualization

# Gene Network 
String network vizualization
Perturbed genes visualization with differentially expressed genes
