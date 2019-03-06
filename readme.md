# Germline Variant Filter

## Introduction

A Python program for filtering variants based on a variety of features. Produces a TSV file describing the variants found for a particular sample. Avoids parsing VCFs by first converting the file to a CSV and doing filtering using Pandas.


## When should I use this program?

- When you have multisample VCFs containing samples which can be analysed as a trio or samples which should be analysed on their own.

## Install

Install Miniconda from [1].

Install all dependencies using YAML config file.

`conda env create -f env/germline_variant_filter.yaml `


## Test 


`python -m unittest test.tests`

## Run

### Data Preparation

The following programs are required for the data preperation. It is reccomended you use bioconda to install these.

- vt = 0.57721
- ensembl-vep = 94.5
- gatk4 = 4.0.11.0

The following annotation files are also required:

- Refseq VEP Cache: ftp://ftp.ensembl.org/pub/release-94/variation/VEP/homo\_sapiens\_refseq\_vep_94\_GRCh37.tar.gz
- CCR BED Files: Merge both autosomal and X BEDs into a single file. https://github.com/quinlan-lab/ccr
- GnomadGenomes: ftp://ftp.ensembl.org/pub/data\_files/homo\_sapiens/GRCh37/variation_genotype/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz
- GnomadExomes: https://storage.googleapis.com/gnomad-public/release-170228/vcf/exomes/gnomad.exomes.r2.0.1.sites.vcf.gz
- SpliceAI VCF exomes: https://github.com/Illumina/SpliceAI

The program requires a TSV file for input. The TSV file needed can be generated from a generic VCF using the following commands:

```
# split multiallellics and normalise
cat input.vcf | vt decompose -s - | vt normalize -r reference.fasta - > input.norm.vcf

# Annotate with VEP
vep --verbose --format vcf --everything --fork 1 --species homo_sapiens --assembly GRCh37 --input_file input.norm.vcf \
--output_file input.norm.vep.vcf --force_overwrite --cache --dir vep_cache_location \
--fasta reference.fasta --offline --cache_version 94 -no_escape --shift_hgvs 1 --exclude_predicted --vcf --refseq --flag_pick \
--custom gnomad.genomes.vcf.gz,gnomADg,vcf,exact,0,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH,AF_POPMAX  \
--custom gnomad.exomes.vcf.gz,gnomADe,vcf,exact,0,AF_POPMAX \
--custom ccrs.bed.gz,ccrs,bed,overlap,0 \
--custom exome_spliceai_scores.vcf.gz,SpliceAI,vcf,exact,0,DS_AG,DS_AL,DS_DG,DS_DL,SYMBOL

# Convert to CSV using GATK
gatk VariantsToTable -V input.norm.vep.vcf -O input.norm.vep.csv -F CHROM -F POS -F REF -F ALT -F ID -F QUAL -F FILTER -F CSQ -F AC -GF GT -GF GQ -GF DP

```

### Run Program

To get help type:

```
python germline_variant_filter.py -h

```
#### Config File

There is a configuration file that is required to run the program. An example of this is seen in the config/ directory. The comments in this file described what each variable does.

#### Example

```
# Get CSQ string descrbing VEP fields
CSQ=$(grep "^##INFO=<ID=CSQ" input.norm.vep.vcf | awk 'BEGIN { FS = ":" } ; { print $2 }' | tr -d '>" ')

python germline_variant_filter.py --config config.yaml --ped pedigree.ped \
--input input.norm.vep.csv \
--panelapp --local-panel-app-dump panelapp.csv --spliceai --add-ccrs \
--gnomad-constraint-scores  --worksheet my_worksheet --results results --csq $CSQ --smart-synonymous

```

#### Option Description

  - config: Filepath to YAML config file. See config/ directory for an example of what these look like.
  - ped: Filepath to PED file. The PED file should describe the family relationships between samples.
  - input: Filepath to the input CSV file.
  - panelapp: Whether to add PanelApp annotations.
  - local-panel-app-dump: Filepath local panelapp store.
  - csq: The VEP CSQ string. For example Allele|Consequence|IMPACT|SYMBOL...
  - spliceai: Attempt to parse SpliceAI annotations added by VEP. 
  - smart-synonymous: Smart synonymous variant filtering. See readme for details.
  - add-ccrs: Add the CCR annotations from VEP. 
  - gnomad-constraint-scores: Add the per gene gnomad constraint scores.
  - patient-hpos: Filepath to file containing patient HPO terms. See examples/ directory for information on the format of this file.
  - worksheet: The worksheet ID.
  - results-dir: Where to put the results.

## Algorithm

### Stage 1 - Quality Filter

- Keep variants which have a FILTER status of PASS or where the FILTER status is empty.

### Stage 2 - Initial Frequency Filter

- Remove all variants which have a default\_cutoff\_gnomad\_genomes of below 0.01 and a default\_cutoff\_gnomad\_exomes of below 0.01
- These values can be changed in the config file.

### Stage 3 - Consequence Filter

- Calculate the worst consequence for each variant. The order of severity is defined in the config file.
- If the worst consequence is in the to_keep_consequences dictionary then we keep that variant otherwise we discard.
- Optionally apply smart synoymous filtering where we remove synoymous variants unless they are listed as pathogenic in clinvar or have a predicted affect on splicing. See Smart Synonymous Filtering.

### Stage 4 - Per Sample Filter

#### Stage 4 A

- For each sample look at the variants which are present in that sample and where the genotype passes the DP and GQ filters. The cutoffs can be specified in the config file.

#### Stage 4 B - Trio

- For samples which are in a trio annotate the workflow. This can be one of:

- UNIPARENTAL\_ISODISOMY - Variant is on an Autosome and the genotype is Homozygous for the ALT allele and (the mother is Heterozygous and the father is Homozygous Reference OR  mother is Homozygous Reference and father is Heterozygous). Alternatively, Variant is on the X chromosome and the sample sex is Female and the genotype is Homozygous for the ALT allele and (the mother is Heterozygous and the father is Homozygous Reference OR mother is Homozygous Reference and father is homozygous for the ALT allele.) All samples must pass the minimum DP and GQ requirements set in the config file.
- COMPOUND\_HET - The variant is either on an Autosome or (on the X chromosome and sample sex is Female) and the genotype is Heterozygous and there is more than one variant in the transcript.
- MITOCHONDRIAL - Variant is on the MT chromosome.
- RECCESSIVE\_X\_FEMALE - Variant is on the X chromosome and the genotype is Homozygous for the ALT allele and the sample sex is Female.
- RECCESSIVE\_AUTOSOMAL - Variant is on an Autosome and the genotype is Homozygous for the ALT allele.
- X\_LINKED_MALE - Variant is on the X chromosome and sample sex is Male.
- Y\_LINKED_MALE - Variant is on the Y chromosome and sample sex is Male.
- DOMINANT\_AUTOSOMAL - Variant is on an Autosome and genotype is Heterozygous.
- DOMINANT\_X\_FEMALE - Variant is on the X chromosome and the sample sex is Female and the genotype is Heterozygous.
- DENOVO_HC - Variant is in the proband and not in either of the parents. All samples pass mimimum DP and GQ requirements set in the config file.
- DENOVO_LC - Variant is in the proband and not in either of the parents. Only proband sample needs to pass mimimum DP and GQ requirements set in the config file.


#### Stage 4 B - Single

- For samples which are by themselves annotate the workflow as:

- COMPOUND\_HET
- MITOCHONDRIAL,
- RECCESSIVE\_X\_FEMALE
- RECCESSIVE\_AUTOSOMAL
- X\_LINKED\_MALE
- Y\_LINKED\_MALE
- DOMINANT\_AUTOSOMAL
- DOMINANT\_X\_FEMALE


#### Stage 4 C


- Filter each variant on both the frequency in gnomad and the within run allele count (AC) depending on which workflow the variant fits into. The cutoff values can be specified in the config file.
- If a variant matches multiple workflows then use the least restrictive one.


#### Stage 4 D

- Write variants to file

## Additional Features

### Smart Synonymous Filtering

When this option is selected synoymous variants will be filtered out unless they are listed as being pathogenic in Clinvar or have a predicted affect on splicing over the threshold specified in the config file. Synoymous variants with no availible splicing prediction will also be kept.

### Patient HPOs

Give the program a TSV file containing HPO terms for each patient. Annotate each variant with a count of how many HPO terms match the gene the variant is found in. Requires the ALL\_SOURCES\_ALL\_FREQUENCIES\_genes\_to\_phenotype.txt file found at https://hpo.jax.org/app/download/annotation

### Gnomad Constraint Scores

Add the Gnomad Constraint Score to the output. Requires the release\_2.1\_ht\_constraint\_constraint\_refseq.txt which can be found on the cluster. https://gnomad.broadinstitute.org/downloads

## Help

- The config file contains the final\_fields\_trio and final\_fields\_single variables these contain the fields which the final output TSV file should contain.

- Pandas will give the SettingWithCopyWarning - this can be ignored we want to be setting on a copy.

## Known Limitations

- No phasing of compound HETs

## References

[1] https://conda.io/en/latest/miniconda.html









