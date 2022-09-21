# kGWAS_watkins

## Description
kGWAS\_watkins is an extension of [AgRenSeq_GLM](https://github.com/kgaurav1208/AgRenSeq_GLM) to identify resistance genes in Watkins wheat diversity panel using R-gene enrichment sequencing and multiple reference genomes. 

### 1: k-mer presence/absence matrix construction

- Count k-mers from either trimmed or raw read files for each accession using jellyfish.

```
zcat accession1_R?.fastq.gz | jellyfish count -C -m 51 -s 3G -L 10 -o accession1.jf /dev/fd/0
jellyfish dump -ct accession1.jf > accession1.dump.txt
```

- Create a configuration file, `jellies.cfg`, containing the name of the
accession, tab-separated by the path to its jellyfish

```
accession1	path/to/accession1.jf
accession2	path/to/accession2.jf
...
```
- Create a header file `presenceMatrix_header.txt` from the first column of `jellies.cfg`. 

```
accession1
accession2
...
```

- Run the following script for the jellyfish dump file of each accession in parallel. The output is a tab-separated file with k-mers in the first column and a string of `1`s and `0`s (indicating k-mer's presence and absence, respectively, in each accession) in the second column.

```
python create_matrix.py -a accession1 -j accession1.dump.txt -c jellies.cfg -o presenceMatrix_accession1.txt -mc 4
```  



Parameter (long version)| Parameter (short version) | Argument | Description
--- | --- | --- | ---
--accname | -a | string | Mandatory. Name of the accession.
--jfdump | -j | filepath | Mandatory. Jellyfish dump file of the accession.
--config | -c | filepath | Mandatory. Configuration file containing the name of each accession tab-separated from the path to the jellyfish of that line.
--kmersize | -k | int | Default 51. Kmer size specified while making jellyfish.
--mincount | -mc | int | Default 2. Only those k-mers are retained which are present/absent in more than this number of accessions.
--output | -o | filepath | Mandatory. Output k-mer matrix file.
  
- The output of above script for each accession can be concatenated, sorted and re-split into similar-sized chunks followed by compression using `gzip`. 


The matrix chunks and the header file obtained after applying these steps on Watkins RenSeq data are available at [Zenodo](https://zenodo.org/record/5557564).

### 2: Subsampling k-mers to compute PCA and correct for population structure in GWAS 

Subsample a set of _N_ (eg: 5000) k-mers from the entire k-mer matrix :

````
zcat presenceMatrix.txt.gz | shuf -n 5000 > subsampledMatrix.txt
````

Run the following script to transform the set of subsampled k-mers into a tab-separated form for PCA computation:

````
python transform_subsampled_kmers.py subsampledMatrix.txt presenceMatrix_header.txt subsampled_kmers.tsv
````

A table with presence/absence of ~5000 k-mers subsampled from Watkins k-mer matrix, and used in PCA computation and phylogenetics, is available [here](https://github.com/arorasanu/watkins_renseq/blob/master/phylogeny/watkins_nlr_kmers.tsv).

### 3: Preparing NLR assembly for k-mer mapping 

Depending on the NLR assembly used to map k-mers, follow either (a) or (b):

#### (a) Extract NLR sequences from a reference genome
Predict NLRs in a reference genome using [NLR Annotator](https://github.com/steuernb/NLR-Annotator). For each predicted NLR region, extract its sequence along with 3kb sequence (if available) from both upstream and downstream region to create a reference NLR assembly. For example, if an NLR is predicted on Chromosome 1 from 15kb to 20kb, then the relevant sequence can be extracted using:
 
````
samtools faidx reference.fasta Chr1:12000-23000
````


#### (b) Assign coordinates to NLR contigs of a non-reference RenSeq assembly by anchoring to a reference genome 
Predict NLR contigs in a non-reference RenSeq assembly using [NLR Parser](https://github.com/steuernb/NLR-Parser) and extract the NLR contigs to create a non-reference NLR assembly. Map each of these non-reference NLR contigs to a reference genome using minimap2 and retain the longest hit:

````
minimap2 -f 0.01 reference.fasta accession_nlr.fasta > mmap.paf
````

````
awk -v OFS='\t' -F'\t' '{print $1,$4-$3,$6,$8,$9}' mmap.paf | sort -k1,1 -k2,2nr |sort -u -k1,1|sort -k3,3 -k4,4n > mapping.txt
````

Rename the contigs of non-reference NLR assembly to reflect the mapping coordinates with respect to the reference genome:
 
````
python anchor.py -a accession_nlr.fasta -m mapping.txt -o accession_nlr.anchored.fasta
````

Parameter (long version)| Parameter (short version) | Argument | Description
--- | --- | --- | ---
--assembly | -a | filepath | Mandatory. Path to assembly fasta file being anchored.
--mapping | -m | filepath | Mandatory. File containing the unique reference mapping coordiantes of each scaffold.
--output | -o | filepath | Mandatory. Output fasta file.

### 4: Create phenotype file

The phenotype file should be in a tab-separated two column format with the accession names in the first column and the corresponding numerical phenotype scores in the second column. When there are multiple reps for an accession, the average phenotype score is used. 

```
accession1	average_phenotype_score_accession1
accession2	average_phenotype_score_accession2
...
```


### 5: Generate association scores of _k_-mers and project onto the chopped assembly

Run the following command for each chunk of the _k_-mer matrix:

```
python RunAssociation_GLM.py -i presenceMatrix.txt.gz -hd presenceMatrix_header.txt -a nlr_assembly.fasta -p phenotype.txt -s subsampled_kmers.tsv  -u usable.txt -o output.txt -c 0.1 -pv 2
```
Parameter (long version)| Parameter (short version) | Argument | Description
--- | --- | --- | ---
--inputmatrix | -i | filepath | Mandatory. The path to file containing the gzipped version of presence/absence matrix of k-mers.
--header | -hd | filepath | Mandatory. The path to file containing the list of accessions in the order that their presence is scored in the presence/absence matrix.
--assembly | -a | filepath | Mandatory. The path to assembly onto which k-mers are mapped for plotting.
--phenotype | -p | filepath | Mandatory. The path to phenotype file.
--usable | -u | filepath | Optional. The path to file containing the list of usable accessions.
--snp | -s | filepath | Mandatory. The path to file containing tab-separated table with either SNP markers or subsampled k-mers to compute PCA and correct for population structure.
--subsample | -sub | integer | Optional. Run association analysis by taking a random subsample of _this size_ from the given accessions.
--permute | -per | | Optional. Permute the phenotype scores.
--pcadimensions | -dim | integer | Default 3. The Number of significant PCA dimensions used as covariates for regression analysis.
--correlationthreshold  | -c | float | Default 0.2. Only those k-mers whose correlation with phenotype is greater than _this value_ are retained for  regression analysis
--mincount  | -mc | int | Default 4. Only those k-mers are retained for regression analysis which are present/absent in more than this number of accessions.
--pvalthreshold  | -pv | float | Default 6.0. Only those k-mers with log10 of pvalue greater than _this value_ are retained
--output | -o | filepath | Mandatory. The path to output file to store the association values corresponding to nlr contigs in the given assembly.

Put all the resultant outputs in a directory and merge them using the following command:

```
python merge_output.py gwas_output_directory gwas_merged_results.txt
```

### 6: Bonferroni correction

The following script can be used to estimate the number of multiple tests for Bonferroni correction:

```
python knum_bonf.py -i presenceMatrix.txt.gz -hd presenceMatrix_header.txt -u usable.txt -o n_kmers.txt
```

Parameter (long version)| Parameter (short version) | Argument | Description
--- | --- | --- | ---
--inputmatrix | -i | filepath | Mandatory. The path to file containing the gzipped version of presence/absence matrix of k-mers.
--header | -hd | filepath | Mandatory. The path to file containing the list of accessions in the order that their presence is scored in the presence/absence matrix.
--usable | -u | filepath | Optional. The path to file containing the list of usable accessions.
--mincount  | -mc | int | Default 4. Only those k-mers are counted which are present/absent in more than this number of accessions.
--output | -o | filepath | Mandatory. The path to output file to store the number of k-mers present in at least `-mc` number of "usable" accessions.


### 7: Plotting
Prepare a file, `chromosome_lengths.txt`, containing the name of each chromosome in the reference genome tab-separated by the corresponding chromosome length

```
chromosome1	520000000
chromosome2	610000000
...
```

The results can then be plotted using the following script:

```
python gwas_plot.py gwas_merged_results.txt chromosome_lengths.txt gwas_plot.eps
```


## Pre-requisites

### Jellyfish
Install [Jellyfish](https://github.com/gmarcais/Jellyfish) version 2.2.6 or above with Python binding.

### Python 3 and above

The code has been tested in Python 3.6.13. 

The following Python modules are required for GWAS:

* `numpy` (tested with v1.17.0)
* `pandas` (tested with v0.23.0)
* `Biopython` (tested with v1.78): to parse assembly file
* `scikit-learn` (tested with v0.24.2): to compute PCA from SNP markers matrix
* `statsmodels` (tested with v0.12.2): for regression analysis
* `bitarray` (tested with v2.3.0)
* `matplotlib` (tested with v3.3.0): for plotting



