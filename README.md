# euTyper

**Work In Progress**

[Chewie](https://github.com/B-UMMI/chewBBACA) steps into the Eukarya Domain

euTyper is a tool to create gene-by-gene typing schemas and perform allele calling for eukaryotes.

## Dependencies

Python >= 3.7  
[BioPython](https://github.com/biopython/biopython)  
[Augustus](https://github.com/Gaius-Augustus/Augustus)  
[BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279671/)
[Parallel](https://www.howtoinstall.me/ubuntu/18-04/parallel/)

## Usage

```
-h, --help      show this help message and exit
-i INPUT_FILES  Path to a FASTA file with contig sequences for a single
                genome.
-o OUTPUT_PATH  Path to the directory that will store output files created
                by the process.
-s sample_metadata path to samplemetadata (TSV file containing fastafile name, species, codon number.
--bsr BSR       BLAST Score Ratio value used as threshold.
--t THREADS     Number of threads used by BLASTp to align translated genes.
--p probabbility probability change of a gene found by Augsustus, used as an cutoff


```
### exstra functionality for flexibility and not needing to rerun euTyper witht the smalles changes
#### add loci
```
--add_loci path path to fasta file containing the CDS of each gene in it. This parameters enables to directly add genes to a folder containing all genes
--codon Required codon table number to calculate the AA sequence of the CDS
--force Forces the genes from the add_loci parameter to be directly added to the clustering dir (genes will have a forced prefix)
-o OUTPUT_PATH  Path to the directory that will store output files created
                by the process.
```                
#### add genomes
```
--add_gene path to directory containing the fastafiles of genomes to be added to processed genomes
-s sample_metadata path to samplemetadata (TSV file containing fastafile name, species, codon number.
-o OUTPUT_PATH  Path to the directory that will store output files created
                by the process.
-recluster If recluster = True, after adding the gene, all known genes will be directly evaluated by Clustering and blasting to retrieve the 
--bsr BSR       BLAST Score Ratio value used as threshold.
--t THREADS     Number of threads used by BLASTp to align translated genes.
--p probabbility probability change of a gene found by Augsustus, used as an cutoff

```

## example commant with test files

### normal run
```
python3 euTyper.py -i /path/to/your/favourite/directory/containing/genomes -o /output/folder/path/of/directory/be/created/and/files/stored -s /path/to/sample_metadata --t 30 --p 0.9
```
### incase euTyper already ran, and some adjustments are required
### add loci directly to scheme seed
```
python3 euTyper.py -add_loci /path/to/your/favourite/directory/containing/genes -o --t 30 --codon 12 --forced True 
```
### If genome to list of genes
```
python3 euTyper.py -add_gene /path/to/your/favourite/directory/containing/genomes -o /output/folder/path/of/directory/be/created/and/files/stored -s /path/to/sample_metadata --t 30 --p 0.7 -recluster True

```
### generating the scheme seed again from already processed genes
```
python3 euTyper.py --reclyster -o  /output/folder/path/of/directory/be/created/and/files/stored --t 20
