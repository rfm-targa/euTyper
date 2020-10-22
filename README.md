# euTyper

**Work In Progress**

[Chewie](https://github.com/B-UMMI/chewBBACA) steps into the Eukarya Domain

euTyper is a tool to create gene-by-gene typing schemas and perform allele calling for eukaryotes.

## Dependencies

Python >= 3.7  
[BioPython](https://github.com/biopython/biopython)  
[Augustus](https://github.com/Gaius-Augustus/Augustus)  
[BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279671/)  

## Usage

```
-h, --help      show this help message and exit
-i INPUT_FILES  Path to a FASTA file with contig sequences for a single
                genome.
-o OUTPUT_PATH  Path to the directory that will store output files created
                by the process.
-s SPECIES      Species identifier passed to AUGUSTUS.
--bsr BSR       BLAST Score Ratio value used as threshold.
--t THREADS     Number of threads used by BLASTp to align translated genes.
```
