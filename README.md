

# Pathway tools / mpwt wrapper

Wrapper for mpwt (https://github.com/AuReMe/mpwt) to take a number of gff annotated bacterial genomes, with representative protein sequences from a pan-genome analysis annotated using eggnog-mapper (http://eggnog-mapper.embl.de/), and run them through pathway tools (http://bioinformatics.ai.sri.com/ptools/) to predict presence/absence of metabolic pathways.

## Dependencies

* Pathway tools - http://bioinformatics.ai.sri.com/ptools/
* mpwt - https://github.com/AuReMe/mpwt

### Python packages
* Biopython
* tqdm

### R packages
* stringr
* tidyr
* dplyr


## Example

### Overview

Given a number of bacterial genomes of known species, we want to get a nice table of presence/absence of metabolic pathways.

We have three _S. marcescens_ annotated genomes in gff format that we will use as input.

In order to to do this, the following steps will need to be performed:
1.  Get representative sequences of all genes across all strains (using `panaroo`)
2.  Get functional annotation of these genes to plug into metabolic reaction prediction (e.g. EC numbers, GO terms), using `eggnog-mapper`
3.  Apply functional annotation to each annotated genome and rewrite out as genbank files in the directory hierarchy needed for `mpwt` to run
4.  Create a tsv file with the strain name and the NCBI taxon number for the strain - for `mpwt`
5.  Run pathway tools by using `mpwt`
6.  Collate all presence/absence of pathways within each genome


### Worked example

1.  Get representative sequences of all genes across all strains (using `panaroo`)

Firstly we will take the _S. marcescens_ genomes that are in the example folder, and then perform a pan genome analysis using `panaroo`.

So, running from a terminal in the example folder, we will run panaroo, using four threads and `--clean-mode moderate`

```
panaroo -i 1_assemblies/*.gff -o 2_panaroo -t 4 --clean-mode moderate
```


2.  Get functional annotation of these genes to plug into metabolic reaction prediction (e.g. EC numbers, GO terms), using `eggnog-mapper`

Panaroo will output a fasta file of representative genes from the pan genome: `pan_genome_reference.fa`. This will be used as input to eggnog-mapper, which can either be run locally or using the webservice http://eggnog-mapper.embl.de/. `pan_genome_reference.fa` can either be used directly as input as nucleotide gene sequences (which will then be translated by eggnog-mapper) or translated manually beforehand, e.g. using `fastaq translate`.

Eggnog mapper will produce a number of files, of which we need `out.emapper.annotations`.



3.  Apply functional annotation to each annotated genome and rewrite out as genbank files in the directory hierarchy needed for `mpwt` to run

The directory structure (and files within) should now look something like this (using `tree` to show the recursive directory lists)

```
.
├── 1_assemblies
│   ├── SJC1033.gff
│   ├── SJC1036.gff
│   └── SJC1039.gff
├── 2_panaroo
│   ├── combined_DNA_CDS.fasta
│   ├── combined_protein_cdhit_out.txt
│   ├── combined_protein_cdhit_out.txt.clstr
│   ├── combined_protein_CDS.fasta
│   ├── final_graph.gml
│   ├── gene_data.csv
│   ├── gene_presence_absence.csv
│   ├── gene_presence_absence_roary.csv
│   ├── gene_presence_absence.Rtab
│   ├── pan_genome_reference.fa
│   ├── pan_genome_reference.faa
│   ├── pre_filt_graph.gml
│   ├── struct_presence_absence.Rtab
│   └── summary_statistics.txt
├── 3_eggnog
│   └── out.emapper.annotations
```

From this base directory, we will now run `gffs2gbks.py` which will take these gff assembly files, along with the `eggnog` annotation data and re-write the gff files into pathologic format (.pf) and fasta sequences of each contig in the assembly, creating a new directory `4_mpwt_input`, where these files will be put. We will use 1 thread as there are only three assemblies, but if processing hundreds/thousands of strains, then more may be better to use.

```
python3 ../gffs2gbks.py -g 1_assemblies/*.gff -p 2_panaroo/gene_presence_absence.csv -e 3_eggnog/out.emapper.annotations -o 4_mpwt_input -t 1
```


4.  Create a tsv file with the strain name and the NCBI taxon number for the strain - for `mpwt`

We now need to put a file named `taxon_id.tsv` into `4_mpwt_input` with two tab delimited columns: 1 - the strain name in the gff file, and 2 - the NCBI taxon id (https://www.ncbi.nlm.nih.gov/taxonomy/). For this example, the contents of the file is shown below, and can be copied into the appropriately named file as above.

```
species	taxon_id
SJC1033	615
SJC1036	615
SJC1039	615
```


5.  Run pathway tools by using `mpwt`

We can now run pathway tools by using this directory as input to `mpwt`.

```
mpwt -f 4_mpwt_input/ -o 5_mpwt_output --patho --cpu 4 -v
```


6.  Collate all presence/absence of pathways within each genome

Providing step 5 has run without any errors (errors for individual genomes will be put into the respective genome folder within `4_mpwt_input`), we can now compile all the predicted presence/absence of pathways into one file using the Rscript `compare_pgdbs.R`


```
Rscript ../compare_pgdbs.R 5_mpwt_output
```

This will output a tsv file of all metabolic pathways present in these assemblies where there is _only_ gene evidence (not species level evidence)
