Nextstrain_build for H5Nx pipeline
A Snakemake-based pipeline for automating phylogenetic analysis and visualization via the Nextstrain platform.

About this workflow
Phylogenetic trees are created using augur tools by aligning sequences to the selected reference sequence. Clade information are automatically retrieved during pipeline execution using nextclade, ensuring up-to-date clade assignment. Tree building is performed using MAFFT and IQ-TREE, followed by time-resolved reconstruction using timetree. Internal node traits are inferred using the traits command in augur, allowing annotation of ancestral nodes with clade information.

Dependencies
Ensure the following tools and libraries are installed:

biopython=1.85
iqtree=2.4.0
mafft=7.525
nextclade=3.12.0
nextstrain-augur=29.0.0
phylo-treetime=0.11.4
pandas=2.2.3
Installation and Setup
Follow these steps to clone and install the repository (Linux-based systems):

```
git clone ssh https://git@github.com:aradahir/H5_auspice.git
cd .\H5_auspice\
conda env create -f requirements.yml

```
Activate the snakemake environment for running the pipeline.

```
conda activate nextstrian_build

```
Usage
Prepare Input Place your .fasta and metadata files into the ./data/ directory: ./data/ ├── filename.fasta └── meta_filename.tsv

Open the environment

 ```
 conda activate nextstrain_build
 ```
Run the pipeline - thread can be adjusted by changing from 1 into the specific number of thread

 ```
 snakemake -j 1
 ```
Output Structure
├── data/
│    # Input sequences and metadata
├── nextclade/
│   └── H5Nx/
│   	└── filename/
│           # tsv with clades and metadata together
├── results/
    └── H5Nx/
       ├── auspice/
       │   ├── filename_auspice.json : json files using for visualize the data in https://auspice.us/
       │   └── filename_auspice_root-sequence.json
       │   └── filename_auspice_tip-frequencies.json
       │   └── filename_auspice_frequencies.json
       │   └── filename_auspice_tree.json
       └── tree/
           ├── filename_aa_muts.json
           ├── filename_aligned.fasta
           ├── filename_aligned.fasta.insertion.csv
           ├── filename_branch_length.json
           ├── filename_nt_muts.json
           ├── filename_traits.json
           ├── filename_traits.json
           ├── filename_tree.nwk
           └── filename_tree_raw.nwk


Important Note

This pipeline is designed for assembled the concatenated segmented influenza sequences with known subtypes.

Before rerunning the pipeline, please remove the following: results/, data/, formatted_metadata.tsv
