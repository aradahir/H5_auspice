# Nextstrain Build for H5Nx Pipeline

A Snakemake-based pipeline for automating phylogenetic analysis and visualization of H5Nx influenza sequences via the Nextstrain platform.

## About This Workflow

This pipeline creates phylogenetic trees using Augur tools by aligning sequences to a selected reference sequence. Key features include:

- **Automated clade assignment**: Clade information is automatically retrieved during pipeline execution using Nextclade, ensuring up-to-date classification
- **Phylogenetic inference**: Tree building is performed using MAFFT for alignment and IQ-TREE for phylogenetic reconstruction
- **Time-resolved analysis**: Temporal reconstruction using TreeTime to infer divergence dates
- **Ancestral trait inference**: Internal node traits are inferred using Augur's `traits` command, allowing annotation of ancestral nodes with clade information
- **Interactive visualization**: Outputs are formatted for visualization on [auspice.us](https://auspice.us/)

## Dependencies

Ensure the following tools and libraries are installed:

- `biopython=1.85`
- `iqtree=2.4.0`
- `mafft=7.525`
- `nextclade=3.12.0`
- `nextstrain-augur=29.0.0`
- `phylo-treetime=0.11.4`
- `pandas=2.2.3`
- `snakemake`
- `csvtk`

## Installation and Setup

Follow these steps to clone and install the repository (Linux/macOS):
```bash
# Clone the repository
git clone git@github.com:aradahir/H5_auspice.git
cd H5_auspice/

# Create conda environment
conda env create -f requirements.yml

# Activate the environment
conda activate nextstrain_build
```

## Usage

### 1. Prepare Input Files

Place your input files into the `./data/` directory:
```
./data/
├── filename.fasta           # Aligned or unaligned sequences
└── meta_filename.tsv        # Metadata file with strain information
```

**Metadata format**: The metadata file should be tab-separated with at minimum the following columns:
- `strain`: Sequence identifier (must match FASTA headers)
- `date`: Collection date in format `DD/MM/YYYY`
- `Host`: Host species
- `Location`: Geographic location
- Additional optional columns: `clade`, `NA_deletion`, `PB2_D701_mutation`, `Submitting_Lab`

### 2. Configure the Pipeline

Edit `config/config.yaml` to adjust pipeline parameters:
```yaml
# Example configuration
clade_calling: "H5Nx (all clades)"
alignment_option: "full_seq"  # or "trimmed"
reference_name: "A/bobcat/Wisconsin/22-016051-001/2022"
clock_rate: 0.00249
clock_std_dev: 0.000329
```

### 3. Run the Pipeline
```bash
# Activate environment
conda activate nextstrain_build

# Dry run to check workflow
snakemake -n

# Run pipeline (adjust -j for number of threads)
snakemake -j 4
```

### 4. Visualize Results

View your phylogenetic tree using Nextstrain:
```bash
# Start local Auspice server
nextstrain view auspice/

# Or upload JSON files to https://auspice.us/
```

## Output Structure
```
├── data/
│   ├── filename.fasta
│   ├── meta_filename.tsv
│   └── filename_metadata.tsv        # Processed metadata
│
├── nextclade/
│   └── H5Nx/
│       └── filename/
│           └── clade.tsv            # Nextclade output with clade assignments
│
├── results/
│   └── H5Nx/
│       ├── auspice/
│       │   ├── filename_auspice.json                  # Main Auspice visualization file
│       │   ├── filename_auspice_root-sequence.json
│       │   ├── filename_auspice_tip-frequencies.json  # Frequency dynamics
│       │   └── filename_frequencies.json
│       │
│       └── tree/
│           ├── filename_aligned.fasta                  # Aligned sequences
│           ├── filename_tree.nwk                       # Time-resolved tree
│           ├── filename_tree_raw.nwk                   # Raw phylogenetic tree
│           ├── filename_aa_muts.json                   # Amino acid mutations
│           ├── filename_nt_muts.json                   # Nucleotide mutations
│           ├── filename_branch_lengths.json            # Branch length data
│           ├── filename_traits.json                    # Ancestral trait reconstruction
│           └── filename_reference.gb                   # Reference GenBank file
│
└── auspice/
    ├── filename_auspice.json                          # Copy for visualization
    └── filename_auspice_tip-frequencies.json
```

## Configuration Files

### `config/config.yaml`
Main configuration file containing:
- Clade calling options
- Reference sequence name
- Molecular clock parameters
- Alignment options

### `config/auspice_config.json`
Auspice visualization settings:
- Color schemes
- Geographic resolutions
- Display defaults
- Metadata filters

## Important Notes

- **Input requirements**: This pipeline is designed for influenza sequences with known subtypes (concatenated segments)
- **File naming**: Ensure FASTA filenames match the prefix in metadata filenames (e.g., `sample.fasta` and `meta_sample.tsv`)
- **Before rerunning**: Clean previous outputs to avoid conflicts:
```bash
  rm -rf results/ nextclade/ auspice/
  rm data/*_metadata.tsv data/*_trimmed.fasta
```

## Workflow Steps

1. **Clade assignment**: Nextclade assigns clades to sequences
2. **Metadata processing**: Combines metadata with clade information
3. **Sequence alignment**: MAFFT aligns sequences to reference
4. **Tree inference**: IQ-TREE builds maximum likelihood tree
5. **Time-tree reconstruction**: TreeTime infers divergence dates
6. **Ancestral reconstruction**: Augur infers mutations and traits
7. **Frequency estimation**: Calculates clade frequencies over time
8. **Export**: Formats output for Auspice visualization

## Troubleshooting

**Error: "No records found in handle"**
- Check that FASTA files are properly formatted
- Verify reference sequence name matches exactly

**Error: Port already in use**
- Use a different port: `nextstrain view --port 4001 auspice/`

**Empty frequency plots**
- Check date formatting in metadata (should be DD/MM/YYYY)
- Ensure sufficient temporal sampling

## Citation

If you use this pipeline, please cite:
- [Nextstrain](https://nextstrain.org/)
- [Augur](https://docs.nextstrain.org/projects/augur/)
- [Nextclade](https://clades.nextstrain.org/)

## Contributors

- Arada Hirankitti ([@aradahir](https://github.com/aradahir))
- Clyde Dapat ([@clyde-dapat](https://github.com/clyde-dapat))
- Michelle Wille ([@michellewille2](https://github.com/michellewille2))

## Contact

For questions or issues, please open an issue on GitHub or contact the maintainers.
