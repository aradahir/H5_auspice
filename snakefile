#nextstrain build protocol
import os
import datetime
import pandas as pd
import pathlib
import glob
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment

configfile: "config/config.yaml"

# Extract config values
clades_dict = config["clades"]
clade_calling = clades_dict[config["clade_calling"]]
alignment_option = config["alignment_option"]
reference_name = config["reference_name"]
header_name = config["header_name"]
clock_rate = config["clock_rate"]
clock_std_dev = config["clock_std_dev"]
gene_name = config.get("gene_name", "all")

def set_alignment(alignment_option, filename):
    """
    Prepare data for processing: 
    - option1 ('full_seq'): use all sequences as-is
    - option2 (other): only overlapped positions within all sequences
    
    Returns: filename to use for downstream processing
    """
    if alignment_option == 'full_seq':
        # Use original file without modification
        print(f"Using full sequences from {filename}")
        return filename
    else:
        # Create output filename
        output_file = str(filename).replace('.fasta', '_trimmed.fasta')
        
        # Check if trimmed file already exists
        if os.path.exists(output_file):
            print(f"Trimmed file already exists: {output_file}")
            return output_file
        
        # Check if file exists and is not empty
        if not os.path.exists(filename):
            print(f"ERROR: File {filename} does not exist!")
            return filename
        
        if os.path.getsize(filename) == 0:
            print(f"ERROR: File {filename} is empty!")
            return filename
        
        try:
            # Only use the overlapped nucleotides
            alignment = AlignIO.read(filename, "fasta")
            
            # Get alignment length and number of sequences
            aln_len = alignment.get_alignment_length()
            n_seqs = len(alignment)
        
            # Find positions present in ALL sequences (no gaps)
            shared_positions = []
            for i in range(aln_len):
                column = [alignment[j][i] for j in range(n_seqs)]
                if all(c != '-' and c != 'N' for c in column):  # No gaps or Ns
                    shared_positions.append(i)
            
            print(f"File: {filename}")
            print(f"Original length: {aln_len}, Shared positions: {len(shared_positions)}")
            
            # Check if we have any shared positions
            if len(shared_positions) == 0:
                print(f"WARNING: No shared positions found in {filename}. Using original file.")
                return filename
            
            # Trim each sequence to shared positions
            trimmed_seqs = []
            for seq_record in alignment:
                trimmed_seq = ''.join(seq_record.seq[i] for i in shared_positions)
                trimmed_seqs.append(SeqRecord(
                    Seq(trimmed_seq), 
                    id=seq_record.id, 
                    description=seq_record.description
                ))

            trimmed_alignment = MultipleSeqAlignment(trimmed_seqs)
            
            # Write trimmed alignment
            AlignIO.write(trimmed_alignment, output_file, "fasta")
            print(f"Trimmed alignment saved to {output_file}\n")
            
            return output_file
            
        except Exception as e:
            print(f"ERROR processing {filename}: {str(e)}")
            print(f"Using original file instead.")
            return filename

# Get list of FASTA files - EXCLUDE already trimmed files
filenames_paths = [f for f in pathlib.Path('./data').glob('*.fasta') 
                   if '_trimmed' not in f.stem]

# Only process files if trimming is requested AND files exist
if alignment_option == 'trimmed' and len(filenames_paths) > 0:
    print(f"Found {len(filenames_paths)} original FASTA files to process")
    processed_files = []
    for filename in filenames_paths:
        output_filename = set_alignment(alignment_option, filename)
        processed_files.append(output_filename)
    # Get stems from processed files
    filenames = [pathlib.Path(f).stem for f in processed_files]
else:
    # Use original filenames
    if len(filenames_paths) == 0:
        print("WARNING: No FASTA files found in ./data/ directory")
        filenames = []
    else:
        print(f"Using {len(filenames_paths)} original FASTA files (no trimming)")
        filenames = [f.stem for f in filenames_paths]

print(f"\nProcessing files: {filenames}\n")

rule all:
    input:
        auspice_json  = expand("results/H5Nx/auspice/{filename}_auspice.json", filename = filenames), 
        nexus_out = expand("results/H5Nx/tree/{filename}_tree.nwk", filename = filenames),
        nextclade_out = expand("nextclade/H5Nx/{filename}/clade.tsv", filename = filenames),
        metadata_out = expand("data/{filename}_metadata.tsv", filename = filenames)

rule clade:
    input: 
        sequence = 'data/{filename}.fasta'
    output:
        "nextclade/H5Nx/{filename}/clade.tsv"
    params:
        clade = clade_calling
    shell:
        """
        nextclade3 dataset get -n {params.clade} --output-dir nextclade/H5Nx 
        nextclade3 run -j 5 -D nextclade/H5Nx \
                  {input.sequence} --quiet --output-tsv {output}
        """

rule fix_encoding_metadata:
    input: 
        meta = lambda wildcards: f"data/meta_{wildcards.filename.replace('_trimmed', '')}.tsv"
    output: 
        meta_utf8 = "data/meta_{filename}_utf8.tsv"  # Fixed: removed lambda, let Snakemake handle the wildcard
    run:
        import pandas as pd
        # Need to reconstruct the actual input filename for trimmed cases
        base_filename = wildcards.filename.replace('_trimmed', '')
        input_file = f"data/meta_{base_filename}.tsv"
        
        df = pd.read_csv(input_file, sep='\t', encoding='latin-1')  # or 'iso-8859-1', 'macroman'
        df.to_csv(output.meta_utf8, sep='\t', index=False, encoding='utf-8')
        print(f"Converted {input_file} to UTF-8")

rule combined_metadata:
    input:
        meta = "data/meta_{filename}_utf8.tsv",  # This will now match the output of fix_encoding_metadata
        nextclade = "nextclade/H5Nx/{filename}/clade.tsv"
    output:
        meta_tree = "data/{filename}_metadata.tsv",
        meta_freq = "data/{filename}_metadata_formatted.tsv"
    shell:
        """
        # Then join with nextclade data
        csvtk join -t -H --fields "strain;seqName" {input.meta} {input.nextclade} > temp_{wildcards.filename}_joined.tsv
        # Remove duplicate lines from the joined output
        csvtk uniq -t -H temp_{wildcards.filename}_joined.tsv > {output.meta_tree}
        
        # Change format of the date into YYYY-MM-DD
        csvtk replace -t -H -f 6 -p '^(\\d{{1,2}})/(\\d{{1,2}})/(\\d{{4}})$' -r '$3-$2-$1' temp_{wildcards.filename}_joined.tsv > {output.meta_freq}
        
        # Clean up temp files
        rm temp_{wildcards.filename}_joined.tsv
        """

rule augur_index:
    input: 
        sequence = 'data/{filename}.fasta'
    output:
        "results/H5Nx/tree/{filename}_index.tsv"
    shell:
        """ 
        augur index --sequences {input.sequence} --output {output} 2>&1
        """

rule augur_align:
    input:
        sequence = 'data/{filename}.fasta'
    output:
        "results/H5Nx/tree/{filename}_aligned.fasta"
    params:
        ref = reference_name
    log:
        "logs/{filename}_align.log"    
    shell:
        """
        augur align \
          --sequences {input.sequence} \
          --reference-name {params.ref} \
          --output {output} \
          --fill-gaps &> {log}
        """    
    
rule augur_tree:
    input: 
        'results/H5Nx/tree/{filename}_aligned.fasta'
    output:
        "results/H5Nx/tree/{filename}_tree_raw.nwk"
    log:
        "logs/{filename}_tree.log"
    shell:
        """
        augur tree --alignment {input} --output {output} &> {log}
        """


rule augur_refine:
    input:
        tree = "results/H5Nx/tree/{filename}_tree_raw.nwk",
        alignment = "results/H5Nx/tree/{filename}_aligned.fasta",
        metadata = "data/{filename}_metadata.tsv"
    output:
        tree = "results/H5Nx/tree/{filename}_tree.nwk",
        node = "results/H5Nx/tree/{filename}_branch_lengths.json",
        formatted_metadata = "results/H5Nx/tree/{filename}_metadata_formatted.tsv"
    params:
        clock_rate = clock_rate,
        clock_std_dev = clock_std_dev,
        root = reference_name
    log:
        "logs/{filename}_refine.log"
    shell:
        """
        augur curate format-dates \
          --metadata {input.metadata} \
          --date-fields date \
          --expected-date-formats "%d/%m/%Y" --output-metadata {output.formatted_metadata}

        augur refine \
          --tree {input.tree} \
          --alignment {input.alignment} \
          --metadata {output.formatted_metadata} \
          --output-tree {output.tree} \
          --output-node-data {output.node} \
          --timetree \
          --clock-rate {params.clock_rate} \
          --clock-std-dev {params.clock_std_dev} \
          --root {params.root} \
          --date-inference marginal \
          --keep-polytomies  &> {log}
        """

rule traits:  
    input:  
        tree =  "results/H5Nx/tree/{filename}_tree.nwk",  
        metadata = "data/{filename}_metadata.tsv" 
    output:  
        node_data = "results/H5Nx/tree/{filename}_traits.json"
    log:
        "logs/{filename}_traits.log"    
    shell:  
        """  
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output-node-data {output.node_data} \
            --columns clade \
            --confidence \
            --sampling-bias-correction 2.0 &> {log}
        """  

rule augur_ancestral:
    input:
        tree =  "results/H5Nx/tree/{filename}_tree.nwk",
        alignment = "results/H5Nx/tree/{filename}_aligned.fasta"
    log:
        "logs/{filename}_ancestral.log"
    output:
        "results/H5Nx/tree/{filename}_nt_muts.json"
    shell:
        """
        augur ancestral \
          --tree {input.tree} \
          --alignment {input.alignment} \
          --output-node-data {output} \
          --inference joint &> {log}
        """

rule create_reference_genbank:
    input:
        sequence = 'data/{filename}.fasta'
    output:
        ref_gb = "results/H5Nx/tree/{filename}_reference.gb"
    params:
        ref_name = reference_name,
        gene_name = "HA",  # Change based on your segment
    run:
        from Bio import SeqIO
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        
        # Find the reference sequence
        ref_record = None
        for record in SeqIO.parse(input.sequence, "fasta"):
            if params.ref_name in record.id or params.ref_name in record.description:
                ref_record = record
                break
        
        if ref_record is None:
            raise ValueError(f"Reference '{params.ref_name}' not found in {input.sequence}")
        
        # Set required annotations
        ref_record.annotations["molecule_type"] = "DNA"
        ref_record.annotations["organism"] = "Influenza A virus"
        ref_record.annotations["date"] = "2022"
        
        # Add mandatory SOURCE feature (spans entire sequence)
        source_feature = SeqFeature(
            FeatureLocation(0, len(ref_record.seq)),
            type="source",
            qualifiers={
                "organism": ["Influenza A virus"],
                "mol_type": ["viral cRNA"]
            }
        )
        
        # Add CDS feature for the gene
        cds_feature = SeqFeature(
            FeatureLocation(0, len(ref_record.seq)),
            type="CDS",
            qualifiers={
                "gene": [params.gene_name],
                "product": [params.gene_name]
            }
        )
        
        # Add features to record (source MUST come first)
        ref_record.features = [source_feature, cds_feature]
        
        # Write GenBank file
        SeqIO.write(ref_record, output.ref_gb, "genbank")
        print(f"Created GenBank: {ref_record.id}, gene: {params.gene_name}, length: {len(ref_record.seq)}")

rule augur_translate:
    input:
        tree = "results/H5Nx/tree/{filename}_tree.nwk",
        ancestral_sequence = "results/H5Nx/tree/{filename}_nt_muts.json",
        reference = "results/H5Nx/tree/{filename}_reference.gb" 
    log:
        "logs/{filename}_translate.log"    
    output:
        "results/H5Nx/tree/{filename}_aa_muts.json"
    shell:
        """
        augur translate \
          --tree {input.tree} \
          --ancestral-sequences {input.ancestral_sequence} \
          --reference-sequence {input.reference} \
          --output-node-data {output} &> {log}
        """

rule frequency:
    input:
        metadata = "results/H5Nx/tree/{filename}_metadata_formatted.tsv",
        tree =  "results/H5Nx/tree/{filename}_tree.nwk",
        traits = "results/H5Nx/tree/{filename}_traits.json",
        script = "script/estimate_freq_annotated.py"
    log:
        "logs/{filename}_frequencies.log"   
    output:
        freq_augur = "results/H5Nx/auspice/{filename}_frequencies.json",  
        freq_sidecar = "results/H5Nx/auspice/{filename}_auspice_tip-frequencies.json"
    params:
        "results/H5Nx/auspice/{filename}_tree.json"
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata}\
            --node-data {input.traits} \
            --output {params}

        python {input.script} \
           --tree-json {params} \
           --metadata {input.metadata} \
           --narrow-bandwidth 0.08333333333333333 \
           --proportion-wide 0.2 \
           --min-date 2025-01-01 \
           --pivot-interval 6 \
           --clade-column clade \
           --output {output.freq_augur} \
           --output-tip-frequencies {output.freq_sidecar} &> {log}
        """


rule export:
    message:
        "Exporting data files for auspice"
    input:
        tree="results/H5Nx/tree/{filename}_tree.nwk",
        metadata="data/{filename}_metadata.tsv",
        branch_lengths = "results/H5Nx/tree/{filename}_branch_lengths.json",
        aa_mut = "results/H5Nx/tree/{filename}_aa_muts.json",
        nt_mut = "results/H5Nx/tree/{filename}_nt_muts.json",
        traits = "results/H5Nx/tree/{filename}_traits.json",
        freq_augur = "results/H5Nx/auspice/{filename}_frequencies.json",
        freq_auspice = "results/H5Nx/auspice/{filename}_auspice_tip-frequencies.json",
        auspice_config= "config/auspice_config.json"
    output:
        auspice_json="results/H5Nx/auspice/{filename}_auspice.json",
    params:
        header = header_name,
        filename = "{filename}"
    log:
        "logs/{filename}_export.log"    
    shell:
        """
        # Select only the columns we need, wrapped in quotes to handle spaces
        tsv-select -H -f 'strain,Host,Location,Broad_Location,PB2_D701_mutation,NA_deletion,Submitting_Lab,clade' \
            {input.metadata} > {input.metadata}_{params.filename}.tmp
        
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata}_{params.filename}.tmp \
            --node-data {input.branch_lengths} {input.aa_mut} {input.nt_mut} {input.traits} {input.freq_augur} \
            --auspice-config {input.auspice_config} \
            --color-by-metadata "Host" "Location" "Broad_Location" "PB2_D701_mutation" "NA_deletion" "Submitting_Lab" "clade"  \
            --minify-json \
            --title "Influenza: {params.header}" \
            --include-root-sequence \
            --output {output.auspice_json} &> {log}
        
        mkdir -p auspice
        cp {output.auspice_json} auspice/{params.filename}_auspice.json
        cp {input.freq_auspice} auspice/{params.filename}_auspice_tip-frequencies.json
        
        rm {input.metadata}_{params.filename}.tmp
        """