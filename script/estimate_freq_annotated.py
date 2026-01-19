#!/usr/bin/env python3
"""
Frequency estimation script that outputs frequency data in dict format.
Outputs:
1. Augur-compatible dict format (for augur export --node-data)
2. Optional tip-frequencies sidecar (same dict format, for Auspice visualization)
"""
import argparse
import numpy as np
import pandas as pd
from augur.dates import get_numerical_dates, numeric_date_type
from augur.frequencies import format_frequencies
from augur.frequency_estimators import get_pivots, KdeFrequencies 
import json
from augur.io import read_metadata
from augur.utils import write_json

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Estimate sequence frequencies from metadata using the KDE method, stratified by clade.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # REQUIRED INPUTS
    parser.add_argument("--metadata", required=True, 
                       help="TSV file of metadata with at least 'strain', 'date', and specified clade columns")
    parser.add_argument("--tree-json", required=True, 
                       help="Annotated tree JSON file (Used only for tip names/clade aggregation)")
    
    # KDE-SPECIFIC PARAMETERS
    parser.add_argument("--narrow-bandwidth", required=True, type=float, 
                       help="Narrow bandwidth (sigma) for KDE calculation")
    parser.add_argument("--proportion-wide", type=float, default=0.0, 
                       help="Proportion of wide bandwidth to use for smoothing")

    # GENERAL PARAMETERS
    parser.add_argument("--pivot-interval", type=int, default=2, 
                       help="Interval between pivots in weeks")
    parser.add_argument("--min-date", type=numeric_date_type, 
                       help="Minimum date to estimate frequencies for")
    parser.add_argument("--max-date", type=numeric_date_type, 
                       help="Maximum date to estimate frequencies for")
    parser.add_argument("--clade-column", nargs='+', default=['clade'], 
                       help="Metadata columns to group strains by (e.g., 'clade subclade')")
    
    # OUTPUT OPTIONS
    parser.add_argument("--output", required=True, 
                       help="JSON file for augur export (frequencies included in node-data)")
    parser.add_argument("--output-tip-frequencies", 
                       help="JSON file for Auspice tip-frequencies sidecar. If not specified, will not create sidecar file.")
    parser.add_argument("--include-tips", action='store_true',
                       help="Include individual tip frequencies in main output (increases file size)")
    
    args = parser.parse_args()

    print("="*70)
    print("Frequency Estimation - KDE Method")
    print("="*70)
    print(f"Metadata: {args.metadata}")
    print(f"Tree JSON: {args.tree_json}")
    print(f"KDE bandwidth: {args.narrow_bandwidth}")
    print(f"Pivot interval: {args.pivot_interval} weeks")
    print(f"Grouping by: {', '.join(args.clade_column)}")
    print(f"Output (for augur export): {args.output}")
    if args.output_tip_frequencies:
        print(f"Output (tip-frequencies sidecar): {args.output_tip_frequencies}")
    print("="*70)

    # --- 1. Load Data ---
    print("\nStep 1: Loading data...")
    with open(args.tree_json, 'r') as f:
        tree_data = json.load(f)
    print(f"  ✓ Loaded tree JSON")
    
    columns_to_load = ["strain", "date"] + args.clade_column
    metadata = read_metadata(
        args.metadata,
        columns=columns_to_load,
        dtype="string",
    )
    print(f"  ✓ Loaded metadata: {len(metadata)} strains")
    
    dates = get_numerical_dates(metadata, fmt='%Y-%m-%d')
    print(f"  ✓ Parsed dates")
    
    # --- 2. Collect Valid Strains and Observations Globally ---
    print("\nStep 2: Filtering valid strains...")
    valid_strains_df = metadata.copy()
    valid_strains_df['num_date'] = pd.Series({strain: np.mean(dates[strain]) for strain in dates})
    
    # Filter on date AND all specified clade columns
    valid_strains_df.dropna(subset=['num_date'] + args.clade_column, inplace=True)
    
    all_strains = valid_strains_df.index.tolist()
    all_observations = valid_strains_df['num_date'].values
    
    print(f"  ✓ Valid strains: {len(all_strains)}")
    print(f"  ✓ Date range: {all_observations.min():.2f} to {all_observations.max():.2f}")
    
    # Print group statistics
    for col in args.clade_column:
        n_groups = valid_strains_df[col].nunique()
        print(f"  ✓ {col}: {n_groups} unique values")
    
    # --- 3. Calculate Global Pivots and Tip Frequency Matrix (KDE) ---
    print("\nStep 3: Computing frequencies with KDE...")
    pivots = get_pivots(
        all_observations,
        args.pivot_interval,
        args.min_date,
        args.max_date,
        "months",
    )
    print(f"  ✓ Generated {len(pivots)} pivot points")
    
    # Initialize the KDE Estimator 
    frequencies_estimator = KdeFrequencies(
        sigma_narrow=args.narrow_bandwidth,
        proportion_wide=args.proportion_wide,
        pivot_frequency=args.pivot_interval,
        start_date=args.min_date,
        end_date=args.max_date,
    )
    
    # This matrix contains the KDE time-series for every sequence (tip)
    print(f"  Computing frequency matrix for {len(all_strains)} strains...")
    global_frequency_matrix = frequencies_estimator.estimate_frequencies(
        all_observations,
        pivots,
    )
    print(f"  ✓ Computed matrix: {global_frequency_matrix.shape}")

    # Normalize
    print("  Normalizing frequencies...")
    pivot_sums = np.sum(global_frequency_matrix, axis=0)
    normalized_matrix = np.divide(
        global_frequency_matrix,
        pivot_sums,
        where=pivot_sums!=0,
        out=np.zeros_like(global_frequency_matrix)
    )
    global_frequency_matrix = normalized_matrix
    print("  ✓ Normalization complete")

    # --- 4. Output Format: DICT FORMAT (dict with pivots) ---
    print(f"\nStep 4: Generating frequency outputs...")
    
    frequency_dict = {"pivots": list(pivots)}
    
    # Optionally include tip frequencies
    if args.include_tips:
        tip_count = 0
        for index, strain in enumerate(all_strains):
            tip_freq_array = global_frequency_matrix[index]
            if tip_freq_array.sum() > 0:
                frequency_dict[strain] = {
                    "frequencies": format_frequencies(tip_freq_array)
                }
                tip_count += 1
        print(f"  ✓ Added {tip_count} tip frequencies")
    else:
        print(f"  ℹ Skipping tip frequencies (use --include-tips to add)")
    
    # Aggregate Group Frequencies
    for group_column in args.clade_column:
        unique_group_values = valid_strains_df[group_column].unique()
        group_count = 0
        
        for group_value in unique_group_values:
            indices_in_matrix = [
                all_strains.index(s) 
                for s in valid_strains_df.index[valid_strains_df[group_column] == group_value].tolist()
            ]
            
            if indices_in_matrix:
                aggregate_freq = np.sum(global_frequency_matrix[indices_in_matrix], axis=0)
                aggregate_freq = np.clip(aggregate_freq, 0.0, 1.0)
                
                node_name = f"{group_column}|{group_value}"
                
                if aggregate_freq.sum() > 0:
                    frequency_dict[node_name] = {
                        "frequencies": format_frequencies(aggregate_freq)
                    }
                    group_count += 1
        
        print(f"  ✓ Added {group_count} '{group_column}' group frequencies")
    
    # Add Global Frequency
    global_freq_array = np.array([1.0] * len(pivots))
    frequency_dict["global"] = {
        "frequencies": format_frequencies(global_freq_array)
    }
    print(f"  ✓ Added global frequency")
    
    # Write Augur format (for node-data in augur export)
    write_json(frequency_dict, args.output)
    print(f"\n✓ Wrote frequency data to: {args.output}")
    print(f"  Total keys: {len(frequency_dict)}")
    print(f"  This file should be included in 'augur export v2 --node-data'")

    # --- 5. Optionally create tip-frequencies sidecar ---
    if args.output_tip_frequencies:
        print(f"\nStep 5: Creating tip-frequencies sidecar...")
        
        # The sidecar file uses the SAME dict format as the augur output
        # It should contain tip frequencies for the dropdown visualization
        tip_freq_dict = {"pivots": list(pivots)}
        
        # Include ALL tips for the sidecar (this is what Auspice uses for the frequency panel)
        tip_count = 0
        for index, strain in enumerate(all_strains):
            tip_freq_array = global_frequency_matrix[index]
            if tip_freq_array.sum() > 0:
                tip_freq_dict[strain] = {
                    "frequencies": format_frequencies(tip_freq_array)
                }
                tip_count += 1
        
        # Write tip-frequencies sidecar
        write_json(tip_freq_dict, args.output_tip_frequencies)
        print(f"✓ Wrote tip-frequencies sidecar to: {args.output_tip_frequencies}")
        print(f"  Total tips: {tip_count}")
        print(f"  This file will be automatically loaded by Auspice if named correctly")

    print("\n" + "="*70)
    print("✓ SUCCESS!")
    print("="*70)