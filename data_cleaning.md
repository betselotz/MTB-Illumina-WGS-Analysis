modify 5.fastp_filtered_reads_plot.tsv in place and remove _1 from the Sample column directly.
```bash
sed -i -E 's/^([^[:space:]]+)_1/\1/' 5.fastp_filtered_reads_plot.tsv
```

```bash
sed -i -E 's/^([^[:space:]]+)_1/\1/' 6.general_stats_table.tsv
```

```bash
sed '/^-/d' checkm.txt | awk '{$1=$1}1' OFS=',' > checkm.csv
```

for aggregate  result 

merge_all_summary_files.sh
```bash

#!/bin/bash
# Robust merge of 1–14 CSV/TSV files using pandas

python3 - <<'EOF'
import pandas as pd

files = [
    "1_fastq_read_counts.csv",
    "2_read_length_summary.csv",
    "3_trimmed_read_counts.csv",
    "4_trimmed_read_length_summary.csv",
    "5_fastp_general_stats_table.tsv",
    "6_fastp_filtered_reads_plot.tsv",
    "7_tbprofiler_collated.csv",
    "8_qualimap_general_stats_table.tsv",
    "9_variant_filter_summary.csv",
    "10_consensus_lengths.csv",
    "11_shovill_assembly_summary.tsv",
    "12_quast_summary_shovill.csv",
    "13_backmap_summary.tsv",
    "14_prokka_annotation_summary.tsv"
]

def read_file(f):
    # First try comma
    try:
        df = pd.read_csv(f)
        if 'Sample' in df.columns:
            return df
    except:
        pass
    # Fallback to tab
    df = pd.read_csv(f, sep='\t')
    # Strip whitespace from column names
    df.columns = df.columns.str.strip()
    return df

# Load first file
merged = read_file(files[0])

# Merge sequentially on 'Sample'
for f in files[1:]:
    df = read_file(f)
    if 'Sample' not in df.columns:
        print(f"❌ Warning: {f} does not have 'Sample' column, skipping")
        continue
    merged = pd.merge(merged, df, on='Sample', how='outer')
    print(f"Merged {f}")

# Save merged file
merged.to_csv("merged_all_summary_files.csv", index=False)
print("\n✅ All files merged successfully into merged_all_summary_files.csv")
EOF

```
