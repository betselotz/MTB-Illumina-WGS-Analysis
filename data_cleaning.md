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
set -euo pipefail

# Directory containing your files
INPUT_DIR="aggregate_result"
OUTPUT_FILE="$INPUT_DIR/merged_summary.csv"

# List of files to merge in order
FILES=(
  "1_fastq_read_counts.csv"
  "2_read_length_summary.csv"
  "3_trimmed_read_counts.csv"
  "4_trimmed_read_length_summary.csv"
  "5_fastp_general_stats_table.tsv"
  "6_fastp_filtered_reads_plot.tsv"
  "7_tbprofiler_collated.csv"
  "8_qualimap_general_stats_table.tsv"
  "9_variant_filter_summary.csv"
  "10_consensus_lengths.csv"
  "11_shovill_assembly_summary.tsv"
  "12_quast_summary_shovill.csv"
  "13_backmap_summary.tsv"
  "14_prokka_annotation_summary.tsv"
)

# Temporary working directory
TMP_DIR=$(mktemp -d)
cd "$INPUT_DIR"

# Convert all TSVs to CSV for consistency
for f in "${FILES[@]}"; do
  if [[ "$f" == *.tsv ]]; then
    tr '\t' ',' < "$f" > "$TMP_DIR/$f.csv"
  else
    cp "$f" "$TMP_DIR/$f.csv"
  fi
done

# Merge step by step on "Sample" column
cd "$TMP_DIR"
FIRST_FILE="${FILES[0]}.csv"
cp "$FIRST_FILE" merged.csv

for f in "${FILES[@]:1}"; do
  csvjoin --no-inference -c "Sample" merged.csv "$f.csv" > tmp.csv
  mv tmp.csv merged.csv
done

# Clean column duplicates (optional)
csvcut -n merged.csv | awk -F: '{print $2}' | awk '{$1=$1};1' | sort | uniq -d | while read -r col; do
  echo "Removing duplicate column: $col"
  csvcut -C "$col" merged.csv > tmp.csv && mv tmp.csv merged.csv
done

# Move merged file to output location
mv merged.csv "$OUTPUT_FILE"

echo "✅ Merged summary file created at: $OUTPUT_FILE"
```
