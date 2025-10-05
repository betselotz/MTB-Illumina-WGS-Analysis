# 1️⃣ Shovill

<details>
<summary>🏗️ Shovill: Bacterial Genome Assembler</summary>

**Shovill** is a fast and easy-to-use **bacterial genome assembler** designed for Illumina short-read data. It wraps around popular assemblers like **SPAdes** or **SKESA**, streamlining genome assembly from paired-end reads.

### Key points for TB genomes
- Optimized for **small bacterial genomes** (~4–5 Mb).  
- Supports multithreading (`--cpus`) and efficient RAM usage (`--ram`) for faster assemblies.  
- Allows customization of **minimum contig length** (`--minlen`) and **coverage depth** (`--mincov` / `--depth`).  
- Automatically renames output contigs for clarity and downstream analyses.  
- Works well with **preprocessed FASTQ files** from `fastp`.  

> ⚠ **Note:** Best used with high-quality, paired-end Illumina reads. Low-quality or fragmented data may require additional QC before assembly.

</details>

##### Step 1: Create or edit the script
```bash
nano run_shovill.sh
```
#####  Step 2: Paste the following into `run_shovill.sh`
```bash
#!/bin/bash

INPUT_DIR="fastp_results_min_50"
OUTDIR="shovill_results"
GSIZE=4411532
MIN_COVERAGE=30
MIN_READS=10000
SAMPLE_SEQ_LINES=1000
REPORT="$OUTDIR/assembly_summary.tsv"

mkdir -p "$OUTDIR"
echo -e "Sample\tTotal_Reads\tAvg_Read_Length\tEstimated_Coverage\tStatus\tReason" > "$REPORT"

shopt -s nullglob
for R1 in "$INPUT_DIR"/*_1.trim.fastq.gz; do
  [[ -e "$R1" ]] || continue
  R2="${R1/_1.trim.fastq.gz/_2.trim.fastq.gz}"
  sample=$(basename "$R1" _1.trim.fastq.gz)
  sample_out="$OUTDIR/$sample"
  mkdir -p "$sample_out"

  # check if paired file exists
  if [[ ! -f "$R2" ]]; then
    echo -e "${sample}\tNA\tNA\tNA\tDropped\tMissing R2" >> "$REPORT"
    continue
  fi

  # skip if already assembled
  if [[ -f "$sample_out/${sample}_contigs.fa" ]]; then
    echo -e "${sample}\tNA\tNA\tNA\tSkipped\tAlready assembled" >> "$REPORT"
    continue
  fi

  # calculate total reads
  total_lines=$(zcat "$R1" | wc -l)
  total_reads=$(( total_lines / 4 ))
  if [[ $total_reads -le 0 ]]; then
    echo -e "${sample}\t0\t0\t0\tDropped\tNo reads" >> "$REPORT"
    continue
  fi

  # calculate average read length from first SAMPLE_SEQ_LINES sequences
  avg_len=$(zcat "$R1" | sed -n '2~4p' | head -n "$SAMPLE_SEQ_LINES" | awk '{total+=length($0); count++} END {if(count==0) print 0; else printf "%.0f", total/count}')
  if [[ "$avg_len" -le 0 ]]; then
    echo -e "${sample}\t${total_reads}\t0\t0\tDropped\tInvalid read length" >> "$REPORT"
    continue
  fi

  # estimated coverage
  estimated_coverage=$(awk -v r="$total_reads" -v l="$avg_len" -v g="$GSIZE" 'BEGIN{cov=(r*l*2)/g; printf "%.2f", cov}')
  if [[ $total_reads -lt $MIN_READS ]]; then
    echo -e "${sample}\t${total_reads}\t${avg_len}\t${estimated_coverage}\tDropped\tLow read count" >> "$REPORT"
    continue
  elif (( $(echo "$estimated_coverage < $MIN_COVERAGE" | bc -l) )); then
    echo -e "${sample}\t${total_reads}\t${avg_len}\t${estimated_coverage}\tDropped\tLow coverage (<${MIN_COVERAGE}x)" >> "$REPORT"
    continue
  fi

  # run SPAdes via Shovill
  echo "==> Running Shovill (paired-end, SPAdes) on $sample (reads=${total_reads}, avg_len=${avg_len}, est_cov=${estimated_coverage}x)"
  shovill \
    --R1 "$R1" \
    --R2 "$R2" \
    --gsize "$GSIZE" \
    --outdir "$sample_out" \
    --assembler spades \
    --minlen 500 \
    --mincov 30 \
    --depth 100 \
    --namefmt "${sample}_%05d" \
    --cpus 16 \
    --ram 120 \
    --tmpdir "${TMPDIR:-/tmp}" \
    --force

  # check if assembly succeeded
  if [[ ! -s "$sample_out/contigs.fa" && ! -s "$sample_out/${sample}_contigs.fa" ]]; then
    echo -e "${sample}\t${total_reads}\t${avg_len}\t${estimated_coverage}\tFailed\tAssembly error" >> "$REPORT"
    continue
  fi

  echo -e "${sample}\t${total_reads}\t${avg_len}\t${estimated_coverage}\tAssembled\tOK" >> "$REPORT"

  # rename output files consistently
  for f in "$sample_out"/*; do
    base=$(basename "$f")
    mv "$f" "$sample_out/${sample}_$base"
  done
done

echo "All samples processed. See $REPORT"

```





```bash
#!/bin/bash

INPUT_DIR="fastp_results_min_50"
OUTDIR="shovill_results"
mkdir -p "$OUTDIR"

GSIZE=4411532
MIN_COVERAGE=30
MIN_READS=10000
SAMPLE_SEQ_LINES=1000
REPORT="$OUTDIR/assembly_summary.tsv"

echo -e "Sample\tTotal_Reads\tAvg_Read_Length\tEstimated_Coverage\tStatus\tReason" > "$REPORT"

shopt -s nullglob
for R1 in "$INPUT_DIR"/*.fastq.gz; do
  [[ -e "$R1" ]] || continue

  # Remove .fastq.gz and optional .trim suffix
  sample=$(basename "$R1" .fastq.gz)
  sample=${sample%.trim}
  sample_out="$OUTDIR/$sample"
  mkdir -p "$sample_out"

  # skip if already assembled
  if [[ -f "$sample_out/contigs.fa" ]]; then
    echo ">> Skipping $sample (already assembled)"
    echo -e "${sample}\tNA\tNA\tNA\tSkipped\tAlready assembled" >> "$REPORT"
    continue
  fi

  # count total reads
  total_lines=$(zcat "$R1" | wc -l)
  total_reads=$(( total_lines / 4 ))
  if [[ $total_reads -le 0 ]]; then
    echo -e "${sample}\t0\t0\t0\tDropped\tNo reads" >> "$REPORT"
    continue
  fi

  # calculate average read length from first SAMPLE_SEQ_LINES sequences
  avg_len=$(zcat "$R1" | sed -n '2~4p' | head -n "$SAMPLE_SEQ_LINES" | awk '{total+=length($0); count++} END {if(count==0) print 0; else printf "%.0f", total/count}')
  if [[ "$avg_len" -le 0 ]]; then
    echo -e "${sample}\t${total_reads}\t0\t0\tDropped\tInvalid read length" >> "$REPORT"
    continue
  fi

  # estimate coverage
  estimated_coverage=$(awk -v r="$total_reads" -v l="$avg_len" -v g="$GSIZE" 'BEGIN{cov=(r*l)/g; printf "%.2f", cov}')
  if [[ $total_reads -lt $MIN_READS ]]; then
    echo -e "${sample}\t${total_reads}\t${avg_len}\t${estimated_coverage}\tDropped\tLow read count" >> "$REPORT"
    continue
  elif (( $(echo "$estimated_coverage < $MIN_COVERAGE" | bc -l) )); then
    echo -e "${sample}\t${total_reads}\t${avg_len}\t${estimated_coverage}\tDropped\tLow coverage (<${MIN_COVERAGE}x)" >> "$REPORT"
    continue
  fi

  echo "==> Running Shovill (single-read workaround, SPAdes) on $sample (reads=${total_reads}, avg_len=${avg_len}, est_cov=${estimated_coverage}x)"
  shovill \
    --R1 "$R1" \
    --R2 "$R1" \
    --gsize "$GSIZE" \
    --outdir "$sample_out" \
    --assembler spades \
    --minlen 500 \
    --mincov 30 \
    --depth 100 \
    --namefmt "${sample}_%05d" \
    --cpus 16 \
    --ram 120 \
    --tmpdir "${TMPDIR:-/tmp}" \
    --force

  # check if assembly succeeded
  if [[ ! -s "$sample_out/contigs.fa" && ! -s "$sample_out/${sample}_contigs.fa" ]]; then
    echo -e "${sample}\t${total_reads}\t${avg_len}\t${estimated_coverage}\tFailed\tAssembly error" >> "$REPORT"
    continue
  fi

  echo -e "${sample}\t${total_reads}\t${avg_len}\t${estimated_coverage}\tAssembled\tOK" >> "$REPORT"

  # rename output files consistently
  for f in "$sample_out"/*; do
    base=$(basename "$f")
    mv "$f" "$sample_out/${sample}_$base"
  done
done

echo "All samples processed. See $REPORT"

```

##### Step 3: Save and exit nano
Press Ctrl + O → Enter (to write the file)
Press Ctrl + X → Exit nano
###### Step 4: Make the script executable
``` bash
chmod +x run_shovill.sh
```
###### Step 5: Activate environment and run
``` bash
conda activate shovill_env
./run_shovill.sh
```
### Checking the output
View large files page by page
```bash
less shovill_results/ET1135_S12/ET1135_S12_contigs.fa
```
View the firt ten lines
```bash
head -n 20 shovill_results/ET1135_S12/ET1135_S12_contigs.fa
```
Counted how many contigs we have
```bash
grep -c ">" ./shovill_results/ET1135_S12/ET1135_S12_contigs.fa
```

To get contig name, length, coverage from the FASTA headers:
```bash
grep ">" ./shovill_results/ET1135_S12/ET1135_S12_contigs.fa | \
awk -F'[ =]' '{print $1, $2, $3, $4, $5}'
```
#  Assembly Evaluation

Before downstream analyses, it is important to verify the quality of the assembled genome. This ensures reliable results in variant calling, consensus generation, and phylogenetic analyses.
<details>
<summary>🧪 Why Do We Evaluate Genome Assemblies in <i>Mycobacterium tuberculosis</i>? (Click to Expand)</summary>

Evaluating genome assemblies is a **critical step** before using them for downstream analyses.  
For *M. tuberculosis* (MTB), this step is especially important due to its **clinical, evolutionary, and genomic characteristics**.

---

### 🔹 1. Ensuring Genome Completeness
- MTB genomes are ~4.4 Mb in size.  
- A high-quality assembly should **cover nearly the full genome** with minimal gaps.  
- Missing regions can result in:  
  - Incorrect phylogenetic placement.  
  - Loss of resistance-associated genes.  
  - Misleading annotation results.  

✅ **Goal:** Verify that assemblies are close to expected genome size and contain all core genes.

---

### 🔹 2. Checking for Contamination
- Clinical MTB isolates can be contaminated with:  
  - Host DNA (human reads).  
  - Other bacteria from sputum samples.  
  - Laboratory contaminants.  
- Contaminants can skew:  
  - Variant calling results.  
  - Lineage classification.  
  - Drug-resistance predictions.  

✅ **Goal:** Confirm GC content (~65%) and check species purity.
---

### 🔹 3. Evaluating Assembly Accuracy
- MTB genomes are **highly conserved** but contain repetitive regions (e.g., PE/PPE gene families).  
- Assemblers can misplace or split these repeats, causing **chimeric contigs** or false duplications.  
- Misassemblies lead to:  
  - Incorrect SNP calls.  
  - Faulty gene presence/absence analyses.  

✅ **Goal:** Use metrics like **N50, contig counts, and read mapping** to verify assembly integrity.

### 🔹 4. Enabling Reliable Phylogenetics
- MTB phylogenetic studies depend on **precise SNP alignments**.  
- Low-quality assemblies introduce:  
  - Missing loci.  
  - False SNPs from sequencing/assembly errors.  
- This can distort lineage/sub-lineage assignments and outbreak tracking.  

✅ **Goal:** Confirm assemblies are suitable for inclusion in phylogenomic pipelines.

---

### 🔹 5. Interpretation Guidelines (MTB Assemblies)

When evaluating *Mycobacterium tuberculosis* assemblies, keep these quality thresholds in mind:

- **Genome size:** ~4.4 Mb (±0.2 Mb).  
- **Number of contigs:** Ideally <200 (draft assemblies). Fewer contigs indicate better assembly continuity.  
- **N50:** >50,000 bp (higher = better).  
  - *(N50 = length of the shortest contig such that 50% of the total assembly length is contained in contigs of this size or longer)*  
- **L50:** Typically <50 for good assemblies.  
  - *(L50 = the minimum number of largest contigs needed to cover 50% of the assembly)*  
- **GC content:** ~65% (should be consistent across MTB isolates; deviations may indicate contamination).  
- **Max contig length:** >200 kb is desirable; indicates assembler successfully reconstructed long genomic regions.  
- **Mean/median contig length:** Larger average lengths suggest fewer fragmented contigs.  
- **% of genome in small contigs (<500 bp):** Should be minimal; many short contigs usually indicate poor assembly or contamination.  
- **Coverage depth (if reads are mapped back):** 30×–100× is typical for MTB. Low coverage (<20×) risks missing variants; extremely high coverage (>200×) may cause assembler artifacts.  

---

</details>

### 1.run QUAST on all Shovill assemblies
Collect the key statistics in a single CSV file
running QUAST on all SPAdes assemblies and collecting key statistics into a single report.tsv. 

##### Step 1: Create or edit the script
``` bash
nano run_quast_shovill.sh
``` 
#####  Step 2: Paste the following into `run_quast_shovill.sh`
``` bash
#!/bin/bash
set -euo pipefail

SHOVILL_DIR="shovill_results"
QUAST_PARENT="quast_results"
QUAST_DIR="$QUAST_PARENT/quast_results_shovill"
CSV_OUTDIR="csv_output"

mkdir -p "$QUAST_DIR" "$CSV_OUTDIR"

CSV_FILE="$CSV_OUTDIR/quast_summary_shovill.csv"
echo "Sample,NumContigs,TotalLength,MinLen,MaxLen,AverageLen,N50,N75,GC%" > "$CSV_FILE"

for sample_out in "$SHOVILL_DIR"/*; do
  [[ -d "$sample_out" ]] || continue
  sample=$(basename "$sample_out")
  
  contigs_file=("$sample_out"/*_contigs.fa)
  [[ -f "${contigs_file[0]}" ]] || continue
  contigs="${contigs_file[0]}"

  outdir="$QUAST_DIR/$sample"
  mkdir -p "$outdir"

  quast "$contigs" -o "$outdir" > /dev/null 2>&1

  stats_file="$outdir/report.tsv"
  if [[ -f "$stats_file" ]]; then
    num_contigs=$(awk -F'\t' '$1=="# contigs (>= 0 bp)"{print $2}' "$stats_file")
    total_len=$(awk -F'\t' '$1=="Total length (>= 0 bp)"{print $2}' "$stats_file")
    min_len=$(awk -F'\t' '$1=="Shortest contig"{print $2}' "$stats_file")
    max_len=$(awk -F'\t' '$1=="Largest contig"{print $2}' "$stats_file")
    avg_len=$(awk -F'\t' '$1=="Average contig length"{print $2}' "$stats_file")
    n50=$(awk -F'\t' '$1=="N50"{print $2}' "$stats_file")
    n75=$(awk -F'\t' '$1=="N75"{print $2}' "$stats_file")
    gc=$(awk -F'\t' '$1=="GC (%)"{print $2}' "$stats_file")

    echo "$sample,$num_contigs,$total_len,$min_len,$max_len,$avg_len,$n50,$n75,$gc" >> "$CSV_FILE"
  fi
done
```

##### Step 3: Save and exit nano
Press Ctrl + O → Enter (to write the file)
Press Ctrl + X → Exit nano

###### Step 4: Make the script executable
``` bash
chmod +x run_quast_shovill.sh
```
###### Step 5: Activate environment and run
``` bash
conda activate quast_env
./run_quast_shovill.sh
```

##### Step 1: Create or edit the script
``` bash
nano run_assembly_scan_shovill.sh
``` 
#####  Step 2: Paste the following into `run_assembly_scan_shovill.sh`

``` bash
#!/bin/bash
set -euo pipefail

SHOVILL_DIR="shovill_results"
CSV_OUTDIR="csv_output"
mkdir -p "$CSV_OUTDIR"
CSV_FILE="$CSV_OUTDIR/spades_assembly_scan.csv"

echo "Sample,total_contig,total_contig_length,max_contig_length,mean_contig_length,median_contig_length,min_contig_length,n50_contig_length,l50_contig_count,contig_percent_a,contig_percent_c,contig_percent_g,contig_percent_t,contigs_greater_1k,contigs_greater_10k,contigs_greater_100k" > "$CSV_FILE"

for sample_dir in "$SHOVILL_DIR"/*/; do
    sample=$(basename "$sample_dir")
    contig_file="${sample_dir}${sample}_contigs.fa"

    if [[ -f "$contig_file" ]]; then
        tmp=$(mktemp)
        assembly-scan "$contig_file" --transpose | tail -n +2 > "$tmp"

        declare -A data
        while IFS=$'\t' read -r contig metric value; do
            data["$metric"]="$value"
        done < "$tmp"

        rm "$tmp"

        printf "%s" "$sample" >> "$CSV_FILE"
        for col in total_contig total_contig_length max_contig_length mean_contig_length median_contig_length min_contig_length n50_contig_length l50_contig_count contig_percent_a contig_percent_c contig_percent_g contig_percent_t contigs_greater_1k contigs_greater_10k contigs_greater_100k; do
            printf ",%s" "${data[$col]:-0}" >> "$CSV_FILE"
        done
        echo "" >> "$CSV_FILE"
    fi
done
``` 


##### Step 3: Save and exit nano
Press Ctrl + O → Enter (to write the file)
Press Ctrl + X → Exit nano

###### Step 4: Make the script executable
``` bash
chmod +x run_assembly_scan_shovill.sh
```
###### Step 5: Activate environment and run
``` bash
conda activate assembly_scan_env
./run_assembly_scan_shovill.sh
```

