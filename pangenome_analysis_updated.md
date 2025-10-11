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
set -euo pipefail
IFS=$'\n\t'

INPUT_DIR="fastp_results_min_50"
OUTDIR="shovill_results"
GSIZE=4411532
MIN_COVERAGE=30
MIN_READS=10000
SAMPLE_SEQ_LINES=1000
CPUS=${CPUS:-8}          # per-job CPU usage
RAM=8                    # per-job RAM (GB)
PARALLEL_JOBS=8          # 8 jobs × 8 GB = 64 GB total
KRAKEN_DB="/path/to/kraken_db"
REPORT="$OUTDIR/assembly_summary.tsv"
REPORT_LOCK="$OUTDIR/report.lock"

# Safety check: if kraken2 isn't found, disable contamination check
if ! command -v kraken2 &>/dev/null; then
  echo "[WARN] Kraken2 not found in environment — skipping contamination checks."
  KRAKEN_DB=""
fi

mkdir -p "$OUTDIR"
printf "Sample\tTotal_Reads\tAvg_Read_Length\tEstimated_Coverage\tAssembly_Size(bp)\tNum_Contigs\tN50\tRuntime_min\tStatus\tReason\n" > "$REPORT"

run_sample() {
  R1="$1"
  R2="${R1/_1.trim.fastq.gz/_2.trim.fastq.gz}"
  sample=$(basename "$R1" _1.trim.fastq.gz)
  sample_out="$OUTDIR/$sample"
  rm -rf "$sample_out"
  mkdir -p "$sample_out"

  total_lines=$(zcat "$R1" | wc -l)
  total_reads=$(( total_lines / 4 ))
  avg_len=$(zcat "$R1" | sed -n '2~4p' | head -n "$SAMPLE_SEQ_LINES" | awk '{s+=length($0)} END{if(NR)printf "%.0f",s/NR;else print 0}')
  estimated_coverage=$(awk -v r="$total_reads" -v l="$avg_len" -v g="$GSIZE" 'BEGIN{printf "%.2f", (r*l*2)/g}')

  status="Ready"; reason="OK"
  if [[ ! -f "$R2" ]]; then status="Dropped"; reason="Missing R2"; fi
  if [[ $total_reads -lt 1 ]]; then status="Dropped"; reason="No reads"; fi
  if [[ $total_reads -lt $MIN_READS ]]; then status="Skipped"; reason="Low read count (<$MIN_READS)"; fi
  if (( $(echo "$estimated_coverage < $MIN_COVERAGE" | bc -l) )); then status="Skipped"; reason="Low coverage (<${MIN_COVERAGE}x)"; fi

  if [[ "$status" != "Ready" ]]; then
    { flock -x 200; printf "%s\t%s\t%s\t%s\t0\t0\t0\t0\t%s\t%s\n" "$sample" "$total_reads" "$avg_len" "$estimated_coverage" "$status" "$reason" >> "$REPORT"; } 200>"$REPORT_LOCK"
    rm -rf "$sample_out"
    echo "[WARN] [$sample] $reason"
    return
  fi

  echo "[INFO] [$sample] Running Shovill..."
  start_time=$(date +%s)
  shovill --R1 "$R1" --R2 "$R2" --gsize "$GSIZE" --outdir "$sample_out" \
          --assembler spades --minlen 200 --mincov 5 --depth 150 --trim yes \
          --cpus "$CPUS" --ram "$RAM" --tmpdir "${TMPDIR:-/tmp}" --force \
          > "$sample_out/shovill.log" 2>&1 || true
  runtime=$(( ( $(date +%s) - start_time ) / 60 ))

  # --- Rename contigs file to include sample name ---
  if [[ -f "$sample_out/contigs.fa" ]]; then
      mv "$sample_out/contigs.fa" "$sample_out/${sample}_contigs.fa"
  fi
  contigs_file="$sample_out/${sample}_contigs.fa"
  # --------------------------------------------------

  total_len=0; n_contigs=0; n50=0
  if [[ -s "$contigs_file" ]]; then
      mapfile -t contig_lengths < <(awk '/^>/{if(seq){print length(seq)}; seq=""} !/^>/{seq=seq $0} END{if(seq) print length(seq)}' "$contigs_file" | sort -nr)
      n_contigs=${#contig_lengths[@]}
      for len in "${contig_lengths[@]}"; do ((total_len+=len)); done
      half=$(( total_len / 2 )); sum=0
      for len in "${contig_lengths[@]}"; do ((sum+=len)); if (( sum >= half )); then n50=$len; break; fi; done
      if [[ "$total_len" -lt 1000000 ]]; then status="Failed"; reason="Very small assembly (${total_len} bp)"; else status="Assembled"; reason="OK"; fi
  else
      status="Failed"; reason="Assembly error (no contigs.fa)"
  fi

  # ---- Fixed and safe Kraken2 command ----
  if [[ -n "$KRAKEN_DB" && -s "$contigs_file" ]]; then
      echo "[INFO] [$sample] Running Kraken2 contamination check..."
      kraken2 --db "$KRAKEN_DB" \
              --threads 4 \
              --use-names \
              --report "$sample_out/kraken_report.txt" \
              --output "$sample_out/kraken_output.txt" \
              "$contigs_file" > "$sample_out/kraken.log" 2>&1 || true

      if [[ -s "$sample_out/kraken_report.txt" ]]; then
          top_species=$(awk -F'\t' 'NR==1{gsub(/^[ \t]+|[ \t]+$/, "", $6); print $6}' "$sample_out/kraken_report.txt" | head -n1)
          if [[ "$top_species" != *"Mycobacterium tuberculosis"* ]]; then
              status="Flagged"
              reason="Possible contamination ($top_species)"
          fi
      else
          echo "[WARN] [$sample] Kraken report empty or no classification results."
      fi
  fi
  # ----------------------------------------

  { flock -x 200; printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$sample" "$total_reads" "$avg_len" "$estimated_coverage" "${total_len:-0}" "${n_contigs:-0}" "${n50:-0}" "$runtime" "$status" "$reason" >> "$REPORT"; } 200>"$REPORT_LOCK"

  echo "[INFO] [$sample] Finished. Status: $status"
}

export -f run_sample
export INPUT_DIR OUTDIR GSIZE MIN_COVERAGE MIN_READS SAMPLE_SEQ_LINES CPUS RAM KRAKEN_DB REPORT REPORT_LOCK

echo "[INFO] Starting parallel Shovill (max 64 GB total)..."
find "$INPUT_DIR" -name "*_1.trim.fastq.gz" | sort | \
  parallel -j "$PARALLEL_JOBS" run_sample {}

echo "[INFO] All samples processed. Final report: $REPORT"

```
Single end  read


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
CSV_OUTDIR="quast_results"
LOG_FILE="$QUAST_DIR/quast_log.txt"

mkdir -p "$QUAST_DIR" "$CSV_OUTDIR"

CSV_FILE="$CSV_OUTDIR/quast_summary_shovill.csv"
echo "Sample,NumContigs,TotalLength,MinLen,MaxLen,AverageLen,N50,N75,GC%" > "$CSV_FILE"

for sample_out in "$SHOVILL_DIR"/*; do
    [[ -d "$sample_out" ]] || continue
    sample=$(basename "$sample_out")

    # Reset variables per sample
    unset num_contigs total_len min_len max_len avg_len n50 n75 gc lengths

    contigs_file=("$sample_out"/*_contigs.fa)
    [[ -f "${contigs_file[0]}" ]] || continue
    contigs="${contigs_file[0]}"

    outdir="$QUAST_DIR/$sample"
    mkdir -p "$outdir"

    # Run QUAST, log errors
    if ! quast "$contigs" -o "$outdir" > "$LOG_FILE" 2>&1; then
        echo "QUAST failed for $sample, see $LOG_FILE" >&2
    fi

    # Extract lengths robustly from headers using len=
    readarray -t lengths < <(grep "^>" "$contigs" | sed -E 's/.*len=([0-9]+).*/\1/')

    if (( ${#lengths[@]} > 0 )); then
        # Basic stats
        num_contigs=${#lengths[@]}
        total_len=$(awk '{sum+=$1} END{print sum}' <<<"${lengths[*]}")
        min_len=$(printf "%s\n" "${lengths[@]}" | sort -n | head -n1)
        max_len=$(printf "%s\n" "${lengths[@]}" | sort -nr | head -n1)
        avg_len=$(awk -v t="$total_len" -v n="$num_contigs" 'BEGIN{if(n>0) print t/n; else print "NA"}')

        # Sort lengths descending for N50/N75
        IFS=$'\n' sorted=($(sort -nr <<<"${lengths[*]}"))
        unset IFS
        cum=0
        n50="NA"
        n75="NA"
        n50_target=$((total_len / 2))
        n75_target=$((total_len * 3 / 4))
        for l in "${sorted[@]}"; do
            cum=$((cum + l))
            [[ "$n50" == "NA" && $cum -ge $n50_target ]] && n50=$l
            [[ "$n75" == "NA" && $cum -ge $n75_target ]] && n75=$l
        done
    else
        # Fallback if lengths extraction failed
        num_contigs=0
        total_len=0
        min_len=NA
        max_len=NA
        avg_len=NA
        n50=NA
        n75=NA
    fi

    # GC% calculation from sequences
    gc=$(awk 'BEGIN{gc=0; total=0} /^[^>]/ {seq=toupper($0); for(i=1;i<=length(seq);i++){b=substr(seq,i,1); if(b=="G"||b=="C")gc++}; total+=length(seq)} END{if(total>0) print (gc/total)*100; else print "NA"}' "$contigs")

    # Write to CSV
    echo "$sample,$num_contigs,$total_len,$min_len,$max_len,$avg_len,$n50,$n75,$gc" >> "$CSV_FILE"

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


#### 4. estimate genome completeness/contamination for all assemblies with CheckM
######   Step 1: Create the Shovill CheckM script
``` bash
nano run_checkm_shovill.sh
``` 
######   Step 2: Paste this into the file:
``` bash
#!/bin/bash
set -euo pipefail

# Enable nullglob so globs that don't match anything expand to nothing
shopt -s nullglob

# Directories
SHOVILL_DIR="shovill_results"
CHECKM_PARENT="checkm_results"
CHECKM_DIR="$CHECKM_PARENT/checkm_results_shovill"
INPUT_DIR="$CHECKM_DIR/input"
LOG_FILE="$CHECKM_DIR/skipped_samples.log"

# Create directories
mkdir -p "$CHECKM_DIR" "$INPUT_DIR"

# Clear previous input files and logs
rm -f "$INPUT_DIR"/*.fasta "$LOG_FILE"

# Check if CheckM is installed
if ! command -v checkm &> /dev/null; then
    echo "Error: CheckM is not installed or not in PATH."
    exit 1
fi

echo "Preparing contigs files for CheckM..."
for sample_out in "$SHOVILL_DIR"/*; do
    [[ -d "$sample_out" ]] || continue
    sample=$(basename "$sample_out")
    
    # Collect contigs files
    contigs_files=("$sample_out"/*_contigs.fa)
    
    # Skip if no contigs found
    if [[ ${#contigs_files[@]} -eq 0 ]]; then
        echo "Warning: No contigs file found for sample $sample"
        echo "$sample : No contigs file found" >> "$LOG_FILE"
        continue
    fi

    # Skip if contigs are empty
    total_size=0
    for f in "${contigs_files[@]}"; do
        if [[ ! -s "$f" ]]; then
            echo "Warning: Contigs file $f is empty"
        else
            (( total_size += $(stat -c%s "$f") ))
        fi
    done

    if [[ $total_size -eq 0 ]]; then
        echo "Warning: All contigs files for sample $sample are empty"
        echo "$sample : All contigs empty" >> "$LOG_FILE"
        continue
    fi

    # Concatenate contigs into a single file for CheckM
    cat "${contigs_files[@]}" > "$INPUT_DIR/${sample}.fasta"
done

# Check if there is any input for CheckM
if [[ $(ls "$INPUT_DIR"/*.fasta 2>/dev/null | wc -l) -eq 0 ]]; then
    echo "No valid contigs found. Exiting."
    exit 0
fi

echo "Running CheckM lineage workflow..."
checkm lineage_wf -x fasta "$INPUT_DIR" "$CHECKM_DIR" -t 8

# Verify lineage.ms exists
if [[ ! -f "$CHECKM_DIR/lineage.ms" ]]; then
    echo "Error: lineage.ms file not found after CheckM run."
    exit 1
fi

# Generate CSV summary
CSV_FILE="$CHECKM_DIR/checkm_summary_shovill.csv"
checkm qa "$CHECKM_DIR/lineage.ms" "$CHECKM_DIR" -o 2 -t 8 > "$CSV_FILE"

echo "CheckM summary saved to $CSV_FILE"

# Summary of skipped samples
if [[ -f "$LOG_FILE" ]]; then
    echo "Skipped samples logged in $LOG_FILE"
fi


``` 
###### Step 3: Make scripts executable
``` bash
chmod +x run_checkm_spades.sh run_checkm_shovill.sh
``` 
###### Step 4: Run the scripts
``` bash
conda activate checkm_env
./run_checkm_shovill.sh
``` 

Backmapping
Shovill

Open a new file in nano
```bash
nano backmap_shovill.sh
```
Paste the script

Copy this into nano:
```bash
#!/bin/bash
set -euo pipefail
shopt -s nullglob

# -----------------------------
# CONFIG
# -----------------------------
ASSEMBLY_DIR="shovill_results"
READ_DIR="fastp_results_min_50"
OUTDIR="backmap_stats"
THREADS=8
MAX_AVG_COVERAGE=1000
MIN_AVG_COVERAGE=1

mkdir -p "$OUTDIR"/{bams,depths,logs}

SUMMARY="$OUTDIR/backmap_summary.tsv"
echo -e "Sample\tTotal_Reads\tProper_Pairs\tPercent_Mapped\tAvg_Coverage\tMedian_Depth\tMin_Depth\tMax_Depth\tMean_Mapping_Quality\tWarnings" > "$SUMMARY"

echo "🔍 Searching for FASTQ files in: $READ_DIR"

# -----------------------------
# MAIN LOOP
# -----------------------------
for R1 in "$READ_DIR"/*_1.trim.fastq.gz; do
    [[ -e "$R1" ]] || { echo "⚠️ No R1 FASTQs found."; break; }

    sample=$(basename "$R1" | sed 's/_1\.trim\.fastq\.gz//')
    R2="$READ_DIR/${sample}_2.trim.fastq.gz"
    [[ -f "$R2" ]] || { echo "⚠️ Missing R2 for $sample"; continue; }

    asm=$(ls "$ASSEMBLY_DIR/${sample}"/*contigs*.fa* 2>/dev/null | head -n1 || true)
    [[ -f "$asm" ]] || { echo "⚠️ No assembly for $sample"; continue; }

    log="$OUTDIR/logs/${sample}.log"
    bam="$OUTDIR/bams/${sample}.bam"
    depth="$OUTDIR/depths/${sample}.depth"
    warnings=""

    echo "🧬 Processing $sample..." | tee "$log"

    # Build BWA index if not found
    [[ ! -f "${asm}.bwt" ]] && bwa index "$asm" >>"$log" 2>&1

    # Align reads and create sorted BAM
    bwa mem -t "$THREADS" "$asm" "$R1" "$R2" 2>>"$log" \
        | samtools view -bS - \
        | samtools sort -@ "$THREADS" -o "$bam"
    samtools index "$bam"

    # -----------------------------
    # Total reads (R1 only)
    # -----------------------------
    total_r1=$(zcat "$R1" | wc -l)
    total_r1=$(( total_r1 / 4 ))
    total=$total_r1

    # -----------------------------
    # Proper paired reads (divide by 2 to avoid double counting)
    # -----------------------------
    proper_pairs_both_mates=$(samtools view -c -f 2 "$bam")
    proper_pairs=$(( proper_pairs_both_mates / 2 ))
    percent_mapped=$(awk -v p=$proper_pairs -v t=$total 'BEGIN{if(t>0) printf "%.2f", (p/t)*100; else print 0}')

    # -----------------------------
    # Mean mapping quality
    # -----------------------------
    mean_mapq=$(samtools view "$bam" | awk '{sum+=$5; n++} END{if(n>0) printf "%.4f", sum/n; else print 0}')

    # -----------------------------
    # Depth calculation (all positions, including zeros)
    # -----------------------------
    samtools depth -a "$bam" > "$depth"
    total_positions=$(wc -l < "$depth" || echo 0)
    if [[ "$total_positions" -eq 0 ]]; then
        echo -e "$sample\t$total\t0\t0\tNA\tNA\tNA\tNA\t$mean_mapq\tNo depth data" >>"$SUMMARY"
        continue
    fi

    # Average coverage
    avgcov=$(awk '{sum+=$3} END{if(NR>0) printf "%.3f", sum/NR; else print 0}' "$depth")

    # Median depth (including zeros)
    mediancov=$(awk '{a[NR]=$3} END{
        n=NR
        if(n==0){print 0}
        else if(n%2==1){print a[(n+1)/2]}
        else{print (a[n/2]+a[n/2+1])/2}
    }' "$depth" | sort -n)

    # Min and max depth
    mindepth=$(awk 'NR==1{min=$3} $3<min{min=$3} END{print min+0}' "$depth")
    maxdepth=$(awk 'NR==1{max=$3} $3>max{max=$3} END{print max+0}' "$depth")

    # -----------------------------
    # Warnings
    # -----------------------------
    if (( $(echo "$percent_mapped > 95" | bc -l) )); then
        warnings="Excellent: >95% properly paired — expected if same strain."
    elif (( $(echo "$percent_mapped >= 80" | bc -l) )); then
        warnings="Good: 80–95% properly paired — same species, slight divergence."
    elif (( $(echo "$percent_mapped >= 50" | bc -l) )); then
        warnings="Moderate/Suspicious: 50–80% properly paired — possible contamination or mismatch."
    elif (( $(echo "$percent_mapped >= 10" | bc -l) )); then
        warnings="Poor: 10–50% properly paired — likely reference mismatch."
    elif (( $(echo "$percent_mapped >= 1" | bc -l) )); then
        warnings="Very Poor: <10% properly paired — likely wrong reference or severe contamination."
    else
        warnings="Critical: <1% properly paired — almost certainly wrong reference or contamination."
    fi

    if (( $(echo "$avgcov > $MAX_AVG_COVERAGE" | bc -l) )); then
        warnings="${warnings}; High coverage (> ${MAX_AVG_COVERAGE}×)"
    elif (( $(echo "$avgcov < $MIN_AVG_COVERAGE" | bc -l) )); then
        warnings="${warnings}; Low coverage (< ${MIN_AVG_COVERAGE}×)"
    fi

    # -----------------------------
    # Output summary
    # -----------------------------
    echo -e "$sample\t$total\t$proper_pairs\t$percent_mapped\t$avgcov\t$mediancov\t$mindepth\t$maxdepth\t$mean_mapq\t$warnings" >> "$SUMMARY"

done

echo "✅ Backmapping summary completed: $SUMMARY"

```

Save and exit nano

Press Ctrl + O → Enter to save

Press Ctrl + X → exit

Make the script executable
``` bash
chmod +x backmap_shovill.sh
``` 
Run the script
``` bash
conda activate backmap_env
./backmap_shovill.sh
``` 


# 1️⃣4️⃣ Prokka
Prokka is a rapid **prokaryotic genome annotation tool** that predicts genes, coding sequences (CDS), rRNAs, tRNAs, and other genomic features from assembled contigs or genomes.  

Key points for TB genomes:

- Annotates **Mycobacterium tuberculosis** genomes with correct taxonomy using `--genus` and `--species`.  
- Produces multiple output files, including **GFF3**, **FASTA of proteins**, and **GenBank format**, which are useful for downstream analysis.  
- Supports **multi-threading** (`--cpus`) to speed up processing of multiple genomes.  
- Works seamlessly with **Shovill-assembled contigs** and   **spades-assembled contigs**
- Output files are organized per sample directory with a consistent naming prefix for easy pipeline integration.  

> ⚠ Note: Prokka relies on the quality of the assembly; fragmented or low-coverage assemblies may result in incomplete annotations.

Prokka for Shovill assemblies
##### Step 1: Create or edit the script
```bash
nano run_prokka_shovill.sh
```
##### Step 2: Paste the following into the script

```bash
#!/bin/bash
set -euo pipefail

# === CONFIGURATION ===
SHOVILL_DIR="shovill_results"
PROKKA_DIR="prokka_results/shovill"
LOG="$PROKKA_DIR/prokka_summary.log"
SUMMARY_TSV="$PROKKA_DIR/prokka_annotation_summary.tsv"
GBK_COLLECT="$PROKKA_DIR/gbk_collection"
FAA_COLLECT="$PROKKA_DIR/faa_collection"

mkdir -p "$PROKKA_DIR" "$GBK_COLLECT" "$FAA_COLLECT"
: > "$LOG"

# === CHECK DEPENDENCIES ===
if ! command -v prokka &>/dev/null; then
  echo "Error: Prokka not found in PATH."
  exit 1
fi

# === FUNCTION ===
run_prokka() {
  sample_out="$1"
  sample=$(basename "$sample_out")

  contigs=("$sample_out"/*contigs*.fa*)
  contigs="${contigs[0]:-}"
  [[ -z "$contigs" ]] && { echo "[WARN] No contigs for $sample" | tee -a "$LOG"; return; }

  outdir="$PROKKA_DIR/$sample"
  [[ -d "$outdir" && -s "$outdir/$sample.gff" ]] && {
    echo "[INFO] Skipping $sample (already annotated)" | tee -a "$LOG"
    return
  }

  echo "[INFO] Annotating $sample..." | tee -a "$LOG"

  if prokka --outdir "$outdir" \
            --prefix "$sample" \
            --locustag "$sample" \
            --kingdom Bacteria \
            --genus Mycobacterium \
            --species tuberculosis \
            --strain "$sample" \
            --cpus 4 \
            --evalue 1e-9 \
            --centre TBLab \
            --compliant \
            --addgenes \
            --force \
            --quiet \
            "$contigs"; then
    echo "[OK] $sample annotated successfully" | tee -a "$LOG"
  else
    echo "[ERROR] Prokka failed for $sample" | tee -a "$LOG"
  fi
}

export -f run_prokka
export PROKKA_DIR LOG

# === RUN ANNOTATIONS IN PARALLEL ===
echo "[INFO] Starting Prokka annotation on Shovill results..."
find "$SHOVILL_DIR" -maxdepth 1 -mindepth 1 -type d | parallel -j "$(nproc)" run_prokka {}

# === SUMMARY TABLE ===
echo "[INFO] Generating annotation summary..."
echo -e "Sample\tGenes\tCDS\trRNA\ttRNA\ttmRNA\tContigs\tBases" > "$SUMMARY_TSV"

for folder in "$PROKKA_DIR"/*; do
  [[ -d "$folder" ]] || continue
  sample=$(basename "$folder")
  stats_file="$folder/${sample}.txt"
  [[ -s "$stats_file" ]] || { echo "[WARN] Missing stats file for $sample" | tee -a "$LOG"; continue; }

  genes=$(grep -m1 '^Genes:' "$stats_file" | awk '{print $2}')
  cds=$(grep -m1 '^CDS:' "$stats_file" | awk '{print $2}')
  rrna=$(grep -m1 '^rRNA:' "$stats_file" | awk '{print $2}')
  trna=$(grep -m1 '^tRNA:' "$stats_file" | awk '{print $2}')
  tmrna=$(grep -m1 '^tmRNA:' "$stats_file" | awk '{print $2}')
  contigs=$(grep -m1 '^Contigs:' "$stats_file" | awk '{print $2}')
  bases=$(grep -m1 '^Bases:' "$stats_file" | awk '{print $2}')

  echo -e "${sample}\t${genes:-0}\t${cds:-0}\t${rrna:-0}\t${trna:-0}\t${tmrna:-0}\t${contigs:-0}\t${bases:-0}" >> "$SUMMARY_TSV"

  # === COLLECT FILES FOR GET_HOMOLOGUES ===
  gbk_src="$folder/${sample}.gbk"
  faa_src="$folder/${sample}.faa"

  if [[ -s "$gbk_src" ]]; then
    cp "$gbk_src" "$GBK_COLLECT/${sample}.gbk"
  fi

  if [[ -s "$faa_src" ]]; then
    cp "$faa_src" "$FAA_COLLECT/${sample}.faa"
  fi
done

# === FINAL MESSAGES ===
echo "[INFO] Annotation summary written to: $SUMMARY_TSV"
echo "[INFO] GenBank files collected in: $GBK_COLLECT"
echo "[INFO] Protein FASTA files collected in: $FAA_COLLECT"
echo "[INFO] You can now run GET_HOMOLOGUES with:"
echo ""
echo "   get_homologues.pl -d $GBK_COLLECT -M -n 10 -t 0"
echo ""
echo "[INFO] Log file saved to: $LOG"


```
##### Step 3: Save and exit nano
Press Ctrl + O → Enter (to write the file)
Press Ctrl + X → Exit nano

###### Step 4: Make the script executable
``` bash
chmod +x run_prokka_shovill.sh
```
###### Step 5: Activate environment and run
``` bash
conda activate prokka_env
./run_prokka_shovill.sh
```

Shovill summary → CSV

``` bash
#!/bin/bash
set -euo pipefail

mkdir -p csv_output
output_file="csv_output/prokka_shovill_summary.csv"
echo "Sample,CDS,rRNA,tRNA,Genome_size,GC_content" > "$output_file"

for sample_dir in prokka_results/prokka_results_shovill/*; do
    sample=$(basename "$sample_dir")
    stats_file="$sample_dir/$sample.txt"
    if [[ -f "$stats_fimerplot visualile" ]]; then
        cds=$(grep "CDS:" "$stats_file" | awk '{print $2}')
        rrna=$(grep "rRNA:" "$stats_file" | awk '{print $2}')
        trna=$(grep "tRNA:" "$stats_file" | awk '{print $2}')
        genome=$(grep "Bases:" "$stats_file" | awk '{print $2}')
        gc=$(grep "GC:" "$stats_file" | awk '{print $2}')
        echo "$sample,$cds,$rrna,$trna,$genome,$gc" >> "$output_file"
    fi
done
``` 

Shovill function → CSV

``` bash
#!/bin/bash

BASE_DIR="./prokka_results/shovill"

for SAMPLE_DIR in "$BASE_DIR"/*/; do
    TSV_FILE=$(find "$SAMPLE_DIR" -maxdepth 1 -type f -name "*.tsv" | head -n 1)
    [ -z "$TSV_FILE" ] && continue
    BASENAME=$(basename "$TSV_FILE" .tsv)
    OUTPUT_FILE="$SAMPLE_DIR/${BASENAME}_product_stats.csv"
    TOTAL=$(tail -n +2 "$TSV_FILE" | wc -l)
    [ "$TOTAL" -eq 0 ] && continue
    echo "Product,Count,Percentage" > "$OUTPUT_FILE"
    tail -n +2 "$TSV_FILE" | cut -f7 | sed '/^$/d' | sort | uniq -c | sort -nr | \
    awk -v total="$TOTAL" '{count=$1; $1=""; sub(/^ /,""); perc=(count/total)*100; printf "\"%s\",%d,%.2f%%\n",$0,count,perc}' >> "$OUTPUT_FILE"
    echo "✅ Saved product stats for $(basename "$SAMPLE_DIR") to $OUTPUT_FILE"
done
``` 

