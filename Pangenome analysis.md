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

INPUT_DIR="fastp_results_min_50"
OUTDIR="shovill_results"

GSIZE=4411532

shopt -s nullglob
for R1 in "$INPUT_DIR"/*_1.trim.fastq.gz; do
  [[ -e "$R1" ]] || continue

  R2="${R1/_1.trim.fastq.gz/_2.trim.fastq.gz}"

  if [[ ! -f "$R2" ]]; then
    echo ">> Skipping $(basename "$R1") (no matching R2 found)" >&2
    continue
  fi

  sample=$(basename "$R1" _1.trim.fastq.gz)
  sample_out="$OUTDIR/$sample"

  if [[ -f "$sample_out/${sample}_contigs.fa" ]]; then
    echo ">> Skipping $sample (already assembled)"
    continue
  fi

  echo "==> Running Shovill (paired-end) on: $sample"
  mkdir -p "$sample_out"

  shovill \
    --R1 "$R1" \
    --R2 "$R2" \
    --gsize "$GSIZE" \
    --outdir "$sample_out" \
    --assembler skesa \
    --minlen 500 \
    --mincov 30 \
    --depth 100 \
    --namefmt "${sample}_%05d" \
    --cpus 16 \
    --ram 120 \
    --tmpdir "${TMPDIR:-/tmp}" \
    --force

  for f in "$sample_out"/*; do
    base=$(basename "$f")
    mv "$f" "$sample_out/${sample}_$base"
  done
done

```
<details>
<summary>📖 Explanation of Shovill Pipeline Script</summary>

- `INPUT_DIR="fastp_results_min_50"` → directory with preprocessed FASTQ files.  
- `OUTDIR="shovill_results"` → directory to store Shovill assemblies.  
- `mkdir -p "$OUTDIR"` → ensures output directory exists.  
- `GSIZE=4411532` → approximate genome size for M. tuberculosis (~4.41 Mb).  
- `shopt -s nullglob` → makes the loop skip if no matching files exist.  

**Loop over samples:**  
- `for R1 in "$INPUT_DIR"/*_1.trim.fastq.gz; do ... done` → loops over all R1 FASTQ files.  
- `R2="${R1/_1.trim.fastq.gz/_2.trim.fastq.gz}"` → guesses corresponding R2 file.  
- `if [[ ! -f "$R2" ]]; then ... fi` → skips sample if R2 is missing.  
- `sample=$(basename "$R1" _1.trim.fastq.gz)` → extracts sample name.  
- `sample_out="$OUTDIR/$sample"` → defines sample-specific output folder.  
- `if [[ -f "$sample_out/${sample}_contigs.fa" ]]; then ... fi` → skips assembly if contigs already exist.  

**Shovill command parameters:**  
- `--R1/--R2` → input paired-end reads  
- `--gsize` → genome size  
- `--outdir` → output directory  
- `--assembler skesa` → use SKESA assembler  
- `--minlen 500` → minimum contig length  
- `--mincov 5` → minimum coverage  
- `--depth 100` → target depth  
- `--namefmt "${sample}_%05d"` → output naming format  
- `--cpus 4` → CPU threads  
- `--ram 16` → RAM in GB  
- `--tmpdir` → temporary directory  
- `--force` → overwrite existing files  

**Post-processing:**  
- `for f in "$sample_out"/*; do ... done` → renames all files to include sample prefix.

</details>
Shovill single-read script

```bash
#!/bin/bash
set -euo pipefail

INPUT_DIR="fastp_results_min_50"
OUTDIR="shovill_results"
mkdir -p "$OUTDIR"

GSIZE=4411532

shopt -s nullglob
for R1 in "$INPUT_DIR"/*.fastq.gz; do
  [[ -e "$R1" ]] || continue

  # Remove .fastq.gz and optional .trim suffix
  sample=$(basename "$R1" .fastq.gz)
  sample=${sample%.trim}

  sample_out="$OUTDIR/$sample"

  if [[ -f "$sample_out/contigs.fa" ]]; then
    echo ">> Skipping $sample (already assembled)"
    continue
  fi

  echo "==> Running Shovill (single-read workaround) on: $sample"
  mkdir -p "$sample_out"

  shovill \
    --R1 "$R1" \
    --R2 "$R1" \
    --gsize "$GSIZE" \
    --outdir "$sample_out" \
    --assembler skesa \
    --minlen 500 \
    --mincov 30 \
    --depth 100 \
    --namefmt "%05d" \
    --cpus 16 \
    --ram 120 \
    --tmpdir "${TMPDIR:-/tmp}" \
    --force
done
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


# 2️⃣ Spades

SPAdes (St. Petersburg genome assembler) is a popular **de novo assembler** used for bacterial genomes, including *Mycobacterium tuberculosis*. It works well with **paired-end Illumina reads**, especially after trimming and quality control.

---

#### Why SPAdes for TB?

- **Accurate short-read assembly**: Handles Illumina reads effectively.
- **Supports multiple library types**: Paired-end, single-end, and mate-pair.
- **Produces high-quality contigs**: Useful for downstream analyses like variant calling, resistance profiling, and phylogenetics.
- **Customizable parameters**: Memory, threads, k-mer sizes, coverage cutoffs.
---

  
###### Step 1: Activate the environment
``` bash
nano run_spades.sh
```

###### Step 2: Copy–paste the following script into the nano editor:
```bash
#!/bin/bash
set -euo pipefail

INPUT_DIR="fastp_results_min_50"
OUTDIR="spades_results"

THREADS=16
MEMORY=120
MIN_CONTIG=500
COV_CUTOFF=30

shopt -s nullglob
for R1 in "$INPUT_DIR"/*_1.trim.fastq.gz; do
  [[ -e "$R1" ]] || continue

  R2="${R1/_1.trim.fastq.gz/_2.trim.fastq.gz}"
  if [[ ! -f "$R2" ]]; then
    echo ">> Skipping $(basename "$R1") (no matching R2 found)" >&2
    continue
  fi

  sample=$(basename "$R1" _1.trim.fastq.gz)
  sample_out="$OUTDIR/$sample"

  if [[ -f "$sample_out/${sample}_contigs.fasta" ]]; then
    echo ">> Skipping $sample (already assembled)"
    continue
  fi

  echo "==> Running SPAdes (paired-end) on: $sample"
  mkdir -p "$sample_out"

  spades.py \
    -1 "$R1" \
    -2 "$R2" \
    -o "$sample_out" \
    -t "$THREADS" \
    -m "$MEMORY" \
    --only-assembler \
    --cov-cutoff "$COV_CUTOFF" \
    > "$sample_out/${sample}_spades.log" 2>&1

  if [[ -f "$sample_out/contigs.fasta" ]]; then
    awk -v minlen="$MIN_CONTIG" 'BEGIN{RS=">"; ORS=""} length($0)>minlen+1 {print ">"$0}' \
      "$sample_out/contigs.fasta" > "$sample_out/${sample}_contigs.fasta"
  fi
done

```
single-end version
```bash
#!/bin/bash
set -euo pipefail

INPUT_DIR="fastp_results_min_50"
OUTDIR="spades_results"
mkdir -p "$OUTDIR"

THREADS=16
MEMORY=120
MIN_CONTIG=500
COV_CUTOFF=30

shopt -s nullglob
for R1 in "$INPUT_DIR"/*.fastq.gz; do
  [[ -e "$R1" ]] || continue
  sample=$(basename "$R1" .fastq.gz)
  sample=${sample%.trim}
  sample_out="$OUTDIR/$sample"
  if [[ -f "$sample_out/${sample}_contigs.fasta" ]]; then
    echo ">> Skipping $sample (already assembled)"
    continue
  fi
  echo "==> Running SPAdes (single-end) on: $sample"
  mkdir -p "$sample_out"
  spades.py -s "$R1" -o "$sample_out" -t "$THREADS" -m "$MEMORY" --only-assembler --cov-cutoff "$COV_CUTOFF" > "$sample_out/spades.log" 2>&1
  if [[ -f "$sample_out/contigs.fasta" ]]; then
    awk -v minlen="$MIN_CONTIG" 'BEGIN{RS=">"; ORS=""} length($0)>minlen+1 {print ">"$0}' \
      "$sample_out/contigs.fasta" > "$sample_out/${sample}_contigs.filtered.fasta"
    mv "$sample_out/contigs.fasta" "$sample_out/${sample}_contigs.fasta"
  fi
done

```


###### Step 3: Save & exit in nano:

Press CTRL + O → hit ENTER (to save).
Press CTRL + X (to exit).

###### Step 4: Make the script executable:
```bash
chmod +x run_spades.sh
```
###### Step 5: Activate and Run it:
```bash
conda activate spades_env
./run_spades.sh
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
######   Step 1: Open nano and create the script
``` bash
nano quast_compare_parallel.sh
``` 
Paste this into the file:

``` bash
#!/bin/bash
set -euo pipefail

SPADES_DIR="spades_results"
SHOVILL_DIR="shovill_results"
QUAST_PARENT="quast_results"
CSV_OUTDIR="csv_output"

mkdir -p "$QUAST_PARENT" "$CSV_OUTDIR"

SPADES_QUAST="$QUAST_PARENT/quast_results_spades"
SHOVILL_QUAST="$QUAST_PARENT/quast_results_shovill"

SPADES_CSV="$CSV_OUTDIR/quast_summary_spades.csv"
SHOVILL_CSV="$CSV_OUTDIR/quast_summary_shovill.csv"
COMBINED_CSV="$CSV_OUTDIR/quast_summary_combined.csv"

rm -f "$SPADES_CSV" "$SHOVILL_CSV" "$COMBINED_CSV"

echo "Sample,NumContigs,TotalLength,MinLen,MaxLen,AverageLen,N50,N75,GC%" > "$SPADES_CSV"
echo "Sample,NumContigs,TotalLength,MinLen,MaxLen,AverageLen,N50,N75,GC%" > "$SHOVILL_CSV"

compute_metrics() {
    CONTIG_FILE="$1"
    lengths=$(awk '/^>/{if(seq) print length(seq); seq=""} !/^>/{seq=seq $0} END{if(seq) print length(seq)}' "$CONTIG_FILE")
    NUMCONTIGS=$(echo "$lengths" | wc -l)
    TOTALLEN=$(echo "$lengths" | awk '{s+=$1} END{print s}')
    MINLEN=$(echo "$lengths" | sort -n | head -1)
    MAXLEN=$(echo "$lengths" | sort -nr | head -1)
    AVG=$(echo "$lengths" | awk -v n="$NUMCONTIGS" '{s+=$1} END{print s/n}')
    lengths_desc=$(echo "$lengths" | sort -nr)
    N50=$(awk -v total="$TOTALLEN" 'BEGIN{c=0} {c+=$1; if(c>=total/2){print $1; exit}}' <<< "$lengths_desc")
    N75=$(awk -v total="$TOTALLEN" 'BEGIN{c=0} {c+=$1; if(c>=total*0.75){print $1; exit}}' <<< "$lengths_desc")
    GC=$(awk '/^>/{next} {g+=gsub(/[Gg]/,"")+gsub(/[Cc]/,""); n+=length($0)} END{if(n>0) print g/n*100}' "$CONTIG_FILE")
    echo "$NUMCONTIGS,$TOTALLEN,$MINLEN,$MAXLEN,$AVG,$N50,$N75,$GC"
}

run_quast_parallel() {
    ASSEMBLER_DIR="$1"
    CONTIG_PATTERN="$2"
    QUAST_DIR="$3"
    CSV_FILE="$4"

    mkdir -p "$QUAST_DIR"

    for sample_out in "$ASSEMBLER_DIR"/*; do
        [[ -d "$sample_out" ]] || continue
        sample=$(basename "$sample_out")
        contigs_file=("$sample_out"/$CONTIG_PATTERN)
        [[ -f "${contigs_file[0]}" ]] || continue
        contigs="${contigs_file[0]}"
        outdir="$QUAST_DIR/$sample"
        mkdir -p "$outdir"

        (
            quast "$contigs" -o "$outdir" --min-contig 0 > /dev/null 2>&1
            stats_file="$outdir/report.tsv"
            if [[ -f "$stats_file" ]]; then
                NUMCONTIGS=$(awk -F'\t' '$1~/contigs \(>= 0 bp\)/{print $2}' "$stats_file")
                TOTALLEN=$(awk -F'\t' '$1~/Total length \(>= 0 bp\)/{print $2}' "$stats_file")
                MINLEN=$(awk -F'\t' '$1~/Shortest contig|Min contig/{print $2}' "$stats_file")
                MAXLEN=$(awk -F'\t' '$1~/Largest contig|max contig/{print $2}' "$stats_file")
                AVG=$(awk -F'\t' '$1~/Average contig length|Mean contig length/{print $2}' "$stats_file")
                N50=$(awk -F'\t' '$1=="N50"{print $2}' "$stats_file")
                N75=$(awk -F'\t' '$1=="N75"{print $2}' "$stats_file")
                GC=$(awk -F'\t' '$1~/GC/ {print $2}' "$stats_file")
                if [[ -z "$MINLEN" || "$MINLEN" == "NA" ]]; then
                    metrics=$(compute_metrics "$contigs")
                    MINLEN=$(echo $metrics | cut -d, -f3)
                    AVG=$(echo $metrics | cut -d, -f5)
                    N75=$(echo $metrics | cut -d, -f7)
                fi
                echo "$sample,${NUMCONTIGS:-NA},${TOTALLEN:-NA},${MINLEN:-NA},${MAXLEN:-NA},${AVG:-NA},${N50:-NA},${N75:-NA},${GC:-NA}" >> "$CSV_FILE"
            fi
        ) &
    done
    wait
}

run_quast_parallel "$SPADES_DIR" "*_contigs.fasta" "$SPADES_QUAST" "$SPADES_CSV" &
run_quast_parallel "$SHOVILL_DIR" "*contigs*.fa*" "$SHOVILL_QUAST" "$SHOVILL_CSV" &
wait

join -t, -1 1 -2 1 \
    <(tail -n +2 "$SPADES_CSV" | sort -t, -k1,1) \
    <(tail -n +2 "$SHOVILL_CSV" | sort -t, -k1,1) | \
awk -F, 'BEGIN {OFS=","} {print $1,$2,$10,$3,$11,$4,$12,$5,$13,$6,$14,$7,$15,$8,$16,$9,$17}' > "$COMBINED_CSV.tmp"

echo "Sample,NumContigs_SPADES,NumContigs_SHOVILL,TotalLength_SPADES,TotalLength_SHOVILL,MinLen_SPADES,MinLen_SHOVILL,MaxLen_SPADES,MaxLen_SHOVILL,AverageLen_SPADES,AverageLen_SHOVILL,N50_SPADES,N50_SHOVILL,N75_SPADES,N75_SHOVILL,GC_SPADES,GC_SHOVILL" | cat - "$COMBINED_CSV.tmp" > "$COMBINED_CSV"

rm -f "$SPADES_CSV" "$SHOVILL_CSV" "$COMBINED_CSV.tmp"

``` 
Delete any leftover temp files from previous runs:
``` bash
rm -f csv_output/spades_tmp.csv csv_output/shovill_tmp.csv csv_output/combined_raw.csv
``` 
###### Step 3: Make scripts executable
``` bash
chmod +x quast_compare_parallel.sh
``` 
###### Step 4: Run the scripts
``` bash
conda activate quast_env
./quast_compare_parallel.sh
``` 

### 3. Assembly summary with assembly-scan
Collect the key statistics in a single CSV file
running QUAST on all SPAdes assemblies and collecting key statistics into a single report.tsv. 
######   Step 1: Create the SPAdes assembly-scan  script
``` bash
nano run_assembly_scan_spades.sh
``` 
Paste this into the file:

``` bash
#!/bin/bash
set -euo pipefail

SPADES_DIR="spades_results"
CSV_OUTDIR="csv_output"
mkdir -p "$CSV_OUTDIR"
CSV_FILE="$CSV_OUTDIR/spades_assembly_scan.csv"

echo "Sample,total_contig,total_contig_length,max_contig_length,mean_contig_length,median_contig_length,min_contig_length,n50_contig_length,l50_contig_count,contig_percent_a,contig_percent_c,contig_percent_g,contig_percent_t,contigs_greater_1k,contigs_greater_10k,contigs_greater_100k" > "$CSV_FILE"

for sample_dir in "$SPADES_DIR"/*/; do
    sample=$(basename "$sample_dir")
    contig_file="${sample_dir}${sample}_contigs.fasta"

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
######   Step 2:Create the Shovill assembly-scan script
``` bash
nano run_assembly_scan_shovill.sh
``` 
Paste this into the file:
``` bash
#!/bin/bash
set -euo pipefail

SHOVILL_DIR="shovill_results"
CSV_OUTDIR="csv_output"
mkdir -p "$CSV_OUTDIR"
CSV_FILE="$CSV_OUTDIR/shovill_assembly_scan.csv"

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
###### Step 3: Make scripts executable
``` bash
chmod +x run_assembly_scan_spades.sh run_assembly_scan_shovill.sh
``` 
###### Step 4: Run the scripts
``` bash
conda activate assembly_scan_env
./run_assembly_scan_spades.sh
./run_assembly_scan_shovill.sh
``` 

#### 4. estimate genome completeness/contamination for all assemblies with CheckM
######   Step 1: Create the SPAdes CheckM script
``` bash
nano run_checkm_spades.sh
``` 
Paste this into the file:

``` bash
#!/bin/bash
set -euo pipefail

SPADES_DIR="spades_results"
CHECKM_PARENT="checkm_results"
CHECKM_DIR="$CHECKM_PARENT/checkm_results_spades"
CSV_OUTDIR="csv_output"

mkdir -p "$CHECKM_DIR" "$CSV_OUTDIR"
INPUT_DIR="$CHECKM_DIR/input"
mkdir -p "$INPUT_DIR"

# Copy all contigs to a single input folder with standardized names
for sample_out in "$SPADES_DIR"/*; do
  [[ -d "$sample_out" ]] || continue
  sample=$(basename "$sample_out")
  contigs_file=("$sample_out"/*_contigs.fasta)
  [[ -f "${contigs_file[0]}" ]] || continue
  cp "${contigs_file[0]}" "$INPUT_DIR/${sample}.fasta"
done

# Run CheckM lineage workflow (completeness, contamination, heterogeneity)
checkm lineage_wf -x fasta "$INPUT_DIR" "$CHECKM_DIR" -t 8

# Generate CSV summary
CSV_FILE="$CSV_OUTDIR/checkm_summary_spades.csv"
checkm qa "$CHECKM_DIR/lineage.ms" "$CHECKM_DIR" -o 2 -t 8 > "$CSV_FILE"

``` 
######   Step 2:Create the Shovill CheckM script
``` bash
nano run_checkm_shovill.sh
``` 
Paste this into the file:
``` bash
#!/bin/bash
set -euo pipefail

SHOVILL_DIR="shovill_results"
CHECKM_PARENT="checkm_results"
CHECKM_DIR="$CHECKM_PARENT/checkm_results_shovill"
CSV_OUTDIR="csv_output"

mkdir -p "$CHECKM_DIR" "$CSV_OUTDIR"
INPUT_DIR="$CHECKM_DIR/input"
mkdir -p "$INPUT_DIR"

# Copy all contigs to a single input folder
for sample_out in "$SHOVILL_DIR"/*; do
  [[ -d "$sample_out" ]] || continue
  sample=$(basename "$sample_out")
  contigs_file=("$sample_out"/*_contigs.fa)
  [[ -f "${contigs_file[0]}" ]] || continue
  cp "${contigs_file[0]}" "$INPUT_DIR/${sample}.fasta"
done

# Run CheckM lineage workflow
checkm lineage_wf -x fasta "$INPUT_DIR" "$CHECKM_DIR" -t 8

# Generate CSV summary
CSV_FILE="$CSV_OUTDIR/checkm_summary_shovill.csv"
checkm qa "$CHECKM_DIR/lineage.ms" "$CHECKM_DIR" -o 2 -t 8 > "$CSV_FILE"

``` 
###### Step 3: Make scripts executable
``` bash
chmod +x run_checkm_spades.sh run_checkm_shovill.sh
``` 
###### Step 4: Run the scripts
``` bash
conda activate checkm_env

./run_checkm_spades.sh
./run_checkm_shovill.sh
``` 
5. 


###### Step 1 — Open a new script
``` bash
nano run_mummer_comparison.sh
``` 
###### Step 2 — Paste this clean script
``` bash
#!/bin/bash

SHOVILL_DIR="shovill_results"
SPADES_DIR="spades_results"
OUTDIR="mummer_results"
mkdir -p "$OUTDIR"

ERROR_LOG="$OUTDIR/errors.log"
> "$ERROR_LOG"

if ! command -v gnuplot &>/dev/null; then
    echo "⚠️ gnuplot not found. PNG plots will be skipped."
fi

for SHOVILL_SAMPLE in "$SHOVILL_DIR"/*; do
    SAMPLE_ID=$(basename "$SHOVILL_SAMPLE")
    SHOVILL="$SHOVILL_SAMPLE/${SAMPLE_ID}_contigs.fa"
    SPADES="$SPADES_DIR/$SAMPLE_ID/contigs.fasta"

    if [[ ! -f "$SHOVILL" || ! -f "$SPADES" ]]; then
        echo "Skipping $SAMPLE_ID (missing files)" | tee -a "$ERROR_LOG"
        continue
    fi

    echo "🔍 Comparing $SAMPLE_ID..."
    SAMPLE_OUTDIR="$OUTDIR/$SAMPLE_ID"
    mkdir -p "$SAMPLE_OUTDIR"
    PREFIX="$SAMPLE_OUTDIR/$SAMPLE_ID"

    {
        nucmer --mincluster=500 --prefix="$PREFIX" "$SHOVILL" "$SPADES" || { echo "⚠️ nucmer failed for $SAMPLE_ID" | tee -a "$ERROR_LOG"; continue; }
        delta-filter -1 -l 1000 "$PREFIX.delta" > "$PREFIX.filtered.delta" || { echo "⚠️ delta-filter failed for $SAMPLE_ID" | tee -a "$ERROR_LOG"; continue; }

        if command -v gnuplot &>/dev/null; then
            mummerplot --png --large --layout --color -p "$PREFIX" "$PREFIX.filtered.delta" || echo "⚠️ mummerplot failed for $SAMPLE_ID" | tee -a "$ERROR_LOG"
        fi

        dnadiff -p "$PREFIX" "$SHOVILL" "$SPADES" || { echo "⚠️ dnadiff failed for $SAMPLE_ID" | tee -a "$ERROR_LOG"; continue; }

        SNP_FILE="$PREFIX.snps"
        show-snps -ClrT "$PREFIX.filtered.delta" > "$SNP_FILE"

        if [[ -s "$SNP_FILE" ]]; then
            NUM_SNPS=$(($(wc -l < "$SNP_FILE") - 1))
            echo "🧬 $NUM_SNPS SNPs detected for $SAMPLE_ID"
        else
            NUM_SNPS=0
            echo "⚠️ No SNPs detected for $SAMPLE_ID" | tee -a "$ERROR_LOG"
        fi

        echo "✅ Completed $SAMPLE_ID"
    } || {
        echo "⚠️ Unexpected error for $SAMPLE_ID" | tee -a "$ERROR_LOG"
        continue
    }
done

echo "📊 MUMmer processing done. Check $OUTDIR for results."
``` 

###### Step 3 — Save and exit

Press CTRL + O → Enter → CTRL + X

###### Step 4 — Make it executable

``` bash
chmod +x run_mummer_comparison.sh
``` 
###### Step 5 — Run it
``` bash
conda activate mummer_env
./run_mummer_comparison.sh
``` 
Step 1: Open nano
``` bash
nano generate_summary.sh
``` 
Step 2: Paste the script

Inside nano, paste the following:
``` bash
#!/bin/bash

OUTDIR="mummer_results"
SUMMARY="$OUTDIR/summary.tsv"
ERROR_LOG="$OUTDIR/errors_summary.log"

echo -e "Sample\tShovill_Length\tSPAdes_Length\tTotal_Aligned_Bases\tPercent_Aligned\tAvg_Identity\tTotal_SNPs\tLength_Difference\tPercent_Length_Diff" > "$SUMMARY"
> "$ERROR_LOG"

for REPORT in "$OUTDIR"/*/*.report; do
    SAMPLE_ID=$(basename "$REPORT" .report)

    if [[ ! -f "$REPORT" ]]; then
        echo "⚠️ Report missing for $SAMPLE_ID" | tee -a "$ERROR_LOG"
        continue
    fi

    Shovill_Len=$(grep 'TotalBases' "$REPORT" | head -1 | awk '{print $2}')
    Spades_Len=$(grep 'TotalBases' "$REPORT" | head -1 | awk '{print $3}')

    Total_Aligned_Bases=$(grep 'AlignedBases' "$REPORT" | head -1 | awk '{print $2}' | sed 's/([^)]*)//')
    Percent_Aligned=$(grep 'AlignedBases' "$REPORT" | head -1 | awk '{print $2}' | sed 's/.*(\([0-9.]*\)%).*/\1/')

    Avg_Identity=$(grep -A1 'M-to-M' "$REPORT" | tail -1 | awk '{print $3}')

    SNP_FILE="${REPORT%.report}.snps"
    if [[ -f "$SNP_FILE" ]]; then
        Total_SNPs=$(awk 'NR>3 && $2 != "." && $3 != "." {count++} END {print count}' "$SNP_FILE")
    else
        Total_SNPs="NA"
        echo "⚠️ SNP file missing for $SAMPLE_ID" | tee -a "$ERROR_LOG"
    fi

    if [[ -z "$Shovill_Len" || -z "$Spades_Len" ]]; then
        echo "⚠️ Could not extract lengths for $SAMPLE_ID" | tee -a "$ERROR_LOG"
        continue
    fi

    Length_Diff=$((Spades_Len - Shovill_Len))
    Abs_Diff=${Length_Diff#-}
    Percent_Diff=$(awk -v d=$Abs_Diff -v s=$Shovill_Len 'BEGIN{printf "%.2f", (d/s)*100}')

    echo -e "$SAMPLE_ID\t$Shovill_Len\t$Spades_Len\t$Total_Aligned_Bases\t$Percent_Aligned\t$Avg_Identity\t$Total_SNPs\t$Length_Diff\t$Percent_Diff" >> "$SUMMARY"
done

echo "📊 Summary generated: $SUMMARY"
echo "Errors logged (if any): $ERROR_LOG"

``` 
Step 3: Save and exit nano

Press Ctrl + O → then Enter to save.

Press Ctrl + X → to exit nano.
Step 4: Make it executable
```bash
chmod +x generate_summary.sh
```
Step 5: Run the script
```bash
./generate_summary.sh
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

ASSEMBLY_DIR="shovill_results"
READS_DIR="raw_data"
OUTDIR="backmapping_results/shovill"
THREADS=8

mkdir -p "$OUTDIR"/{bams,depths,logs}

SUMMARY_FILE="$OUTDIR/summary.tsv"
echo -e "Sample\tTotal_Reads\tMapped_Reads\tPercent_Mapped\tAverage_Coverage" > "$SUMMARY_FILE"

for asm_dir in ${ASSEMBLY_DIR}/SRR*; do
    sample=$(basename "$asm_dir")
    asm=$(ls "$asm_dir"/*contigs*.fa* 2>/dev/null | head -n1)

    if [[ ! -f "$asm" ]]; then
        echo "Assembly not found for $sample" | tee -a "$OUTDIR/logs/$sample.log"
        continue
    fi

    r1="$READS_DIR/${sample}_1.fastq.gz"
    r2="$READS_DIR/${sample}_2.fastq.gz"

    if [[ ! -f "$r1" || ! -f "$r2" ]]; then
        echo "Reads not found for $sample" | tee -a "$OUTDIR/logs/$sample.log"
        continue
    fi

    echo "[$sample] Mapping reads..." | tee "$OUTDIR/logs/$sample.log"

    if [[ ! -f "${asm}.bwt" ]]; then
        bwa index "$asm"
    fi

    bwa mem -t "$THREADS" "$asm" "$r1" "$r2" \
      | samtools view -bS - \
      | samtools sort -@ "$THREADS" -o "$OUTDIR/bams/${sample}.bam"

    samtools index "$OUTDIR/bams/${sample}.bam"

    samtools depth -a "$OUTDIR/bams/${sample}.bam" > "$OUTDIR/depths/${sample}.depth"

    total_r1=$(zcat "$r1" | wc -l)
    total_r2=$(zcat "$r2" | wc -l)
    total=$(( (total_r1 + total_r2)/4 ))

    mapped=$(samtools view -c -F 4 "$OUTDIR/bams/${sample}.bam")

    percent=$(awk -v m=$mapped -v t=$total 'BEGIN{printf "%.2f", (m/t)*100}')

    avgcov=$(awk '{sum+=$3} END{if(NR>0) print sum/NR; else print 0}' "$OUTDIR/depths/${sample}.depth")

    echo -e "$sample\t$total\t$mapped\t$percent\t$avgcov" >> "$SUMMARY_FILE"
done

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

SPAdes back-mapping script
Step 1: Open a new file
``` bash
nano backmap_spades.sh
``` 
Step 2: Paste this script
``` bash
#!/bin/bash
set -euo pipefail

ASSEMBLY_DIR="spades_results"
READS_DIR="raw_data"
OUTDIR="backmapping_results/spades"
THREADS=8

mkdir -p "$OUTDIR"/{bams,depths,logs}

SUMMARY_FILE="$OUTDIR/summary.tsv"
echo -e "Sample\tTotal_Reads\tMapped_Reads\tPercent_Mapped\tAverage_Coverage" > "$SUMMARY_FILE"

for asm_dir in ${ASSEMBLY_DIR}/SRR*; do
    sample=$(basename "$asm_dir")
    asm=$(ls "$asm_dir"/*contigs*.fa* 2>/dev/null | head -n1)

    if [[ ! -f "$asm" ]]; then
        echo "Assembly not found for $sample" | tee -a "$OUTDIR/logs/$sample.log"
        continue
    fi

    r1="$READS_DIR/${sample}_1.fastq.gz"
    r2="$READS_DIR/${sample}_2.fastq.gz"

    if [[ ! -f "$r1" || ! -f "$r2" ]]; then
        echo "Reads not found for $sample" | tee -a "$OUTDIR/logs/$sample.log"
        continue
    fi

    echo "[$sample] Mapping reads..." | tee "$OUTDIR/logs/$sample.log"

    if [[ ! -f "${asm}.bwt" ]]; then
        bwa index "$asm"
    fi

    bwa mem -t "$THREADS" "$asm" "$r1" "$r2" \
      | samtools view -bS - \
      | samtools sort -@ "$THREADS" -o "$OUTDIR/bams/${sample}.bam"

    samtools index "$OUTDIR/bams/${sample}.bam"

    samtools depth -a "$OUTDIR/bams/${sample}.bam" > "$OUTDIR/depths/${sample}.depth"

    total_r1=$(zcat "$r1" | wc -l)
    total_r2=$(zcat "$r2" | wc -l)
    total=$(( (total_r1 + total_r2)/4 ))

    mapped=$(samtools view -c -F 4 "$OUTDIR/bams/${sample}.bam")

    percent=$(awk -v m=$mapped -v t=$total 'BEGIN{printf "%.2f", (m/t)*100}')

    avgcov=$(awk '{sum+=$3} END{if(NR>0) print sum/NR; else print 0}' "$OUTDIR/depths/${sample}.depth")

    echo -e "$sample\t$total\t$mapped\t$percent\t$avgcov" >> "$SUMMARY_FILE"
done

``` 
Step 3: Save & exit

Ctrl + O → Enter → Ctrl + X
Step 4: Make executable
``` bash
chmod +x backmap_spades.sh
``` 
Run the script
``` bash
conda activate backmap_env
./backmap_spades.sh
```
ready-to-run Python script that will:

Load our summary.tsv files for Shovill and SPAdes.

Compare Percent_Mapped and Average_Coverage per sample.

Highlight which assembler performs better per sample.

Generate a simple bar plot for visual comparison.

Save the script

In your terminal:
``` bash
nano compare_assemblers.py
``` 
paste
``` bash
#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

shovill = pd.read_csv("backmapping_results/shovill/summary.tsv", sep="\t")
spades = pd.read_csv("backmapping_results/spades/summary.tsv", sep="\t")

merged = pd.merge(shovill, spades, on="Sample", suffixes=("_shovill", "_spades"))

merged['Best_Percent_Mapped'] = np.where(
    merged['Percent_Mapped_shovill'] >= merged['Percent_Mapped_spades'], "Shovill", "SPAdes"
)
merged['Best_Avg_Coverage'] = np.where(
    merged['Average_Coverage_shovill'] >= merged['Average_Coverage_spades'], "Shovill", "SPAdes"
)

merged.to_csv("backmapping_results/assembly_comparison_summary.tsv", sep="\t", index=False)

output_dir = "backmapping_results"
os.makedirs(output_dir, exist_ok=True)

x = np.arange(len(merged))
width = 0.35

plt.figure(figsize=(10,6))
plt.bar(x - width/2, merged['Percent_Mapped_shovill'], width, label='Shovill')
plt.bar(x + width/2, merged['Percent_Mapped_spades'], width, label='SPAdes')
plt.xticks(x, merged['Sample'], rotation=90)
plt.ylabel('Percent Mapped Reads')
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "percent_mapped_comparison.png"), dpi=300)
plt.close()

plt.figure(figsize=(10,6))
plt.bar(x - width/2, merged['Average_Coverage_shovill'], width, label='Shovill')
plt.bar(x + width/2, merged['Average_Coverage_spades'], width, label='SPAdes')
plt.xticks(x, merged['Sample'], rotation=90)
plt.ylabel('Average Coverage')
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "average_coverage_comparison.png"), dpi=300)
plt.close()

```
Run the script
``` bash
python3 compare_assemblers.py
``` 

MUMMER
create nano script 
``` bash
run_mummer_comparison.sh
``` 

paste
``` bash
#!/bin/bash

SHOVILL_DIR="shovill_results"
SPADES_DIR="spades_results"
OUTDIR="mummer_results"
mkdir -p "$OUTDIR"

ERROR_LOG="$OUTDIR/errors.log"
> "$ERROR_LOG"

if ! command -v gnuplot &>/dev/null; then
    echo "⚠️ gnuplot not found. PNG plots will be skipped."
fi

for SHOVILL_SAMPLE in "$SHOVILL_DIR"/*; do
    SAMPLE_ID=$(basename "$SHOVILL_SAMPLE")
    SHOVILL="$SHOVILL_SAMPLE/${SAMPLE_ID}_contigs.fa"
    SPADES="$SPADES_DIR/$SAMPLE_ID/contigs.fasta"

    if [[ ! -f "$SHOVILL" || ! -f "$SPADES" ]]; then
        echo "Skipping $SAMPLE_ID (missing files)" | tee -a "$ERROR_LOG"
        continue
    fi

    echo "🔍 Comparing $SAMPLE_ID..."
    SAMPLE_OUTDIR="$OUTDIR/$SAMPLE_ID"
    mkdir -p "$SAMPLE_OUTDIR"
    PREFIX="$SAMPLE_OUTDIR/$SAMPLE_ID"

    {
        nucmer --mincluster=500 --prefix="$PREFIX" "$SHOVILL" "$SPADES" || { echo "⚠️ nucmer failed for $SAMPLE_ID" | tee -a "$ERROR_LOG"; continue; }
        delta-filter -1 -l 1000 "$PREFIX.delta" > "$PREFIX.filtered.delta" || { echo "⚠️ delta-filter failed for $SAMPLE_ID" | tee -a "$ERROR_LOG"; continue; }

        if command -v gnuplot &>/dev/null; then
            mummerplot --png --large --layout --color -p "$PREFIX" "$PREFIX.filtered.delta" || echo "⚠️ mummerplot failed for $SAMPLE_ID" | tee -a "$ERROR_LOG"
        fi

        dnadiff -p "$PREFIX" "$SHOVILL" "$SPADES" || { echo "⚠️ dnadiff failed for $SAMPLE_ID" | tee -a "$ERROR_LOG"; continue; }

        SNP_FILE="$PREFIX.snps"
        show-snps -ClrT "$PREFIX.filtered.delta" > "$SNP_FILE"

        if [[ -s "$SNP_FILE" ]]; then
            NUM_SNPS=$(($(wc -l < "$SNP_FILE") - 1))
            echo "🧬 $NUM_SNPS SNPs detected for $SAMPLE_ID"
        else
            NUM_SNPS=0
            echo "⚠️ No SNPs detected for $SAMPLE_ID" | tee -a "$ERROR_LOG"
        fi

        echo "✅ Completed $SAMPLE_ID"
    } || {
        echo "⚠️ Unexpected error for $SAMPLE_ID" | tee -a "$ERROR_LOG"
        continue
    }
done

echo "📊 MUMmer processing done. Check $OUTDIR for results."

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

SPADES_DIR="spades_results"
PROKKA_DIR="prokka_results/spades"

mkdir -p "$PROKKA_DIR"

if ! command -v prokka &>/dev/null; then
  echo "Prokka not found"
  exit 1
fi

run_prokka() {
  sample_out="$1"
  sample=$(basename "$sample_out")

  contigs=("$sample_out"/*_contigs.fasta)
  contigs="${contigs[0]:-}"
  [[ -z "$contigs" ]] && return

  echo "Processing $sample..."

  outdir="$PROKKA_DIR/$sample"
  mkdir -p "$outdir"

  prokka --outdir "$outdir" \
         --prefix "$sample" \
         --kingdom Bacteria \
         --genus Mycobacterium \
         --species tuberculosis \
         --cpus 4 \
         --evalue 1e-9 \
         --coverage 90 \
         --force "$contigs"
}

export -f run_prokka
export PROKKA_DIR

find "$SPADES_DIR" -maxdepth 1 -mindepth 1 -type d | parallel -j 8 run_prokka {}
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

Prokka for SPAdes assemblies
##### Step 1: Create or edit the script
```bash
nano run_prokka_spades.sh
```
##### Step 2: Paste the following into the script

```bash
#!/bin/bash
set -euo pipefail

SPADES_DIR="spades_results"
PROKKA_DIR="prokka_results/spades"

mkdir -p "$PROKKA_DIR"

if ! command -v prokka &>/dev/null; then
  echo "Prokka not found"
  exit 1
fi

run_prokka() {
  sample_out="$1"
  sample=$(basename "$sample_out")

  contigs=("$sample_out"/*_contigs.fasta)
  contigs="${contigs[0]:-}"
  [[ -z "$contigs" ]] && return

  echo "Processing $sample..."

  outdir="$PROKKA_DIR/$sample"
  mkdir -p "$outdir"

  prokka --outdir "$outdir" \
         --prefix "$sample" \
         --kingdom Bacteria \
         --genus Mycobacterium \
         --species tuberculosis \
         --cpus 4 \
         --evalue 1e-9 \
         --coverage 90 \
         --force "$contigs"
}

export -f run_prokka
export PROKKA_DIR

find "$SPADES_DIR" -maxdepth 1 -mindepth 1 -type d | parallel -j 8 run_prokka {}
```
##### Step 3: Save and exit nano
Press Ctrl + O → Enter (to write the file)
Press Ctrl + X → Exit nano

###### Step 4: Make the script executable
``` bash
chmod +x run_prokka_spades.sh
```
###### Step 5: Activate environment and run
``` bash
conda activate prokka_env
./run_prokka_spades.sh
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

Spades function → CSV
``` bash
#!/bin/bash

BASE_DIR="./prokka_results/spades"

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

GET_HOMOLOGUES

Prepare a GET_HOMOLOGUES workspace
``` bash
mkdir -p gethomologues_results
cd gethomologues_results
``` 
Copy all .faa files into the workspace:
``` bash
for f in ../prokka_results/shovill/*/*.faa; do
    cp "$f" .
done
``` 
Verify the files:
``` bash
ls -lh *.faa
``` 
Instead of letting GET_HOMOLOGUES automatically run BLAST, do it manually to optimize RAM usage:
``` bash
makeblastdb -in *.faa -dbtype prot
blastp -query *.faa -db *.faa -outfmt 6 -evalue 1e-5 -num_threads 32 -out all_vs_all.blast
``` 

Step — Cluster homologous genes with COGtriangles
``` bash
conda activate gethomologues_env
get_homologues.pl -d ./ -c -t 32 -n 2
``` 

Step 4 — Generate pan-genome matrix
``` bash
compare_clusters.pl -d ./ -o pan_genome_matrix
``` 
Step 5 — Parse matrix and classify genes
``` bash
parse_pangenome_matrix.pl -m pan_genome_matrix_cluster.tab -o summary
```
Step 6 — Plot pan-genome and core-genome curves
``` bash
plot_pancore_matrix.pl -m pan_genome_matrix_cluster.tab -o plots -p 75 -t 1000
```
Step 7 — Visualization in R
Use the R scripts generated by parse_pangenome_matrix.pl -r for cluster size distribution.

Plot:

Core genome decay curve

Pan-genome growth curve

Shell/cloud gene distributions

With 128 GB RAM, you can load all  genomes’ matrices into R simultaneously for plotting without memory issues.



create the script.

``` bash
nano run_gethomologues.sh
``` 
Paste the following:
``` bash
#!/bin/bash
set -euo pipefail

WORKDIR=$(pwd)
PROKKA_DIR="$WORKDIR/prokka_results/shovill"
THREADS=32
MIN_CLUSTER_SIZE=2
PAN_OUTPUT="pan_genome_matrix"
PLOT_OUTPUT="plots"

mkdir -p "$WORKDIR/gethomologues_results"
cd "$WORKDIR/gethomologues_results"

# Copy and rename FAA files
for f in $PROKKA_DIR/*/*.faa; do
    sample=$(basename "$(dirname "$f")")
    cp "$f" "${sample}.faa"
done

# Run GET_HOMOLOGUES with 90% identity and 75% coverage (MCL clustering)
get_homologues.pl -d ./ -M -S 90 -C 75 -t $THREADS -n $MIN_CLUSTER_SIZE

# Build pan-genome matrix
compare_clusters.pl -d ./ -o $PAN_OUTPUT

# Categorize core, soft-core, shell, cloud genes
parse_pangenome_matrix.pl -m ${PAN_OUTPUT}_cluster.tab -o summary

# Plot core/pangenome curves with Tettelin model
plot_pancore_matrix.pl -m ${PAN_OUTPUT}_cluster.tab \
    -o $PLOT_OUTPUT \
    -p $(ls $PROKKA_DIR | wc -l) \
    -t 1000 -r

echo "✅ Pan-genome construction complete"
echo "📊 Results: summary/ (core/soft-core/shell/cloud)"
echo "📈 Plots: $PLOT_OUTPUT/"

``` 
Save and exit

Ctrl+O → Enter → Ctrl+X
Make it executable and run
``` bash
conda activate gethomologues_env
chmod +x run_gethomologues.sh
./run_gethomologues.sh
``` 

``` bash
nano run_pangenome_analysis.R
``` 

``` bash
library(ggplot2)
library(reshape2)

workdir <- getwd()
results_dir <- file.path(workdir, "gethomologues_results")
tab_file <- list.files(results_dir, pattern = "pan_genome_matrix.*\\.tab$", full.names = TRUE)

if(length(tab_file) == 0){
  stop("No pan-genome matrix .tab file found in gethomologues_results/")
}

pangenome <- read.table(tab_file[1], header = TRUE, sep = "\t", check.names = FALSE)
pangenome_long <- melt(pangenome, id.vars = "Cluster")

presence_absence <- pangenome_long
presence_absence$value[presence_absence$value != "" & presence_absence$value != "0"] <- 1
presence_absence$value[presence_absence$value == "" | presence_absence$value == "0"] <- 0
presence_absence$value <- as.numeric(presence_absence$value)

cluster_summary <- aggregate(value ~ Cluster, data = presence_absence, sum)
num_genomes <- ncol(pangenome) - 1
cluster_summary$type <- cut(cluster_summary$value,
                            breaks = c(-1, 2, floor(0.95*num_genomes)-1, num_genomes-1, num_genomes),
                            labels = c("Cloud", "Shell", "Soft-core", "Core"))

write.table(cluster_summary, file = file.path(results_dir, "cluster_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

pdf(file.path(results_dir, "pan_core_plot.pdf"), width = 8, height = 6)
ggplot(cluster_summary, aes(x = type, fill = type)) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..), vjust=-0.5) +
  labs(title = "Pan-Core Genome Plot", x = "Gene Type", y = "Number of Clusters") +
  theme_minimal() +
  scale_fill_manual(values = c("Cloud"="#FF9999","Shell"="#FFCC66","Soft-core"="#66CC99","Core"="#3399FF"))
dev.off()
```
``` bash
chmod +x run_pangenome_analysis.R
``` 

``` bash
conda activate r_env
Rscript run_pangenome_analysis.R
``` 

``` bash
nano run_summarize_clusters.R
``` 

Step 3 – Summarize clusters and assign types (summarize_clusters.R)
``` bash
library(reshape2)

workdir <- getwd()
results_dir <- file.path(workdir, "gethomologues_results")
tab_file <- list.files(results_dir, pattern = "pan_genome_matrix.*\\.tab$", full.names = TRUE)

if(length(tab_file) == 0){
  stop("No pan-genome matrix .tab file found in gethomologues_results/")
}

pangenome <- read.table(tab_file[1], header = TRUE, sep = "\t", check.names = FALSE)
pangenome_long <- melt(pangenome, id.vars = "Cluster")

presence_absence <- pangenome_long
presence_absence$value[presence_absence$value != "" & presence_absence$value != "0"] <- 1
presence_absence$value[presence_absence$value == "" | presence_absence$value == "0"] <- 0
presence_absence$value <- as.numeric(presence_absence$value)

cluster_summary <- aggregate(value ~ Cluster, data = presence_absence, sum)
num_genomes <- ncol(pangenome) - 1
cluster_summary$type <- cut(cluster_summary$value,
                            breaks = c(-1, 2, floor(0.95*num_genomes)-1, num_genomes-1, num_genomes),
                            labels = c("Cloud", "Shell", "Soft-core", "Core"))

write.table(cluster_summary, file = file.path(results_dir, "cluster_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
``` 
``` bash
conda activate r_env
Rscript run_summarize_clusters.R
``` 


Step 4 – Plot cluster distribution (plot_clusters.R)
``` bash
nano run_plot_clusters.R
``` 

``` bash
library(ggplot2)

workdir <- getwd()
results_dir <- file.path(workdir, "gethomologues_results")
cluster_file <- file.path(results_dir, "cluster_summary.tsv")

if(!file.exists(cluster_file)){
  stop("cluster_summary.tsv not found. Run summarize_clusters.R first.")
}

cluster_summary <- read.table(cluster_file, header = TRUE, sep = "\t", check.names = FALSE)

pdf(file.path(results_dir, "pan_core_plot.pdf"), width = 8, height = 6)
ggplot(cluster_summary, aes(x = type, fill = type)) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..), vjust=-0.5) +
  labs(title = "Pan-Core Genome Plot", x = "Gene Type", y = "Number of Clusters") +
  theme_minimal() +
  scale_fill_manual(values = c("Cloud"="#FF9999","Shell"="#FFCC66","Soft-core"="#66CC99","Core"="#3399FF"))
dev.off()

# -------------------------
# Save Tettelin pan/core curve data as CSV
# -------------------------
plots_dir <- file.path(results_dir, "plots")
curve_file <- list.files(plots_dir, pattern = "_pancore_curve.txt$", full.names = TRUE)[1]

if(length(curve_file) > 0){
  curve_data <- read.table(curve_file, header = TRUE, sep = "\t")
  write.csv(curve_data, file=file.path(plots_dir, "pan_core_curve_data.csv"), row.names=FALSE)
}
``` 
Step 5: pan-genome openness/closeness numerically so you have both the plots and values for reporting

Using GET_HOMOLOGUES output

plot_pancore_matrix.pl generates:

Core gene decay curve (exponential decay)

Pan-genome growth curve (exponential growth)

These are usually saved in a PDF, but the script also outputs the underlying curve parameters in .txt files inside $PLOT_OUTPUT.

You can extract the parameters (Tettelin alpha for core decay, gamma for pan-genome growth) with this approach:

``` bash
NUM_GENOMES=$(ls *.faa | wc -l)
plot_pancore_matrix.pl -m ${PAN_OUTPUT}_cluster.tab \
    -o $PLOT_OUTPUT \
    -p $NUM_GENOMES \
    -t 1000 -r

# Extract Tettelin model parameters
grep -i "alpha\|gamma" $PLOT_OUTPUT/*_pancore_curve.txt > $PLOT_OUTPUT/tettelin_parameters.txt
echo "Tettelin parameters saved to $PLOT_OUTPUT/tettelin_parameters.txt"
``` 
Optional: R script to read and visualize curves

If you want the curves numerically and plotted in R:

``` 
library(ggplot2)

results_dir <- file.path(getwd(), "gethomologues_results", "plots")
curve_file <- list.files(results_dir, pattern = "_pancore_curve.txt$", full.names = TRUE)[1]

curve_data <- read.table(curve_file, header = TRUE, sep = "\t")
# Usually columns: NumGenomes, Core, PanGenome, MeanCore, MeanPan

pdf(file.path(results_dir, "pan_core_growth.pdf"), width = 8, height = 6)
ggplot(curve_data, aes(x = NumGenomes)) +
  geom_line(aes(y = Core, color="Core")) +
  geom_line(aes(y = PanGenome, color="Pan-genome")) +
  labs(title="Pan-genome and Core Genome Curves", x="Number of Genomes", y="Genes") +
  scale_color_manual(values=c("Core"="blue","Pan-genome"="red")) +
  theme_minimal()
dev.off()
``` 


