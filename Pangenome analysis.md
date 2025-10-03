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
set -euo pipefail

mkdir -p mummer_results

if ! command -v gnuplot &>/dev/null; then
    echo "⚠️ gnuplot not found. PNG plots will be skipped."
fi

for SHOVILL_DIR in shovill_results/*; do
    SAMPLE_ID=$(basename "$SHOVILL_DIR")
    SHOVILL=$(ls "$SHOVILL_DIR"/${SAMPLE_ID}_contigs.fa 2>/dev/null || true)
    SPADES=$(ls spades_results/"$SAMPLE_ID"/contigs.fasta 2>/dev/null || true)

    if [[ -z "$SHOVILL" || -z "$SPADES" ]]; then
        echo "Skipping $SAMPLE_ID (missing Shovill or SPAdes assembly)"
        continue
    fi

    echo "🔍 Comparing $SAMPLE_ID (Shovill vs SPAdes)..."

    SAMPLE_OUTDIR="mummer_results/$SAMPLE_ID"
    mkdir -p "$SAMPLE_OUTDIR"
    OUT_PREFIX="$SAMPLE_OUTDIR/$SAMPLE_ID"

    nucmer --mincluster=500 --prefix="$OUT_PREFIX" "$SHOVILL" "$SPADES"
    delta-filter -1 -l 1000 "$OUT_PREFIX.delta" > "$OUT_PREFIX.filtered.delta"
    show-coords -rcl "$OUT_PREFIX.filtered.delta" > "$OUT_PREFIX.coords"
    dnadiff -p "$OUT_PREFIX" "$SHOVILL" "$SPADES"

    if command -v gnuplot &>/dev/null; then
        mummerplot --png --large --layout --color \
            -p "$OUT_PREFIX" \
            "$OUT_PREFIX.filtered.delta"

        echo "✅ Plot generated for $SAMPLE_ID: $OUT_PREFIX.png"
    else
        echo "⚠️ Plot skipped for $SAMPLE_ID (gnuplot missing)"
    fi
done

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
###### Step 6 — Dotplot Visualization

MUMmer provides mummerplot to visualize genome assembly alignments. This shows synteny, rearrangements, inversions, and gaps between your Shovill vs SPAdes assemblies.
Command for a single sample:
``` bash
mummerplot --png --large --layout --color \
  -p mummer_results/SRR32861807/SRR32861807 \
  mummer_results/SRR32861807/SRR32861807.filtered.delta
``` 

`--layout` separates forward and reverse matches into two distinct panels, making inversions and rearrangements easier to spot.

Still lightweight; no extra software needed beyond gnuplot.

``` bash
mummerplot --png --large --layout -p mummer_results/SRR32861807 mummer_results/SRR32861807.filtered.delta
``` 
Zoom in on regions

If you’re interested in specific regions:
``` bash
mummerplot --png --large --xrange=100000-500000 --yrange=100000-500000 -p zoom_plot mummer_results/SRR32861807.filtered.delta
``` 
Use external visualization tools
dnadiff + mummerplot

Run dnadiff on your assemblies:
``` bash
dnadiff reference.fasta query.fasta -p dnadiff_results
``` 

Loop to generate dotplots for all samples

If you want dotplots for all your samples automatically:
``` bash
for DELTA in mummer_results/*.filtered.delta; do
    SAMPLE=$(basename "$DELTA" .filtered.delta)
    echo "Generating dotplot for $SAMPLE..."
    mummerplot --png --large -p mummer_results/"$SAMPLE" "$DELTA"
done
``` 


comparison script
###### Step 1 — Open a new script
``` bash
nano run_quast_comparison.sh
``` 

###### Step 2 — Paste this clean script
``` bash
#!/bin/bash
set -euo pipefail

mkdir -p mummer_results
mkdir -p csv_output

TMP_CSV="csv_output/mummer_alignment_comparison_unsorted.csv"
FINAL_CSV="csv_output/mummer_alignment_comparison.csv"

echo "Sample,Shovill_Length,SPAdes_Length,Total_Aligned_Bases,Avg_Identity,Num_Aligned_Blocks,Length_Difference,Percent_Sh_vs_SP,Percent_SP_vs_Sh,Percent_Length_Diff" > "$TMP_CSV"

for SAMPLE_ID in $(ls shovill_results); do
    SHOVILL="shovill_results/${SAMPLE_ID}/${SAMPLE_ID}_contigs.fa"
    SPADES="spades_results/${SAMPLE_ID}/contigs.fasta"

    if [[ ! -f "$SHOVILL" || ! -f "$SPADES" ]]; then
        echo "Skipping $SAMPLE_ID (missing assembly)"
        continue
    fi

    SAMPLE_OUT="mummer_results/${SAMPLE_ID}"
    mkdir -p "$SAMPLE_OUT"

    nucmer --maxmatch --prefix="$SAMPLE_OUT/$SAMPLE_ID" "$SHOVILL" "$SPADES"
    dnadiff -p "$SAMPLE_OUT/$SAMPLE_ID" "$SHOVILL" "$SPADES"

    SHOVILL_LEN=$(awk '/^>/ {if(seq){total+=length(seq)}; seq=""} /^[ACGTN]+$/ {seq=seq $0} END{total+=length(seq); print total}' "$SHOVILL")
    SPADES_LEN=$(awk '/^>/ {if(seq){total+=length(seq)}; seq=""} /^[ACGTN]+$/ {seq=seq $0} END{total+=length(seq); print total}' "$SPADES")

    TOTAL_ALIGNED=$(awk '/^AlignedBases/{gsub(",","",$2); print $2}' "$SAMPLE_OUT/$SAMPLE_ID.report")
    AVG_IDENTITY=$(awk '/^AvgIdentity/{gsub("%","",$2); print $2}' "$SAMPLE_OUT/$SAMPLE_ID.report")
    NUM_BLOCKS=$(awk '/^TotalRefMatches/{gsub(",","",$2); print $2}' "$SAMPLE_OUT/$SAMPLE_ID.report")

    LENGTH_DIFF=$(( SPADES_LEN - SHOVILL_LEN ))
    PERCENT_SH_vs_SP=$(awk -v aligned="$TOTAL_ALIGNED" -v len="$SHOVILL_LEN" 'BEGIN{if(len>0) printf "%.2f", (aligned/len)*100; else print 0}')
    PERCENT_SP_vs_SH=$(awk -v aligned="$TOTAL_ALIGNED" -v len="$SPADES_LEN" 'BEGIN{if(len>0) printf "%.2f", (aligned/len)*100; else print 0}')
    PERCENT_LENGTH_DIFF=$(awk -v diff="$LENGTH_DIFF" -v len="$SHOVILL_LEN" 'BEGIN{if(len>0) printf "%.2f", (diff/len)*100; else print 0}')

    echo "${SAMPLE_ID},${SHOVILL_LEN},${SPADES_LEN},${TOTAL_ALIGNED},${AVG_IDENTITY},${NUM_BLOCKS},${LENGTH_DIFF},${PERCENT_SH_vs_SP},${PERCENT_SP_vs_SH},${PERCENT_LENGTH_DIFF}" >> "$TMP_CSV"

    mummerplot --png --large -p "$SAMPLE_OUT/$SAMPLE_ID" "$SAMPLE_OUT/$SAMPLE_ID.delta"
    rm -f "$SAMPLE_OUT/$SAMPLE_ID.gp"
done

(head -n 1 "$TMP_CSV" && tail -n +2 "$TMP_CSV" | sort -t',' -k8,8nr) > "$FINAL_CSV"
rm -f "$TMP_CSV"

``` 
###### Step 3 — Save and exit

Press CTRL + O → Enter → CTRL + X
###### Step 4 — Make the script executable
```bash
chmod +x run_quast_comparison.sh
```
###### Step 5 — Run the script
```bash
conda activate quast_env
./run_quast_comparison.sh
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




Contigutor

``` bash
nano run_shovill_contigutor.sh
``` 
``` bash
#!/bin/bash
set -euo pipefail

INPUT_DIR="fastp_results_min_50"
OUTDIR="spades_results"
THREADS=16
MEMORY=120
MIN_CONTIG=500
COV_CUTOFF=30
CONTIGUTOR_DIR="CONTIGutor_results/spades_CONTIGutor_results"
REF_GENOME="H37Rv.gbk"

mkdir -p "$CONTIGUTOR_DIR"

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
    mkdir -p "$sample_out"

    echo "==> Running SPAdes (paired-end) on: $sample"
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

    sample_contig="$sample_out/${sample}_contigs.fasta"
    if [[ -f "$sample_contig" ]]; then
        sample_cont_out="$CONTIGUTOR_DIR/$sample"
        mkdir -p "$sample_cont_out"

        echo "==> Running CONTIGutor on $sample"
        perl CONTIGutor.pl \
            -c "$sample_contig" \
            -r "$REF_GENOME" \
            -o "$sample_cont_out" \
            > "$sample_cont_out/CONTIGutor.log" 2>&1
    fi
done

echo ">> All SPAdes and CONTIGutor runs completed. Results are in $CONTIGUTOR_DIR"

``` 
``` bash
chmod +x run_shovill_contigutor.sh
```

``` bash
./run_shovill_contigutor.sh
``` 
``` bash
nano run_spades_contigutor.sh
``` 

``` bash
#!/bin/bash
set -euo pipefail

INPUT_DIR="fastp_results_min_50"
OUTDIR="spades_results"
THREADS=16
MEMORY=120
MIN_CONTIG=500
COV_CUTOFF=30
CONTIGUTOR_DIR="CONTIGutor_results/spades_CONTIGutor_results"
REF_GENOME="H37Rv.gbk"

mkdir -p "$CONTIGUTOR_DIR"

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
    else
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
    fi

    sample_contig="$sample_out/${sample}_contigs.fasta"
    if [[ -f "$sample_contig" ]]; then
        sample_cont_out="$CONTIGUTOR_DIR/$sample"
        mkdir -p "$sample_cont_out"

        if [[ -f "$sample_cont_out/CONTIGutor.log" ]]; then
            echo ">> Skipping $sample (CONTIGutor already run)"
            continue
        fi

        echo "==> Running CONTIGutor on $sample"
        perl CONTIGutor.pl \
            -c "$sample_contig" \
            -r "$REF_GENOME" \
            -o "$sample_cont_out" \
            > "$sample_cont_out/CONTIGutor.log" 2>&1
    fi
done

echo ">> All SPAdes and CONTIGutor runs completed. Results are in $CONTIGUTOR_DIR"
``` 
``` bash
chmod +x run_spades_contigutor.sh
./run_spades_contigutor.sh
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

SHOVILL_DIR="shovill_results"
PROKKA_DIR="prokka_results/shovill"

mkdir -p "$PROKKA_DIR"

if ! command -v prokka &>/dev/null; then
  echo "Prokka not found"
  exit 1
fi

run_prokka() {
  sample_out="$1"
  sample=$(basename "$sample_out")

  contigs=("$sample_out"/*_contigs.fa)
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

find "$SHOVILL_DIR" -maxdepth 1 -mindepth 1 -type d | parallel -j 8 run_prokka {}

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

SPAdes summary → CSV
``` bash
#!/bin/bash
set -euo pipefail

mkdir -p csv_output
output_file="csv_output/prokka_spades_summary.csv"
echo "Sample,CDS,rRNA,tRNA,Genome_size,GC_content" > "$output_file"

for sample_dir in prokka_results/prokka_results_spades/*; do
    sample=$(basename "$sample_dir")
    stats_file="$sample_dir/$sample.txt"
    if [[ -f "$stats_file" ]]; then
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
WORKDIR=$(pwd)
PROKKA_DIR="$WORKDIR/prokka_results/shovill"
THREADS=32
MIN_CLUSTER_SIZE=2
PAN_OUTPUT="pan_genome_matrix"
PLOT_OUTPUT="plots"

mkdir -p "$WORKDIR/gethomologues_results"
cd "$WORKDIR/gethomologues_results"

for f in $PROKKA_DIR/*/*.faa; do
    cp "$f" .
done

cat *.faa > all_genomes.faa

makeblastdb -in all_genomes.faa -dbtype prot -parse_seqids -out all_genomes_db

blastp -query all_genomes.faa -db all_genomes_db -outfmt 6 -evalue 1e-5 -num_threads $THREADS -out all_vs_all.blast

get_homologues.pl -d ./ -c -t $THREADS -n $MIN_CLUSTER_SIZE

compare_clusters.pl -d ./ -o $PAN_OUTPUT

parse_pangenome_matrix.pl -m ${PAN_OUTPUT}_cluster.tab -o summary

plot_pancore_matrix.pl -m ${PAN_OUTPUT}_cluster.tab -o $PLOT_OUTPUT -p $(ls $PROKKA_DIR | wc -l) -t 1000 -r

echo "Pan-genome analysis complete"
echo "Core/Soft-core/Shell/Cloud genes: summary/"
echo "Pan/Core genome plots: $PLOT_OUTPUT/"

``` 
Save and exit

Ctrl+O → Enter → Ctrl+X
Make it executable and run
``` bash
conda activate gethomologues_env
chmod +x run_gethomologues.sh
./run_gethomologues.sh
``` 















panoroo 

``` bash
cond activate panaroo_env
``` 


``` bash
mkdir -p panaroo_results && \
panaroo -i $(find prokka_results -name "*.gff") \
        -o panaroo_results \
        --clean-mode strict \
        --remove-invalid-genes \
        --core_threshold 99 \
        -t 8
```




Pan-genome COG Functional Analysis Workflow
Objective: Assign COG categories to core, accessory, and singleton genes from Panaroo output, and visualize their distribution.
Setup Project and Output Folder
Navigated to our project directory:
``` bash
cd /media/betselot_z/DATADRIVE0/betselot/TB/TB_project/Illumina/WGS/paired_end/PRJNA1247743
``` 

Create a folder for COG functional annotation results:
``` bash
mkdir -p cog_results
```
Install COGclassifier (if not already installed)

COGclassifier simplifies functional annotation by downloading and using the COG database automatically.
``` bash
pip install cogclassifier
``` 
Download / Prepare COG Database (optional if you want local copy)

If you prefer to manually keep a local copy of the COG database, you can create a folder:
``` bash
mkdir -p ~/databases/COG2020
cd ~/databases/COG2020
``` 


Run COGclassifier on Panaroo Representative Proteins
Input: panaroo_results/pan_genome_reference.fa
Output folder: cog_results

``` bash
#!/bin/bash
set -euo pipefail

PROJECT_DIR="/media/betselot_z/DATADRIVE0/betselot/TB/TB_project/Illumina/WGS/paired_end/PRJNA1247743"
PANAROO_PROTEINS="$PROJECT_DIR/panaroo_results/pan_genome_reference.fa"
COG_OUTPUT="$PROJECT_DIR/cog_results"
THREADS=16

mkdir -p "$COG_OUTPUT"

echo "Input file: $PANAROO_PROTEINS"
echo "Output folder: $COG_OUTPUT"
echo "Using $THREADS threads"

COGclassifier -i "$PANAROO_PROTEINS" \
              -o "$COG_OUTPUT" \
              --thread_num "$THREADS"

echo "COGclassifier finished. Key output files:"
echo " - cog_classify.tsv (per-gene COG assignments)"
echo " - cog_count.tsv (summary counts per COG category)"
echo " - rpsblast.tsv (raw RPS-BLAST results)"
echo " - cog_count_barchart.html (interactive bar chart)"
echo " - cog_count_piechart.html (interactive pie chart)"
echo "Note: PNG plots may fail on Linux. Use HTML charts or install vl-convert to fix."

``` 

