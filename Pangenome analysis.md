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

  # Remove .fastq.gz and optional .trim suffix
  sample=$(basename "$R1" .fastq.gz)
  sample=${sample%.trim}

  sample_out="$OUTDIR/$sample"

  if [[ -f "$sample_out/contigs.fasta" ]]; then
    echo ">> Skipping $sample (already assembled)"
    continue
  fi

  echo "==> Running SPAdes (single-end) on: $sample"
  mkdir -p "$sample_out"

  spades.py \
    -s "$R1" \
    -o "$sample_out" \
    -t "$THREADS" \
    -m "$MEMORY" \
    --only-assembler \
    --cov-cutoff "$COV_CUTOFF" \
    > "$sample_out/spades.log" 2>&1

  if [[ -f "$sample_out/contigs.fasta" ]]; then
    awk -v minlen="$MIN_CONTIG" 'BEGIN{RS=">"; ORS=""} length($0)>minlen+1 {print ">"$0}' \
      "$sample_out/contigs.fasta" > "$sample_out/contigs.filtered.fasta"
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



### 1. Using BBMap to explore assembly

We run `stats.sh` on all Shovill and spades assemblies in your shovill_results directory and save the summary results into a CSV file

script that will loop over all contig FASTA files, run stats.sh from BBMap, and save the results to a CSV file:

######   Step 1: Create the SPAdes CheckM script
``` bash
nano run_spades_stats.sh
``` 
Paste this into the file:

``` bash
#!/bin/bash
set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate bbmap_env

SPADES_DIR="spades_results"
TSV_OUTDIR="csv_output"
mkdir -p "$TSV_OUTDIR"

first_sample=true
OUTFILE="$TSV_OUTDIR/spades_assembly_stats.tsv"

for sample_dir in "$SPADES_DIR"/*; do
    sample=$(basename "$sample_dir")
    contig_file="$sample_dir/${sample}_contigs.fasta"

    if [[ -f "$contig_file" ]]; then
        echo "Processing $sample..."
        if $first_sample; then
            echo -e "Sample\t$(stats.sh in="$contig_file" format=3 | head -n1)" > "$OUTFILE"
            first_sample=false
        fi
        stats.sh in="$contig_file" format=3 | tail -n +2 | awk -v s="$sample" 'BEGIN{OFS="\t"} {print s,$0}' >> "$OUTFILE"
    else
        echo ">> Contig file not found for $sample, skipping."
    fi
done

echo "All SPAdes assembly stats saved to $OUTFILE"
``` 
######   Step 2:Create the Shovill BBMap script
``` bash
nano run_shovill_stats.sh
``` 
paste
``` bash
#!/bin/bash
set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate bbmap_env

SHOVILL_DIR="shovill_results"
TSV_OUTDIR="csv_output"
mkdir -p "$TSV_OUTDIR"

first_sample=true
OUTFILE="$TSV_OUTDIR/shovill_assembly_stats.tsv"

for sample_dir in "$SHOVILL_DIR"/*; do
    sample=$(basename "$sample_dir")
    contig_file="$sample_dir/${sample}_contigs.fa"

    if [[ -f "$contig_file" ]]; then
        echo "Processing $sample..."
        if $first_sample; then
            echo -e "Sample\t$(stats.sh in="$contig_file" format=3 | head -n1)" > "$OUTFILE"
            first_sample=false
        fi
        stats.sh in="$contig_file" format=3 | tail -n +2 | awk -v s="$sample" 'BEGIN{OFS="\t"} {print s,$0}' >> "$OUTFILE"
    else
        echo ">> Contig file not found for $sample, skipping."
    fi
done

echo "All Shovill assembly stats saved to $OUTFILE"
``` 
###### Step 3: Make scripts executable
``` bash
chmod +x run_spades_stats.sh run_shovill_stats.sh
``` 
###### Step 4: Run the scripts
``` bash
conda activate bbmap_env
./run_spades_stats.sh
./run_shovill_stats.sh
``` 

### 2.run QUAST on all Shovill assemblies
Collect the key statistics in a single CSV file
running QUAST on all SPAdes assemblies and collecting key statistics into a single report.tsv. 
######   Step 1: Create the SPAdes quast script
``` bash
nano run_quast_spades.sh
``` 
Paste this into the file:

``` bash
#!/bin/bash
set -euo pipefail

SPADES_DIR="spades_results"
QUAST_PARENT="quast_results"
QUAST_DIR="$QUAST_PARENT/quast_results_spades"
CSV_OUTDIR="csv_output"

mkdir -p "$QUAST_DIR" "$CSV_OUTDIR"

CSV_FILE="$CSV_OUTDIR/quast_summary_spades.csv"
echo "Sample,NumContigs,TotalLength,MinLen,MaxLen,AverageLen,N50,N75,GC%" > "$CSV_FILE"

for sample_out in "$SPADES_DIR"/*; do
  [[ -d "$sample_out" ]] || continue
  sample=$(basename "$sample_out")
  
  contigs_file=("$sample_out"/*_contigs.fasta)
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
######   Step 2:Create the Shovill assembly-scan script
``` bash
nano run_quast_shovill.sh
``` 
Paste this into the file:
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
###### Step 3: Make scripts executable
``` bash
chmod +x run_quast_spades.sh run_quast_shovill.sh
``` 
###### Step 4: Run the scripts
``` bash
conda activate quast_env
./run_quast_spades.sh
./run_quast_shovill.sh
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

BASE_DIR="./prokka_results/prokka_results_shovill"

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

BASE_DIR="./prokka_results/prokka_results_spades"

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






Prokka annotation results for both Shovill and SPAdes, we can visualize them in multiple ways. The simplest approach is to summarize key annotation metrics per sample and plot them using Python.

Here’s a clean workflow using matplotlib/seaborn.

##### Step 1: Install necessary Python packages
``` bash
conda activate your_env
conda install pandas matplotlib seaborn -y
``` 
##### Step 2: Create a Python script
``` bash
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import os

# Function to parse Prokka summary.txt
def parse_prokka_summary(prokka_dir):
    data = []
    for sample_dir in glob.glob(f"{prokka_dir}/*"):
        sample = os.path.basename(sample_dir)
        summary_file = os.path.join(sample_dir, f"{sample}.txt")
        if os.path.isfile(summary_file):
            summary = {line.split(':')[0].strip(): line.split(':')[1].strip() 
                       for line in f if ':' in line} if (f:=open(summary_file)) else {}
            summary['Sample'] = sample
            data.append(summary)
            f.close()
    return pd.DataFrame(data)

# Load Shovill and SPAdes results
shovill_df = parse_prokka_summary("prokka_results_shovill")
spades_df = parse_prokka_summary("prokka_results_spades")

shovill_df['Assembler'] = 'Shovill'
spades_df['Assembler'] = 'SPAdes'

df = pd.concat([shovill_df, spades_df], ignore_index=True)

# List of numeric columns to include
features = ['CDS','tRNAs','rRNAs','tmRNAs','Repeat_regions','pseudogenes','rRNA operons','CRISPRs']
features = [f for f in features if f in df.columns]

# Convert columns to numeric
for col in features:
    df[col] = pd.to_numeric(df[col], errors='coerce')

# Melt for plotting
df_melt = df.melt(id_vars=['Sample','Assembler'], value_vars=features, 
                  var_name='Feature', value_name='Count')

# Grouped barplot
plt.figure(figsize=(14,8))
sns.barplot(x='Sample', y='Count', hue='Feature', data=df_melt, ci=None)
plt.xticks(rotation=90)
plt.ylabel('Count')
plt.title('Prokka Annotation Features per Sample (Shovill vs SPAdes)')
plt.legend(bbox_to_anchor=(1.05,1), loc='upper left')
plt.tight_layout()
plt.savefig('csv_output/prokka_all_features_comparison.png')
plt.show()

# Optional: heatmap of features per sample and assembler
heatmap_df = df.set_index(['Sample','Assembler'])[features]
plt.figure(figsize=(12,8))
sns.heatmap(heatmap_df.fillna(0), annot=True, fmt=".0f", cmap="YlGnBu")
plt.title('Heatmap of Prokka Annotation Features')
plt.tight_layout()
plt.savefig('csv_output/prokka_heatmap_features.png')
plt.show()
``` 

##### Step 3: Run the script
``` bash
python visualize_prokka.py
```
merge this with the assembly statistics (N50, contigs, GC%) so that both assembly and annotation metrics are in one visualization

##### Step 1: Prepare the Python script
``` bash
nano visualize_combined_metrics.py
``` 
##### Step 2: Create a Python script
``` bash
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load assembly metrics (from assembly-scan CSVs)
shovill_assembly = pd.read_csv("csv_output/shovill_assembly_scan.csv")
spades_assembly = pd.read_csv("csv_output/spades_assembly_scan.csv")

shovill_assembly['Assembler'] = 'Shovill'
spades_assembly['Assembler'] = 'SPAdes'

assembly_df = pd.concat([shovill_assembly, spades_assembly], ignore_index=True)

# Aggregate assembly metrics per sample
assembly_summary = assembly_df.groupby(['Sample','Assembler']).agg(
    NumContigs=('Contig','count'),
    TotalLength=('Length','sum'),
    AvgGC=('GC%','mean')
).reset_index()

# Load Prokka annotation metrics
def parse_prokka_summary(prokka_dir):
    import glob, os
    data = []
    for sample_dir in glob.glob(f"{prokka_dir}/*"):
        sample = os.path.basename(sample_dir)
        summary_file = os.path.join(sample_dir, f"{sample}.txt")
        if os.path.isfile(summary_file):
            with open(summary_file) as f:
                summary = {line.split(':')[0].strip(): line.split(':')[1].strip() 
                           for line in f if ':' in line}
            summary['Sample'] = sample
            data.append(summary)
    return pd.DataFrame(data)

shovill_prokka = parse_prokka_summary("prokka_results_shovill")
spades_prokka = parse_prokka_summary("prokka_results_spades")

shovill_prokka['Assembler'] = 'Shovill'
spades_prokka['Assembler'] = 'SPAdes'

prokka_df = pd.concat([shovill_prokka, spades_prokka], ignore_index=True)

# List of annotation features to include
features = ['CDS','tRNAs','rRNAs','tmRNAs','Repeat_regions','pseudogenes','CRISPRs']
features = [f for f in features if f in prokka_df.columns]

for col in features:
    prokka_df[col] = pd.to_numeric(prokka_df[col], errors='coerce')

prokka_summary = prokka_df[['Sample','Assembler'] + features]

# Merge assembly and annotation metrics
combined_df = assembly_summary.merge(prokka_summary, on=['Sample','Assembler'], how='outer')

# Save combined CSV
combined_df.to_csv("csv_output/combined_assembly_annotation_metrics.csv", index=False)

# -------------------------
# Visualization
# -------------------------

# Melt for grouped barplot
melt_features = ['NumContigs','TotalLength','AvgGC'] + features
df_melt = combined_df.melt(id_vars=['Sample','Assembler'], value_vars=melt_features,
                           var_name='Feature', value_name='Value')

plt.figure(figsize=(16,8))
sns.barplot(x='Sample', y='Value', hue='Feature', data=df_melt)
plt.xticks(rotation=90)
plt.title('Combined Assembly and Annotation Metrics per Sample')
plt.tight_layout()
plt.savefig('csv_output/combined_metrics_barplot.png')
plt.show()

# Heatmap version
heatmap_df = combined_df.set_index(['Sample','Assembler'])[melt_features]
plt.figure(figsize=(14,10))
sns.heatmap(heatmap_df.fillna(0), annot=True, fmt=".1f", cmap='YlOrRd')
plt.title('Heatmap of Combined Assembly and Annotation Metrics')
plt.tight_layout()
plt.savefig('csv_output/combined_metrics_heatmap.png')
plt.show()

``` 

##### Step 3: Run the script
``` bash
python visualize_combined_metrics.py
```
