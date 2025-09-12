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
mkdir -p "$OUTDIR"

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

  echo "==> Running Shovill on: $sample"
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
    --cpus 4 \
    --ram 16 \
    --tmpdir "${TMPDIR:-/tmp}" \
    --force

  echo "==> Renaming output files in: $sample_out"
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


### 1. Quick assembly stats using `stats.sh`

The `stats.sh` script from BBMap provides basic assembly statistics such as total length, N50, number of contigs, and GC content.
###### Activate the environment
```bash
conda activate bbmap_env
```
###### Run stats.sh to Verify Installation
```bash
stats.sh
```

######  run `stats.sh`
```bash
stats.sh in=./shovill_results/ET1135_S12/ET1135_S12_contigs.fa
```
> **Tip:** This command will display key statistics including:
> 
> - Total bases
> - Number of contigs
> - Minimum, maximum, and N50 contig lengths
> - GC content

We run `stats.sh` on all Shovill assemblies in your shovill_results directory and save the summary results into a CSV file

loop over multiple assemblies:
```bash
for f in ./shovill_results/*/*_contigs.fa; do
    stats.sh in="$f"
done
```
script that will loop over all contig FASTA files, run stats.sh from BBMap, and save the results to a CSV file:
##### Step 1: Create or edit the script
```bash
nano run_assembly_stats.sh
```
##### Step 2: Paste the following into the script
```bash
#!/bin/bash

CONTIG_DIR="./shovill_results"
OUTPUT_CSV="assembly_stats.csv"

echo "Sample,Total_Bases,Num_Contigs,Min_Contig,Max_Contig,N50,GC_Content" > "$OUTPUT_CSV"

for f in "$CONTIG_DIR"/*/*_contigs.fa; do
    sample=$(basename "$f" _contigs.fa)
    
    stats_output=$(stats.sh in="$f" format=tsv 2>/dev/null | tail -n 1)
    
    total=$(echo "$stats_output" | cut -f1)
    num=$(echo "$stats_output" | cut -f3)
    min=$(echo "$stats_output" | cut -f4)
    max=$(echo "$stats_output" | cut -f5)
    n50=$(echo "$stats_output" | cut -f6)
    gc=$(echo "$stats_output" | cut -f8)
    
    echo "$sample,$total,$num,$min,$max,$n50,$gc" >> "$OUTPUT_CSV"
done

echo "✅ Assembly stats saved to $OUTPUT_CSV"
```

##### Step 3: Save and exit nano
Press Ctrl + O → Enter (to write the file)
Press Ctrl + X → Exit nano

###### Step 4: Make the script executable
``` bash
chmod +x run_assembly_stats.sh
```
###### Step 5: Activate environment and run
``` bash
conda activate bbmap_env
./run_assembly_stats.sh
```

### 2. Using seqkit to explore assembly
###### Activate the environment
``` bash
conda activate seqkit_env
```
###### Display help
``` bash
seqkit -h
```
###### Convert FASTA to tab-delimited table (sequence length and name)
``` bash
seqkit fx2tab -nl ./shovill_results/ET1135_S12/ET1135_S12_contigs.fa
```
### run QUAST on all Shovill assemblies
Collect the key statistics in a single CSV file
##### Step 1: Create the script
```bash
nano run_seqkit_on_shovill.sh
```
#####  Step 2: Paste the following into `run_seqkit_on_shovill.sh`

``` bash
#!/bin/bash
set -euo pipefail

# Directories
SHOVILL_DIR="shovill_results"
QUAST_DIR="quast_results"
mkdir -p "$QUAST_DIR"

# CSV output
CSV_FILE="quast_summary.csv"
echo "Sample,NumContigs,TotalLength,MinLen,MaxLen,AverageLen,N50,N75,GC%" > "$CSV_FILE"

# Loop over all samples
for sample_out in "$SHOVILL_DIR"/*; do
  [[ -d "$sample_out" ]] || continue

  sample=$(basename "$sample_out")
  contigs=$(ls "$sample_out"/*_contigs.fa 2>/dev/null | head -n 1)
  if [[ -z "$contigs" ]]; then
    echo ">> Skipping $sample (no contigs found)" >&2
    continue
  fi

  # Output folder for QUAST
  outdir="$QUAST_DIR/$sample"
  mkdir -p "$outdir"

  echo "==> Running QUAST on sample: $sample"
  quast "$contigs" -o "$outdir" --csv

  # Extract statistics from QUAST CSV
  stats_file="$outdir/report.tsv"
  if [[ -f "$stats_file" ]]; then
    # QUAST tsv has header, take the second line
    stats=$(sed -n '2p' "$stats_file" | tr '\t' ',')
    echo "${sample},${stats}" >> "$CSV_FILE"
  else
    echo ">> Warning: QUAST report missing for $sample" >&2
  fi
done

echo "✅ All QUAST stats saved in $CSV_FILE"
```
<details>
<summary>📖 Explanation of Assembly Statistics Script</summary>

- `CONTIG_DIR="./shovill_results"` → directory containing contig FASTA files.  
- `OUTPUT_CSV="assembly_stats.csv"` → CSV file to save assembly statistics.  
- `echo "Sample,Total_Bases,...,GC_Content" > "$OUTPUT_CSV"` → writes CSV header.  

**Loop over contigs:**  
- `for f in "$CONTIG_DIR"/*/*_contigs.fa; do ... done` → loops over all contig FASTA files in subdirectories.  
- `sample=$(basename "$f" _contigs.fa)` → extracts sample name from filename.  
- `stats_output=$(stats.sh in="$f" format=tsv 2>/dev/null | tail -n 1)` → runs `stats.sh` to get assembly metrics in TSV format and takes the last line.  

**Extract metrics:**  
- `total=$(echo "$stats_output" | cut -f1)` → total bases.  
- `num=$(echo "$stats_output" | cut -f3)` → number of contigs.  
- `min=$(echo "$stats_output" | cut -f4)` → minimum contig length.  
- `max=$(echo "$stats_output" | cut -f5)` → maximum contig length.  
- `n50=$(echo "$stats_output" | cut -f6)` → N50 statistic.  
- `gc=$(echo "$stats_output" | cut -f8)` → GC content percentage.  

**Append to CSV:**  
- `echo "$sample,$total,$num,$min,$max,$n50,$gc" >> "$OUTPUT_CSV"` → adds the stats as a new row.  
- `echo "✅ Assembly stats saved to $OUTPUT_CSV"` → prints completion message.

</details>

##### Step 3: Save and exit nano
Press Ctrl + O → Enter (to write the file)
Press Ctrl + X → Exit nano

##### Step 4: Make the script executable
```bash
chmod +x run_seqkit_on_shovill.sh
```
##### Step 5: Activate environment and run
```bash
conda activate tbprofiler_env
./run_seqkit_on_shovill.sh
```

#### 3. Assembly summary with assembly-scan
We can use another tool assembly-scan to generate summary statistics of the assembly.
###### Activate the environment
``` bash
conda activate assembly_scan_env
```
###### Verify the Installation
``` bash
assembly-scan --version
```
######  Generate summary statistics
``` bash
assembly-scan ./shovill_results/ET1135_S12/ET1135_S12_contigs.fa \
  --transpose \
  | tee ./shovill_results/ET1135_S12/ET1135_S12-assembly-scan.tsv
```
##### 4. Compute GC content from assembly-scan output
``` bash
grep 'contig_percent_[cg]' \
  ./shovill_results/ET1135_S12/ET1135_S12-assembly-scan.tsv \
  | awk -F '\t' '{sum+=$3} END {print "GC%=",sum}'
```

#1️⃣4️⃣ Spades

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
mkdir -p "$OUTDIR"

THREADS=4
MEMORY=16  # in GB
MIN_CONTIG=500
COV_CUTOFF=30  # minimum coverage for contigs

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

  echo "==> Running SPAdes on: $sample"
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

  # Filter contigs by minimum length
  if [[ -f "$sample_out/contigs.fasta" ]]; then
    awk -v minlen="$MIN_CONTIG" 'BEGIN{RS=">"; ORS=""} length($0)>minlen+1 {print ">"$0}' "$sample_out/contigs.fasta" > "$sample_out/${sample}_contigs.fasta"
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

#1️⃣4️⃣ Prokka
Prokka is a rapid **prokaryotic genome annotation tool** that predicts genes, coding sequences (CDS), rRNAs, tRNAs, and other genomic features from assembled contigs or genomes.  

Key points for TB genomes:

- Annotates **Mycobacterium tuberculosis** genomes with correct taxonomy using `--genus` and `--species`.  
- Produces multiple output files, including **GFF3**, **FASTA of proteins**, and **GenBank format**, which are useful for downstream analysis.  
- Supports **multi-threading** (`--cpus`) to speed up processing of multiple genomes.  
- Works seamlessly with **Shovill-assembled contigs**.  
- Output files are organized per sample directory with a consistent naming prefix for easy pipeline integration.  

> ⚠ Note: Prokka relies on the quality of the assembly; fragmented or low-coverage assemblies may result in incomplete annotations.
##### Step 1: Create or edit the script
```bash
nano run_prokka.sh
```
##### Step 2: Paste the following into the script

```bash
#!/bin/bash
set -euo pipefail

SHOVILL_DIR="shovill_results"
PROKKA_DIR="prokka_results"
mkdir -p "$PROKKA_DIR"

for sample_out in "$SHOVILL_DIR"/*; do
  [[ -d "$sample_out" ]] || continue

  sample=$(basename "$sample_out")

  contigs=$(ls "$sample_out"/*_contigs.fa 2>/dev/null | head -n 1)
  if [[ -z "$contigs" ]]; then
    echo ">> Skipping $sample (no contigs.fa found)" >&2
    continue
  fi

  echo "==> Running Prokka on sample: $sample"

  outdir="$PROKKA_DIR/$sample"
  mkdir -p "$outdir"

  prokka \
    --outdir "$outdir" \
    --prefix "$sample" \
    --kingdom Bacteria \
    --genus Mycobacterium \
    --species tuberculosis \
    --cpus 4 \
    --force \
    "$contigs"
done


```
<details><summary>🧬 Prokka Pipeline Overview (Click to Expand)</summary>

**📂 Input:** `SHOVILL_DIR="shovill_results"` → Shovill assemblies  
**📁 Output:** `PROKKA_DIR="prokka_results"` → Prokka annotations  

**⚙️ Steps:**  
- `mkdir -p "$PROKKA_DIR"` → ensures output directory exists  
- `for sample_out in "$SHOVILL_DIR"/*; do ... done` → loop over each sample folder  
- `[[ -d "$sample_out" ]] || continue` → skip non-directories  
- `sample=$(basename "$sample_out")` → extract sample name  
- `contigs=$(ls "$sample_out"/*_contigs.fa 2>/dev/null | head -n 1)` → locate contigs FASTA  
- `if [[ -z "$contigs" ]]; then ... fi` → skip sample if no contigs found  
- `outdir="$PROKKA_DIR/$sample"` → define sample-specific output folder  
- `mkdir -p "$outdir"` → ensure folder exists  
- `prokka --outdir "$outdir" --prefix "$sample" --kingdom Bacteria --genus Mycobacterium --species tuberculosis --cpus 4 --force "$contigs"` → run Prokka with TB-specific annotation, **overwriting old results if present**  

**✅ Result:** Each sample gets a fully annotated genome folder ready for downstream analysis, with old results replaced automatically when re-running.  

</details>

##### Step 3: Save and exit nano
Press Ctrl + O → Enter (to write the file)
Press Ctrl + X → Exit nano


###### Step 4: Make the script executable
``` bash
chmod +x run_prokka.sh
```
###### Step 5: Activate environment and run
``` bash
conda activate prokka_env
./run_prokka.sh
```

