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

  # Extract metrics if available
  if [[ -f "$stats_file" ]]; then
    num_contigs=$(awk -F'\t' '$1 ~ /# contigs/ {print $2; exit}' "$stats_file")
    total_len=$(awk -F'\t' '$1 ~ /Total length/ {print $2; exit}' "$stats_file")
    min_len=$(awk -F'\t' '$1 ~ /(Shortest|Smallest) contig/ {print $2; exit}' "$stats_file")
    max_len=$(awk -F'\t' '$1 ~ /(Largest|Longest) contig/ {print $2; exit}' "$stats_file")
    avg_len=$(awk -F'\t' '$1 ~ /(Average|Mean) contig length/ {print $2; exit}' "$stats_file")
    n50=$(awk -F'\t' '$1 ~ /^N50/ {print $2; exit}' "$stats_file")
    n75=$(awk -F'\t' '$1 ~ /^N75/ {print $2; exit}' "$stats_file")
    gc=$(awk -F'\t' '$1 ~ /GC/ {print $2; exit}' "$stats_file")
  fi

  # Compute missing values from FASTA
  if [[ -z "${min_len:-}" || -z "${avg_len:-}" || -z "${n75:-}" ]]; then
    readarray -t lengths < <(awk '/^>/{if (seqlen){print seqlen}; seqlen=0; next} {seqlen+=length($0)} END{print seqlen}' "$contigs")
    if (( ${#lengths[@]} > 0 )); then
      IFS=$'\n' sorted=($(sort -nr <<<"${lengths[*]}"))
      unset IFS
      total_len_calc=$(awk '{sum+=$1} END{print sum}' <<<"${sorted[*]}")
      num_contigs_calc=${#sorted[@]}
      min_len_calc=$(awk 'END{print $1}' <<<"${sorted[*]}")
      max_len_calc=$(awk 'NR==1{print $1}' <<<"${sorted[*]}")
      avg_len_calc=$(awk -v t="$total_len_calc" -v n="$num_contigs_calc" 'BEGIN{if(n>0) print t/n; else print "NA"}')
      cum=0; n50_calc=NA; n75_calc=NA
      n50_target=$((total_len_calc / 2))
      n75_target=$((total_len_calc * 3 / 4))
      for len in "${sorted[@]}"; do
        cum=$((cum + len))
        if [[ "$n50_calc" == "NA" && $cum -ge $n50_target ]]; then n50_calc=$len; fi
        if [[ "$n75_calc" == "NA" && $cum -ge $n75_target ]]; then n75_calc=$len; fi
      done
      gc_calc=$(awk '/^>/{next}{seq=toupper($0); for(i=1;i<=length(seq);i++){b=substr(seq,i,1); if(b=="G"||b=="C")gc++}; total+=length(seq)} END{if(total>0) print (gc/total)*100; else print "NA"}' "$contigs")

      num_contigs=${num_contigs:-$num_contigs_calc}
      total_len=${total_len:-$total_len_calc}
      min_len=${min_len:-$min_len_calc}
      max_len=${max_len:-$max_len_calc}
      avg_len=${avg_len:-$avg_len_calc}
      n50=${n50:-$n50_calc}
      n75=${n75:-$n75_calc}
      gc=${gc:-$gc_calc}
    fi
  fi

  echo "$sample,${num_contigs:-NA},${total_len:-NA},${min_len:-NA},${max_len:-NA},${avg_len:-NA},${n50:-NA},${n75:-NA},${gc:-NA}" >> "$CSV_FILE"
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
CSV_OUTDIR="quast_results"
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



#### 4. estimate genome completeness/contamination for all assemblies with CheckM
######   Step 1: Create the Shovill CheckM script
``` bash
nano run_checkm_shovill.sh
``` 
######   Step 2: Paste this into the file:
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

ASSEMBLY_DIR="shovill_results"
READS_DIR="fastp_results_min_50"
OUTDIR="backmap_results"
THREADS=8
MIN_MAPPED_PERCENT=1        # Minimum percent mapped reads to consider sample valid
MAX_AVG_COVERAGE=1000       # Threshold to warn about extremely high coverage

mkdir -p "$OUTDIR"/{bams,depths,logs}

SUMMARY_FILE="$OUTDIR/shovill_assembly_summary.tsv"
echo -e "Sample\tTotal_Reads\tMapped_Reads\tPercent_Mapped\tAvg_Coverage\tMedian_Depth\tMin_Depth\tMax_Depth\tMean_Mapping_Quality\tPercent_Proper_Pairs\tCoverage_Breadth\tWarnings" > "$SUMMARY_FILE"

for asm_dir in "${ASSEMBLY_DIR}"/SRR*; do
    sample=$(basename "$asm_dir")
    asm=$(ls "$asm_dir"/*contigs*.fa* 2>/dev/null | head -n1 || true)

    warnings=""

    if [[ -z "$asm" || ! -f "$asm" ]]; then
        echo "[$sample] No assembly FASTA found." | tee -a "$OUTDIR/logs/$sample.log"
        continue
    fi

    # Flexible read file matching
    r1=$(ls "$READS_DIR/${sample}"*_R1*.fastq.gz "$READS_DIR/${sample}"*_1*.fastq.gz 2>/dev/null | head -n1 || true)
    r2=$(ls "$READS_DIR/${sample}"*_R2*.fastq.gz "$READS_DIR/${sample}"*_2*.fastq.gz 2>/dev/null | head -n1 || true)

    if [[ -z "$r1" || -z "$r2" ]]; then
        echo "[$sample] Missing read files." | tee -a "$OUTDIR/logs/$sample.log"
        continue
    fi

    echo "[$sample] Mapping reads..." | tee "$OUTDIR/logs/$sample.log"

    # Build BWA index if missing
    if [[ ! -f "${asm}.bwt" ]]; then
        bwa index "$asm" >> "$OUTDIR/logs/$sample.log" 2>&1
    fi

    # Map reads and create BAM
    bwa mem -t "$THREADS" "$asm" "$r1" "$r2" \
        | samtools view -bS - \
        | samtools sort -@ "$THREADS" -o "$OUTDIR/bams/${sample}.bam"

    samtools index "$OUTDIR/bams/${sample}.bam"

    # Compute coverage depth
    depth_file="$OUTDIR/depths/${sample}.depth"
    samtools depth -a "$OUTDIR/bams/${sample}.bam" > "$depth_file"

    # Total reads
    total_r1=$(zcat "$r1" | wc -l)
    total_r2=$(zcat "$r2" | wc -l)
    total=$(( (total_r1 + total_r2)/4 ))

    # Mapped reads
    mapped=$(samtools view -c -F 4 "$OUTDIR/bams/${sample}.bam")
    percent_mapped=$(awk -v m=$mapped -v t=$total 'BEGIN{if(t>0) printf "%.2f", (m/t)*100; else print 0}')

    # Check if mapping is too low
    if (( $(echo "$percent_mapped < $MIN_MAPPED_PERCENT" | bc -l) )); then
        warnings="Low mapping (<${MIN_MAPPED_PERCENT}%)"
        echo "[$sample] Warning: $warnings" | tee -a "$OUTDIR/logs/$sample.log"
        echo -e "$sample\t$total\t$mapped\t$percent_mapped\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t$warnings" >> "$SUMMARY_FILE"
        continue
    fi

    # Proper pairs
    proper_pairs=$(samtools view -c -f 2 "$OUTDIR/bams/${sample}.bam")
    percent_proper=$(awk -v p=$proper_pairs -v t=$total 'BEGIN{if(t>0) printf "%.2f", (p/t)*100; else print 0}')

    # Depth statistics
    avgcov=$(awk '{sum+=$3} END{if(NR>0) print sum/NR; else print 0}' "$depth_file")
    mediancov=$(awk '{print $3}' "$depth_file" | sort -n | awk '{a[NR]=$1} END{if(NR==0) print 0; else if(NR%2==1) print a[(NR+1)/2]; else print (a[NR/2]+a[NR/2+1])/2}')
    mindepth=$(awk 'NR==1 || $3<min{min=$3} END{print min+0}' "$depth_file")
    maxdepth=$(awk 'NR==1 || $3>max{max=$3} END{print max+0}' "$depth_file")

    # Average mapping quality
    mean_qual=$(samtools view "$OUTDIR/bams/${sample}.bam" | awk '{sum+=$5; n++} END{if(n>0) print sum/n; else print 0}')

    # Coverage breadth (% of assembly positions with depth >=1)
    total_positions=$(wc -l < "$depth_file")
    covered_positions=$(awk '$3>=1{count++} END{print count+0}' "$depth_file")
    breadth=$(awk -v cov=$covered_positions -v tot=$total_positions 'BEGIN{if(tot>0) printf "%.2f", (cov/tot)*100; else print 0}')

    # Check for extremely high coverage
    if (( $(echo "$avgcov > $MAX_AVG_COVERAGE" | bc -l) )); then
        warnings="${warnings} High coverage (> ${MAX_AVG_COVERAGE}×)"
        echo "[$sample] Warning: High coverage (> ${MAX_AVG_COVERAGE}×)" | tee -a "$OUTDIR/logs/$sample.log"
    fi

    # Write to summary
    echo -e "$sample\t$total\t$mapped\t$percent_mapped\t$avgcov\t$mediancov\t$mindepth\t$maxdepth\t$mean_qual\t$percent_proper\t$breadth\t$warnings" >> "$SUMMARY_FILE"

done

echo "✅ Assembly quality summary saved to: $SUMMARY_FILE"

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

SHOVILL_DIR="shovill_results"
PROKKA_DIR="prokka_results/shovill"

mkdir -p "$PROKKA_DIR"

if ! command -v prokka &>/dev/null; then
  echo "Error: Prokka not found in PATH."
  exit 1
fi

run_prokka() {
  sample_out="$1"
  sample=$(basename "$sample_out")

  contigs=("$sample_out"/*contigs*.fa*)
  contigs="${contigs[0]:-}"
  [[ -z "$contigs" ]] && return

  echo "Annotating $sample..."

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
         --centre TBLab \
         --compliant \
         --force --quiet \
         "$contigs"
}

export -f run_prokka
export PROKKA_DIR

find "$SHOVILL_DIR" -maxdepth 1 -mindepth 1 -type d | parallel -j "$(nproc)" run_prokka {}

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

