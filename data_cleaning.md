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
