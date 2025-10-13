modify 5.fastp_filtered_reads_plot.tsv in place and remove _1 from the Sample column directly.
```bash
sed -i -E 's/^([^[:space:]]+)_1/\1/' 5.fastp_filtered_reads_plot.tsv

```
