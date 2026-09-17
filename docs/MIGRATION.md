# Migration from the original scripts

| Original script or group | CHARTseq_tools command |
|---|---|
| `run_QRseq.sh` | `chartseq run` |
| `spilt_fastq.py` | `chartseq split-fastq` |
| `CHART-seqSpiltBam.py` | `chartseq split-bam` |
| `sort_tsv.py` | `chartseq pivot-counts` |
| `summarytsv.py` | `chartseq summarize-counts` |
| `CHART-seqDownsample.py` + mapping shell script | `chartseq saturation-run` |
| Downsampling summary scripts | `chartseq summarize-downsample` |
| Three saturation plotting scripts | `chartseq plot-saturation` |
| SE/PE RPKM/FPKM scripts | `chartseq normalize-featurecounts [--paired]` |
| SE/PE isoform scripts | `chartseq isoform-summary` |
| SE/PE AS event summaries | `chartseq summarize-splicing` |
| Fixed-prefix UMI removal | `chartseq trim-prefix` |

The old absolute paths are intentionally not preserved. Executables are found
on `PATH` or set in YAML; reference files are YAML values or explicit command
arguments. The old `nohup.out` file is runtime output and is not part of the
package.
