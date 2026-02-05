下面是一份可以直接保存为 `nipah_ml_pipeline.md` 的 Markdown 记录模板，涵盖你这次从 NCBI Virus 下载到 ML 树的整个流程、关键参数与实际命令（含你最终确认可用的 `-T 64` 版本）。

# Nipah virus (Henipavirus nipahense) ML phylogeny pipeline (1999–2026)

> Goal: Build a Maximum Likelihood (ML) phylogenetic tree for complete Nipah virus genomes (1999–2026).  
> Current dataset: 80 GenBank records downloaded from NCBI Virus → deduplicated to 78 unique sequences.

## 0. Environment / tools

- OS: Linux
- CPU: AMD EPYC (128 cores available)
- conda env: `(base)`  
- Tools:
  - `seqkit` (for FASTA manipulation & dedup)
  - `mafft` (multiple sequence alignment)
  - `iqtree2` v2.0.7 (ML tree inference)
  - `bioawk` (optional QC)
  - `trimal` (not installed; trimming step skipped)

> Notes:
> - `iqtree2 -T AUTO` did not effectively use multiple threads in our run; explicit `-T 64` was used for speed.
> - We initially attempted `trimal` trimming but `trimal` was not available on the system.

## 1. Data acquisition (NCBI Virus)

### 1.1 Query / filters (GUI)
NCBI Virus → Explore Virus Data → Advanced Filters for GenBank Sequences

Applied filters:
- Virus/Taxonomy: `Henipavirus nipahense (taxid:3052225)`
- Sequence Quality:
  - Nucleotide Completeness: `Complete`
  - Assembly Completeness: `Complete`
  - Sequence Length: `17000–20000`
- Dates:
  - Collection date: `1999–2026`

Result: 80 records

### 1.2 Download
Download All Results → **Sequence Data (FASTA format)** → **Nucleotide**
- Saved as: `nipah_80_nt.fasta`

(Recommended) also download metadata:
- Download All Results → **Results Table** → TSV/CSV
- Include fields such as accession, isolate/strain, country, host, collection date, length.



## 2. FASTA normalization & cleaning

### 2.1 Convert sequences to single-line FASTA
This avoids downstream tools mis-reading wrapped sequences.

```bash
seqkit seq -w 0 nipah_80_nt.fasta > nipah_80_nt.oneline.fasta
```

Sanity check:

```bash
seqkit stats nipah_80_nt.oneline.fasta
```

```bash
file                       format  type  num_seqs    sum_len  min_len   avg_len  max_len
nipah_80_nt.oneline.fasta  FASTA   DNA         80  1,454,822   18,027  18,185.3   18,252
```

### 2.2 Clean non-IUPAC characters (replace with N)
To ensure parsers/QC tools work reliably and to avoid unexpected characters.

```bash
seqkit seq -w 0 nipah_80_nt.oneline.fasta \
| awk '{
  if($0 ~ /^>/){print; next}
  gsub(/[^ACGTRYSWKMBDHVNacgtryswkmbdhvn]/,"N");
  print
}' > nipah_80_nt.clean.fasta
```

Optional check: should return nothing if clean.

```bash
grep -v "^>" nipah_80_nt.clean.fasta | grep -n "[^ACGTRYSWKMBDHVNacgtryswkmbdhvn]" | head
```

## 3. Deduplication (remove identical sequences)

We removed exact duplicate sequences at the nucleotide level.

```bash
seqkit rmdup -s nipah_80_nt.clean.fasta -o nipah_80_dedup.fasta
```

Count sequences:

```bash
grep -c "^>" nipah_80_dedup.fasta
```

Observed:
- `[INFO] 2 duplicated records removed`
- Resulting unique sequences: `78`



## 4. Multiple sequence alignment (MAFFT)

```bash
mafft --auto nipah_80_dedup.fasta > nipah_80.aln.fasta
```

Notes:
- MAFFT strategy selected by `--auto` (observed FFT-NS-2 fast strategy).



## 5. Alignment trimming (optional; skipped)

We attempted to use trimAl but it was not installed:

```bash
trimal -in nipah_80.aln.fasta -out nipah_80.aln.trim.fasta -automated1
```


## 6. Maximum Likelihood tree inference (IQ-TREE 2)

Final command (confirmed working + faster with explicit threads):

```bash
iqtree2 -s nipah_80.aln.trim.fasta -m MFP -B 1000 --alrt 1000 -T 128
```

Parameters:
- `-s nipah_80.aln.trim.fasta` : alignment
- `-m MFP` : ModelFinder Plus (automatic best-fit model selection)
- `-B 1000` : ultrafast bootstrap (1000 replicates)
- `--alrt 1000` : SH-aLRT support (1000 replicates)
- `-T 128` : use 128 CPU threads (explicitly set for speed)

Key outputs:
- `nipah_80.aln.trim.fasta.treefile` : best ML tree with support values
- `nipah_80.aln.trim.fasta.iqtree` : summary (selected model, likelihood, etc.)
- `nipah_80.aln.trim.fasta.log` : run log
- `nipah_80.aln.trim.fasta.ufboot` : bootstrap trees (UFBoot)
- `nipah_80.aln.trim.fasta.contree` : consensus tree (if generated)
