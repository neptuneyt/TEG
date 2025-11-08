# Mapping the microbiome on the Tibetan Plateau for bioprospecting
## Overview
# Metagenomic quality control, assembly, binning and downstream analyses — Code examples

Below are reproducible code examples (command-line, R scripts and small Python helpers) translated to use English comments. Each section includes key commands, example parameters and notes to help you run the pipeline. I kept the same tools and parameters you provided. Confirm software is installed and available in PATH (or run inside appropriate containers/conda environments).

---

## 1. Raw read quality control (Trim Galore)
```bash
# Example: trim adapters and low-quality bases from paired-end reads,
# trim 3' end bases with Phred < 30 and keep paired reads >= 100 bp.
trim_galore --paired \
  --quality 30 \            # Trim bases with Q < 30
  --length 100 \            # Keep reads >= 100 bp after trimming
  --fastqc \                # Optional: run FastQC for QC reports
  -o trimmed_reads \
  sample_R1.fastq.gz sample_R2.fastq.gz
```

---

## 2. Assembly (MEGAHIT)
```bash
# Assemble cleaned paired-end reads with MEGAHIT using specified k-list
# and retain contigs > 500 bp.
megahit -1 trimmed_reads/sample_R1_val_1.fq.gz -2 trimmed_reads/sample_R2_val_2.fq.gz \
  -o megahit_out \
  --k-list 21,29,39,59,79,99,119,141 \
  --min-contig-len 500 \
  -t 32
# Output: megahit_out/final.contigs.fa
```

---

## 3. Metagenome binning (MetaWRAP + VAMB + SemiBin2) and refinement (MAGScoT)
MetaWRAP binning (runs multiple binning tools and produces consolidated bins):
```bash
# Use metaWRAP binning to run multiple binning tools and produce bin sets.
metaWRAP binning -o metaWRAP_binning_out -t 32 -a megahit_out/final.contigs.fa \
  samples/sample_R1.fastq.gz samples/sample_R2.fastq.gz
```

Run MetaWRAP with explicit MetaBAT2 and MaxBin2:
```bash
# Run only MetaBAT2 and MaxBin2 via metaWRAP
metaWRAP binning -o metaWRAP_binning_out -t 32 -a megahit_out/final.contigs.fa \
  --metabat2 --maxbin2 \
  samples/sample_R1.fastq.gz samples/sample_R2.fastq.gz
```

VAMB (deep-learning based binning) example:
```bash
# Index contigs and map reads to get BAM, then create coverage and run vamb
bwa index megahit_out/final.contigs.fa
bwa mem -t 16 megahit_out/final.contigs.fa sample_R1.fastq.gz sample_R2.fastq.gz | \
  samtools view -bS - | samtools sort -o sample.sorted.bam
samtools index sample.sorted.bam

# Generate depth table (using jgi_summarize_bam_contig_depths or coverm)
jgi_summarize_bam_contig_depths --outputDepth depth.txt sample.sorted.bam

# Run VAMB
vamb --outdir vamb_out --fasta megahit_out/final.contigs.fa --bamfiles sample.sorted.bam \
  --minfasta 200000 --nbit 256 -t 32
```

SemiBin2 example:
```bash
# Run SemiBin2 in automatic mode using contigs and BAM
semibin auto --contigs megahit_out/final.contigs.fa --bamfiles sample.sorted.bam \
  --output semibin_out --threads 32
```

MAGScoT scoring and refinement:
```bash
# Score and refine bins with MAGScoT
magscot score --bins-dir combined_bins --out magscot_scores.tsv --threads 16
magscot refine --scores magscot_scores.tsv --bins-dir combined_bins --out refined_bins
```

---

## 4. Genome quality assessment with CheckM2
```bash
# Predict completeness and contamination for bins with CheckM2
checkm2 predict --input_dir refined_bins --output checkm2_results.tsv --threads 32
# Filter MAGs that meet medium quality or above according to MIMAG:
# e.g., completeness >= 50% and contamination <= 10%
awk -F"\t" 'NR==1 || ($3>=50 && $4<=10){print $0}' checkm2_results.tsv > mags_medium_plus.tsv
```

---

## 5. Dereplication and species-level clustering with dRep
Strain-level dereplication (99.99% ANI):
```bash
# dRep dereplicate using strict ANI for strain-level clustering
dRep dereplicate drep_strain_out -g genomes_dir/*.fa \
  -p 32 -comp 50 -con 10 -sa 0.9999 --multiround_primary_clustering
# Output: dereplicated genomes in drep_strain_out
```

Species-level clustering (95% ANI):
```bash
# Cluster representatives to species-level mOTUs (95% ANI)
dRep dereplicate drep_species_out -g drep_strain_out/dereplicated_genomes/*.fa \
  -p 32 -comp 50 -con 10 -sa 0.95
# Output: species-level representative genomes
```

---

## 6. Compare mOTUs to public genomes using dRep
```bash
# Compare your mOTUs to a public reference set (GTDB/SMAG/GEM)
dRep compare --genomes1 motus/*.fa --genomes2 refs/gtdb_r220/*.fa -o compare_motus_gtdb -p 16
# Mark mOTU as novel if the best ANI to all references < 95%
```

---

## 7. Taxonomic annotation and phylogenetic tree (GTDB-Tk + IQ-TREE)
```bash
# Classify genomes with GTDB-Tk using GTDB release R220
gtdbtk classify_wf --genome_dir motus/ --out_dir gtdbtk_out --cpus 32 --db_dir /path/to/gtdb_r220

# Build concatenated-marker phylogeny using GTDB-Tk alignments
ALIGN=gtdbtk_out/marker_alignments/concat.aln

# Trim alignment gaps with trimAl
trimal -in $ALIGN -out concat.trim.aln -gt 0.1

# Build tree with IQ-TREE using LG+G and 1000 bootstraps
iqtree2 -s concat.trim.aln -m LG+G -B 1000 -wbtl -nt AUTO
# Use resulting tree in iTOL for visualization
```

---

## 8. BGC prediction and classification (antiSMASH + BiG-SCAPE)
```bash
# Run antiSMASH on contigs >= 5kb
antismash --output-dir antismash_out --taxon bacteria --minimal --cpus 16 \
  --genefinding-tool prodigal megahit_out/final.contigs.fa

# Group BGCs into GCFs/GCCs with BiG-SCAPE
bigscape.py -i antismash_out/ -o bigscape_out -c 16 --cutoffs 0.3,0.7
# Clusters are annotated (Terpene, NRPS, RiPP, PKS, hybrids, etc.)
```

---

## 9. BGC novelty assessment (BiG-SLiCE with MIBiG and BGCAtlas references)
```bash
# Prepare BiG-SLiCE database with MIBiG GenBank files beforehand, then query mode 2
bigslice prepare --mibig_dir mibig_gb --out bgc_reference_db
bigslice query --db bgc_reference_db --input antismash_out/bgc_files/ --mode 2 --out bigslice_query.tsv

# Mark BGCs as novel if membership score > 0.4 (L2-normalization-based cosine-like distance)
python3 - <<'PY'
import csv
novel=[]
with open('bigslice_query.tsv') as f:
    reader=csv.DictReader(f,delimiter='\t')
    for r in reader:
        score=float(r['membership_score'])
        if score>0.4:
            novel.append(r['bgc_id'])
print("Novel BGC count:",len(novel))
PY
```

---

## 10. Correlation of BGC density vs corrected genome size (R)
```r
# R script: bgc_density_vs_genome_size.R
library(tidyverse)
library(ggpubr)
library(ggpmisc)

# Read metadata: observed_size (bp), completeness (%), contamination (%), bgc_count, phylum
df <- read_tsv("genomes_metadata.tsv")

# Correct genome size and calculate BGC density (BGCs per bp)
df <- df %>%
  mutate(corrected_genome_size = (observed_size * 100 / completeness) - (observed_size * contamination / 100),
         bgc_density = bgc_count / corrected_genome_size)

# Keep phyla with >= 50 genomes
df_sub <- df %>% group_by(phylum) %>% filter(n() >= 50) %>% ungroup()

# Plot BGC density vs corrected genome size with LOESS smoothing and Spearman correlation
p <- ggplot(df_sub, aes(x = corrected_genome_size/1e6, y = bgc_density)) +
  geom_point(alpha=0.4) +
  stat_smooth(method="loess", se=TRUE) +
  facet_wrap(~phylum, scales="free") +
  stat_cor(method = "spearman")

ggsave("bgc_density_vs_genome_size.png", p, width=12, height=8)
```

Run:
```bash
Rscript bgc_density_vs_genome_size.R
```

---

## 11. cAMP prediction and clustering (Macrel + Prodigal smORF + CD-HIT with reduced alphabet)
Run Macrel contigs mode:
```bash
# Predict candidate antimicrobial peptides (cAMPs) from genomes
macrel contigs -i motus_representatives.fna -o macrel_out -t 16
# macrel_out includes predicted peptide fasta files
```

Predict small ORFs using a custom Prodigal that allows short ORFs (example):
```bash
# Prodigal meta-mode with small ORF support; replace with custom prodigal if required
prodigal -i motus_representatives.fna -a prodigal_smORFs.faa -p meta -m -g 11
```

Map sequences to reduced 8-letter alphabet then cluster with CD-HIT:
```bash
# Convert to reduced alphabet before clustering
python3 reduce_alphabet.py prodigal_smORFs.faa prodigal_smORFs_reduced.faa

# CD-HIT clustering with identity >= 0.75 and coverage >= 0.80
cd-hit -i prodigal_smORFs_reduced.faa -o cdhit_clusters -c 0.75 -s 0.80 -T 16 -M 16000
```

reduce_alphabet.py:
```python
#!/usr/bin/env python3
# reduce_alphabet.py
# Map standard amino acids into an 8-letter reduced alphabet
from Bio import SeqIO
import sys

# Group mapping:
# [LVIMC] -> A
# [AG]    -> B
# [ST]    -> C
# [FYW]   -> D
# [EDNQ]  -> E
# [KR]    -> F
# [P]     -> G
# [H]     -> H
mapping = {}
for aa in "LVIMC": mapping[aa] = "A"
for aa in "AG":   mapping[aa] = "B"
for aa in "ST":   mapping[aa] = "C"
for aa in "FYW":  mapping[aa] = "D"
for aa in "EDNQ": mapping[aa] = "E"
for aa in "KR":   mapping[aa] = "F"
mapping['P'] = "G"
mapping['H'] = "H"

def reduce_seq(s):
    # Map each residue to its reduced symbol; unknowns -> X
    return ''.join(mapping.get(ch.upper(),'X') for ch in s)

in_f, out_f = sys.argv[1], sys.argv[2]
with open(out_f, 'w') as out:
    for rec in SeqIO.parse(in_f, 'fasta'):
        rec.seq = type(rec.seq)(reduce_seq(str(rec.seq)))
        SeqIO.write(rec, out, 'fasta')
```

Select representative sequences per cluster (longest sequence wins; alphabetical tie-breaker).

Compute AMP density per genome:
- AMP density = number of AMPs / assembled base pairs (per genome). Map clusters back to genome origin to compute.

---

## 12. cAMP physicochemical properties (R + Peptides package)
```r
# R script: camp_physchem.R
library(tidyverse)
library(Peptides)
library(Biostrings)

# Read peptide fasta using Biostrings
fasta <- readAAStringSet("camp_representatives.faa")
df <- tibble(id = names(fasta), seq = as.character(fasta))

# Compute properties: molecular weight, amino acid composition, pI, net charge at pH 7,
# Boman index, instability index, hydrophobicity
amps_features <- df %>% rowwise() %>% mutate(
  mw = mw(seq),
  pI = pI(seq),
  net_charge = charge(seq, pH=7.0),
  boman = boman(seq),
  instability = instability(seq),
  hydrophobicity = hydrophobicity(seq)
) %>% ungroup()

write_tsv(amps_features, "amps_physchem.tsv")
```

Predict peptide structures with ColabFold (example command when installed locally or via container):
```bash
# Use colabfold_batch to predict peptide structures (adjust options as needed)
colabfold_batch --fasta camp_representatives.faa --outdir colabfold_out --num_relax 0
```

---

## 13. Compare TEG cAMPs to public AMP databases (MMseqs2)
```bash
# Prepare MMseqs2 databases and perform easy-search
mmseqs createdb camp_representatives.faa camp_db
mmseqs createdb public_amps_combined.faa public_db
mmseqs easy-search camp_db public_db mmseqs_search_res tmpdir --min-seq-id 0.75 -e 1e-5 -s 7.5

# Parse results to retain hits with identity >= 75% and E-value <= 1e-5
```

---

## 14. CRISPR–Cas system mining (CRISPRCasTyper, align Cas9, build tree)
```bash
# Detect Cas operons and CRISPR arrays with CRISPRCasTyper
CRISPRCasTyper -i motus_representatives.fna -o crispr_out -t 32

# Compare Cas proteins to CasPEDIA Cas9 sequences using DIAMOND or BLAST
diamond makedb --in caspedia_cas9.faa -d caspedia_db
diamond blastp -q crispr_out/cas_proteins.faa -d caspedia_db -o cas9_matches.m8 -e 1e-5 -k 5 -p 16

# Align complete Cas9 proteins with MUSCLE and trim gaps
muscle -in cas9_complete.faa -out cas9_align.fasta
trimal -in cas9_align.fasta -out cas9_align.trim.fasta -gt 0.95 -cons 50

# Build tree using IQ-TREE with automated model selection
iqtree2 -s cas9_align.trim.fasta -MFP -B 1000 -nt AUTO
```

---

## 15. Screening for plastic-degrading homologues (DIAMOND + PlasticDB)
```bash
# Create DIAMOND DB from PlasticDB verified proteins
diamond makedb --in PlasticDB_verified_proteins.faa -d plasticdb

# Search predicted proteins from MAGs against PlasticDB
diamond blastp -q proteins_from_MAGs.faa -d plasticdb -o diamond_plastic_hits.m8 -e 1e-5 -k 1 -p 16

# Filter high-confidence hits with identity >= 75% and query coverage >= 80%
awk '$3>=75 && $13>=80 {print $0}' diamond_plastic_hits.m8 > plastic_high_confidence_hits.m8
```

Normalize counts per genome (example Python aggregation):
```python
# aggregate_plastic_hits.py
# Aggregate plastic-degradation hits per MAG and normalize by genome size (bp)
import pandas as pd
hits = pd.read_csv('plastic_high_confidence_hits.m8', sep='\t', header=None)
# Extract genome id from query id (adjust parsing as needed)
hits['genome'] = hits[0].str.split('|').str[0]
counts = hits.groupby('genome').size().reset_index(name='plastic_gene_count')
sizes = pd.read_csv('genome_sizes.tsv', sep='\t')  # columns: genome, basepairs
merged = counts.merge(sizes, on='genome', how='left')
merged['plastic_genes_per_mb'] = merged['plastic_gene_count'] / (merged['basepairs']/1e6)
merged.to_csv('plastic_genes_per_genome.tsv', sep='\t', index=False)
```

Select highest-confidence candidates with >= 75% identity and >= 80% coverage (812 sequences in your reported example).

---

## 16. Identification and analysis of PETase candidates
1. Use DIAMOND results to extract candidate PET hydrolase sequences.
2. Align candidates to reference IsPETase (GAP38373.1) with MUSCLE.
3. Detect catalytic triad (Ser-Asp-His) and motif G158-X-S160-X-G162 using a script.
4. Run SignalP v6.0 to predict signal peptides and keep sequences with a signal peptide.
5. Predict structures via ESMFold and compare to engineered CaPETaseM9 (PDB: 7YME) using TM-align.
6. Predict optimal temperature and melting temperature via Seq2Topt; predict MAG growth temperature with Tome.

Motif and catalytic triad detection script:
```python
#!/usr/bin/env python3
# detect_petase_motifs.py
# Check sequences for PETase-like motif GxSxG and presence of catalytic residues S, D, H.
from Bio import SeqIO
import re
import sys

infile = sys.argv[1]
for rec in SeqIO.parse(infile, 'fasta'):
    seq = str(rec.seq)
    motif = bool(re.search(r'G.{1}S.{1}G', seq))  # matches GxSxG
    triad = all(res in seq for res in ('S','D','H'))  # crude check; refine with alignment
    if motif and triad:
        print(rec.id)
```

SignalP:
```bash
# Run SignalP to keep sequences with predicted signal peptide
signalp -fasta petase_candidates.faa -organism other -format short -gff3 -verbose -o signalp_out
```

ESMFold prediction and TM-align comparison (examples):
```bash
# Predict structures with ESMFold (adjust to your installed invocation)
python -m esmfold.run_prediction --input petase_candidates_filtered.faa --output_dir esmfold_out

# Compare predicted PDBs to reference 7YME using TM-align
TMalign pred.pdb 7YME.pdb > tmalign_report.txt
```

Predict biochemical properties:
```bash
# Predict protein optimal temperature etc. with Seq2Topt and predict MAG growth opt temp with Tome
seq2topt predict --input petase_candidates_filtered.faa --output seq2topt_out.tsv
tome predict --genomes motus_representatives.fna --out tome_out.tsv
```

Visualize 3D structures with PyMOL:
```bash
# Load predicted PDBs into PyMOL for manual inspection
pymol pred.pdb 7YME.pdb
```

---

## 17. Statistics and reproducibility (R)
- Use R v4.3.2 and RStudio for analysis.
- Plotting libraries: ggplot2, ggpubr, ggpmisc, maps for geographic visualization, ggtukey for compact letter displays.
- Save sessionInfo to capture package versions and reproducibility metadata.

Example of saving session info:
```r
# Save R session info to ensure reproducibility
writeLines(capture.output(sessionInfo()), "R_sessionInfo.txt")
```

Notes on plotting:
- Use scatter, box, violin, heatmap, and raincloud plots as needed.
- For pairwise tests: use Kruskal-Wallis + Dunn's test with Bonferroni correction.
- Overlay individual data points on violin/boxplots to display raw data distributions.
- Use ggtukey to add compact letter displays to plots for multiple comparisons at alpha = 0.05.

---

## Included helper scripts (names and brief descriptions)
- reduce_alphabet.py — Map amino acids to 8-letter reduced alphabet (used prior to CD-HIT).
- detect_petase_motifs.py — Screen candidate PETase sequences for motif and catalytic residues.
- bgc_density_vs_genome_size.R — Compute corrected genome sizes, BGC density, Spearman correlation and plot.
- camp_physchem.R — Compute peptide physicochemical properties using Peptides package.

---

If you want these scripts provided as separate downloadable files (with executable headers), or if you prefer a full Snakemake/Nextflow pipeline that wires these steps together for HPC execution, tell me which format you prefer and I will produce those files next.
