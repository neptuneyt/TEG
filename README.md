# Mapping the microbiome on the Tibetan Plateau for bioprospecting
## Overview
# Metagenomic quality control, assembly, binning and downstream analyses — Code examples

Below are reproducible code examples (command-line, R scripts and small Python helpers). Each section includes key commands, example parameters and notes to help you run the pipeline. Confirm software is installed and available in PATH (or run inside appropriate containers/conda environments).

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
MAGScoT_folder="/path/to/MAGScoT"
cd $MAGScoT_folder/example

Rscript $MAGScoT_folder/MAGScoT.R -i example.contigs_to_bin.tsv --hmm example.hmm
```

---

## 4. Genome quality assessment with CheckM2
```bash
# Predict completeness and contamination for bins with CheckM2
# mags.path: One MAG each line
time checkm2 predict -i $(cat in.path) -o cm2_in2 -t 30 --tmpdir . --force

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
dRep dereplicate drep_comapre_species_out -g merged_genomes/*.fa \
  -p 32 -comp 50 -con 10 -sa 0.95
# Mark mOTU as novel if the best ANI to all references < 95%
```

---

## 7. Taxonomic annotation and phylogenetic tree (GTDB-Tk + IQ-TREE)
```bash
# Classify genomes with GTDB-Tk using GTDB release R220
gtdbtk classify_wf --genome_dir motus/ -x fa --out_dir gtdbtk_out --cpus 32 --db_dir /path/to/gtdb_r220

# Build concatenated-marker phylogeny using GTDB-Tk alignments
ALIGN=gtdbtk_out/align/

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





## 13. Compare TEG cAMPs to public AMP databases (MMseqs2)
```bash
# Prepare MMseqs2 databases and perform easy-search
mmseqs createdb camp_representatives.faa camp_db
mmseqs createdb public_amps_combined.faa public_db
mmseqs easy-search camp_db public_db mmseqs_search_res tmpdir --min-seq-id 0.75 -e 1e-5 -s 7.5

# Parse results to retain hits with identity >= 75% and E-value <= 1e-5
```


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
