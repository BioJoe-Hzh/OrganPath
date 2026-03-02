# OrganPath

`organpath` is a command-line package for building phylogenetic trees from a `vcf.gz` file:

1. sample coverage/depth statistics and removal of highly-missing samples
2. uncertain genotype masking to missing (`./.`)
3. per-sample consensus generation from reference (`N` for missing)
4. MSA with `mafft` and trimming with `trimal`
5. PhyView stage:
   - ML tree inference with `iqtree2`/`iqtree` (UFBoot)
   - pairwise distance matrix and NJ tree
   - PopART haplotype-network input preparation

## Install

```bash
cd OrganPath
python -m pip install -e .
```

## Server Dependencies

`OrganPath` requires external tools in `PATH`:

- `mafft`
- `trimal`
- `iqtree2` (preferred; fallback to `iqtree`)

Python dependency:
- `biopython` (for NJ tree construction)

Recommended installation on Linux server:

```bash
conda install -c bioconda mafft trimal iqtree
```

Then validate:

```bash
OrganPath check
# or
organpath check
```

## Usage

Step 1: Batch organelle assembly with GetOrganelle:

```bash
OrganPath getOrgan \
  -i reads_dir \
  -o getorgan_out \
  -s seed.fa \
  --r1-suffix _1.fastq.gz \
  --r2-suffix _2.fastq.gz \
  --organelle-type embplant_pt \
  --jobs 5 \
  --threads 16
```

Step 2: Sort/orient assembled contigs and build per-sample fasta:

```bash
OrganPath sortOrgan \
  -i getorgan_out \
  -o sortorgan_out \
  -s seed.fa \
  --min-identity 0.95 \
  --min-len 1000
```

Outputs:
- `sortorgan_summary.tsv` (per-sample summary)
- `assembled_samples.fasta` (merged multifasta)
- per-sample `*.organellar.fasta`

Three organelle channels:

```bash
# Plant chloroplast channel
OrganPath plant_pt -i reads_dir -o out_pt -s seed_pt.fa --jobs 5 --threads 16

# Plant mitochondrial channel (with panmanUtils block route)
OrganPath plant_mt -i reads_dir -o out_mt -s seed_mt.fa --jobs 5 --threads 16 --run-panman --panman-bin panmanUtils --panman-args ...

# Animal mitochondrial channel (simpler direct alignment route)
OrganPath animal_mt -i reads_dir -o out_animal_mt -s seed_animal_mt.fa --jobs 5 --threads 16 --run-ml
```

```bash
OrganPath run \
  -v input.vcf.gz \
  -r ref.fa \
  -o organpath_out \
  --sample-threads 8 \
  --run-phyview
```

Default filtering thresholds:
- `--min-coverage 0.5`
- `--min-mean-depth 10`
- `--min-dp 8`
- `--min-gq 30`

Sample rename options (`--id-map`):
- Rename mode (`oldID newID`, second column unique):
  `A01 Sample_01`
- Population-prefix mode (`oldID POP`, second column repeated):
  `A01 POP1` -> output name becomes `POP1_A01`

Outputs include:
- `name_map.tsv` (`old_id -> output_id`)
- `sample_stats.tsv` (both old and output IDs)
- renamed IDs in `filtered.kept.masked.vcf` and fasta headers

Or run `PhyView` as explicit single-mode jobs (must choose one):

```bash
# 1) ML tree only (IQ-TREE, default safe mode):
OrganPath PhyView -i organpath_out/trimmed.fasta -o organpath_out/phyview_ml --run_ml

# 2) NJ tree only:
OrganPath PhyView -i organpath_out/trimmed.fasta -o organpath_out/phyview_nj --run_nj

# 3) PopART input only:
OrganPath PhyView -i organpath_out/trimmed.fasta -o organpath_out/phyview_popart --run_popart

# Optional: execute PopART CLI after preparing input files:
OrganPath PhyView -i organpath_out/trimmed.fasta -o organpath_out/phyview_popart --run_popart --exec-popart --popart-bin popart
```

Rename IDs in an existing tree only:

```bash
OrganPath RenameTree \
  -i input.treefile \
  -m id_map.txt \
  -o renamed.treefile
```
