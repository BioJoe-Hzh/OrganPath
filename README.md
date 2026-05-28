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
- `snp-sites` (for MSA -> VCF)
- `bcftools` (for atomize/split/filter SNP VCF)
- `beast` + `treeannotator` (optional, for BEAST tree mode)

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
  --organelle-mode plant_mt \
  --cp-exclude-ref cp_ref.fa \
  --cp-exclude-min-identity 0.98 \
  --cp-exclude-min-len 5000 \
  --cp-exclude-min-query-cov 0.8 \
  --orient-min-query-cov 0.6 \
  --min-identity 0.90 \
  --min-len 3000 \
  --min-non-n-len 8000
```

`sortOrgan -i` accepts three input layouts:
- an `OrganPath getOrgan` / GetOrganelle output directory with one subdirectory per sample
- a directory containing FASTA files directly; each FASTA filename prefix becomes the sample ID
- a single FASTA file; its filename prefix becomes the sample ID

`--organelle-mode` strategies:
- `plant_pt`: chloroplast-oriented sorting; output is rotated to seed start/orientation.
- `plant_mt`: mitochondrial homolog filtering against seed (reduce nuclear contamination).
  You can also provide `--cp-exclude-ref` to remove large chloroplast-like contigs while keeping short plastid-derived transfers.
  Output keeps selected contigs as separate FASTA records per sample; contigs are not joined by `N`.
  In this mode, `--gap-n` is ignored.
- `animal_mt`: mitochondrial sorting with seed-start rotation for consistent coordinates.
- `generic`: original contig-order behavior.

Built-in default profile by type (`min_identity / min_len / gap_n`):
- `plant_pt`: `0.95 / 1000 / 100`
- `plant_mt`: `0.90 / 3000 / 100`
- `animal_mt`: `0.95 / 1000 / 100`
- `generic`: `0.95 / 1000 / 100`

Plant chloroplast single-IR mode (LSC + IR + SSC) with `cpstools`:

```bash
OrganPath sortOrgan \
  -i getorgan_out \
  -o sortorgan_out \
  -s seed_pt.fa \
  --organelle-mode plant_pt \
  --pt-keep-ir auto \
  --cpstools-bin cpstools \
  --min-identity 0.95 \
  --min-len 1000
```

Notes:
- In `--organelle-mode plant_pt`, OrganPath now enables single-IR logic by default and requires seed-based cp regions to be resolvable; if seed cannot be split into chloroplast regions, command fails.
- If `--cpstools-args` is not provided, OrganPath now runs a built-in default workflow:
  `cpstools IR -> cpstools Seq -m LSC -> cpstools IR`
  and auto-generates `cp_regions.tsv`.
- You can skip running `cpstools` inside OrganPath by passing `--cp-regions cp_regions.tsv`.
- `cp_regions.tsv` must include labels for `LSC`, `SSC`, and one of `IRB/IRA/IR`.
- In `plant_pt --pt-single-ir`, OrganPath now uses two routes:
  - complete-like (single long contig): sample-level cpstools IR + single-IR reorder
  - fragmented: assign contigs to `LSC/IR/SSC`, stitch in order with `N` gaps
- `sortorgan_summary.tsv` records per-sample type and missing estimate in `message`
  (`type:complete` or `type:fragmented`, plus `missing_bp` and `expected_len`).
- In `plant_pt` mode, `message` also includes `part_orient:LSC:...,SSC:...` to show
  per-partition orientation normalization against seed.
  Orientation is now scored by `blastn` (FWD vs RC) per partition.

Outputs:
- `sortorgan_summary.tsv` (per-sample summary)
- `assembled_samples.fasta` (merged multifasta)
- per-sample `*.organellar.fasta`
  In `plant_mt`, this is a multi-record FASTA with one record per retained contig.
- for `plant_pt` single-IR: per-sample `sample.LSC.fasta`, `sample.IR.fasta`, `sample.SSC.fasta`
- for `plant_pt` single-IR: merged `partitions/LSC_samples.fasta`, `partitions/IR_samples.fasta`, `partitions/SSC_samples.fasta`

`--min-non-n-len` only controls whether a sample is included in `assembled_samples.fasta`.
Per-sample fasta is still written, and summary marks filtered rows as `FILTERED` with `non_n_len`.
`--orient-min-query-cov` requires a seed/reference hit to cover at least this fraction of the query contig before it is used for sorting/orientation and `ref_covered_frac` summary statistics.

For `plant_mt`, chloroplast exclusion is optional and conservative:
- `--cp-exclude-ref`: chloroplast reference fasta used only for exclusion checking
- `--cp-exclude-min-identity`: default `0.98`
- `--cp-exclude-min-len`: default `5000`, so only larger cp-like alignments are considered for exclusion
- `--cp-exclude-min-query-cov`: default `0.8`, so the cp hit must cover most of the contig
- short cp-derived inserts are therefore kept by default, reducing the chance of removing intracellular plastid transfer fragments

Step 3 (Panel): Align and trim multifasta:

```bash
OrganPath align \
  -i sortorgan_out/assembled_samples.fasta \
  -o align_out \
  --threads 24 \
  --trim
```

For plant chloroplast partition-aware alignment (recommended when SSC inversion/isomer may affect whole-genome alignment):

```bash
OrganPath align \
  -i sortorgan_out \
  -o align_out \
  --plant-pt-partition \
  --threads 24 \
  --trim
```

In partition mode, OrganPath now does:
- MAFFT on `LSC/IR/SSC` separately
- optionally drop samples with high per-partition missing (`--partition-max-missing-frac`, default `1.0`)
- trim/filter each partition separately
- concatenate partitions to final `aligned.fasta` and `trimmed.fasta`
- write `partition_sample_stats.tsv` with per-sample/per-partition missing rate and non-missing length
- write `partition_sample_filter.tsv` with KEEP/DROP decision before trim
- if MAFFT adds `_R_` prefixes during auto-reverse, OrganPath normalizes sample IDs back before filtering/concatenation

Use ClipKIT instead of trimAl when trimming:

```bash
OrganPath align \
  -i sortorgan_out \
  -o align_clipkit_out \
  --plant-pt-partition \
  --auto-reverse \
  --trim-tool clipkit \
  --clipkit-mode smart-gap \
  --clipkit-log \
  --max-missing-frac 0.5 \
  --threads 48
```

Step 4 (Panel): ML tree from trimmed alignment:

```bash
OrganPath PhyView \
  -i align_out/trimmed.fasta \
  -o phyview_ml \
  --run_ml \
  --threads AUTO \
  --ufboot 1000 \
  --model MFP
```

Specify one or more outgroups for IQ-TREE rooted output:

```bash
OrganPath PhyView \
  -i align_out/trimmed.fasta \
  -o phyview_ml \
  --run_ml \
  --outgroup OUTGROUP_1,OUTGROUP_2 \
  --threads AUTO \
  --ufboot 3000 \
  --model MFP
```

Prepare Pathphynder panel files from VCF + tree:

```bash
OrganPath Pathphynder \
  --prepare \
  -v input.vcf.gz \
  -r ref.fa \
  -t tree.nwk \
  -o pathphynder_prepare_out \
  --prefix panel_name
```

If your panel VCF has multiple chromosomes/contigs and `phynder -B` fails, flatten to one synthetic coordinate system:

```bash
OrganPath Pathphynder \
  --prepare \
  -v input.vcf.gz \
  -r ref.fa \
  -t tree.nwk \
  -o pathphynder_prepare_out \
  --prefix panel_name \
  --concat-multi-chrom \
  --concat-chrom-name mt_concat
```

This runs:
- `bcftools norm -m -any -a -f` to atomize/split multi-allelic/complex variants (no filtering)
- optional `--concat-multi-chrom`: rewrite normalized VCF to one synthetic `CHROM`, shifting each contig by cumulative reference length and writing `panel_name.concat_pos_map.tsv` to map new positions back to original `CHROM:POS`
- optional `--concat-multi-chrom`: also writes `panel_name.concat_ref.fa` (single synthetic reference sequence) used for coordinate-consistent `--findpath`
- `phynder -B` on normalized VCF to create `panel_name.snp`
- `pathPhynder -s prepare` to build prepare files
- writes `pathphynder_prepare_manifest.tsv` for downstream usage

When `--findpath` is pointed to this prepare output (`--tree_data .../tree_data`), OrganPath auto-detects concat-mode from the prepare manifest and automatically uses `panel_name.concat_ref.fa` for mapping/placement coordinates.

Run Pathphynder placement directly from FASTQ (map to ref -> MAPQ filter -> remove duplicates -> mapDamage rescale -> `pathPhynder -s all`):

```bash
OrganPath Pathphynder \
  --findpath \
  --tree_data pathphynder_prepare_out/tree_data \
  -t tree.nwk \
  -r ref.fa \
  -o pathphynder_findpath_out \
  --fastq sample.fastq.gz \
  --sample-id SAMPLE1 \
  --bwa-aln-seedlen 1024 \
  --bwa-aln-mismatch 0.01 \
  --min-mapq 20 \
  --min-baseq 20 \
  --pathphynder-filtering-mode no-filter \
  --pathphynder-mismatch-threshold 0.5 \
  --pathphynder-max-tolerance 9999
```

To only check existing intermediates (without rerunning):

```bash
OrganPath Pathphynder \
  --findpath \
  --tree_data pathphynder_prepare_out/tree_data \
  -t tree.nwk \
  -r ref.fa \
  --fastq sample.fastq.gz \
  -o pathphynder_findpath_out \
  --sample-id SAMPLE1 \
  --check-only
```

For GetOrganelle outputs, `sortOrgan` now prefers `*path_sequence*.fasta` candidates.
When multiple path candidates exist, it de-duplicates equivalent sequences and selects the best candidate by mapping score to the seed.

After `sortOrgan`, recommended default is PanGraph + user-provided guide tree + panmanUtils (no DIPPER/TWILIGHT):

```bash
OrganPath panman \
  -i out_mt/sortOrgan/assembled_samples.fasta \
  -o out_mt/panman \
  --run-pangraph \
  --guide-tree your_tree.nwk \
  --run-panman
```

`OrganPath panman` supports two routes:
- Default route: `pangraph -> panmanUtils` (with `--guide-tree`)
- Optional full route: `pangraph -> dipper -> twilight -> panmanUtils` (enable with `--run-dipper --run-twilight`)

Three organelle channels:

```bash
# Plant chloroplast channel
OrganPath plant_pt -i reads_dir -o out_pt -s seed_pt.fa --jobs 5 --threads 16

# Plant mitochondrial channel (sortOrgan -> Pangraph LCBs -> mtBlocks tree)
OrganPath plant_mt \
  -i reads_dir \
  -o out_mt \
  -s seed_mt.fa \
  --jobs 5 \
  --threads 16 \
  --ufboot 1000

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

# Or directly from multifasta (e.g. sortOrgan output):
OrganPath run \
  -i out_pt/sortOrgan/assembled_samples.fasta \
  -o Pt_ml_tree

# More robust for highly divergent genomes (extra strict post-trim filtering):
OrganPath run \
  -i out_pt/sortOrgan/assembled_samples.fasta \
  -o Pt_ml_tree_strict \
  --auto-reverse \
  --max-missing-frac 0.2 \
  --snp-only
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

# 4) BEAST tree mode (from aligned/trimmed fasta):
OrganPath PhyView \
  -i organpath_out/trimmed.fasta \
  -o organpath_out/phyview_beast \
  --run_beast \
  --beast-template default \
  --beast-prefix mito_panel \
  --beast-chain-length 20000000 \
  --beast-store-every 5000
```

`--beast-template` accepts either:
- `default` (built-in template path in OrganPath)
- a custom XML template path

MSA to strict biallelic SNP VCF (for `assembled_samples.fasta`):

```bash
OrganPath MSA2VCF \
  -i out_pt/sortOrgan/assembled_samples.fasta \
  -o out_pt/msa2vcf \
  --auto-reverse \
  --trim
```

This command performs:
1. MAFFT alignment with auto reverse-complement detection (`--adjustdirection`)
2. optional trimAl
3. `snp-sites -v` to generate raw VCF from MSA
4. `bcftools norm -a` to atomize complex/MNV records
5. `bcftools norm -m -any` to split multi-allelic records
6. keep only biallelic SNPs (`bcftools view -v snps -m2 -M2`)

Rename IDs in an existing tree only:

```bash
OrganPath RenameTree \
  -i input.treefile \
  -m id_map.txt \
  -o renamed.treefile
```

Concatenate aligned multi-FASTA blocks by sample ID:

```bash
OrganPath concatFasta \
  -i gene1.fasta gene2.fasta gene3.fasta \
  -o concatenated.fasta \
  --samples sample_order.txt \
  --missing-report concatenated.missing.tsv \
  --partitions concatenated.partitions.tsv
```

Inputs are concatenated in the exact order provided. If a sample is missing from
one block, OrganPath fills that block with `-` characters matching the block
length and records the event in the missing report. If `--samples` is omitted,
sample order is inferred from first appearance across the input files.

Build ngsLCA-ready taxonomy from a tree + MSA (+ optional tip rename map):

```bash
OrganPath tree2taxonomy \
  -t Salicaceae_cp_final.treefile \
  -i Salicaceae_cp_final_MSA.fasta \
  --tip-map tip_map.tsv \
  --leaf-taxid-map leaf_taxid_map.tsv \
  --base-taxonomy-tsv Populus16.taxonomy.tree.tsv \
  -o tree2taxonomy_out
```

This writes:
- `renamed.treefile`
- `renamed.msa.fasta`
- `renamed.ungapped.fasta`
- `renamed.acc2taxid`
- `nodes.dmp`
- `names.dmp`
- `clade_taxonomy.tsv`
- `tip_name_map.audit.tsv`

Key behavior:
- if `--tip-map` is given, tree tips and MSA headers must both be fully covered by the map
- internal node labels like `99.7/100` are kept only when all numeric parts exceed `--support-cutoff`
- accepted internal nodes are added as `clade` taxa starting from taxid `100000000`
- leaf taxa are re-parented to the nearest accepted clade
- acc2tax uses the renamed tip headers and gap-free sequences

## Plant mtBlocks Workflow

`mtBlocks` now supports this order:
1. Read `sortOrgan` output directory or mt multifasta input
2. Run `pangraph build`
3. Run `pangraph export block-sequences`
4. Flag duplicated samples in each block
5. If one sample has duplicate identical sequences in a block, keep one copy; if duplicates differ, keep the copy most similar to the other samples
6. initial per-block `mafft --adjustdirection`
7. drop low-quality samples within each block by consensus identity/coverage from the initial alignment, then count sample presence
8. rerun `mafft --adjustdirection` on the remaining samples for that block
9. run `trimal` on the final per-block alignment
10. By default, keep only blocks still present in at least 50% of samples
11. optional site filtering
12. concatenate kept blocks from long to short into `mt_supermatrix.fasta` (`mt_partitions.txt`)
13. build an ML tree from the concatenated supermatrix with IQ-TREE

This is the recommended plant mitochondrial route:
- keep `sortOrgan` `plant_mt` outputs as multi-contig FASTA per sample
- let Pangraph derive homologous blocks directly from those contigs
- align each block separately and concatenate only after block detection

Example:

```bash
OrganPath mtBlocks \
  -i out_mt/sortOrgan \
  -o out_mt/mtBlocks \
  --block-jobs 8 \
  --block-min-sample-frac 0.5 \
  --block-min-sites 300 \
  --run-ml
```

`mtBlocks -i` accepts:
- a `sortOrgan` output directory
- a directory containing FASTA files directly; each FASTA filename prefix becomes the sample ID
- a legacy multifasta input

Directory input is preferred because it preserves each sample's separate mt contigs.

`mtBlocks` now handles Pangraph internally:
- `pangraph build --circular --output-json pangraph/pangraph_output.json`
- `pangraph export block-sequences --output pangraph_blocks pangraph/pangraph_output.json`
- by default this internal Pangraph step is on; you only need `--no-run-pangraph` if you already have `--pangraph-json` and/or `--blocks-dir`
- in sortOrgan directory mode, samples are pre-filtered by `sortorgan_summary.tsv`; default keeps only `ref_covered_frac >= 0.5`
- the final concatenated alignment keeps exactly that filtered sample set; if a sample is missing in a block it is padded with `-`

Useful `mtBlocks` block filters:
- `--sample-min-ref-cover-frac`: default `0.5`; pre-filter samples before Pangraph using `sortorgan_summary.tsv`
- `--threads`: one unified thread setting; Pangraph uses it, mtBlocks block-parallelism inherits it by default, and IQ-TREE also uses it
- `--block-jobs`: optional override for per-block MAFFT/trim/filter parallelism; default inherits from `--threads`
- `--block-min-sample-frac`: default `0.5`; keep only blocks present in at least half the samples
- `--block-sample-min-identity`: default `0.85`; within each block, drop samples below this post-trim consensus identity before block presence filtering
- `--block-sample-min-cover-frac`: default `0.5`; within each block, drop samples whose non-missing coverage is below this fraction before block presence filtering
- `--block-max-missing-frac`: optional extra post-trim site filter; default is disabled so trimAl output is used directly
- `--block-snp-only`: keep SNP columns only within each block
- `--block-min-sites`: default `300`; skip blocks that are too short after trim/filter

`mtblocks_summary.tsv` now records raw record count, unique sample count, sample presence fraction,
duplicate-sample flags, conflicting multi-copy samples resolved per block, short-sample removals, quality-filtered samples,
aligned/trimmed/final block length, final missing fraction, kept output fasta, and skip/fail reason for each block.
Final primary outputs are:
- `mt_supermatrix.fasta`
- `mtblocks.treefile` (copied from IQ-TREE output when `--run-ml` is enabled)
