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
  --min-identity 0.95 \
  --min-len 1000 \
  --min-non-n-len 8000
```

`--organelle-mode` strategies:
- `plant_pt`: chloroplast-oriented sorting; output is rotated to seed start/orientation.
- `plant_mt`: mitochondrial homolog filtering against seed (reduce nuclear contamination).
- `animal_mt`: mitochondrial sorting with seed-start rotation for consistent coordinates.
- `generic`: original contig-order behavior.

Built-in default profile by type (`min_identity / min_len / gap_n`):
- `plant_pt`: `0.95 / 1000 / 100`
- `plant_mt`: `0.95 / 3000 / 100`
- `animal_mt`: `0.95 / 1000 / 100`
- `generic`: `0.95 / 1000 / 100`

Plant chloroplast single-IR mode (LSC + IR + SSC) with `cpstools`:

```bash
OrganPath sortOrgan \
  -i getorgan_out \
  -o sortorgan_out \
  -s seed_pt.fa \
  --organelle-mode plant_pt \
  --pt-single-ir \
  --pt-keep-ir auto \
  --cpstools-bin cpstools \
  --min-identity 0.95 \
  --min-len 1000
```

Notes:
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

Outputs:
- `sortorgan_summary.tsv` (per-sample summary)
- `assembled_samples.fasta` (merged multifasta)
- per-sample `*.organellar.fasta`

`--min-non-n-len` only controls whether a sample is included in `assembled_samples.fasta`.
Per-sample fasta is still written, and summary marks filtered rows as `FILTERED` with `non_n_len`.

Step 3 (Panel): Align and trim multifasta:

```bash
OrganPath align \
  -i sortorgan_out/assembled_samples.fasta \
  -o align_out \
  --threads 24 \
  --auto-reverse \
  --trim
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

For GetOrganelle outputs, `sortOrgan` now prefers `*path_sequence*.fasta` candidates.
When multiple path candidates exist, it de-duplicates equivalent sequences and selects the best candidate by mapping score to the seed.

After `sortOrgan`, you can run PanGraph -> DIPPER -> TWILIGHT -> panmanUtils directly:

```bash
OrganPath panman \
  -i out_mt/sortOrgan/assembled_samples.fasta \
  -o out_mt/panman \
  --run-pangraph \
  --run-dipper \
  --run-twilight \
  --run-panman
```

`OrganPath panman` supports two routes:
- Full route: `pangraph -> dipper -> twilight -> panmanUtils`
- User-tree route: provide `--guide-tree your_tree.nwk` and skip TWILIGHT (`--no-run-twilight`)

Three organelle channels:

```bash
# Plant chloroplast channel
OrganPath plant_pt -i reads_dir -o out_pt -s seed_pt.fa --jobs 5 --threads 16

# Plant mitochondrial channel (PanGraph -> DIPPER/TWILIGHT -> panmanUtils -> blocks)
OrganPath plant_mt \
  -i reads_dir \
  -o out_mt \
  -s seed_mt.fa \
  --jobs 5 \
  --threads 16 \
  --run-pangraph \
  --pangraph-args <PANGRAPH_ARGS...> \
  --run-dipper \
  --dipper-args <DIPPER_ARGS...> \
  --run-twilight \
  --twilight-args <TWILIGHT_ARGS...> \
  --run-panman \
  --panman-args <PANMAN_ARGS...>

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
1. MAFFT alignment with auto reverse-complement detection (`--adjustdirectionaccurately`)
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

## Plant mtBlocks Workflow

`mtBlocks` now supports this order:
1. PanGraph: build pangenome graph (`pangraph_output.json`)
2. DIPPER and/or TWILIGHT: build alignment/guide tree (user-defined args)
3. panmanUtils: derive homologous local blocks
4. per-block `mafft` + `trimal`
5. concatenate all kept blocks to `mt_supermatrix.fasta` (`mt_partitions.txt`)

Example:

```bash
OrganPath mtBlocks \
  -i out_mt/sortOrgan/assembled_samples.fasta \
  -o out_mt/mtBlocks \
  --run-pangraph \
  --pangraph-args <PANGRAPH_ARGS...> \
  --run-dipper \
  --dipper-args <DIPPER_ARGS...> \
  --run-twilight \
  --twilight-args <TWILIGHT_ARGS...> \
  --run-panman \
  --panman-args <PANMAN_ARGS...>
```

Argument placeholders supported in `--pangraph-args/--dipper-args/--twilight-args/--panman-args`:
- `{input_fasta}`: input multifasta
- `{pangraph_out}`: PanGraph working directory
- `{pangraph_json}`: PanGraph JSON path (default: `pangraph/pangraph_output.json`)
- `{dipper_out}`: DIPPER working directory
- `{dipper_graph}`: expected graph path (default: `dipper/graph.gfa`)
- `{aln_fasta}`: alignment file used downstream
- `{twilight_out}`: TWILIGHT working directory
- `{guide_tree}`: guide tree path (default: `twilight/guide_tree.nwk`)
- `{panman_out}`: panman working directory
- `{panman_file}`: default panman output path (`panman_run/result.panman`)
