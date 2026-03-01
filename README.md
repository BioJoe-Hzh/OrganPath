# OrganPath

`organpath` is a command-line package for building phylogenetic trees from a `vcf.gz` file:

1. sample coverage/depth statistics and removal of highly-missing samples
2. uncertain genotype masking to missing (`./.`)
3. per-sample consensus generation from reference (`N` for missing)
4. MSA with `mafft` and trimming with `trimal`
5. tree inference with `iqtree`/`iqtree2` (UFBoot), via `OrganPath PhyView`

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

Or run tree inference separately:

```bash
OrganPath PhyView -i organpath_out/trimmed.fasta -o organpath_out/phyview
# or call PhyView directly:
PhyView -i organpath_out/trimmed.fasta -o organpath_out/phyview

# default uses IQ-TREE -safe; disable only if needed:
OrganPath PhyView -i organpath_out/trimmed.fasta -o organpath_out/phyview --unsafe
```

Rename IDs in an existing tree only:

```bash
OrganPath RenameTree \
  -i input.treefile \
  -m id_map.txt \
  -o renamed.treefile
```
