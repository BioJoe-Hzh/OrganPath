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

Or run tree inference separately:

```bash
OrganPath PhyView -i organpath_out/trimmed.fasta -o organpath_out/phyview
# or call PhyView directly:
PhyView -i organpath_out/trimmed.fasta -o organpath_out/phyview

# default uses IQ-TREE -safe; disable only if needed:
OrganPath PhyView -i organpath_out/trimmed.fasta -o organpath_out/phyview --unsafe
```
