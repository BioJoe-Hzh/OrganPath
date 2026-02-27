# OrganPath

`OrganPath` is a command-line package for building phylogenetic trees from a `vcf.gz` file:

1. sample coverage/depth statistics and removal of highly-missing samples
2. uncertain genotype masking to missing (`./.`)
3. per-sample consensus generation from reference (`N` for missing)
4. MSA with `mafft` and trimming with `trimal`
5. tree inference with `iqtree`/`iqtree2` (UFBoot), via `OrganPath PhyView`

## Install

```bash
cd vcf2tree
python -m pip install -e .
```

## Server Dependencies

`OrganPath` requires external tools in `PATH`:

- `mafft`
- `trimal`
- `iqtree2` (or `iqtree`)

Recommended installation on Linux server:

```bash
conda install -c bioconda mafft trimal iqtree
```

Then validate:

```bash
OrganPath check
```

## Usage

```bash
OrganPath run \
  -v input.vcf.gz \
  -r ref.fa \
  -o organpath_out \
  --run-phyview
```

Or run tree inference separately:

```bash
OrganPath PhyView -i organpath_out/trimmed.fasta -o organpath_out/phyview
```
