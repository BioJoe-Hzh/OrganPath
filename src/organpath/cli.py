import argparse
import concurrent.futures
import gzip
import logging
import shutil
import subprocess
import sys
import importlib.util
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from . import __version__

logger = logging.getLogger(__name__)


MISSING_GT = {".", "./.", ".|.", "././.", ".|.|."}
DEFAULT_MIN_COVERAGE = 0.5
DEFAULT_MIN_MEAN_DEPTH = 10.0
DEFAULT_MIN_DP = 8
DEFAULT_MIN_GQ = 30.0


@dataclass
class SampleStats:
    total_sites: int = 0
    called_sites: int = 0
    dp_sum: float = 0.0

    @property
    def coverage(self) -> float:
        if self.total_sites == 0:
            return 0.0
        return self.called_sites / self.total_sites

    @property
    def mean_depth(self) -> float:
        if self.total_sites == 0:
            return 0.0
        return self.dp_sum / self.total_sites


def open_maybe_gzip(path: Path, mode: str = "rt"):
    if path.suffix == ".gz":
        return gzip.open(path, mode)
    return path.open(mode)


def parse_ref_fasta(path: Path) -> Dict[str, List[str]]:
    refs: Dict[str, List[str]] = {}
    curr = None
    chunks: List[str] = []
    with path.open("rt") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if curr is not None:
                    refs[curr] = list("".join(chunks).upper())
                curr = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if curr is not None:
        refs[curr] = list("".join(chunks).upper())
    if not refs:
        raise ValueError(f"Reference fasta is empty: {path}")
    return refs


def split_gt(gt: str) -> List[str]:
    if gt in MISSING_GT:
        return []
    if "/" in gt:
        return gt.split("/")
    if "|" in gt:
        return gt.split("|")
    if gt == ".":
        return []
    return [gt]


def is_gt_missing(gt: str) -> bool:
    return gt in MISSING_GT or gt == "." or "." in gt


def is_uncertain_genotype(
    sample_field: str,
    format_keys: List[str],
    min_dp: int,
    min_gq: float,
    het_as_missing: bool,
) -> bool:
    vals = sample_field.split(":")
    f = {k: vals[i] if i < len(vals) else "." for i, k in enumerate(format_keys)}
    gt = f.get("GT", ".")
    if is_gt_missing(gt):
        return True

    alleles = split_gt(gt)
    if het_as_missing and len(alleles) >= 2 and len(set(alleles)) > 1:
        return True

    dp_txt = f.get("DP", ".")
    if dp_txt not in {".", ""}:
        try:
            if float(dp_txt) < min_dp:
                return True
        except ValueError:
            return True

    gq_txt = f.get("GQ", ".")
    if gq_txt not in {".", ""}:
        try:
            if float(gq_txt) < min_gq:
                return True
        except ValueError:
            return True

    return False


def extract_depth(sample_field: str, format_keys: List[str]) -> float:
    vals = sample_field.split(":")
    f = {k: vals[i] if i < len(vals) else "." for i, k in enumerate(format_keys)}
    dp_txt = f.get("DP", ".")
    if dp_txt in {".", ""}:
        return 0.0
    try:
        return float(dp_txt)
    except ValueError:
        return 0.0


def choose_base_from_gt(gt: str, ref: str, alts: List[str]) -> str:
    if is_gt_missing(gt):
        return "N"
    alleles = split_gt(gt)
    if not alleles:
        return "N"

    idx_txt = alleles[0]
    if idx_txt == ".":
        return "N"
    try:
        idx = int(idx_txt)
    except ValueError:
        return "N"

    if idx == 0:
        return ref.upper()
    alt_i = idx - 1
    if 0 <= alt_i < len(alts):
        alt = alts[alt_i]
        if len(ref) == 1 and len(alt) == 1:
            return alt.upper()
    return "N"


def read_header_and_samples(vcf_path: Path) -> Tuple[List[str], List[str]]:
    headers: List[str] = []
    samples: List[str] = []
    with open_maybe_gzip(vcf_path, "rt") as fh:
        for raw in fh:
            if raw.startswith("##"):
                headers.append(raw.rstrip("\n"))
                continue
            if raw.startswith("#CHROM"):
                cols = raw.rstrip("\n").split("\t")
                samples = cols[9:]
                headers.append(raw.rstrip("\n"))
                break
    if not samples:
        raise ValueError(f"Failed to read samples from VCF: {vcf_path}")
    return headers, samples


def compute_sample_stats(vcf_path: Path, samples: List[str]) -> Dict[str, SampleStats]:
    stats = {s: SampleStats() for s in samples}
    with open_maybe_gzip(vcf_path, "rt") as fh:
        for raw in fh:
            if raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            fmt = fields[8].split(":")
            sample_fields = fields[9:]
            for i, sf in enumerate(sample_fields):
                s = samples[i]
                ss = stats[s]
                ss.total_sites += 1
                vals = sf.split(":")
                gt = vals[fmt.index("GT")] if "GT" in fmt and fmt.index("GT") < len(vals) else "."
                if not is_gt_missing(gt):
                    ss.called_sites += 1
                ss.dp_sum += extract_depth(sf, fmt)
    return stats


def write_stats_table_with_name_map(
    stats: Dict[str, SampleStats],
    name_map: Dict[str, str],
    path: Path,
) -> None:
    with path.open("wt") as out:
        out.write("old_id\toutput_id\ttotal_sites\tcalled_sites\tcoverage\tmean_depth\n")
        for sample, s in stats.items():
            out.write(
                f"{sample}\t{name_map.get(sample, sample)}\t"
                f"{s.total_sites}\t{s.called_sites}\t{s.coverage:.4f}\t{s.mean_depth:.4f}\n"
            )


def load_sample_name_map(samples: List[str], id_map_path: Optional[Path]) -> Tuple[Dict[str, str], str]:
    name_map = {s: s for s in samples}
    if id_map_path is None:
        return name_map, "none"
    if not id_map_path.exists():
        raise FileNotFoundError(f"id-map file not found: {id_map_path}")

    raw_map: Dict[str, str] = {}
    with id_map_path.open("rt") as fh:
        for ln, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split()
            if len(cols) < 2:
                raise ValueError(f"Invalid id-map line {ln}: expected at least 2 columns")
            raw_map[cols[0]] = cols[1]

    mapped_vals = [raw_map[s] for s in samples if s in raw_map]
    mode = "rename"
    if len(mapped_vals) != len(set(mapped_vals)):
        mode = "pop"

    for s in samples:
        if s not in raw_map:
            continue
        name_map[s] = raw_map[s] if mode == "rename" else f"{raw_map[s]}_{s}"

    new_names = list(name_map.values())
    if len(new_names) != len(set(new_names)):
        raise ValueError(
            "Output sample names are not unique after id-map conversion; "
            "please adjust id_map/pop assignments."
        )
    return name_map, mode


def safe_filename(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name)


def write_name_map(path: Path, samples: List[str], name_map: Dict[str, str], mode: str) -> None:
    with path.open("wt") as out:
        out.write(f"# mode={mode}\n")
        out.write("old_id\toutput_id\n")
        for s in samples:
            out.write(f"{s}\t{name_map.get(s, s)}\n")


def extract_tree_tip_names(newick_text: str) -> List[str]:
    names: List[str] = []

    for m in re.finditer(r"(?<=\(|,)'([^']*)'(?=\s*[:),;])", newick_text):
        names.append(m.group(1))
    for m in re.finditer(r"(?<=\(|,)\s*([^'():;,][^():;,]*)\s*(?=\s*[:),;])", newick_text):
        cand = m.group(1).strip()
        if cand:
            names.append(cand)

    seen = set()
    ordered: List[str] = []
    for n in names:
        if n not in seen:
            seen.add(n)
            ordered.append(n)
    return ordered


def rename_newick_tips(newick_text: str, name_map: Dict[str, str]) -> str:
    def repl_quoted(match: re.Match) -> str:
        old = match.group(1)
        return "'" + name_map.get(old, old) + "'"

    def repl_unquoted(match: re.Match) -> str:
        old = match.group(1).strip()
        return name_map.get(old, old)

    txt = re.sub(r"(?<=\(|,)'([^']*)'(?=\s*[:),;])", repl_quoted, newick_text)
    txt = re.sub(r"(?<=\(|,)\s*([^'():;,][^():;,]*)\s*(?=\s*[:),;])", repl_unquoted, txt)
    return txt


def keep_samples_by_threshold(
    stats: Dict[str, SampleStats],
    min_coverage: float,
    min_depth: float,
) -> Tuple[List[str], List[str]]:
    keep: List[str] = []
    drop: List[str] = []
    for sample, s in stats.items():
        low_cov = s.coverage < min_coverage
        low_dp = s.mean_depth <= min_depth
        if low_cov and low_dp:
            drop.append(sample)
        else:
            keep.append(sample)
    return keep, drop


def mask_gt_to_missing(sample_field: str, format_keys: List[str]) -> str:
    vals = sample_field.split(":")
    if "GT" in format_keys:
        i = format_keys.index("GT")
        if i < len(vals):
            vals[i] = "./."
        else:
            while len(vals) <= i:
                vals.append(".")
            vals[i] = "./."
    return ":".join(vals)


def filter_vcf_and_build_consensus(
    vcf_path: Path,
    ref_fasta: Path,
    out_vcf: Path,
    sample_fastas_dir: Path,
    multifasta_path: Path,
    keep_samples: List[str],
    sample_name_map: Dict[str, str],
    min_dp: int,
    min_gq: float,
    het_as_missing: bool,
    sample_threads: int = 1,
) -> None:
    refs = parse_ref_fasta(ref_fasta)
    sample_seqs: Dict[str, Dict[str, List[str]]] = {
        s: {ctg: seq.copy() for ctg, seq in refs.items()} for s in keep_samples
    }

    sample_index: Dict[str, int] = {}

    executor: Optional[concurrent.futures.ThreadPoolExecutor] = None
    if sample_threads > 1:
        executor = concurrent.futures.ThreadPoolExecutor(max_workers=sample_threads)

    with open_maybe_gzip(vcf_path, "rt") as inp, out_vcf.open("wt") as out:
        for raw in inp:
            if raw.startswith("##"):
                out.write(raw)
                continue

            if raw.startswith("#CHROM"):
                cols = raw.rstrip("\n").split("\t")
                all_samples = cols[9:]
                sample_index = {s: i for i, s in enumerate(all_samples)}
                keep_set = set(keep_samples)
                selected = [sample_name_map.get(s, s) for s in all_samples if s in keep_set]
                out.write("\t".join(cols[:9] + selected) + "\n")
                continue

            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue

            chrom, pos_txt, _vid, ref, alt_txt = fields[:5]
            if chrom not in refs:
                continue
            pos = int(pos_txt)
            if pos < 1 or pos > len(refs[chrom]):
                continue

            alts = alt_txt.split(",") if alt_txt else []
            fmt_keys = fields[8].split(":")
            all_sample_fields = fields[9:]

            def process_one_sample(s: str) -> Tuple[str, str, str]:
                i = sample_index[s]
                sf = all_sample_fields[i] if i < len(all_sample_fields) else "."
                uncertain = is_uncertain_genotype(
                    sf,
                    fmt_keys,
                    min_dp=min_dp,
                    min_gq=min_gq,
                    het_as_missing=het_as_missing,
                )
                if uncertain:
                    sf = mask_gt_to_missing(sf, fmt_keys)
                vals = sf.split(":")
                gt = vals[fmt_keys.index("GT")] if "GT" in fmt_keys and fmt_keys.index("GT") < len(vals) else "."
                base = "N"
                if len(ref) == 1:
                    base = choose_base_from_gt(gt, ref, alts)
                return s, sf, base

            if executor is not None:
                results = list(executor.map(process_one_sample, keep_samples))
            else:
                results = [process_one_sample(s) for s in keep_samples]

            selected_fields = []
            for s, sf, base in results:
                if len(ref) == 1:
                    sample_seqs[s][chrom][pos - 1] = base
                selected_fields.append(sf)

            out.write("\t".join(fields[:9] + selected_fields) + "\n")

    if executor is not None:
        executor.shutdown(wait=True)

    sample_fastas_dir.mkdir(parents=True, exist_ok=True)
    with multifasta_path.open("wt") as multi:
        for s in keep_samples:
            out_name = sample_name_map.get(s, s)
            seq = "".join(
                "".join(sample_seqs[s][ctg]) for ctg in sorted(sample_seqs[s].keys())
            )
            indiv = sample_fastas_dir / f"{safe_filename(out_name)}.fa"
            with indiv.open("wt") as out:
                out.write(f">{out_name}\n")
                out.write(f"{seq}\n")
            multi.write(f">{out_name}\n{seq}\n")


def run_command(cmd: List[str], stdout_path: Optional[Path] = None) -> None:
    logger.info("Running: %s", " ".join(cmd))
    if stdout_path is None:
        subprocess.run(cmd, check=True)
        return
    with stdout_path.open("wt") as out:
        subprocess.run(cmd, stdout=out, check=True)


def ensure_tool(name: str) -> str:
    p = shutil.which(name)
    if not p:
        install_hint = {
            "mafft": "Please install MAFFT and ensure `mafft` is in PATH.",
            "trimal": "Please install trimAl and ensure `trimal` is in PATH.",
            "iqtree": "Please install IQ-TREE (iqtree/iqtree2) and ensure it is in PATH.",
            "iqtree2": "Please install IQ-TREE (iqtree/iqtree2) and ensure it is in PATH.",
        }.get(name, f"Please install `{name}` and add it to PATH.")
        raise RuntimeError(f"Required tool not found in PATH: {name}. {install_hint}")
    return p


def run_alignment_and_trimming(
    multifasta: Path, aligned: Path, trimmed: Path, mafft_bin: str, trimal_bin: str
) -> None:
    run_command([mafft_bin, "--auto", str(multifasta)], stdout_path=aligned)
    run_command([trimal_bin, "-automated1", "-in", str(aligned), "-out", str(trimmed)])


def run_phyview(
    trimmed_fasta: Path,
    out_dir: Path,
    ufboot: int = 1000,
    threads: str = "AUTO",
    model: str = "MFP",
    safe: bool = True,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    iqtree_bin = shutil.which("iqtree2") or shutil.which("iqtree")
    if not iqtree_bin:
        raise RuntimeError("Required tool not found: iqtree2 or iqtree")
    prefix = out_dir / "organpath_tree"
    cmd = [
        iqtree_bin,
        "-s",
        str(trimmed_fasta),
        "-m",
        model,
        "-B",
        str(ufboot),
        "-T",
        threads,
        "--prefix",
        str(prefix),
    ]
    if safe:
        cmd.append("-safe")
    run_command(cmd)


def read_fasta_sequences(path: Path) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    cur = None
    chunks: List[str] = []
    with path.open("rt") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur is not None:
                    seqs[cur] = "".join(chunks).upper()
                cur = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if cur is not None:
        seqs[cur] = "".join(chunks).upper()
    if not seqs:
        raise ValueError(f"No sequences found in fasta: {path}")
    return seqs


def pairwise_distance(seq1: str, seq2: str) -> float:
    valid = 0
    diff = 0
    for a, b in zip(seq1, seq2):
        if a in {"N", "-", "?"} or b in {"N", "-", "?"}:
            continue
        valid += 1
        if a != b:
            diff += 1
    if valid == 0:
        return 1.0
    return diff / valid


def build_distance_matrix(seqs: Dict[str, str]) -> Tuple[List[str], List[List[float]]]:
    names = list(seqs.keys())
    if len(set(len(v) for v in seqs.values())) != 1:
        raise ValueError("Sequences in trimmed fasta do not have equal lengths.")
    n = len(names)
    mat = [[0.0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            d = pairwise_distance(seqs[names[i]], seqs[names[j]])
            mat[i][j] = d
            mat[j][i] = d
    return names, mat


def write_distance_matrix_tsv(names: List[str], mat: List[List[float]], out_path: Path) -> None:
    with out_path.open("wt") as out:
        out.write("sample\t" + "\t".join(names) + "\n")
        for i, name in enumerate(names):
            out.write(name + "\t" + "\t".join(f"{mat[i][j]:.6f}" for j in range(len(names))) + "\n")


def build_nj_tree(names: List[str], mat: List[List[float]], out_path: Path) -> None:
    try:
        from Bio import Phylo
        from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
    except ImportError as exc:
        raise RuntimeError(
            "Biopython is required for NJ tree construction. Install with: pip install biopython"
        ) from exc

    lower_triangle: List[List[float]] = []
    for i in range(len(names)):
        lower_triangle.append([mat[i][j] for j in range(i + 1)])

    dm = DistanceMatrix(names=names, matrix=lower_triangle)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)

    # NJ may yield small negative branch lengths due to non-additive distances
    # and numerical noise. Clamp to 0 for stable downstream visualization.
    neg_count = 0
    for clade in tree.find_clades():
        if clade.branch_length is not None and clade.branch_length < 0:
            clade.branch_length = 0.0
            neg_count += 1
    if neg_count > 0:
        logger.info("NJ tree: clamped %d negative branch lengths to 0.", neg_count)

    Phylo.write(tree, str(out_path), "newick")


def prepare_popart_inputs(trimmed_fasta: Path, out_dir: Path) -> Tuple[Path, Path]:
    seqs = read_fasta_sequences(trimmed_fasta)
    hap_to_samples: Dict[str, List[str]] = defaultdict(list)
    for sample, seq in seqs.items():
        hap_to_samples[seq].append(sample)

    hap_ids: Dict[str, str] = {}
    for i, hap_seq in enumerate(hap_to_samples.keys(), start=1):
        hap_ids[hap_seq] = f"HAP{i}"

    hap_tsv = out_dir / "haplotypes.tsv"
    with hap_tsv.open("wt") as out:
        out.write("haplotype_id\tsize\tsamples\n")
        for seq, samples in hap_to_samples.items():
            out.write(f"{hap_ids[seq]}\t{len(samples)}\t{','.join(samples)}\n")

    nchar = len(next(iter(seqs.values())))
    popart_nex = out_dir / "popart_input.nex"
    with popart_nex.open("wt") as out:
        out.write("#NEXUS\n\n")
        out.write("BEGIN TAXA;\n")
        out.write(f"  DIMENSIONS NTAX={len(seqs)};\n")
        out.write("  TAXLABELS\n")
        for s in seqs:
            out.write(f"    {s}\n")
        out.write("  ;\nEND;\n\n")
        out.write("BEGIN CHARACTERS;\n")
        out.write(f"  DIMENSIONS NCHAR={nchar};\n")
        out.write("  FORMAT DATATYPE=DNA MISSING=? GAP=-;\n")
        out.write("  MATRIX\n")
        for s, seq in seqs.items():
            out.write(f"    {s} {seq}\n")
        out.write("  ;\nEND;\n")
    return hap_tsv, popart_nex


def run_popart_cli(popart_bin: str, popart_input: Path, out_dir: Path, extra_args: List[str]) -> None:
    p = shutil.which(popart_bin)
    if not p:
        raise RuntimeError(
            f"PopART executable not found: {popart_bin}. "
            "Install PopART or set --popart-bin correctly."
        )
    cmd = [p] + extra_args + [str(popart_input)]
    logger.info("Running PopART command: %s", " ".join(cmd))
    run_command(cmd)


def discover_paired_reads(reads_dir: Path, r1_suffix: str, r2_suffix: str) -> Tuple[List[Tuple[str, Path, Path]], List[str]]:
    pairs: List[Tuple[str, Path, Path]] = []
    missing: List[str] = []
    for r1 in sorted(reads_dir.glob(f"*{r1_suffix}")):
        name = r1.name
        if not name.endswith(r1_suffix):
            continue
        sample = name[: -len(r1_suffix)]
        r2 = reads_dir / f"{sample}{r2_suffix}"
        if not r2.exists():
            missing.append(sample)
            continue
        pairs.append((sample, r1, r2))
    return pairs, missing


def run_getorgan_one_sample(
    go_bin: str,
    sample: str,
    r1: Path,
    r2: Path,
    seed: Path,
    out_dir: Path,
    organelle_type: str,
    threads: int,
    rounds: int,
    kmer: str,
    max_reads: int,
    extra_args: List[str],
) -> Tuple[str, str, str]:
    sample_out = out_dir / sample
    cmd = [
        go_bin,
        "-1",
        str(r1),
        "-2",
        str(r2),
        "-o",
        str(sample_out),
        "-F",
        organelle_type,
        "-s",
        str(seed),
        "-t",
        str(threads),
        "-R",
        str(rounds),
    ]
    if kmer:
        cmd += ["-k", kmer]
    if max_reads is not None and max_reads > 0:
        cmd += ["-n", str(max_reads)]
    if extra_args:
        cmd += extra_args

    try:
        run_command(cmd)
        # Validate that GetOrganelle actually produced outputs.
        if not sample_out.exists():
            return sample, "FAIL", "no_output_dir_created"
        expected = (
            list(sample_out.rglob("*contigs*.fasta"))
            + list(sample_out.rglob("*.fastg"))
            + list(sample_out.rglob("*.gfa"))
        )
        if not expected:
            return sample, "FAIL", "no_assembly_outputs_found"
        return sample, "OK", "-"
    except subprocess.CalledProcessError as exc:
        return sample, "FAIL", str(exc.returncode)
    except Exception as exc:
        return sample, "FAIL", str(exc)


def cmd_get_organ(args: argparse.Namespace) -> int:
    reads_dir = Path(args.reads_dir).resolve()
    out_dir = Path(args.outdir).resolve()
    seed = Path(args.seed).resolve()
    if not reads_dir.exists():
        raise FileNotFoundError(f"reads-dir not found: {reads_dir}")
    if not seed.exists():
        raise FileNotFoundError(f"seed fasta not found: {seed}")

    go_bin = shutil.which(args.getorganelle_bin)
    if not go_bin:
        raise RuntimeError(f"GetOrganelle executable not found: {args.getorganelle_bin}")

    out_dir.mkdir(parents=True, exist_ok=True)
    pairs, missing = discover_paired_reads(reads_dir, args.r1_suffix, args.r2_suffix)
    if not pairs:
        raise ValueError("No paired-end samples discovered with given suffixes.")
    if missing:
        logger.warning("Missing R2 for %d samples; skipped.", len(missing))

    summary = out_dir / "getorgan_summary.tsv"
    pair_map = {sample: (r1, r2) for sample, r1, r2 in pairs}
    results: List[Tuple[str, str, str]] = []

    jobs = max(1, int(args.jobs))
    workers = min(jobs, len(pairs))
    logger.info("Running GetOrganelle with %d parallel samples; %d threads per sample.", workers, args.threads)

    if workers == 1:
        for sample, r1, r2 in pairs:
            results.append(
                run_getorgan_one_sample(
                    go_bin=go_bin,
                    sample=sample,
                    r1=r1,
                    r2=r2,
                    seed=seed,
                    out_dir=out_dir,
                    organelle_type=args.organelle_type,
                    threads=args.threads,
                    rounds=args.rounds,
                    kmer=args.kmer,
                    max_reads=args.max_reads,
                    extra_args=args.extra_args,
                )
            )
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:
            futs = []
            for sample, r1, r2 in pairs:
                futs.append(
                    ex.submit(
                        run_getorgan_one_sample,
                        go_bin,
                        sample,
                        r1,
                        r2,
                        seed,
                        out_dir,
                        args.organelle_type,
                        args.threads,
                        args.rounds,
                        args.kmer,
                        args.max_reads,
                        args.extra_args,
                    )
                )
            for fut in concurrent.futures.as_completed(futs):
                results.append(fut.result())

    with summary.open("wt") as out:
        out.write("sample\tr1\tr2\tstatus\tmessage\n")
        for sample, status, msg in sorted(results, key=lambda x: x[0]):
            r1, r2 = pair_map[sample]
            out.write(f"{sample}\t{r1}\t{r2}\t{status}\t{msg}\n")

    if missing:
        with (out_dir / "missing_pairs.txt").open("wt") as f:
            f.write("\n".join(missing) + "\n")
    logger.info("getOrgan completed. Summary: %s", summary)
    return 0


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def pick_sample_contig_and_fastg(sample_dir: Path) -> Tuple[Optional[Path], Optional[Path]]:
    contig = None
    for pat in ("*contigs.fasta", "*contigs.fa", "contigs.fasta", "contigs.fa"):
        cands = sorted(sample_dir.rglob(pat))
        if cands:
            contig = cands[0]
            break
    fastg_cands = sorted(sample_dir.rglob("*.fastg"))
    fastg = fastg_cands[0] if fastg_cands else None
    return contig, fastg


def parse_minimap2_paf(paf_path: Path, min_identity: float, min_len: int) -> List[Tuple[str, int, int, str, float, int]]:
    hits: List[Tuple[str, int, int, str, float, int]] = []
    with paf_path.open("rt") as fh:
        for raw in fh:
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            qname = parts[0]
            strand = parts[4]
            tstart = int(parts[7])
            tend = int(parts[8])
            nmatch = int(parts[9])
            alen = int(parts[10])
            if alen <= 0:
                continue
            ident = nmatch / alen
            if ident < min_identity or alen < min_len:
                continue
            hits.append((qname, tstart, tend, strand, ident, alen))
    return hits


def parse_blast_tab(tab_path: Path, min_identity: float, min_len: int) -> List[Tuple[str, int, int, str, float, int]]:
    hits: List[Tuple[str, int, int, str, float, int]] = []
    with tab_path.open("rt") as fh:
        for raw in fh:
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            qname = parts[0]
            sstart = int(parts[1])
            send = int(parts[2])
            alen = int(parts[3])
            pident = float(parts[4]) / 100.0
            strand = "+" if sstart <= send else "-"
            if pident < min_identity or alen < min_len:
                continue
            hits.append((qname, min(sstart, send), max(sstart, send), strand, pident, alen))
    return hits


def choose_best_hits(hits: List[Tuple[str, int, int, str, float, int]]) -> List[Tuple[str, int, int, str, float, int]]:
    best: Dict[str, Tuple[str, int, int, str, float, int]] = {}
    for h in hits:
        qname = h[0]
        if qname not in best:
            best[qname] = h
            continue
        old = best[qname]
        if (h[4], h[5]) > (old[4], old[5]):
            best[qname] = h
    selected = list(best.values())
    selected.sort(key=lambda x: (x[1], x[2]))
    return selected


def build_sample_assembly_from_contigs(
    contig_fa: Path,
    seed_fa: Path,
    sample_name: str,
    out_dir: Path,
    min_identity: float,
    min_len: int,
    gap_n: int,
    aligner: str,
) -> Tuple[str, int, int]:
    out_dir.mkdir(parents=True, exist_ok=True)
    contigs = read_fasta_sequences(contig_fa)
    if not contigs:
        raise ValueError(f"No contigs found in {contig_fa}")

    tool = aligner
    if tool == "auto":
        tool = "minimap2" if shutil.which("minimap2") else "blastn"

    tmp = out_dir / f"{sample_name}.align.tmp"
    if tool == "minimap2":
        if not shutil.which("minimap2"):
            raise RuntimeError("minimap2 not found in PATH.")
        with tmp.open("wt") as out:
            subprocess.run(
                ["minimap2", "-x", "asm5", str(seed_fa), str(contig_fa)],
                stdout=out,
                check=True,
            )
        hits = parse_minimap2_paf(tmp, min_identity=min_identity, min_len=min_len)
    elif tool == "blastn":
        if not shutil.which("blastn"):
            raise RuntimeError("blastn not found in PATH.")
        with tmp.open("wt") as out:
            subprocess.run(
                [
                    "blastn",
                    "-query",
                    str(contig_fa),
                    "-subject",
                    str(seed_fa),
                    "-outfmt",
                    "6 qseqid sstart send length pident bitscore qstart qend",
                ],
                stdout=out,
                check=True,
            )
        hits = parse_blast_tab(tmp, min_identity=min_identity, min_len=min_len)
    else:
        raise ValueError(f"Unsupported aligner: {aligner}")

    try:
        tmp.unlink(missing_ok=True)
    except Exception:
        pass

    chosen = choose_best_hits(hits)
    if not chosen:
        raise ValueError("No contigs passed identity/length filter against seed.")

    n_gap = "N" * gap_n
    seq_parts: List[str] = []
    kept = 0
    with (out_dir / f"{sample_name}.selected_contigs.tsv").open("wt") as tab:
        tab.write("contig\tseed_start\tseed_end\tstrand\tidentity\taln_len\n")
        for contig_id, sstart, send, strand, ident, alen in chosen:
            if contig_id not in contigs:
                continue
            cseq = contigs[contig_id]
            if strand == "-":
                cseq = reverse_complement(cseq)
            seq_parts.append(cseq)
            tab.write(f"{contig_id}\t{sstart}\t{send}\t{strand}\t{ident:.4f}\t{alen}\n")
            kept += 1

    merged = n_gap.join(seq_parts)
    out_fa = out_dir / f"{sample_name}.organellar.fasta"
    with out_fa.open("wt") as out:
        out.write(f">{sample_name}\n{merged}\n")
    return str(out_fa), kept, len(merged)


def cmd_sort_organ(args: argparse.Namespace) -> int:
    in_dir = Path(args.input_dir).resolve()
    out_dir = Path(args.outdir).resolve()
    seed = Path(args.seed).resolve()
    if not in_dir.exists():
        raise FileNotFoundError(f"input-dir not found: {in_dir}")
    if not seed.exists():
        raise FileNotFoundError(f"seed fasta not found: {seed}")
    out_dir.mkdir(parents=True, exist_ok=True)

    sample_dirs = sorted([p for p in in_dir.iterdir() if p.is_dir()])
    if not sample_dirs:
        raise ValueError("No sample folders found in input-dir.")

    summary = out_dir / "sortorgan_summary.tsv"
    all_multi = out_dir / "assembled_samples.fasta"
    with summary.open("wt") as sumf, all_multi.open("wt") as mf:
        sumf.write(
            "sample\tstatus\tcontigs_fasta\tfastg\tselected_contigs\tassembled_len\tassembled_fasta\tmessage\n"
        )
        for sdir in sample_dirs:
            sample = sdir.name
            contig, fastg = pick_sample_contig_and_fastg(sdir)
            if contig is None:
                sumf.write(f"{sample}\tFAIL\t-\t{fastg or '-'}\t0\t0\t-\tcontigs not found\n")
                continue

            sample_out = out_dir / sample
            try:
                out_fa, nsel, alen = build_sample_assembly_from_contigs(
                    contig_fa=contig,
                    seed_fa=seed,
                    sample_name=sample,
                    out_dir=sample_out,
                    min_identity=args.min_identity,
                    min_len=args.min_len,
                    gap_n=args.gap_n,
                    aligner=args.aligner,
                )
                seqs = read_fasta_sequences(Path(out_fa))
                for sid, seq in seqs.items():
                    mf.write(f">{sid}\n{seq}\n")
                sumf.write(
                    f"{sample}\tOK\t{contig}\t{fastg or '-'}\t{nsel}\t{alen}\t{out_fa}\t-\n"
                )
            except Exception as exc:
                sumf.write(
                    f"{sample}\tFAIL\t{contig}\t{fastg or '-'}\t0\t0\t-\t{str(exc).replace(chr(9), ' ')}\n"
                )

    logger.info("sortOrgan completed. Summary: %s", summary)
    logger.info("Merged multi-fasta: %s", all_multi)
    return 0


def list_block_fastas(blocks_dir: Path) -> List[Path]:
    pats = ["*.fa", "*.fasta", "*.fas", "*.fna"]
    files: List[Path] = []
    for p in pats:
        files.extend(sorted(blocks_dir.glob(p)))
    return files


def align_trim_one_block(block_fa: Path, out_dir: Path, mafft_bin: str, trimal_bin: str) -> Optional[Path]:
    seqs = read_fasta_sequences(block_fa)
    if len(seqs) < 2:
        return None
    aln = out_dir / f"{block_fa.stem}.aln.fasta"
    trimmed = out_dir / f"{block_fa.stem}.trim.fasta"
    run_command([mafft_bin, "--auto", str(block_fa)], stdout_path=aln)
    run_command([trimal_bin, "-automated1", "-in", str(aln), "-out", str(trimmed)])
    return trimmed


def concatenate_blocks(trimmed_blocks: List[Path], out_fa: Path, out_partitions: Path) -> None:
    block_data: List[Tuple[str, Dict[str, str], int]] = []
    all_samples = set()
    for bf in trimmed_blocks:
        seqs = read_fasta_sequences(bf)
        if not seqs:
            continue
        lengths = {len(v) for v in seqs.values()}
        if len(lengths) != 1:
            raise ValueError(f"Block has inconsistent lengths: {bf}")
        blen = next(iter(lengths))
        block_data.append((bf.stem, seqs, blen))
        all_samples.update(seqs.keys())

    all_samples = sorted(all_samples)
    if not block_data:
        raise ValueError("No trimmed blocks available for concatenation.")

    with out_partitions.open("wt") as part:
        start = 1
        for bid, _seqs, blen in block_data:
            end = start + blen - 1
            part.write(f"DNA,{bid}={start}-{end}\n")
            start = end + 1

    with out_fa.open("wt") as out:
        for s in all_samples:
            concat = []
            for _bid, seqs, blen in block_data:
                concat.append(seqs.get(s, "N" * blen))
            out.write(f">{s}\n{''.join(concat)}\n")


def cmd_mt_blocks(args: argparse.Namespace) -> int:
    input_fa = Path(args.input).resolve()
    out_dir = Path(args.outdir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    if not input_fa.exists():
        raise FileNotFoundError(f"Input multifasta not found: {input_fa}")

    dipper_out = out_dir / "dipper"
    twilight_out = out_dir / "twilight"
    panman_out = out_dir / "panman_run"
    dipper_out.mkdir(parents=True, exist_ok=True)
    twilight_out.mkdir(parents=True, exist_ok=True)
    panman_out.mkdir(parents=True, exist_ok=True)

    aln_fa = Path(args.aln_file).resolve() if args.aln_file else input_fa
    guide_tree = Path(args.guide_tree).resolve() if args.guide_tree else (twilight_out / "guide_tree.nwk")
    dipper_graph = Path(args.dipper_graph).resolve() if args.dipper_graph else (dipper_out / "graph.gfa")

    def _render_tokens(items: List[str]) -> List[str]:
        rendered: List[str] = []
        for x in items:
            rendered.append(
                x.replace("{input_fasta}", str(input_fa))
                .replace("{dipper_out}", str(dipper_out))
                .replace("{twilight_out}", str(twilight_out))
                .replace("{dipper_graph}", str(dipper_graph))
                .replace("{aln_fasta}", str(aln_fa))
                .replace("{guide_tree}", str(guide_tree))
                .replace("{panman_out}", str(panman_out))
            )
        return rendered

    if args.run_dipper:
        dipper_bin = shutil.which(args.dipper_bin)
        if not dipper_bin:
            raise RuntimeError(f"DIPPER executable not found: {args.dipper_bin}")
        if not args.dipper_args:
            raise ValueError(
                "--run-dipper requires --dipper-args. "
                "Use placeholders like {input_fasta}, {dipper_out}, {dipper_graph}, {aln_fasta}."
            )
        run_command([dipper_bin] + _render_tokens(list(args.dipper_args)))

    if args.run_twilight:
        twilight_bin = shutil.which(args.twilight_bin)
        if not twilight_bin:
            raise RuntimeError(f"TWILIGHT executable not found: {args.twilight_bin}")
        if not args.twilight_args:
            raise ValueError(
                "--run-twilight requires --twilight-args. "
                "Use placeholders like {dipper_graph}, {aln_fasta}, {twilight_out}, {guide_tree}."
            )
        run_command([twilight_bin] + _render_tokens(list(args.twilight_args)))

    blocks_dir = Path(args.blocks_dir).resolve() if args.blocks_dir else (out_dir / "panman_blocks")
    if args.run_panman:
        panman_bin = shutil.which(args.panman_bin)
        if not panman_bin:
            raise RuntimeError(f"panman executable not found: {args.panman_bin}")
        blocks_dir.mkdir(parents=True, exist_ok=True)
        # panman ecosystem commonly exposes `panmanUtils` with its own CLI syntax.
        # We therefore pass through user-provided args directly.
        if not args.panman_args:
            raise ValueError(
                "--run-panman requires --panman-args for your panmanUtils workflow. "
                "Please provide panmanUtils arguments that write block FASTA files into --blocks-dir."
            )
        cmd = [panman_bin] + _render_tokens(list(args.panman_args))
        run_command(cmd)

    if not blocks_dir.exists():
        raise FileNotFoundError(f"blocks-dir not found: {blocks_dir}")

    block_fastas = list_block_fastas(blocks_dir)
    if not block_fastas:
        raise ValueError(f"No block FASTA files found in {blocks_dir}")

    mafft_bin = ensure_tool(args.mafft_bin)
    trimal_bin = ensure_tool(args.trimal_bin)
    aln_dir = out_dir / "blocks_aln"
    aln_dir.mkdir(parents=True, exist_ok=True)

    trimmed_blocks: List[Path] = []
    with (out_dir / "mtblocks_summary.tsv").open("wt") as sumf:
        sumf.write("block_file\tstatus\ttrimmed_block\n")
        for bf in block_fastas:
            try:
                trimmed = align_trim_one_block(bf, aln_dir, mafft_bin=mafft_bin, trimal_bin=trimal_bin)
                if trimmed is None:
                    sumf.write(f"{bf}\tSKIP\t-\n")
                    continue
                trimmed_blocks.append(trimmed)
                sumf.write(f"{bf}\tOK\t{trimmed}\n")
            except Exception:
                sumf.write(f"{bf}\tFAIL\t-\n")

    supermatrix = out_dir / "mt_supermatrix.fasta"
    partitions = out_dir / "mt_partitions.txt"
    concatenate_blocks(trimmed_blocks, supermatrix, partitions)
    logger.info("Concatenated supermatrix: %s", supermatrix)
    logger.info("Partitions file: %s", partitions)

    if args.run_ml:
        run_phyview(
            trimmed_fasta=supermatrix,
            out_dir=out_dir / "ml",
            ufboot=args.ufboot,
            threads=args.threads,
            model=args.model,
            safe=not args.unsafe,
        )
    return 0


def cmd_channel_plant_pt(args: argparse.Namespace) -> int:
    ns = argparse.Namespace(
        reads_dir=args.reads_dir,
        outdir=str((Path(args.outdir).resolve() / "getOrgan")),
        seed=args.seed,
        r1_suffix=args.r1_suffix,
        r2_suffix=args.r2_suffix,
        organelle_type="embplant_pt",
        jobs=args.jobs,
        threads=args.threads,
        rounds=args.rounds,
        kmer=args.kmer,
        max_reads=args.max_reads,
        getorganelle_bin=args.getorganelle_bin,
        extra_args=args.getorgan_extra_args,
    )
    cmd_get_organ(ns)

    ss = argparse.Namespace(
        input_dir=str((Path(args.outdir).resolve() / "getOrgan")),
        outdir=str((Path(args.outdir).resolve() / "sortOrgan")),
        seed=args.seed,
        aligner=args.aligner,
        min_identity=args.min_identity,
        min_len=args.min_len_pt,
        gap_n=args.gap_n,
    )
    return cmd_sort_organ(ss)


def cmd_channel_plant_mt(args: argparse.Namespace) -> int:
    ns = argparse.Namespace(
        reads_dir=args.reads_dir,
        outdir=str((Path(args.outdir).resolve() / "getOrgan")),
        seed=args.seed,
        r1_suffix=args.r1_suffix,
        r2_suffix=args.r2_suffix,
        organelle_type="embplant_mt",
        jobs=args.jobs,
        threads=args.threads,
        rounds=args.rounds,
        kmer=args.kmer,
        max_reads=args.max_reads,
        getorganelle_bin=args.getorganelle_bin,
        extra_args=args.getorgan_extra_args,
    )
    cmd_get_organ(ns)

    ss = argparse.Namespace(
        input_dir=str((Path(args.outdir).resolve() / "getOrgan")),
        outdir=str((Path(args.outdir).resolve() / "sortOrgan")),
        seed=args.seed,
        aligner=args.aligner,
        min_identity=args.min_identity,
        min_len=args.min_len_mt,
        gap_n=args.gap_n,
    )
    cmd_sort_organ(ss)

    ms = argparse.Namespace(
        input=str((Path(args.outdir).resolve() / "sortOrgan" / "assembled_samples.fasta")),
        outdir=str((Path(args.outdir).resolve() / "mtBlocks")),
        aln_file=args.aln_file,
        guide_tree=args.guide_tree,
        dipper_graph=args.dipper_graph,
        run_dipper=args.run_dipper,
        dipper_bin=args.dipper_bin,
        dipper_args=args.dipper_args,
        run_twilight=args.run_twilight,
        twilight_bin=args.twilight_bin,
        twilight_args=args.twilight_args,
        blocks_dir=args.blocks_dir,
        run_panman=args.run_panman,
        panman_bin=args.panman_bin,
        panman_args=args.panman_args,
        mafft_bin=args.mafft_bin,
        trimal_bin=args.trimal_bin,
        run_ml=args.run_ml,
        ufboot=args.ufboot,
        threads=args.ml_threads,
        model=args.model,
        unsafe=args.unsafe,
    )
    return cmd_mt_blocks(ms)


def cmd_channel_animal_mt(args: argparse.Namespace) -> int:
    ns = argparse.Namespace(
        reads_dir=args.reads_dir,
        outdir=str((Path(args.outdir).resolve() / "getOrgan")),
        seed=args.seed,
        r1_suffix=args.r1_suffix,
        r2_suffix=args.r2_suffix,
        organelle_type="animal_mt",
        jobs=args.jobs,
        threads=args.threads,
        rounds=args.rounds,
        kmer=args.kmer,
        max_reads=args.max_reads,
        getorganelle_bin=args.getorganelle_bin,
        extra_args=args.getorgan_extra_args,
    )
    cmd_get_organ(ns)

    ss = argparse.Namespace(
        input_dir=str((Path(args.outdir).resolve() / "getOrgan")),
        outdir=str((Path(args.outdir).resolve() / "sortOrgan")),
        seed=args.seed,
        aligner=args.aligner,
        min_identity=args.min_identity,
        min_len=args.min_len_mt,
        gap_n=args.gap_n,
    )
    cmd_sort_organ(ss)

    # Animal mt usually linearize/alignment simpler; build ML directly from multifasta.
    sort_multi = Path(args.outdir).resolve() / "sortOrgan" / "assembled_samples.fasta"
    aligned = Path(args.outdir).resolve() / "animal_mt.aligned.fasta"
    trimmed = Path(args.outdir).resolve() / "animal_mt.trimmed.fasta"
    mafft_bin = ensure_tool(args.mafft_bin)
    trimal_bin = ensure_tool(args.trimal_bin)
    run_alignment_and_trimming(sort_multi, aligned, trimmed, mafft_bin=mafft_bin, trimal_bin=trimal_bin)
    if args.run_ml:
        run_phyview(
            trimmed_fasta=trimmed,
            out_dir=Path(args.outdir).resolve() / "ml",
            ufboot=args.ufboot,
            threads=args.ml_threads,
            model=args.model,
            safe=not args.unsafe,
        )
    return 0


def cmd_run(args: argparse.Namespace) -> int:
    out_dir = Path(args.outdir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    vcf_path = Path(args.vcf).resolve()
    ref_fasta = Path(args.ref).resolve()

    _headers, samples = read_header_and_samples(vcf_path)
    logger.info("Detected %d samples in VCF.", len(samples))
    id_map_path = Path(args.id_map).resolve() if args.id_map else None
    sample_name_map, id_map_mode = load_sample_name_map(samples, id_map_path)
    logger.info("Sample naming mode: %s", id_map_mode)
    write_name_map(out_dir / "name_map.tsv", samples, sample_name_map, id_map_mode)

    stats = compute_sample_stats(vcf_path, samples)
    stats_tsv = out_dir / "sample_stats.tsv"
    write_stats_table_with_name_map(stats, sample_name_map, stats_tsv)

    keep, drop = keep_samples_by_threshold(
        stats,
        min_coverage=args.min_coverage,
        min_depth=args.min_mean_depth,
    )
    logger.info("Kept %d samples, dropped %d samples.", len(keep), len(drop))

    with (out_dir / "kept_samples.txt").open("wt") as f:
        f.write("\n".join(sample_name_map[s] for s in keep) + ("\n" if keep else ""))
    with (out_dir / "dropped_samples.txt").open("wt") as f:
        f.write("\n".join(sample_name_map[s] for s in drop) + ("\n" if drop else ""))

    if not keep:
        logger.error("No samples passed filtering; stop.")
        return 2

    filtered_vcf = out_dir / "filtered.kept.masked.vcf"
    sample_fastas = out_dir / "sample_fastas"
    multifasta = out_dir / "all_samples.fasta"
    aligned = out_dir / "aligned.fasta"
    trimmed = out_dir / "trimmed.fasta"

    filter_vcf_and_build_consensus(
        vcf_path=vcf_path,
        ref_fasta=ref_fasta,
        out_vcf=filtered_vcf,
        sample_fastas_dir=sample_fastas,
        multifasta_path=multifasta,
        keep_samples=keep,
        sample_name_map=sample_name_map,
        min_dp=args.min_dp,
        min_gq=args.min_gq,
        het_as_missing=not args.keep_het,
        sample_threads=max(1, args.sample_threads),
    )
    logger.info("Generated filtered VCF and consensus FASTA files.")

    mafft_bin = ensure_tool(args.mafft_bin)
    trimal_bin = ensure_tool(args.trimal_bin)
    run_alignment_and_trimming(multifasta, aligned, trimmed, mafft_bin, trimal_bin)
    logger.info("Alignment and trimming completed.")

    if args.run_phyview:
        run_phyview(
            trimmed_fasta=trimmed,
            out_dir=out_dir / "phyview",
            ufboot=args.ufboot,
            threads=args.threads,
            model=args.model,
            safe=not args.unsafe,
        )
        logger.info("PhyView ML tree inference completed.")

    logger.info("All outputs are in: %s", out_dir)
    return 0


def cmd_phyview(args: argparse.Namespace) -> int:
    trimmed = Path(args.input).resolve()
    out_dir = Path(args.outdir).resolve()
    mode_count = int(args.run_ml) + int(args.run_nj) + int(args.run_popart)
    if mode_count != 1:
        raise ValueError("PhyView requires exactly one mode: --run_ml or --run_nj or --run_popart")
    if args.exec_popart and not args.run_popart:
        raise ValueError("--exec-popart is only valid together with --run_popart")

    if args.run_ml:
        run_phyview(
            trimmed_fasta=trimmed,
            out_dir=out_dir,
            ufboot=args.ufboot,
            threads=args.threads,
            model=args.model,
            safe=not args.unsafe,
        )
        logger.info("PhyView completed (ML): %s", out_dir)
        return 0

    if args.run_nj:
        out_dir.mkdir(parents=True, exist_ok=True)
        seqs = read_fasta_sequences(trimmed)
        names, mat = build_distance_matrix(seqs)
        write_distance_matrix_tsv(names, mat, out_dir / "pairwise_distance.tsv")
        build_nj_tree(names, mat, out_dir / "nj_tree.nwk")
        logger.info("PhyView completed (NJ): %s", out_dir)
        return 0

    if args.run_popart:
        out_dir.mkdir(parents=True, exist_ok=True)
        _hap_tsv, popart_nex = prepare_popart_inputs(trimmed, out_dir)
        if args.exec_popart:
            run_popart_cli(args.popart_bin, popart_nex, out_dir, args.popart_args)
        logger.info("PhyView completed (PopART inputs): %s", out_dir)
        return 0

    logger.info("PhyView completed: %s", out_dir)
    return 0


def cmd_rename_tree(args: argparse.Namespace) -> int:
    tree_in = Path(args.input).resolve()
    tree_out = Path(args.output).resolve()
    id_map_path = Path(args.id_map).resolve()

    if not tree_in.exists():
        raise FileNotFoundError(f"Input tree not found: {tree_in}")

    newick = tree_in.read_text()
    tip_names = extract_tree_tip_names(newick)
    if not tip_names:
        raise ValueError("No tip names were detected from the input tree.")

    name_map, mode = load_sample_name_map(tip_names, id_map_path)
    renamed = rename_newick_tips(newick, name_map)

    tree_out.parent.mkdir(parents=True, exist_ok=True)
    tree_out.write_text(renamed)
    map_path = tree_out.with_suffix(tree_out.suffix + ".name_map.tsv")
    write_name_map(map_path, tip_names, name_map, mode)

    logger.info("Renamed tree written to: %s", tree_out)
    logger.info("Name mapping written to: %s", map_path)
    logger.info("Rename mode: %s", mode)
    return 0


def cmd_check(args: argparse.Namespace) -> int:
    missing_python: List[str] = []
    missing_tools: List[str] = []

    for mod in ["gzip", "argparse", "pathlib"]:
        if importlib.util.find_spec(mod) is None:
            missing_python.append(mod)
    if importlib.util.find_spec("Bio") is None:
        missing_python.append("biopython (Bio)")

    for tool in ["mafft", "trimal"]:
        if shutil.which(tool) is None:
            missing_tools.append(tool)

    if shutil.which("iqtree2") is None and shutil.which("iqtree") is None:
        missing_tools.append("iqtree2/iqtree")

    print("OrganPath environment check")
    print(f"Python: {sys.executable}")
    print(f"Version: {sys.version.split()[0]}")

    if missing_python:
        print("Missing python modules:", ", ".join(missing_python))
    else:
        print("Python modules: OK")

    if missing_tools:
        print("Missing external tools:", ", ".join(missing_tools))
        print("Install hints:")
        print("- mafft: conda install -c bioconda mafft")
        print("- trimal: conda install -c bioconda trimal")
        print("- iqtree2: conda install -c bioconda iqtree")
        return 2

    print("External tools: OK")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="OrganPath",
        description="Organelle VCF to phylogeny pipeline",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")

    subs = parser.add_subparsers(dest="command")

    p_run = subs.add_parser(
        "run",
        help="Run full workflow from VCF to trimmed alignment",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p_run.add_argument("-v", "--vcf", required=True, help="Input VCF or VCF.GZ")
    p_run.add_argument("-r", "--ref", required=True, help="Reference fasta")
    p_run.add_argument("-o", "--outdir", required=True, help="Output directory")
    p_run.add_argument(
        "--id-map",
        help="Two-column file. oldID newID => rename; oldID POP (repeated POP labels) => rename to POP_oldID",
    )
    p_run.add_argument(
        "--min-coverage",
        type=float,
        default=DEFAULT_MIN_COVERAGE,
        help="Min sample coverage",
    )
    p_run.add_argument(
        "--min-mean-depth",
        type=float,
        default=DEFAULT_MIN_MEAN_DEPTH,
        help="Drop sample only when coverage < threshold AND mean depth <= this value",
    )
    p_run.add_argument(
        "--min-dp",
        type=int,
        default=DEFAULT_MIN_DP,
        help="Mark genotype missing if DP < this",
    )
    p_run.add_argument(
        "--min-gq",
        type=float,
        default=DEFAULT_MIN_GQ,
        help="Mark genotype missing if GQ < this",
    )
    p_run.add_argument("--keep-het", action="store_true", help="Keep heterozygous genotype (default: mask)")
    p_run.add_argument(
        "--sample-threads",
        type=int,
        default=1,
        help="Thread number for per-sample genotype/consensus processing",
    )
    p_run.add_argument("--mafft-bin", default="mafft", help="Path or name of mafft executable")
    p_run.add_argument("--trimal-bin", default="trimal", help="Path or name of trimal executable")
    p_run.add_argument("--run-phyview", action="store_true", help="Also run OrganPath PhyView after trimming")
    p_run.add_argument("--ufboot", type=int, default=1000, help="UFBoot replicate number for iqtree")
    p_run.add_argument("--threads", default="AUTO", help="Thread setting passed to iqtree -T")
    p_run.add_argument("--model", default="MFP", help="Model option passed to iqtree -m")
    p_run.add_argument(
        "--unsafe",
        action="store_true",
        help="Disable IQ-TREE safe likelihood kernel (default uses -safe)",
    )
    p_run.set_defaults(func=cmd_run)

    p_phy = subs.add_parser("PhyView", help="Run one selected phylogenetic/relationship task on trimmed multifasta")
    p_phy.add_argument("-i", "--input", required=True, help="Input trimmed multifasta (FASTA)")
    p_phy.add_argument("-o", "--outdir", required=True, help="Output directory for the selected task")
    p_phy.add_argument("--ufboot", type=int, default=1000, help="UFBoot replicate number (used with --run_ml)")
    p_phy.add_argument("--threads", default="AUTO", help="Thread setting for IQ-TREE -T (used with --run_ml)")
    p_phy.add_argument("--model", default="MFP", help="Model option for IQ-TREE -m (used with --run_ml)")
    mode_group = p_phy.add_mutually_exclusive_group(required=True)
    mode_group.add_argument("--run_ml", "--run-ml", action="store_true", help="Run ML tree (IQ-TREE)")
    mode_group.add_argument("--run_nj", "--run-nj", action="store_true", help="Run pairwise distance + NJ tree")
    mode_group.add_argument("--run_popart", "--run-popart", action="store_true", help="Prepare PopART haplotype inputs")
    p_phy.add_argument(
        "--exec-popart",
        action="store_true",
        help="Execute PopART CLI after preparing input files (used with --run_popart)",
    )
    p_phy.add_argument("--popart-bin", default="popart", help="PopART executable name/path (used with --run_popart)")
    p_phy.add_argument(
        "--popart-args",
        nargs="*",
        default=[],
        help="Extra args passed to PopART CLI (used with --run_popart)",
    )
    p_phy.add_argument(
        "--unsafe",
        action="store_true",
        help="Disable IQ-TREE safe likelihood kernel (used with --run_ml; default uses -safe)",
    )
    p_phy.set_defaults(func=cmd_phyview)

    p_rename = subs.add_parser(
        "RenameTree",
        help="Rename tip IDs in Newick tree using id_map/pop file",
    )
    p_rename.add_argument("-i", "--input", required=True, help="Input Newick tree file")
    p_rename.add_argument("-m", "--id-map", required=True, help="Two-column id_map/pop file")
    p_rename.add_argument("-o", "--output", required=True, help="Output Newick tree file")
    p_rename.set_defaults(func=cmd_rename_tree)

    p_get = subs.add_parser(
        "getOrgan",
        help="Batch assemble organellar genomes from paired-end reads with GetOrganelle",
    )
    p_get.add_argument("-i", "--reads-dir", required=True, help="Folder containing raw paired-end reads")
    p_get.add_argument("-o", "--outdir", required=True, help="Output folder for per-sample GetOrganelle results")
    p_get.add_argument("-s", "--seed", required=True, help="Seed fasta sequence for GetOrganelle")
    p_get.add_argument("--r1-suffix", default="_1.fastq.gz", help="R1 filename suffix")
    p_get.add_argument("--r2-suffix", default="_2.fastq.gz", help="R2 filename suffix")
    p_get.add_argument(
        "--organelle-type",
        default="embplant_pt",
        choices=["embplant_pt", "embplant_mt", "animal_mt", "fungus_mt", "other_pt", "other_mt"],
        help="GetOrganelle database type (-F)",
    )
    p_get.add_argument("--jobs", type=int, default=1, help="Number of samples assembled in parallel")
    p_get.add_argument("--threads", type=int, default=8, help="Threads per sample")
    p_get.add_argument("--rounds", type=int, default=15, help="GetOrganelle extension rounds (-R)")
    p_get.add_argument("--kmer", default="", help="Custom kmers for GetOrganelle, e.g. 21,45,65,85,105")
    p_get.add_argument("--max-reads", type=int, default=0, help="Optional max reads for GetOrganelle (-n)")
    p_get.add_argument(
        "--getorganelle-bin",
        default="get_organelle_from_reads.py",
        help="GetOrganelle executable name/path",
    )
    p_get.add_argument("--extra-args", nargs="*", default=[], help="Extra args passed to GetOrganelle")
    p_get.set_defaults(func=cmd_get_organ)

    p_sort = subs.add_parser(
        "sortOrgan",
        help="Sort and orient assembled organellar contigs by seed; output per-sample fasta",
    )
    p_sort.add_argument("-i", "--input-dir", required=True, help="Input folder from OrganPath getOrgan outputs")
    p_sort.add_argument("-o", "--outdir", required=True, help="Output folder for sorted assemblies")
    p_sort.add_argument("-s", "--seed", required=True, help="Seed/reference fasta used for contig ordering")
    p_sort.add_argument(
        "--aligner",
        default="auto",
        choices=["auto", "minimap2", "blastn"],
        help="Aligner for contig-to-seed mapping",
    )
    p_sort.add_argument("--min-identity", type=float, default=0.95, help="Minimum alignment identity")
    p_sort.add_argument("--min-len", type=int, default=1000, help="Minimum aligned length")
    p_sort.add_argument("--gap-n", type=int, default=100, help="Number of Ns inserted between contigs")
    p_sort.set_defaults(func=cmd_sort_organ)

    p_mt = subs.add_parser(
        "mtBlocks",
        help="Plant mitochondrial route: panman blocks -> block MAFFT+trimAl -> concatenated supermatrix",
    )
    p_mt.add_argument("-i", "--input", required=True, help="Input multifasta (per-sample mt assemblies)")
    p_mt.add_argument("-o", "--outdir", required=True, help="Output directory")
    p_mt.add_argument("--aln-file", help="Existing MSA fasta (if already generated)")
    p_mt.add_argument("--guide-tree", help="Existing guide tree nwk (if already generated)")
    p_mt.add_argument("--dipper-graph", help="Expected dipper graph path (for twilight/panman)")
    p_mt.add_argument("--run-dipper", action="store_true", help="Run DIPPER to produce graph/MSA")
    p_mt.add_argument("--dipper-bin", default="dipper", help="DIPPER executable name/path")
    p_mt.add_argument(
        "--dipper-args",
        nargs="*",
        default=[],
        help="Arguments passed to DIPPER. Supports placeholders: {input_fasta} {dipper_out} {dipper_graph} {aln_fasta}",
    )
    p_mt.add_argument("--run-twilight", action="store_true", help="Run TWILIGHT to produce guide tree")
    p_mt.add_argument("--twilight-bin", default="twilight", help="TWILIGHT executable name/path")
    p_mt.add_argument(
        "--twilight-args",
        nargs="*",
        default=[],
        help="Arguments passed to TWILIGHT. Supports placeholders: {dipper_graph} {aln_fasta} {twilight_out} {guide_tree}",
    )
    p_mt.add_argument("--blocks-dir", help="Directory containing block fasta files (if panman already run)")
    p_mt.add_argument("--run-panman", action="store_true", help="Run panman to generate blocks")
    p_mt.add_argument("--panman-bin", default="panmanUtils", help="panman executable name/path")
    p_mt.add_argument(
        "--panman-args",
        nargs="*",
        default=[],
        help="Arguments passed directly to panmanUtils (must generate block FASTA files into --blocks-dir)",
    )
    p_mt.add_argument("--mafft-bin", default="mafft", help="Path or name of mafft executable")
    p_mt.add_argument("--trimal-bin", default="trimal", help="Path or name of trimal executable")
    p_mt.add_argument("--run-ml", action="store_true", help="Run ML tree on concatenated supermatrix")
    p_mt.add_argument("--ufboot", type=int, default=1000, help="UFBoot replicate number")
    p_mt.add_argument("--threads", default="AUTO", help="Thread setting passed to iqtree -T")
    p_mt.add_argument("--model", default="MFP", help="Model option passed to iqtree -m")
    p_mt.add_argument("--unsafe", action="store_true", help="Disable IQ-TREE safe likelihood kernel")
    p_mt.set_defaults(func=cmd_mt_blocks)

    p_ch_pt = subs.add_parser(
        "plant_pt",
        help="Channel: plant chloroplast (GetOrganelle + sortOrgan)",
    )
    p_ch_pt.add_argument("-i", "--reads-dir", required=True, help="Folder containing paired reads")
    p_ch_pt.add_argument("-o", "--outdir", required=True, help="Output directory")
    p_ch_pt.add_argument("-s", "--seed", required=True, help="Seed/reference fasta")
    p_ch_pt.add_argument("--r1-suffix", default="_1.fastq.gz", help="R1 filename suffix")
    p_ch_pt.add_argument("--r2-suffix", default="_2.fastq.gz", help="R2 filename suffix")
    p_ch_pt.add_argument("--jobs", type=int, default=1, help="Parallel sample jobs for GetOrganelle")
    p_ch_pt.add_argument("--threads", type=int, default=8, help="Threads per GetOrganelle sample")
    p_ch_pt.add_argument("--rounds", type=int, default=15, help="GetOrganelle rounds")
    p_ch_pt.add_argument("--kmer", default="", help="GetOrganelle kmer string")
    p_ch_pt.add_argument("--max-reads", type=int, default=0, help="GetOrganelle -n max reads")
    p_ch_pt.add_argument("--getorganelle-bin", default="get_organelle_from_reads.py", help="GetOrganelle executable")
    p_ch_pt.add_argument("--getorgan-extra-args", nargs="*", default=[], help="Extra args for GetOrganelle")
    p_ch_pt.add_argument("--aligner", default="auto", choices=["auto", "minimap2", "blastn"], help="Contig aligner")
    p_ch_pt.add_argument("--min-identity", type=float, default=0.95, help="Minimum identity for contig selection")
    p_ch_pt.add_argument("--min-len-pt", type=int, default=1000, help="Minimum length for cp contig selection")
    p_ch_pt.add_argument("--gap-n", type=int, default=100, help="Ns inserted between selected contigs")
    p_ch_pt.set_defaults(func=cmd_channel_plant_pt)

    p_ch_mt = subs.add_parser(
        "plant_mt",
        help="Channel: plant mitochondria (GetOrganelle + sortOrgan + mtBlocks)",
    )
    p_ch_mt.add_argument("-i", "--reads-dir", required=True, help="Folder containing paired reads")
    p_ch_mt.add_argument("-o", "--outdir", required=True, help="Output directory")
    p_ch_mt.add_argument("-s", "--seed", required=True, help="Seed/reference fasta")
    p_ch_mt.add_argument("--r1-suffix", default="_1.fastq.gz", help="R1 filename suffix")
    p_ch_mt.add_argument("--r2-suffix", default="_2.fastq.gz", help="R2 filename suffix")
    p_ch_mt.add_argument("--jobs", type=int, default=1, help="Parallel sample jobs for GetOrganelle")
    p_ch_mt.add_argument("--threads", type=int, default=8, help="Threads per GetOrganelle sample")
    p_ch_mt.add_argument("--rounds", type=int, default=15, help="GetOrganelle rounds")
    p_ch_mt.add_argument("--kmer", default="", help="GetOrganelle kmer string")
    p_ch_mt.add_argument("--max-reads", type=int, default=0, help="GetOrganelle -n max reads")
    p_ch_mt.add_argument("--getorganelle-bin", default="get_organelle_from_reads.py", help="GetOrganelle executable")
    p_ch_mt.add_argument("--getorgan-extra-args", nargs="*", default=[], help="Extra args for GetOrganelle")
    p_ch_mt.add_argument("--aligner", default="auto", choices=["auto", "minimap2", "blastn"], help="Contig aligner")
    p_ch_mt.add_argument("--min-identity", type=float, default=0.95, help="Minimum identity for contig selection")
    p_ch_mt.add_argument("--min-len-mt", type=int, default=3000, help="Minimum length for mt contig selection")
    p_ch_mt.add_argument("--gap-n", type=int, default=100, help="Ns inserted between selected contigs")
    p_ch_mt.add_argument("--run-panman", action="store_true", help="Run panman to derive conserved blocks")
    p_ch_mt.add_argument("--aln-file", help="Existing MSA fasta for mtBlocks")
    p_ch_mt.add_argument("--guide-tree", help="Existing guide tree nwk for mtBlocks")
    p_ch_mt.add_argument("--dipper-graph", help="Expected dipper graph path")
    p_ch_mt.add_argument("--run-dipper", action="store_true", help="Run DIPPER before TWILIGHT/panman")
    p_ch_mt.add_argument("--dipper-bin", default="dipper", help="DIPPER executable")
    p_ch_mt.add_argument("--dipper-args", nargs="*", default=[], help="Arguments passed to DIPPER")
    p_ch_mt.add_argument("--run-twilight", action="store_true", help="Run TWILIGHT to generate guide tree")
    p_ch_mt.add_argument("--twilight-bin", default="twilight", help="TWILIGHT executable")
    p_ch_mt.add_argument("--twilight-args", nargs="*", default=[], help="Arguments passed to TWILIGHT")
    p_ch_mt.add_argument("--panman-bin", default="panmanUtils", help="panman executable")
    p_ch_mt.add_argument(
        "--panman-args",
        nargs="*",
        default=[],
        help="Arguments passed directly to panmanUtils",
    )
    p_ch_mt.add_argument("--blocks-dir", help="Existing panman blocks dir")
    p_ch_mt.add_argument("--mafft-bin", default="mafft", help="Path or name of mafft executable")
    p_ch_mt.add_argument("--trimal-bin", default="trimal", help="Path or name of trimal executable")
    p_ch_mt.add_argument("--run-ml", action="store_true", help="Run ML tree on concatenated blocks")
    p_ch_mt.add_argument("--ufboot", type=int, default=1000, help="UFBoot replicate number")
    p_ch_mt.add_argument("--ml-threads", default="AUTO", help="Thread setting passed to iqtree -T")
    p_ch_mt.add_argument("--model", default="MFP", help="Model option passed to iqtree -m")
    p_ch_mt.add_argument("--unsafe", action="store_true", help="Disable IQ-TREE safe likelihood kernel")
    p_ch_mt.set_defaults(func=cmd_channel_plant_mt)

    p_ch_amt = subs.add_parser(
        "animal_mt",
        help="Channel: animal mitochondria (GetOrganelle + sortOrgan + direct alignment/tree)",
    )
    p_ch_amt.add_argument("-i", "--reads-dir", required=True, help="Folder containing paired reads")
    p_ch_amt.add_argument("-o", "--outdir", required=True, help="Output directory")
    p_ch_amt.add_argument("-s", "--seed", required=True, help="Seed/reference fasta")
    p_ch_amt.add_argument("--r1-suffix", default="_1.fastq.gz", help="R1 filename suffix")
    p_ch_amt.add_argument("--r2-suffix", default="_2.fastq.gz", help="R2 filename suffix")
    p_ch_amt.add_argument("--jobs", type=int, default=1, help="Parallel sample jobs for GetOrganelle")
    p_ch_amt.add_argument("--threads", type=int, default=8, help="Threads per GetOrganelle sample")
    p_ch_amt.add_argument("--rounds", type=int, default=15, help="GetOrganelle rounds")
    p_ch_amt.add_argument("--kmer", default="", help="GetOrganelle kmer string")
    p_ch_amt.add_argument("--max-reads", type=int, default=0, help="GetOrganelle -n max reads")
    p_ch_amt.add_argument("--getorganelle-bin", default="get_organelle_from_reads.py", help="GetOrganelle executable")
    p_ch_amt.add_argument("--getorgan-extra-args", nargs="*", default=[], help="Extra args for GetOrganelle")
    p_ch_amt.add_argument("--aligner", default="auto", choices=["auto", "minimap2", "blastn"], help="Contig aligner")
    p_ch_amt.add_argument("--min-identity", type=float, default=0.95, help="Minimum identity for contig selection")
    p_ch_amt.add_argument("--min-len-mt", type=int, default=3000, help="Minimum length for mt contig selection")
    p_ch_amt.add_argument("--gap-n", type=int, default=100, help="Ns inserted between selected contigs")
    p_ch_amt.add_argument("--mafft-bin", default="mafft", help="Path or name of mafft executable")
    p_ch_amt.add_argument("--trimal-bin", default="trimal", help="Path or name of trimal executable")
    p_ch_amt.add_argument("--run-ml", action="store_true", help="Run ML tree after alignment")
    p_ch_amt.add_argument("--ufboot", type=int, default=1000, help="UFBoot replicate number")
    p_ch_amt.add_argument("--ml-threads", default="AUTO", help="Thread setting passed to iqtree -T")
    p_ch_amt.add_argument("--model", default="MFP", help="Model option passed to iqtree -m")
    p_ch_amt.add_argument("--unsafe", action="store_true", help="Disable IQ-TREE safe likelihood kernel")
    p_ch_amt.set_defaults(func=cmd_channel_animal_mt)

    p_chk = subs.add_parser("check", help="Check runtime environment and external tools")
    p_chk.set_defaults(func=cmd_check)

    return parser


def main(argv: Optional[Iterable[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(list(argv) if argv is not None else None)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="[%(levelname)s] %(message)s",
    )

    if not hasattr(args, "func"):
        parser.print_help(sys.stderr)
        return 1

    try:
        return args.func(args)
    except Exception as exc:
        logger.error("%s", exc)
        return 2


def phyview_main(argv: Optional[Iterable[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        prog="PhyView",
        description="Run one selected PhyView task on trimmed multifasta",
    )
    parser.add_argument("-i", "--input", required=True, help="Input trimmed multifasta (FASTA)")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory for the selected task")
    parser.add_argument("--ufboot", type=int, default=1000, help="UFBoot replicate number (used with --run_ml)")
    parser.add_argument("--threads", default="AUTO", help="Thread setting for IQ-TREE -T (used with --run_ml)")
    parser.add_argument("--model", default="MFP", help="Model option for IQ-TREE -m (used with --run_ml)")
    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument("--run_ml", "--run-ml", action="store_true", help="Run ML tree (IQ-TREE)")
    mode_group.add_argument("--run_nj", "--run-nj", action="store_true", help="Run pairwise distance + NJ tree")
    mode_group.add_argument("--run_popart", "--run-popart", action="store_true", help="Prepare PopART haplotype inputs")
    parser.add_argument(
        "--exec-popart",
        action="store_true",
        help="Execute PopART CLI after preparing input files (used with --run_popart)",
    )
    parser.add_argument("--popart-bin", default="popart", help="PopART executable name/path (used with --run_popart)")
    parser.add_argument(
        "--popart-args",
        nargs="*",
        default=[],
        help="Extra args passed to PopART CLI (used with --run_popart)",
    )
    parser.add_argument(
        "--unsafe",
        action="store_true",
        help="Disable IQ-TREE safe likelihood kernel (used with --run_ml; default uses -safe)",
    )
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    args = parser.parse_args(list(argv) if argv is not None else None)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="[%(levelname)s] %(message)s",
    )
    try:
        return cmd_phyview(args)
    except Exception as exc:
        logger.error("%s", exc)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
