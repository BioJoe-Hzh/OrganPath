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
        if not args.skip_ml:
            run_phyview(
                trimmed_fasta=trimmed,
                out_dir=out_dir / "phyview",
                ufboot=args.ufboot,
                threads=args.threads,
                model=args.model,
                safe=not args.unsafe,
            )
        else:
            (out_dir / "phyview").mkdir(parents=True, exist_ok=True)
            logger.info("Skipping IQ-TREE ML step (--skip-ml).")
        seqs = read_fasta_sequences(trimmed)
        names, mat = build_distance_matrix(seqs)
        write_distance_matrix_tsv(names, mat, out_dir / "phyview" / "pairwise_distance.tsv")
        if not args.no_nj:
            build_nj_tree(names, mat, out_dir / "phyview" / "nj_tree.nwk")
        if not args.no_popart:
            _hap_tsv, popart_nex = prepare_popart_inputs(trimmed, out_dir / "phyview")
            if args.run_popart:
                run_popart_cli(args.popart_bin, popart_nex, out_dir / "phyview", args.popart_args)
        logger.info("PhyView tree inference completed.")

    logger.info("All outputs are in: %s", out_dir)
    return 0


def cmd_phyview(args: argparse.Namespace) -> int:
    trimmed = Path(args.input).resolve()
    out_dir = Path(args.outdir).resolve()
    if not args.skip_ml:
        run_phyview(
            trimmed_fasta=trimmed,
            out_dir=out_dir,
            ufboot=args.ufboot,
            threads=args.threads,
            model=args.model,
            safe=not args.unsafe,
        )
    else:
        out_dir.mkdir(parents=True, exist_ok=True)
        logger.info("Skipping IQ-TREE ML step (--skip-ml).")
    seqs = read_fasta_sequences(trimmed)
    names, mat = build_distance_matrix(seqs)
    write_distance_matrix_tsv(names, mat, out_dir / "pairwise_distance.tsv")
    if not args.no_nj:
        build_nj_tree(names, mat, out_dir / "nj_tree.nwk")
    if not args.no_popart:
        _hap_tsv, popart_nex = prepare_popart_inputs(trimmed, out_dir)
        if args.run_popart:
            run_popart_cli(args.popart_bin, popart_nex, out_dir, args.popart_args)
    logger.info("PhyView completed: %s", Path(args.outdir).resolve())
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
    p_run.add_argument("--skip-ml", action="store_true", help="Skip IQ-TREE ML tree step in PhyView stage")
    p_run.add_argument("--no-nj", action="store_true", help="Skip NJ tree construction in PhyView stage")
    p_run.add_argument("--no-popart", action="store_true", help="Skip PopART input generation in PhyView stage")
    p_run.add_argument("--run-popart", action="store_true", help="Execute PopART CLI after preparing inputs")
    p_run.add_argument("--popart-bin", default="popart", help="PopART executable name/path")
    p_run.add_argument("--popart-args", nargs="*", default=[], help="Extra args passed to PopART CLI")
    p_run.add_argument(
        "--unsafe",
        action="store_true",
        help="Disable IQ-TREE safe likelihood kernel (default uses -safe)",
    )
    p_run.set_defaults(func=cmd_run)

    p_phy = subs.add_parser("PhyView", help="Build tree with IQ-TREE on trimmed multifasta")
    p_phy.add_argument("-i", "--input", required=True, help="Trimmed multifasta file")
    p_phy.add_argument("-o", "--outdir", required=True, help="Output directory for IQ-TREE")
    p_phy.add_argument("--ufboot", type=int, default=1000, help="UFBoot replicate number")
    p_phy.add_argument("--threads", default="AUTO", help="Thread setting passed to iqtree -T")
    p_phy.add_argument("--model", default="MFP", help="Model option passed to iqtree -m")
    p_phy.add_argument("--skip-ml", action="store_true", help="Skip IQ-TREE ML tree step")
    p_phy.add_argument("--no-nj", action="store_true", help="Skip NJ tree construction")
    p_phy.add_argument("--no-popart", action="store_true", help="Skip PopART input generation")
    p_phy.add_argument("--run-popart", action="store_true", help="Execute PopART CLI after preparing inputs")
    p_phy.add_argument("--popart-bin", default="popart", help="PopART executable name/path")
    p_phy.add_argument("--popart-args", nargs="*", default=[], help="Extra args passed to PopART CLI")
    p_phy.add_argument(
        "--unsafe",
        action="store_true",
        help="Disable IQ-TREE safe likelihood kernel (default uses -safe)",
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
        description="Build tree with IQ-TREE on trimmed multifasta",
    )
    parser.add_argument("-i", "--input", required=True, help="Trimmed multifasta file")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory for IQ-TREE")
    parser.add_argument("--ufboot", type=int, default=1000, help="UFBoot replicate number")
    parser.add_argument("--threads", default="AUTO", help="Thread setting passed to iqtree -T")
    parser.add_argument("--model", default="MFP", help="Model option passed to iqtree -m")
    parser.add_argument("--skip-ml", action="store_true", help="Skip IQ-TREE ML tree step")
    parser.add_argument("--no-nj", action="store_true", help="Skip NJ tree construction")
    parser.add_argument("--no-popart", action="store_true", help="Skip PopART input generation")
    parser.add_argument("--run-popart", action="store_true", help="Execute PopART CLI after preparing inputs")
    parser.add_argument("--popart-bin", default="popart", help="PopART executable name/path")
    parser.add_argument("--popart-args", nargs="*", default=[], help="Extra args passed to PopART CLI")
    parser.add_argument(
        "--unsafe",
        action="store_true",
        help="Disable IQ-TREE safe likelihood kernel (default uses -safe)",
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
