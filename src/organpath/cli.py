import argparse
import gzip
import logging
import shutil
import subprocess
import sys
import importlib.util
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from . import __version__

logger = logging.getLogger(__name__)


MISSING_GT = {".", "./.", ".|.", "././.", ".|.|."}


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


def write_stats_table(stats: Dict[str, SampleStats], path: Path) -> None:
    with path.open("wt") as out:
        out.write("sample\ttotal_sites\tcalled_sites\tcoverage\tmean_depth\n")
        for sample, s in stats.items():
            out.write(
                f"{sample}\t{s.total_sites}\t{s.called_sites}\t"
                f"{s.coverage:.4f}\t{s.mean_depth:.4f}\n"
            )


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
    min_dp: int,
    min_gq: float,
    het_as_missing: bool,
) -> None:
    refs = parse_ref_fasta(ref_fasta)
    sample_seqs: Dict[str, Dict[str, List[str]]] = {
        s: {ctg: seq.copy() for ctg, seq in refs.items()} for s in keep_samples
    }

    sample_index: Dict[str, int] = {}

    with open_maybe_gzip(vcf_path, "rt") as inp, out_vcf.open("wt") as out:
        for raw in inp:
            if raw.startswith("##"):
                out.write(raw)
                continue

            if raw.startswith("#CHROM"):
                cols = raw.rstrip("\n").split("\t")
                all_samples = cols[9:]
                sample_index = {s: i for i, s in enumerate(all_samples)}
                selected = [s for s in all_samples if s in set(keep_samples)]
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

            selected_fields: List[str] = []
            for s in keep_samples:
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
                if len(ref) == 1:
                    base = choose_base_from_gt(gt, ref, alts)
                    sample_seqs[s][chrom][pos - 1] = base

                selected_fields.append(sf)

            out.write("\t".join(fields[:9] + selected_fields) + "\n")

    sample_fastas_dir.mkdir(parents=True, exist_ok=True)
    with multifasta_path.open("wt") as multi:
        for s in keep_samples:
            seq = "".join(
                "".join(sample_seqs[s][ctg]) for ctg in sorted(sample_seqs[s].keys())
            )
            indiv = sample_fastas_dir / f"{s}.fa"
            with indiv.open("wt") as out:
                out.write(f">{s}\n")
                out.write(f"{seq}\n")
            multi.write(f">{s}\n{seq}\n")


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
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    iqtree_bin = shutil.which("iqtree2") or shutil.which("iqtree")
    if not iqtree_bin:
        raise RuntimeError("Required tool not found: iqtree2 or iqtree")
    prefix = out_dir / "organpath_tree"
    run_command(
        [
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
    )


def cmd_run(args: argparse.Namespace) -> int:
    out_dir = Path(args.outdir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    vcf_path = Path(args.vcf).resolve()
    ref_fasta = Path(args.ref).resolve()

    _headers, samples = read_header_and_samples(vcf_path)
    logger.info("Detected %d samples in VCF.", len(samples))

    stats = compute_sample_stats(vcf_path, samples)
    stats_tsv = out_dir / "sample_stats.tsv"
    write_stats_table(stats, stats_tsv)

    keep, drop = keep_samples_by_threshold(
        stats,
        min_coverage=args.min_coverage,
        min_depth=args.min_mean_depth,
    )
    logger.info("Kept %d samples, dropped %d samples.", len(keep), len(drop))

    with (out_dir / "kept_samples.txt").open("wt") as f:
        f.write("\n".join(keep) + ("\n" if keep else ""))
    with (out_dir / "dropped_samples.txt").open("wt") as f:
        f.write("\n".join(drop) + ("\n" if drop else ""))

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
        min_dp=args.min_dp,
        min_gq=args.min_gq,
        het_as_missing=not args.keep_het,
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
        )
        logger.info("PhyView tree inference completed.")

    logger.info("All outputs are in: %s", out_dir)
    return 0


def cmd_phyview(args: argparse.Namespace) -> int:
    run_phyview(
        trimmed_fasta=Path(args.input).resolve(),
        out_dir=Path(args.outdir).resolve(),
        ufboot=args.ufboot,
        threads=args.threads,
        model=args.model,
    )
    logger.info("PhyView completed: %s", Path(args.outdir).resolve())
    return 0


def cmd_check(args: argparse.Namespace) -> int:
    missing_python: List[str] = []
    missing_tools: List[str] = []

    for mod in ["gzip", "argparse", "pathlib"]:
        if importlib.util.find_spec(mod) is None:
            missing_python.append(mod)

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

    p_run = subs.add_parser("run", help="Run full workflow from VCF to trimmed alignment")
    p_run.add_argument("-v", "--vcf", required=True, help="Input VCF or VCF.GZ")
    p_run.add_argument("-r", "--ref", required=True, help="Reference fasta")
    p_run.add_argument("-o", "--outdir", required=True, help="Output directory")
    p_run.add_argument("--min-coverage", type=float, default=0.5, help="Min sample coverage")
    p_run.add_argument(
        "--min-mean-depth",
        type=float,
        default=5.0,
        help="Drop sample only when coverage < threshold AND mean depth <= this value",
    )
    p_run.add_argument("--min-dp", type=int, default=3, help="Mark genotype missing if DP < this")
    p_run.add_argument("--min-gq", type=float, default=20.0, help="Mark genotype missing if GQ < this")
    p_run.add_argument("--keep-het", action="store_true", help="Keep heterozygous genotype (default: mask)")
    p_run.add_argument("--mafft-bin", default="mafft", help="Path or name of mafft executable")
    p_run.add_argument("--trimal-bin", default="trimal", help="Path or name of trimal executable")
    p_run.add_argument("--run-phyview", action="store_true", help="Also run OrganPath PhyView after trimming")
    p_run.add_argument("--ufboot", type=int, default=1000, help="UFBoot replicate number for iqtree")
    p_run.add_argument("--threads", default="AUTO", help="Thread setting passed to iqtree -T")
    p_run.add_argument("--model", default="MFP", help="Model option passed to iqtree -m")
    p_run.set_defaults(func=cmd_run)

    p_phy = subs.add_parser("PhyView", help="Build tree with IQ-TREE on trimmed multifasta")
    p_phy.add_argument("-i", "--input", required=True, help="Trimmed multifasta file")
    p_phy.add_argument("-o", "--outdir", required=True, help="Output directory for IQ-TREE")
    p_phy.add_argument("--ufboot", type=int, default=1000, help="UFBoot replicate number")
    p_phy.add_argument("--threads", default="AUTO", help="Thread setting passed to iqtree -T")
    p_phy.add_argument("--model", default="MFP", help="Model option passed to iqtree -m")
    p_phy.set_defaults(func=cmd_phyview)

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
