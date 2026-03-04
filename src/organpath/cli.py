import argparse
import concurrent.futures
import gzip
import logging
import shutil
import subprocess
import sys
import importlib.util
import re
from xml.etree import ElementTree as ET
import xml.dom.minidom as minidom
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
ORG_SORT_DEFAULTS = {
    "generic": {"min_identity": 0.95, "min_len": 1000, "gap_n": 100},
    "plant_pt": {"min_identity": 0.95, "min_len": 1000, "gap_n": 100},
    "plant_mt": {"min_identity": 0.95, "min_len": 3000, "gap_n": 100},
    "animal_mt": {"min_identity": 0.95, "min_len": 1000, "gap_n": 100},
}
DEFAULT_BEAST_TEMPLATE = (
    "/Users/zh384/Desktop/scripts_dev/pipelines/eDNA_phylogeny/Organelle_pipelineV2/"
    "2_Construct_phylogentic_framework/auxiliary_scripts/ultrametric_tree_template.xml"
)


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


def non_n_length(seq: str) -> int:
    return sum(1 for b in seq.upper() if b != "N")


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
            "snp-sites": "Please install snp-sites and ensure `snp-sites` is in PATH.",
            "bcftools": "Please install bcftools and ensure `bcftools` is in PATH.",
            "beast": "Please install BEAST2 and ensure `beast` is in PATH.",
            "treeannotator": "Please install BEAST2 tools and ensure `treeannotator` is in PATH.",
        }.get(name, f"Please install `{name}` and add it to PATH.")
        raise RuntimeError(f"Required tool not found in PATH: {name}. {install_hint}")
    return p


def run_alignment_and_trimming(
    multifasta: Path, aligned: Path, trimmed: Path, mafft_bin: str, trimal_bin: str
) -> None:
    run_command([mafft_bin, "--auto", str(multifasta)], stdout_path=aligned)
    run_command([trimal_bin, "-automated1", "-in", str(aligned), "-out", str(trimmed)])


def run_alignment_with_direction(
    multifasta: Path, aligned: Path, mafft_bin: str, adjust_direction: bool, threads: str = "AUTO"
) -> None:
    cmd = [mafft_bin, "--auto"]
    t = str(threads).strip().upper()
    if t == "AUTO":
        cmd += ["--thread", "-1"]
    else:
        cmd += ["--thread", str(threads)]
    if adjust_direction:
        cmd.append("--adjustdirection")
    cmd.append(str(multifasta))
    run_command(cmd, stdout_path=aligned)


def write_ref_from_msa(msa_fa: Path, out_ref: Path, ref_id: Optional[str]) -> str:
    seqs = read_fasta_sequences(msa_fa)
    if ref_id:
        if ref_id not in seqs:
            raise ValueError(f"--ref-id not found in MSA: {ref_id}")
        chosen = ref_id
    else:
        chosen = next(iter(seqs.keys()))
    with out_ref.open("wt") as out:
        out.write(f">{chosen}\n{seqs[chosen]}\n")
    return chosen


def _strip_beast_annotations(newick: str, include_posterior: bool) -> str:
    if include_posterior:
        def repl(m: re.Match) -> str:
            payload = m.group(1)
            pm = re.search(r"posterior=([0-9eE.+-]+)", payload)
            return pm.group(1) if pm else ""
        return re.sub(r"\[&([^\]]+)\]", repl, newick)
    return re.sub(r"\[&[^\]]+\]", "", newick)


def nexus_tree_to_newick(nexus_path: Path, out_newick: Path, include_posterior: bool) -> None:
    tree_txt: Optional[str] = None
    pat = re.compile(r"^\s*tree\s+.+?=\s*(.+)$", flags=re.IGNORECASE)
    with nexus_path.open("rt") as fh:
        for raw in fh:
            m = pat.match(raw.strip())
            if m:
                tree_txt = m.group(1).strip()
    if not tree_txt:
        raise ValueError(f"No tree line found in NEXUS: {nexus_path}")
    cleaned = _strip_beast_annotations(tree_txt, include_posterior=include_posterior)
    if not cleaned.endswith(";"):
        cleaned += ";"
    out_newick.write_text(cleaned + "\n")


def build_beast_xml_from_template(
    template_xml: Path,
    aligned_fasta: Path,
    output_xml: Path,
    prefix: str,
    chain_length: int,
    store_every: int,
    old_prefix: str = "ultrametric",
) -> None:
    xml_content = template_xml.read_text(encoding="utf-8")
    m = re.match(r'<\?xml.*?\?>', xml_content)
    xml_decl = m.group(0) if m else '<?xml version="1.0" encoding="UTF-8" standalone="no"?>'
    tree = ET.ElementTree(ET.fromstring(xml_content))
    root = tree.getroot()

    for elem in root.iter():
        for k, v in list(elem.attrib.items()):
            if old_prefix in v:
                elem.set(k, re.sub(rf'([:@.]?){re.escape(old_prefix)}', rf"\1{prefix}", v))

    run_elem = root.find(".//run[@id='mcmc']")
    if run_elem is not None:
        run_elem.set("chainLength", str(chain_length))
        run_elem.set("storeEvery", str(store_every))

    data_block = root.find(f".//data[@id='{prefix}']")
    if data_block is None:
        data_block = root.find(".//data")
    if data_block is None:
        raise ValueError("BEAST template has no <data> block.")

    for child in list(data_block):
        data_block.remove(child)
    seqs = read_fasta_sequences(aligned_fasta)
    for sid, seq in seqs.items():
        seq_elem = ET.SubElement(data_block, "sequence")
        seq_elem.set("id", f"seq_{sid}")
        seq_elem.set("spec", "Sequence")
        seq_elem.set("taxon", sid)
        seq_elem.set("totalcount", "4")
        seq_elem.set("value", seq)

    xml_str = ET.tostring(tree.getroot(), encoding="UTF-8")
    parsed = minidom.parseString(xml_str)
    pretty = parsed.toprettyxml(indent="  ")
    pretty = pretty.split("\n", 1)[1] if pretty.startswith("<?xml") else pretty
    output_xml.write_text(f"{xml_decl}\n{pretty}", encoding="utf-8")


def run_beast_pipeline(
    aligned_fasta: Path,
    out_dir: Path,
    template_xml: Path,
    prefix: str,
    chain_length: int,
    store_every: int,
    beast_bin: str,
    beast_threads: int,
    use_beagle: bool,
    treeannotator_bin: str,
    run_treeannotator: bool,
    burnin: int,
    heights: str,
    include_posterior: bool,
    old_prefix: str = "ultrametric",
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    xml_out = out_dir / f"{prefix}.beast.xml"
    build_beast_xml_from_template(
        template_xml=template_xml,
        aligned_fasta=aligned_fasta,
        output_xml=xml_out,
        prefix=prefix,
        chain_length=chain_length,
        store_every=store_every,
        old_prefix=old_prefix,
    )

    beast = ensure_tool(beast_bin)
    beast_cmd = [beast]
    if use_beagle:
        beast_cmd.append("-beagle")
    beast_cmd += [
        "-threads",
        str(beast_threads),
        "-prefix",
        str(out_dir / prefix),
        "-overwrite",
        str(xml_out),
    ]
    run_command(beast_cmd)

    tree_file = out_dir / f"{prefix}.trees"
    if not tree_file.exists():
        logger.warning("BEAST run finished, but tree file not found: %s", tree_file)
        return

    if run_treeannotator:
        ta = ensure_tool(treeannotator_bin)
        mcc_nexus = out_dir / f"{prefix}.mcc.nexus"
        run_command(
            [
                ta,
                "-burnin",
                str(burnin),
                "-heights",
                heights,
                str(tree_file),
                str(mcc_nexus),
            ]
        )
        if mcc_nexus.exists():
            nexus_tree_to_newick(
                mcc_nexus,
                out_dir / f"{prefix}.mcc.nwk",
                include_posterior=include_posterior,
            )


def resolve_beast_template_arg(template_arg: str) -> Path:
    if template_arg.strip().lower() == "default":
        return Path(DEFAULT_BEAST_TEMPLATE).resolve()
    return Path(template_arg).resolve()


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


def write_fasta_sequences(path: Path, seqs: Dict[str, str]) -> None:
    with path.open("wt") as out:
        for sid, seq in seqs.items():
            out.write(f">{sid}\n{seq}\n")


def filter_alignment_sites(
    in_fa: Path,
    out_fa: Path,
    max_missing_frac: float = 1.0,
    snp_only: bool = False,
) -> Tuple[int, int]:
    seqs = read_fasta_sequences(in_fa)
    names = list(seqs.keys())
    lens = {len(v) for v in seqs.values()}
    if len(lens) != 1:
        raise ValueError(f"Alignment has inconsistent sequence lengths: {in_fa}")
    n = len(names)
    L = next(iter(lens))
    keep_cols: List[int] = []

    for i in range(L):
        col = [seqs[s][i].upper() for s in names]
        called = [b for b in col if b in {"A", "C", "G", "T"}]
        missing_frac = 1.0 - (len(called) / n)
        if missing_frac > max_missing_frac:
            continue
        if snp_only and len(set(called)) < 2:
            continue
        keep_cols.append(i)

    if not keep_cols:
        raise ValueError(
            f"No alignment columns left after filtering (max_missing_frac={max_missing_frac}, snp_only={snp_only})"
        )

    filtered = {s: "".join(seqs[s][i] for i in keep_cols) for s in names}
    write_fasta_sequences(out_fa, filtered)
    return len(keep_cols), L


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


def read_primary_fasta_sequence(path: Path) -> str:
    seqs = read_fasta_sequences(path)
    if not seqs:
        raise ValueError(f"No sequences found in fasta: {path}")
    return next(iter(seqs.values())).upper()


def rotate_sequence_to_seed_start(seq: str, seed_seq: str, k: int = 80) -> Tuple[str, str]:
    if not seq or not seed_seq:
        return seq, "no-seq"
    seed_u = seed_seq.upper()
    s_u = seq.upper()
    rc_u = reverse_complement(s_u)

    k = max(20, min(k, len(seed_u)))
    anchor = seed_u[:k]
    while "N" in anchor and k > 20:
        k -= 10
        anchor = seed_u[:k]

    def _find_and_rotate(q: str) -> Optional[str]:
        idx = q.find(anchor)
        if idx < 0:
            return None
        return q[idx:] + q[:idx]

    fwd = _find_and_rotate(s_u)
    rev = _find_and_rotate(rc_u)
    if fwd is not None and rev is not None:
        return (fwd, "forward") if fwd.count("N") <= rev.count("N") else (rev, "reverse-complement")
    if fwd is not None:
        return fwd, "forward"
    if rev is not None:
        return rev, "reverse-complement"
    return s_u, "anchor-not-found"


def _normalize_region_name(token: str) -> Optional[str]:
    t = token.strip().lower().replace("-", "").replace("_", "")
    if not t:
        return None
    if "lsc" in t:
        return "LSC"
    if "ssc" in t:
        return "SSC"
    if "irb" in t or t == "ir2":
        return "IRB"
    if "ira" in t or t == "ir1":
        return "IRA"
    if t == "ir":
        return "IR"
    return None


def parse_cp_regions_file(path: Path) -> Dict[str, Tuple[int, int]]:
    regions: Dict[str, Tuple[int, int]] = {}
    with path.open("rt") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            cols = re.split(r"\s+", line)
            name = None
            for c in cols:
                nm = _normalize_region_name(c)
                if nm:
                    name = nm
                    break
            if name is None:
                continue
            ints: List[int] = []
            for c in cols:
                try:
                    ints.append(int(float(c)))
                except ValueError:
                    continue
            if len(ints) < 2:
                continue
            regions[name] = (ints[0], ints[1])
    if "IR" in regions:
        regions.setdefault("IRB", regions["IR"])
    need = {"LSC", "SSC"}
    if not need.issubset(set(regions.keys())):
        raise ValueError(f"cp regions file missing required labels (LSC/SSC): {path}")
    if not any(x in regions for x in ("IRA", "IRB", "IR")):
        raise ValueError(f"cp regions file missing IR label: {path}")
    return regions


def parse_cpstools_four_sec(path: Path) -> Dict[str, Tuple[int, int]]:
    txt = path.read_text()
    regions: Dict[str, Tuple[int, int]] = {}
    # Primary parser: tolerate different separators/spaces in cpstools IR outputs.
    for m in re.finditer(r"(?i)\b(LSC|SSC|IRA|IRB|IR1|IR2|IR)\b[^0-9]{0,20}(\d+)[^0-9]+(\d+)", txt):
        raw_name = m.group(1)
        nm = _normalize_region_name(raw_name)
        if nm is None:
            continue
        s = int(m.group(2))
        e = int(m.group(3))
        regions[nm] = (s, e)
    # Fallback token parser for unusual output layouts.
    if "LSC" not in regions or "SSC" not in regions:
        toks = re.split(r"[\s,;]+", txt)
        for t in toks:
            m = re.search(r"(?i)(LSC|SSC|IRA|IRB|IR1|IR2|IR)", t)
            if not m:
                continue
            nm = _normalize_region_name(m.group(1))
            if nm is None:
                continue
            nums = re.findall(r"\d+", t)
            if len(nums) >= 2:
                regions[nm] = (int(nums[0]), int(nums[1]))
    if "IR" in regions:
        regions.setdefault("IRB", regions["IR"])
    if "IR1" in regions:
        regions.setdefault("IRA", regions["IR1"])
    if "IR2" in regions:
        regions.setdefault("IRB", regions["IR2"])
    if "LSC" not in regions or "SSC" not in regions:
        raise ValueError(f"Failed to parse LSC/SSC from cpstools IR output: {path}")
    if not any(x in regions for x in ("IRA", "IRB", "IR")):
        raise ValueError(f"Failed to parse IR from cpstools IR output: {path}")
    return regions


def extract_circular_region(seq: str, start: int, end: int) -> str:
    if not seq:
        return ""
    n = len(seq)
    s = (start - 1) % n
    e = (end - 1) % n
    if s <= e:
        return seq[s : e + 1]
    return seq[s:] + seq[: e + 1]


def reorder_cp_single_ir(seq: str, cp_regions: Dict[str, Tuple[int, int]], keep_ir: str) -> Tuple[str, str]:
    parts, ir_name = reorder_cp_single_ir_parts(seq=seq, cp_regions=cp_regions, keep_ir=keep_ir)
    return parts["LSC"] + parts["IR"] + parts["SSC"], ir_name


def reorder_cp_single_ir_parts(seq: str, cp_regions: Dict[str, Tuple[int, int]], keep_ir: str) -> Tuple[Dict[str, str], str]:
    lsc = extract_circular_region(seq, *cp_regions["LSC"])
    ssc = extract_circular_region(seq, *cp_regions["SSC"])
    keep = keep_ir.lower()
    if keep == "ira" and "IRA" in cp_regions:
        ir_name = "IRA"
    elif keep == "irb" and "IRB" in cp_regions:
        ir_name = "IRB"
    else:
        ir_name = "IRB" if "IRB" in cp_regions else ("IRA" if "IRA" in cp_regions else "IR")
    ir = extract_circular_region(seq, *cp_regions[ir_name])
    return {"LSC": lsc, "IR": ir, "SSC": ssc}, ir_name


def region_len(start: int, end: int, total_len: int) -> int:
    s = (start - 1) % total_len
    e = (end - 1) % total_len
    return (e - s + 1) if s <= e else (total_len - s + e + 1)


def point_in_region(pos1: int, start: int, end: int, total_len: int) -> bool:
    p = (pos1 - 1) % total_len
    s = (start - 1) % total_len
    e = (end - 1) % total_len
    if s <= e:
        return s <= p <= e
    return p >= s or p <= e


def choose_ir_name(cp_regions: Dict[str, Tuple[int, int]], keep_ir: str) -> str:
    keep = keep_ir.lower()
    if keep == "ira" and "IRA" in cp_regions:
        return "IRA"
    if keep == "irb" and "IRB" in cp_regions:
        return "IRB"
    if "IRB" in cp_regions:
        return "IRB"
    if "IRA" in cp_regions:
        return "IRA"
    return "IR"


def _count_kmer_matches(seq: str, ref_kmers: set, k: int) -> int:
    if len(seq) < k:
        return 0
    c = 0
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        if "N" in kmer:
            continue
        if kmer in ref_kmers:
            c += 1
    return c


def orient_partition_to_seed(seq: str, seed_ref: str, k: int = 31) -> Tuple[str, str]:
    seq_u = seq.upper()
    ref_u = seed_ref.upper()
    if not seq_u or not ref_u:
        return seq_u, "NA"
    # Skip orientation for near-empty/gap-only sequences.
    non_n = sum(1 for c in seq_u if c in {"A", "C", "G", "T"})
    if non_n < max(200, k * 3):
        return seq_u, "LOW_INFO"

    ref_kmers = set()
    if len(ref_u) >= k:
        for i in range(len(ref_u) - k + 1):
            km = ref_u[i : i + k]
            if "N" not in km:
                ref_kmers.add(km)
    if not ref_kmers:
        return seq_u, "LOW_INFO"

    fwd = _count_kmer_matches(seq_u, ref_kmers, k)
    rc_seq = reverse_complement(seq_u)
    rev = _count_kmer_matches(rc_seq, ref_kmers, k)
    if rev > fwd:
        return rc_seq, "RC"
    return seq_u, "FWD"


def cp_regions_from_sequence_with_cpstools(seq: str, sample_out: Path, sample_name: str, cpstools_bin: str) -> Dict[str, Tuple[int, int]]:
    cpstools = shutil.which(cpstools_bin)
    if not cpstools:
        raise RuntimeError(f"cpstools executable not found: {cpstools_bin}")
    tmp_fa = sample_out / f"{sample_name}.cpstools_input.fa"
    with tmp_fa.open("wt") as out:
        out.write(f">{sample_name}\n{seq}\n")
    four_sec = sample_out / f"{sample_name}.cpstools_four_sec.txt"
    with four_sec.open("wt") as out:
        subprocess.run([cpstools, "IR", "-i", str(tmp_fa)], stdout=out, check=True)
    return parse_cpstools_four_sec(four_sec)


def resolve_cp_regions(
    seed_fa: Path,
    out_dir: Path,
    cp_regions_path: Optional[str],
    cpstools_bin: str,
    cpstools_args: List[str],
) -> Dict[str, Tuple[int, int]]:
    if cp_regions_path:
        return parse_cp_regions_file(Path(cp_regions_path).resolve())

    cpstools = shutil.which(cpstools_bin)
    if not cpstools:
        raise RuntimeError(f"cpstools executable not found: {cpstools_bin}")

    cp_out = out_dir / "cpstools"
    cp_out.mkdir(parents=True, exist_ok=True)
    cp_regions_tsv = cp_out / "cp_regions.tsv"

    if cpstools_args:
        rendered: List[str] = []
        for x in cpstools_args:
            rendered.append(
                x.replace("{seed_fasta}", str(seed_fa))
                .replace("{cpstools_out}", str(cp_out))
                .replace("{cp_regions_tsv}", str(cp_regions_tsv))
            )
        run_command([cpstools] + rendered)
        if not cp_regions_tsv.exists():
            raise FileNotFoundError(
                f"cpstools finished but expected region file not found: {cp_regions_tsv}. "
                "Please ensure --cpstools-args writes regions to {cp_regions_tsv}."
            )
        return parse_cp_regions_file(cp_regions_tsv)

    # Default cpstools workflow: IR(seed) -> Seq -m LSC -> IR(LSC_adjusted)
    four_sec_1 = cp_out / "seed.four_sec.txt"
    with four_sec_1.open("wt") as out:
        subprocess.run([cpstools, "IR", "-i", str(seed_fa)], stdout=out, check=True)
    seed_regions: Optional[Dict[str, Tuple[int, int]]] = None
    try:
        seed_regions = parse_cpstools_four_sec(four_sec_1)
    except Exception:
        seed_regions = None

    tmp_dir = cp_out / "cpstools_tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    link_path = tmp_dir / seed_fa.name
    try:
        if not link_path.exists():
            link_path.symlink_to(seed_fa)
    except Exception:
        shutil.copy2(seed_fa, link_path)

    # cpstools Seq fails if the target output folder already exists.
    lsc_adj_dir = cp_out / "LSC_adj"
    if lsc_adj_dir.exists():
        shutil.rmtree(lsc_adj_dir, ignore_errors=True)

    regions: Optional[Dict[str, Tuple[int, int]]] = None
    try:
        subprocess.run(
            [cpstools, "Seq", "-d", str(tmp_dir), "-f", str(four_sec_1), "-m", "LSC"],
            check=True,
            cwd=str(cp_out),
        )
        lsc_candidates = sorted(lsc_adj_dir.glob("*.fa")) + sorted(lsc_adj_dir.glob("*.fasta"))
        if not lsc_candidates:
            raise FileNotFoundError(f"cpstools Seq did not produce LSC-adjusted fasta in: {lsc_adj_dir}")
        lsc_seed = lsc_candidates[0]

        four_sec_2 = cp_out / "lsc_start.four_sec.txt"
        with four_sec_2.open("wt") as out:
            subprocess.run([cpstools, "IR", "-i", str(lsc_seed)], stdout=out, check=True)
        regions = parse_cpstools_four_sec(four_sec_2)
    except Exception as exc:
        if seed_regions is None:
            raise RuntimeError(
                "cpstools default workflow failed and fallback parsing of `cpstools IR` output also failed. "
                "Please provide --cp-regions."
            ) from exc
        logger.warning(
            "cpstools Seq workflow failed (%s). Falling back to seed IR coordinates from: %s",
            str(exc),
            four_sec_1,
        )
        regions = seed_regions

    with cp_regions_tsv.open("wt") as out:
        out.write("region\tstart\tend\n")
        for k in ("LSC", "IRB", "IRA", "SSC"):
            if k in regions:
                s, e = regions[k]
                out.write(f"{k}\t{s}\t{e}\n")
    return regions


def discover_sample_candidates(sample_dir: Path) -> Tuple[List[Path], Optional[Path], str]:
    # Prefer GetOrganelle path sequences (best circular candidates).
    path_cands: List[Path] = []
    for pat in (
        "*path_sequence*.fasta",
        "*path_sequence*.fa",
        "*scaffolds.graph*.fasta",
        "*scaffolds.graph*.fa",
    ):
        path_cands.extend(sorted(sample_dir.rglob(pat)))
    # de-dup while preserving order
    seen = set()
    uniq_path_cands = []
    for p in path_cands:
        sp = str(p)
        if sp in seen:
            continue
        seen.add(sp)
        uniq_path_cands.append(p)
    fastg_cands = sorted(sample_dir.rglob("*.fastg"))
    fastg = fastg_cands[0] if fastg_cands else None
    if uniq_path_cands:
        return uniq_path_cands, fastg, "path_sequence"

    # Fallback to SPAdes/GetOrganelle contigs.
    contig_cands: List[Path] = []
    for pat in ("*contigs.fasta", "*contigs.fa", "contigs.fasta", "contigs.fa"):
        contig_cands.extend(sorted(sample_dir.rglob(pat)))
    seen = set()
    uniq_contigs = []
    for p in contig_cands:
        sp = str(p)
        if sp in seen:
            continue
        seen.add(sp)
        uniq_contigs.append(p)
    return uniq_contigs, fastg, "contigs"


def parse_minimap2_paf(
    paf_path: Path, min_identity: float, min_len: int
) -> List[Tuple[str, int, int, str, float, int, int, int]]:
    hits: List[Tuple[str, int, int, str, float, int, int, int]] = []
    with paf_path.open("rt") as fh:
        for raw in fh:
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            qname = parts[0]
            qstart = int(parts[2])
            qend = int(parts[3])
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
            hits.append((qname, tstart, tend, strand, ident, alen, qstart, qend))
    return hits


def parse_blast_tab(
    tab_path: Path, min_identity: float, min_len: int
) -> List[Tuple[str, int, int, str, float, int, int, int]]:
    hits: List[Tuple[str, int, int, str, float, int, int, int]] = []
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
            qstart = int(parts[6])
            qend = int(parts[7])
            strand = "+" if sstart <= send else "-"
            if pident < min_identity or alen < min_len:
                continue
            hits.append((qname, min(sstart, send), max(sstart, send), strand, pident, alen, qstart, qend))
    return hits


def choose_best_hits(
    hits: List[Tuple[str, int, int, str, float, int, int, int]]
) -> List[Tuple[str, int, int, str, float, int, int, int]]:
    best: Dict[str, Tuple[str, int, int, str, float, int, int, int]] = {}
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


def map_hits_for_candidate(
    contig_fa: Path,
    seed_fa: Path,
    min_identity: float,
    min_len: int,
    aligner: str,
    tmp_prefix: Path,
) -> Tuple[List[Tuple[str, int, int, str, float, int, int, int]], str]:
    tool = aligner
    if tool == "auto":
        tool = "minimap2" if shutil.which("minimap2") else "blastn"

    tmp = tmp_prefix.with_suffix(".align.tmp")
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
    return choose_best_hits(hits), tool


def canonical_fasta_key(path: Path) -> str:
    seqs = read_fasta_sequences(path)
    if not seqs:
        return ""
    seq = "".join(seqs.values()).upper()
    rc = reverse_complement(seq)
    return seq if seq <= rc else rc


def select_best_candidate_fasta(
    candidates: List[Path],
    seed_fa: Path,
    sample_name: str,
    sample_out: Path,
    min_identity: float,
    min_len: int,
    aligner: str,
    expected_len: Optional[int],
) -> Tuple[Path, int]:
    unique: List[Path] = []
    seen_keys = set()
    for p in candidates:
        key = canonical_fasta_key(p)
        if key and key in seen_keys:
            continue
        if key:
            seen_keys.add(key)
        unique.append(p)
    if not unique:
        raise ValueError("No usable candidate fasta files found.")
    if len(unique) == 1:
        return unique[0], 1

    best_path: Optional[Path] = None
    best_score: Optional[Tuple[float, float, int, int]] = None
    for i, cand in enumerate(unique, start=1):
        try:
            chosen, _tool = map_hits_for_candidate(
                contig_fa=cand,
                seed_fa=seed_fa,
                min_identity=min_identity,
                min_len=min_len,
                aligner=aligner,
                tmp_prefix=sample_out / f"{sample_name}.cand{i}",
            )
        except Exception:
            chosen = []
        if not chosen:
            score = (0.0, 0.0, 0, -10**9)
        else:
            best_alen = max(h[5] for h in chosen)
            best_ident = max(h[4] for h in chosen)
            total_alen = sum(h[5] for h in chosen)
            cseq_len = sum(len(s) for s in read_fasta_sequences(cand).values())
            closeness = 0
            if expected_len and expected_len > 0:
                closeness = -abs(cseq_len - expected_len)
            score = (float(best_alen), float(best_ident), int(total_alen), int(closeness))
        if best_score is None or score > best_score:
            best_score = score
            best_path = cand
    if best_path is None:
        return unique[0], len(unique)
    return best_path, len(unique)


def slice_query_segment(seq: str, qstart: int, qend: int, one_based: bool) -> str:
    if not seq:
        return seq
    n = len(seq)
    if one_based:
        s = max(1, min(qstart, qend)) - 1
        e = min(n, max(qstart, qend))
    else:
        s = max(0, min(qstart, qend))
        e = min(n, max(qstart, qend))
    if e <= s:
        return seq
    return seq[s:e]


def build_sample_assembly_from_contigs(
    contig_fa: Path,
    seed_fa: Path,
    sample_name: str,
    out_dir: Path,
    min_identity: float,
    min_len: int,
    gap_n: int,
    aligner: str,
    organelle_mode: str = "generic",
    seed_seq: Optional[str] = None,
    pt_single_ir: bool = False,
    cp_regions: Optional[Dict[str, Tuple[int, int]]] = None,
    pt_keep_ir: str = "auto",
    cpstools_bin: str = "cpstools",
    pt_fragment_min_len: int = 1000,
    pt_complete_min_frac: float = 0.85,
    seed_len: int = 0,
    seed_parts: Optional[Dict[str, str]] = None,
) -> Tuple[str, int, int, str, Optional[Dict[str, str]]]:
    out_dir.mkdir(parents=True, exist_ok=True)
    contigs = read_fasta_sequences(contig_fa)
    if not contigs:
        raise ValueError(f"No contigs found in {contig_fa}")

    chosen, tool = map_hits_for_candidate(
        contig_fa=contig_fa,
        seed_fa=seed_fa,
        min_identity=min_identity,
        min_len=min_len,
        aligner=aligner,
        tmp_prefix=out_dir / sample_name,
    )
    # path_sequence fasta often contains multiple alternative graph paths.
    # Keep only the best-scoring path to avoid concatenating alternative assemblies.
    if "path_sequence" in contig_fa.name.lower() and len(chosen) > 1:
        chosen = [max(chosen, key=lambda x: (x[5], x[4]))]
    if not chosen:
        raise ValueError("No contigs passed identity/length filter against seed.")

    # Plant chloroplast specialized flow:
    # 1) complete-like: single long contig, use cpstools on sample itself for LSC/IR/SSC;
    # 2) fragmented: partition contigs by reference regions and stitch with N gaps.
    if organelle_mode == "plant_pt" and pt_single_ir and cp_regions:
        ir_name = choose_ir_name(cp_regions, pt_keep_ir)
        ref_total = max(seed_len, len(seed_seq) if seed_seq else 1)
        expected_lsc = region_len(*cp_regions["LSC"], total_len=ref_total)
        expected_ssc = region_len(*cp_regions["SSC"], total_len=ref_total)
        expected_ir = region_len(*cp_regions[ir_name], total_len=ref_total)
        expected_one_ir = expected_lsc + expected_ir + expected_ssc

        complete_seq: Optional[str] = None
        complete_parts: Optional[Dict[str, str]] = None
        complete_note = ""
        if len(chosen) == 1:
            cid, _sstart, _send, strand, _ident, _alen, _qstart, _qend = chosen[0]
            cseq = contigs.get(cid, "")
            if len(cseq) >= int(expected_one_ir * pt_complete_min_frac):
                try:
                    samp_regions = cp_regions_from_sequence_with_cpstools(
                        cseq, sample_out=out_dir, sample_name=sample_name, cpstools_bin=cpstools_bin
                    )
                    complete_parts, kept_ir = reorder_cp_single_ir_parts(cseq, cp_regions=samp_regions, keep_ir=pt_keep_ir)
                    orient_notes: List[str] = []
                    if seed_parts:
                        for pname in ("LSC", "SSC"):
                            if pname in complete_parts and pname in seed_parts:
                                oseq, odir = orient_partition_to_seed(complete_parts[pname], seed_parts[pname])
                                complete_parts[pname] = oseq
                                orient_notes.append(f"{pname}:{odir}")
                    complete_seq = complete_parts["LSC"] + complete_parts["IR"] + complete_parts["SSC"]
                    complete_note = f"type:complete;single_ir:{kept_ir}"
                    if orient_notes:
                        complete_note += ";part_orient:" + ",".join(orient_notes)
                except Exception as exc:
                    complete_note = f"type:fragmented_fallback;cpstools_sample_fail:{str(exc).replace(chr(9), ' ')}"

        if complete_seq is not None:
            miss_bp = max(expected_one_ir - len(complete_seq), 0)
            note = f"{complete_note};missing_bp:{miss_bp};expected_len:{expected_one_ir}"
            out_fa = out_dir / f"{sample_name}.organellar.fasta"
            with out_fa.open("wt") as out:
                out.write(f">{sample_name}\n{complete_seq}\n")
            with (out_dir / f"{sample_name}.selected_contigs.tsv").open("wt") as tab:
                tab.write("contig\tseed_start\tseed_end\tstrand\tidentity\taln_len\tregion\n")
                h = chosen[0]
                tab.write(f"{h[0]}\t{h[1]}\t{h[2]}\t{h[3]}\t{h[4]:.4f}\t{h[5]}\tcomplete\n")
            return str(out_fa), 1, len(complete_seq), note, complete_parts

        # Fragmented route
        region_order = ["LSC", ir_name, "SSC"]
        region_hits: Dict[str, List[Tuple[str, int, int, str, float, int, int, int]]] = {r: [] for r in region_order}
        for h in chosen:
            cid, sstart, send, strand, ident, alen, qstart, qend = h
            if alen < pt_fragment_min_len:
                continue
            mid = (sstart + send) // 2
            assigned = None
            for r in region_order:
                rs, re = cp_regions[r]
                if point_in_region(mid, rs, re, total_len=ref_total):
                    assigned = r
                    break
            if assigned:
                region_hits[assigned].append(h)

        n_gap = "N" * gap_n
        final_parts: List[str] = []
        part_dict: Dict[str, str] = {"LSC": "", "IR": "", "SSC": ""}
        missing_bp = 0
        kept = 0
        with (out_dir / f"{sample_name}.selected_contigs.tsv").open("wt") as tab:
            tab.write("contig\tseed_start\tseed_end\tstrand\tidentity\taln_len\tregion\n")
            for r in region_order:
                rhs = sorted(region_hits[r], key=lambda x: (x[1], x[2]))
                region_seq_parts: List[str] = []
                covered = 0
                for cid, sstart, send, strand, ident, alen, qstart, qend in rhs:
                    cseq = contigs.get(cid, "")
                    if not cseq:
                        continue
                    # Keep only mapped query segment to reduce non-organelle flanking sequence.
                    cseq = slice_query_segment(cseq, qstart=qstart, qend=qend, one_based=True)
                    if strand == "-":
                        cseq = reverse_complement(cseq)
                    region_seq_parts.append(cseq)
                    covered += max(send - sstart + 1, 0)
                    kept += 1
                    tab.write(f"{cid}\t{sstart}\t{send}\t{strand}\t{ident:.4f}\t{alen}\t{r}\n")
                exp = region_len(*cp_regions[r], total_len=ref_total)
                part_key = "IR" if r in {"IRB", "IRA", "IR"} else r
                if region_seq_parts:
                    region_concat = n_gap.join(region_seq_parts)
                    final_parts.append(region_concat)
                    part_dict[part_key] = region_concat
                    missing_bp += max(exp - covered, 0)
                else:
                    region_gap = "N" * exp
                    final_parts.append(region_gap)
                    part_dict[part_key] = region_gap
                    missing_bp += exp
        orient_notes: List[str] = []
        if seed_parts:
            for pname in ("LSC", "SSC"):
                if pname in part_dict and pname in seed_parts:
                    oseq, odir = orient_partition_to_seed(part_dict[pname], seed_parts[pname])
                    part_dict[pname] = oseq
                    orient_notes.append(f"{pname}:{odir}")
        final_parts = [part_dict["LSC"], part_dict["IR"], part_dict["SSC"]]
        merged = "".join(final_parts)
        note = f"type:fragmented;single_ir:{ir_name};missing_bp:{missing_bp};expected_len:{expected_one_ir}"
        if orient_notes:
            note += ";part_orient:" + ",".join(orient_notes)
        out_fa = out_dir / f"{sample_name}.organellar.fasta"
        with out_fa.open("wt") as out:
            out.write(f">{sample_name}\n{merged}\n")
        return str(out_fa), kept, len(merged), note, part_dict

    n_gap = "N" * gap_n
    seq_parts: List[str] = []
    kept = 0
    with (out_dir / f"{sample_name}.selected_contigs.tsv").open("wt") as tab:
        tab.write("contig\tseed_start\tseed_end\tstrand\tidentity\taln_len\n")
        for contig_id, sstart, send, strand, ident, alen, qstart, qend in chosen:
            if contig_id not in contigs:
                continue
            cseq = contigs[contig_id]
            if tool == "minimap2":
                cseq = slice_query_segment(cseq, qstart=qstart, qend=qend, one_based=False)
            else:
                cseq = slice_query_segment(cseq, qstart=qstart, qend=qend, one_based=True)
            if strand == "-":
                cseq = reverse_complement(cseq)
            seq_parts.append(cseq)
            tab.write(f"{contig_id}\t{sstart}\t{send}\t{strand}\t{ident:.4f}\t{alen}\n")
            kept += 1

    merged = n_gap.join(seq_parts)
    note = "-"
    if organelle_mode in {"plant_pt", "animal_mt"} and seed_seq:
        merged, orient = rotate_sequence_to_seed_start(merged, seed_seq=seed_seq)
        note = f"rotated:{orient}"
    if organelle_mode == "plant_pt" and pt_single_ir and cp_regions:
        merged, ir_name = reorder_cp_single_ir(merged, cp_regions=cp_regions, keep_ir=pt_keep_ir)
        note = f"{note};single_ir:{ir_name}" if note != "-" else f"single_ir:{ir_name}"
    out_fa = out_dir / f"{sample_name}.organellar.fasta"
    with out_fa.open("wt") as out:
        out.write(f">{sample_name}\n{merged}\n")
    return str(out_fa), kept, len(merged), note, None


def cmd_sort_organ(args: argparse.Namespace) -> int:
    in_dir = Path(args.input_dir).resolve()
    out_dir = Path(args.outdir).resolve()
    seed = Path(args.seed).resolve()
    if not in_dir.exists():
        raise FileNotFoundError(f"input-dir not found: {in_dir}")
    if not seed.exists():
        raise FileNotFoundError(f"seed fasta not found: {seed}")
    out_dir.mkdir(parents=True, exist_ok=True)
    seed_seq = read_primary_fasta_sequence(seed)
    prof = ORG_SORT_DEFAULTS.get(args.organelle_mode, ORG_SORT_DEFAULTS["generic"])
    min_identity = args.min_identity if args.min_identity is not None else prof["min_identity"]
    min_len = args.min_len if args.min_len is not None else prof["min_len"]
    gap_n = args.gap_n if args.gap_n is not None else prof["gap_n"]
    logger.info(
        "sortOrgan defaults for mode=%s: min_identity=%.3f min_len=%d gap_n=%d min_non_n_len=%d",
        args.organelle_mode,
        min_identity,
        min_len,
        gap_n,
        args.min_non_n_len,
    )
    cp_regions: Optional[Dict[str, Tuple[int, int]]] = None
    seed_parts: Optional[Dict[str, str]] = None
    expected_one_ir_len: Optional[int] = None
    if args.organelle_mode == "plant_pt" and args.pt_single_ir is None:
        logger.info("plant_pt mode detected: enabling --pt-single-ir by default.")
        args.pt_single_ir = True

    if args.organelle_mode == "plant_pt":
        cp_regions = resolve_cp_regions(
            seed_fa=seed,
            out_dir=out_dir,
            cp_regions_path=args.cp_regions,
            cpstools_bin=args.cpstools_bin,
            cpstools_args=list(args.cpstools_args),
        )
    if args.organelle_mode == "plant_pt" and args.pt_single_ir:
        ir_name = choose_ir_name(cp_regions, args.pt_keep_ir)
        seed_parts, _ = reorder_cp_single_ir_parts(seed_seq, cp_regions=cp_regions, keep_ir=args.pt_keep_ir)
        expected_one_ir_len = (
            region_len(*cp_regions["LSC"], total_len=len(seed_seq))
            + region_len(*cp_regions[ir_name], total_len=len(seed_seq))
            + region_len(*cp_regions["SSC"], total_len=len(seed_seq))
        )

    sample_dirs = sorted([p for p in in_dir.iterdir() if p.is_dir()])
    if not sample_dirs:
        raise ValueError("No sample folders found in input-dir.")

    summary = out_dir / "sortorgan_summary.tsv"
    all_multi = out_dir / "assembled_samples.fasta"
    part_files: Dict[str, Path] = {}
    if args.organelle_mode == "plant_pt" and args.pt_single_ir:
        pdir = out_dir / "partitions"
        pdir.mkdir(parents=True, exist_ok=True)
        part_files = {
            "LSC": pdir / "LSC_samples.fasta",
            "IR": pdir / "IR_samples.fasta",
            "SSC": pdir / "SSC_samples.fasta",
        }

    with summary.open("wt") as sumf, all_multi.open("wt") as mf:
        lsc_fh = part_files.get("LSC").open("wt") if "LSC" in part_files else None
        ir_fh = part_files.get("IR").open("wt") if "IR" in part_files else None
        ssc_fh = part_files.get("SSC").open("wt") if "SSC" in part_files else None
        try:
            part_handles = {"LSC": lsc_fh, "IR": ir_fh, "SSC": ssc_fh}
            sumf.write(
                "sample\tstatus\tmode\tsource\tchosen_fasta\tcandidate_count\tfastg\tselected_contigs\tassembled_len\tnon_n_len\tassembled_fasta\tmessage\n"
            )
            for sdir in sample_dirs:
                sample = sdir.name
                candidates, fastg, source = discover_sample_candidates(sdir)
                if not candidates:
                    sumf.write(
                        f"{sample}\tFAIL\t{args.organelle_mode}\t{source}\t-\t0\t{fastg or '-'}\t0\t0\t0\t-\tno candidate fasta found\n"
                    )
                    continue

                sample_out = out_dir / sample
                try:
                    chosen_fa, cand_n = select_best_candidate_fasta(
                        candidates=candidates,
                        seed_fa=seed,
                        sample_name=sample,
                        sample_out=sample_out,
                        min_identity=min_identity,
                        min_len=min_len,
                        aligner=args.aligner,
                        expected_len=expected_one_ir_len,
                    )
                    out_fa, nsel, alen, note, part_seqs = build_sample_assembly_from_contigs(
                        contig_fa=chosen_fa,
                        seed_fa=seed,
                        sample_name=sample,
                        out_dir=sample_out,
                        min_identity=min_identity,
                        min_len=min_len,
                        gap_n=gap_n,
                        aligner=args.aligner,
                        organelle_mode=args.organelle_mode,
                        seed_seq=seed_seq,
                        pt_single_ir=args.pt_single_ir,
                        cp_regions=cp_regions,
                        pt_keep_ir=args.pt_keep_ir,
                        cpstools_bin=args.cpstools_bin,
                        pt_fragment_min_len=args.pt_fragment_min_len,
                        pt_complete_min_frac=args.pt_complete_min_frac,
                        seed_len=len(seed_seq),
                        seed_parts=seed_parts,
                    )
                    seqs = read_fasta_sequences(Path(out_fa))
                    seq_non_n = max((non_n_length(seq) for seq in seqs.values()), default=0)
                    if seq_non_n >= args.min_non_n_len:
                        for sid, seq in seqs.items():
                            mf.write(f">{sid}\n{seq}\n")
                        if part_seqs:
                            for pname in ("LSC", "IR", "SSC"):
                                pseq = part_seqs.get(pname, "")
                                if pseq:
                                    with (sample_out / f"{sample}.{pname}.fasta").open("wt") as p_out:
                                        p_out.write(f">{sample}\n{pseq}\n")
                            for pname in ("LSC", "IR", "SSC"):
                                fh = part_handles.get(pname)
                                pseq = part_seqs.get(pname, "")
                                if fh and pseq:
                                    fh.write(f">{sample}\n{pseq}\n")
                        status = "OK"
                        note2 = note
                    else:
                        status = "FILTERED"
                        note2 = (
                            f"{note};filtered:min_non_n_len<{args.min_non_n_len};non_n_len:{seq_non_n}"
                            if note != "-"
                            else f"filtered:min_non_n_len<{args.min_non_n_len};non_n_len:{seq_non_n}"
                        )
                    sumf.write(
                        f"{sample}\t{status}\t{args.organelle_mode}\t{source}\t{chosen_fa}\t{cand_n}\t{fastg or '-'}\t{nsel}\t{alen}\t{seq_non_n}\t{out_fa}\t{note2}\n"
                    )
                except Exception as exc:
                    sumf.write(
                        f"{sample}\tFAIL\t{args.organelle_mode}\t{source}\t-\t{len(candidates)}\t{fastg or '-'}\t0\t0\t0\t-\t{str(exc).replace(chr(9), ' ')}\n"
                    )
        finally:
            for fh in (lsc_fh, ir_fh, ssc_fh):
                if fh:
                    fh.close()

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


def run_panman_pipeline(args: argparse.Namespace, input_fa: Path, out_dir: Path, strict_args: bool) -> Dict[str, Path]:
    pangraph_out = out_dir / "pangraph"
    dipper_out = out_dir / "dipper"
    twilight_out = out_dir / "twilight"
    panman_out = out_dir / "panman_run"
    pangraph_out.mkdir(parents=True, exist_ok=True)
    dipper_out.mkdir(parents=True, exist_ok=True)
    twilight_out.mkdir(parents=True, exist_ok=True)
    panman_out.mkdir(parents=True, exist_ok=True)

    pangraph_json = Path(args.pangraph_json).resolve() if args.pangraph_json else (pangraph_out / "pangraph_output.json")
    aln_fa = Path(args.aln_file).resolve() if args.aln_file else input_fa
    guide_tree = Path(args.guide_tree).resolve() if args.guide_tree else (twilight_out / "guide_tree.nwk")
    dipper_graph = Path(args.dipper_graph).resolve() if args.dipper_graph else (dipper_out / "graph.gfa")
    panman_file = panman_out / "result.panman"

    def _render_tokens(items: List[str]) -> List[str]:
        rendered: List[str] = []
        for x in items:
            rendered.append(
                x.replace("{input_fasta}", str(input_fa))
                .replace("{pangraph_out}", str(pangraph_out))
                .replace("{pangraph_json}", str(pangraph_json))
                .replace("{dipper_out}", str(dipper_out))
                .replace("{twilight_out}", str(twilight_out))
                .replace("{dipper_graph}", str(dipper_graph))
                .replace("{aln_fasta}", str(aln_fa))
                .replace("{guide_tree}", str(guide_tree))
                .replace("{panman_out}", str(panman_out))
                .replace("{panman_file}", str(panman_file))
            )
        return rendered

    if args.run_pangraph:
        pangraph_bin = shutil.which(args.pangraph_bin)
        if not pangraph_bin:
            raise RuntimeError(f"PanGraph executable not found: {args.pangraph_bin}")
        if not args.pangraph_args:
            if strict_args:
                raise ValueError(
                    "--run-pangraph requires --pangraph-args. "
                    "Use placeholders like {input_fasta}, {pangraph_out}, {pangraph_json}."
                )
            run_command(
                [pangraph_bin, "build", "--circular", "--len", "1000", str(input_fa)],
                stdout_path=pangraph_json,
            )
        else:
            run_command([pangraph_bin] + _render_tokens(list(args.pangraph_args)))
    elif args.run_panman and not pangraph_json.exists():
        raise FileNotFoundError(
            f"--run-panman requires PanGraph JSON, but not found: {pangraph_json}. "
            "Use --run-pangraph or provide --pangraph-json."
        )

    if args.run_dipper:
        dipper_bin = shutil.which(args.dipper_bin)
        if not dipper_bin:
            raise RuntimeError(f"DIPPER executable not found: {args.dipper_bin}")
        if not args.dipper_args:
            if strict_args:
                raise ValueError(
                    "--run-dipper requires --dipper-args. "
                    "Use placeholders like {input_fasta}, {dipper_out}, {dipper_graph}, {aln_fasta}."
                )
            # Default assumes DIPPER can output an alignment-like file to {aln_fasta}.
            run_command([dipper_bin, "-i", "m", "-I", str(input_fa), "-O", str(aln_fa), "-m", "1", "-d", "4"])
        else:
            run_command([dipper_bin] + _render_tokens(list(args.dipper_args)))
    elif args.run_twilight and not Path(aln_fa).exists():
        raise FileNotFoundError(
            f"--run-twilight requires an alignment file, but not found: {aln_fa}. "
            "Use --run-dipper or provide --aln-file."
        )

    if args.run_twilight:
        twilight_bin = shutil.which(args.twilight_bin)
        if not twilight_bin:
            raise RuntimeError(f"TWILIGHT executable not found: {args.twilight_bin}")
        if not args.twilight_args:
            if strict_args:
                raise ValueError(
                    "--run-twilight requires --twilight-args. "
                    "Use placeholders like {dipper_graph}, {aln_fasta}, {twilight_out}, {guide_tree}."
                )
            # Generic default: alignment in, guide tree out.
            run_command([twilight_bin, "-I", str(aln_fa), "-O", str(guide_tree)])
        else:
            run_command([twilight_bin] + _render_tokens(list(args.twilight_args)))
    elif args.run_panman and not Path(guide_tree).exists():
        raise FileNotFoundError(
            f"--run-panman requires guide tree, but not found: {guide_tree}. "
            "Use --run-twilight or provide --guide-tree."
        )

    blocks_dir = Path(args.blocks_dir).resolve() if args.blocks_dir else (out_dir / "panman_blocks")
    if args.run_panman:
        panman_bin = shutil.which(args.panman_bin)
        if not panman_bin:
            raise RuntimeError(f"panman executable not found: {args.panman_bin}")
        blocks_dir.mkdir(parents=True, exist_ok=True)
        # panman ecosystem commonly exposes `panmanUtils` with its own CLI syntax.
        # We therefore pass through user-provided args directly.
        if not args.panman_args:
            if strict_args:
                raise ValueError(
                    "--run-panman requires --panman-args for your panmanUtils workflow. "
                    "Please provide panmanUtils arguments that write block FASTA files into --blocks-dir."
                )
            cmd = [panman_bin, "-P", str(pangraph_json), "-N", str(guide_tree), "-o", str(panman_file)]
        else:
            cmd = [panman_bin] + _render_tokens(list(args.panman_args))
        run_command(cmd)

    return {
        "pangraph_json": pangraph_json,
        "aln_fa": aln_fa,
        "guide_tree": guide_tree,
        "dipper_graph": dipper_graph,
        "panman_out": panman_out,
        "panman_file": panman_file,
        "blocks_dir": blocks_dir,
    }


def cmd_panman(args: argparse.Namespace) -> int:
    input_fa = Path(args.input).resolve()
    out_dir = Path(args.outdir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    if not input_fa.exists():
        raise FileNotFoundError(f"Input multifasta not found: {input_fa}")
    paths = run_panman_pipeline(args, input_fa=input_fa, out_dir=out_dir, strict_args=False)
    logger.info("PanGraph JSON: %s", paths["pangraph_json"])
    logger.info("Guide tree: %s", paths["guide_tree"])
    logger.info("PanMAN output: %s", paths["panman_file"])
    logger.info("Blocks dir (if generated by panman args): %s", paths["blocks_dir"])
    return 0


def cmd_mt_blocks(args: argparse.Namespace) -> int:
    input_fa = Path(args.input).resolve()
    out_dir = Path(args.outdir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    if not input_fa.exists():
        raise FileNotFoundError(f"Input multifasta not found: {input_fa}")

    paths = run_panman_pipeline(args, input_fa=input_fa, out_dir=out_dir, strict_args=True)
    blocks_dir = paths["blocks_dir"]
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
        organelle_mode="plant_pt",
        pt_single_ir=args.pt_single_ir,
        pt_keep_ir=args.pt_keep_ir,
        cp_regions=args.cp_regions,
        cpstools_bin=args.cpstools_bin,
        cpstools_args=args.cpstools_args,
        pt_fragment_min_len=args.pt_fragment_min_len,
        pt_complete_min_frac=args.pt_complete_min_frac,
        min_identity=args.min_identity,
        min_len=args.min_len_pt,
        gap_n=args.gap_n,
        min_non_n_len=args.min_non_n_len,
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
        organelle_mode="plant_mt",
        pt_single_ir=False,
        pt_keep_ir="auto",
        cp_regions=None,
        cpstools_bin="cpstools",
        cpstools_args=[],
        pt_fragment_min_len=1000,
        pt_complete_min_frac=0.85,
        min_identity=args.min_identity,
        min_len=args.min_len_mt,
        gap_n=args.gap_n,
        min_non_n_len=args.min_non_n_len,
    )
    cmd_sort_organ(ss)

    ms = argparse.Namespace(
        input=str((Path(args.outdir).resolve() / "sortOrgan" / "assembled_samples.fasta")),
        outdir=str((Path(args.outdir).resolve() / "mtBlocks")),
        run_pangraph=args.run_pangraph,
        pangraph_bin=args.pangraph_bin,
        pangraph_args=args.pangraph_args,
        pangraph_json=args.pangraph_json,
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
        organelle_mode="animal_mt",
        pt_single_ir=False,
        pt_keep_ir="auto",
        cp_regions=None,
        cpstools_bin="cpstools",
        cpstools_args=[],
        pt_fragment_min_len=1000,
        pt_complete_min_frac=0.85,
        min_identity=args.min_identity,
        min_len=args.min_len_mt,
        gap_n=args.gap_n,
        min_non_n_len=args.min_non_n_len,
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

    aligned = out_dir / "aligned.fasta"
    trimmed = out_dir / "trimmed.fasta"
    mafft_bin = ensure_tool(args.mafft_bin)
    trimal_bin = ensure_tool(args.trimal_bin)

    has_msa = bool(args.multifasta)
    has_vcf = bool(args.vcf)
    has_ref = bool(args.ref)

    if has_msa and (has_vcf or has_ref):
        raise ValueError("Use either --multifasta OR (--vcf and --ref), not both.")
    if has_msa:
        multifasta_in = Path(args.multifasta).resolve()
        if not multifasta_in.exists():
            raise FileNotFoundError(f"multifasta not found: {multifasta_in}")
        multifasta = out_dir / "all_samples.fasta"
        shutil.copy2(multifasta_in, multifasta)
        run_alignment_with_direction(
            multifasta=multifasta,
            aligned=aligned,
            mafft_bin=mafft_bin,
            adjust_direction=args.auto_reverse,
            threads=args.align_threads,
        )
        pre_trim = out_dir / "trimmed.pre.fasta"
        run_command([trimal_bin, "-automated1", "-in", str(aligned), "-out", str(pre_trim)])
        if args.max_missing_frac < 1.0 or args.snp_only:
            kept, total = filter_alignment_sites(
                in_fa=pre_trim,
                out_fa=trimmed,
                max_missing_frac=args.max_missing_frac,
                snp_only=args.snp_only,
            )
            logger.info("Post-trim site filter kept %d/%d columns (missing<=%.3f, snp_only=%s).", kept, total, args.max_missing_frac, args.snp_only)
        else:
            shutil.move(str(pre_trim), str(trimmed))
        logger.info("Alignment and trimming completed from multifasta input.")
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

    if not (has_vcf and has_ref):
        raise ValueError("For VCF mode, both --vcf and --ref are required.")

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

    run_alignment_with_direction(
        multifasta=multifasta,
        aligned=aligned,
        mafft_bin=mafft_bin,
        adjust_direction=args.auto_reverse,
        threads=args.align_threads,
    )
    pre_trim = out_dir / "trimmed.pre.fasta"
    run_command([trimal_bin, "-automated1", "-in", str(aligned), "-out", str(pre_trim)])
    if args.max_missing_frac < 1.0 or args.snp_only:
        kept, total = filter_alignment_sites(
            in_fa=pre_trim,
            out_fa=trimmed,
            max_missing_frac=args.max_missing_frac,
            snp_only=args.snp_only,
        )
        logger.info("Post-trim site filter kept %d/%d columns (missing<=%.3f, snp_only=%s).", kept, total, args.max_missing_frac, args.snp_only)
    else:
        shutil.move(str(pre_trim), str(trimmed))
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


def cmd_msa2vcf(args: argparse.Namespace) -> int:
    in_fa = Path(args.input).resolve()
    out_dir = Path(args.outdir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    if not in_fa.exists():
        raise FileNotFoundError(f"input fasta not found: {in_fa}")

    mafft_bin = ensure_tool(args.mafft_bin)
    trimal_bin = ensure_tool(args.trimal_bin) if args.trim else args.trimal_bin
    snp_sites = ensure_tool(args.snp_sites_bin)
    bcftools = ensure_tool(args.bcftools_bin)

    aligned = out_dir / "aligned.auto_rc.fasta"
    run_alignment_with_direction(
        multifasta=in_fa,
        aligned=aligned,
        mafft_bin=mafft_bin,
        adjust_direction=args.auto_reverse,
        threads=args.threads,
    )

    msa_for_vcf = aligned
    trimmed = out_dir / "aligned.auto_rc.trimmed.fasta"
    if args.trim:
        run_command([trimal_bin, "-automated1", "-in", str(aligned), "-out", str(trimmed)])
        msa_for_vcf = trimmed

    if args.norm_ref:
        norm_ref = Path(args.norm_ref).resolve()
        if not norm_ref.exists():
            raise FileNotFoundError(f"--norm-ref not found: {norm_ref}")
    else:
        norm_ref = out_dir / "norm_ref.fasta"
        chosen = write_ref_from_msa(msa_for_vcf, norm_ref, ref_id=args.ref_id)
        logger.info("Normalization reference selected from MSA: %s", chosen)

    raw_vcf = out_dir / "raw.snpsites.vcf"
    run_command([snp_sites, "-v", str(msa_for_vcf)], stdout_path=raw_vcf)

    atomized = out_dir / "atomized.vcf.gz"
    split = out_dir / "split.multiallelic.vcf.gz"
    biallelic = out_dir / "biallelic.snps.vcf.gz"

    run_command(
        [
            bcftools,
            "norm",
            "-f",
            str(norm_ref),
            "-a",
            "--atom-overlaps",
            ".",
            str(raw_vcf),
            "-Oz",
            "-o",
            str(atomized),
        ]
    )
    run_command([bcftools, "index", "-t", str(atomized)])

    run_command(
        [
            bcftools,
            "norm",
            "-m",
            "-any",
            str(atomized),
            "-Oz",
            "-o",
            str(split),
        ]
    )
    run_command([bcftools, "index", "-t", str(split)])

    run_command(
        [
            bcftools,
            "view",
            "-v",
            "snps",
            "-m2",
            "-M2",
            str(split),
            "-Oz",
            "-o",
            str(biallelic),
        ]
    )
    run_command([bcftools, "index", "-t", str(biallelic)])

    logger.info("MSA2VCF completed. Outputs: %s", out_dir)
    return 0


def concatenate_three_partition_alignments(
    lsc_aln: Path,
    ir_aln: Path,
    ssc_aln: Path,
    out_fa: Path,
) -> None:
    lsc = read_fasta_sequences(lsc_aln)
    ir = read_fasta_sequences(ir_aln)
    ssc = read_fasta_sequences(ssc_aln)
    if not lsc or not ir or not ssc:
        raise ValueError("Partition alignments must be non-empty for LSC/IR/SSC.")
    l_lens = {len(v) for v in lsc.values()}
    i_lens = {len(v) for v in ir.values()}
    s_lens = {len(v) for v in ssc.values()}
    if len(l_lens) != 1 or len(i_lens) != 1 or len(s_lens) != 1:
        raise ValueError("Inconsistent sequence lengths detected within partition alignments.")
    l_len = next(iter(l_lens))
    i_len = next(iter(i_lens))
    s_len = next(iter(s_lens))
    all_ids = sorted(set(lsc) | set(ir) | set(ssc))
    with out_fa.open("wt") as out:
        for sid in all_ids:
            out.write(
                f">{sid}\n"
                f"{lsc.get(sid, 'N' * l_len)}{ir.get(sid, 'N' * i_len)}{ssc.get(sid, 'N' * s_len)}\n"
            )


def write_partition_sample_stats(partition_fastas: Dict[str, Path], out_tsv: Path) -> None:
    with out_tsv.open("wt") as out:
        out.write("partition\tsample\tlen\tnon_missing\tmissing_frac\n")
        for pname in ("LSC", "IR", "SSC"):
            pfa = partition_fastas[pname]
            seqs = read_fasta_sequences(pfa)
            for sid in sorted(seqs.keys()):
                seq = seqs[sid]
                total = len(seq)
                miss = sum(1 for c in seq if c in "Nn-?")
                non_miss = total - miss
                miss_frac = (miss / total) if total > 0 else 1.0
                out.write(f"{pname}\t{sid}\t{total}\t{non_miss}\t{miss_frac:.6f}\n")


def filter_partition_samples_by_missing(
    partition_fastas: Dict[str, Path],
    out_dir: Path,
    max_missing_frac: float,
) -> Dict[str, Path]:
    raw: Dict[str, Dict[str, str]] = {p: read_fasta_sequences(partition_fastas[p]) for p in ("LSC", "IR", "SSC")}
    all_ids = sorted(set(raw["LSC"]) | set(raw["IR"]) | set(raw["SSC"]))
    miss: Dict[str, Dict[str, float]] = {sid: {} for sid in all_ids}
    keep_ids: List[str] = []
    report = out_dir / "partition_sample_filter.tsv"
    with report.open("wt") as rep:
        rep.write("sample\tLSC_missing\tIR_missing\tSSC_missing\tkeep\n")
        for sid in all_ids:
            ok = True
            for pname in ("LSC", "IR", "SSC"):
                seq = raw[pname].get(sid, "")
                if not seq:
                    m = 1.0
                else:
                    total = len(seq)
                    ms = sum(1 for c in seq if c in "Nn-?")
                    m = ms / total if total > 0 else 1.0
                miss[sid][pname] = m
                if m > max_missing_frac:
                    ok = False
            if ok:
                keep_ids.append(sid)
            rep.write(
                f"{sid}\t{miss[sid]['LSC']:.6f}\t{miss[sid]['IR']:.6f}\t{miss[sid]['SSC']:.6f}\t"
                f"{'KEEP' if ok else 'DROP'}\n"
            )

    if not keep_ids:
        raise ValueError(f"No samples left after partition missing filter (threshold={max_missing_frac}).")

    out_paths: Dict[str, Path] = {}
    for pname in ("LSC", "IR", "SSC"):
        out_p = out_dir / f"{pname}.aligned.filtered.fasta"
        out_paths[pname] = out_p
        kept = {sid: raw[pname][sid] for sid in keep_ids if sid in raw[pname]}
        write_fasta_sequences(out_p, kept)
    logger.info(
        "Partition sample filter kept %d/%d samples at threshold %.3f. Report: %s",
        len(keep_ids),
        len(all_ids),
        max_missing_frac,
        report,
    )
    return out_paths


def cmd_align(args: argparse.Namespace) -> int:
    in_fa = Path(args.input).resolve()
    out_dir = Path(args.outdir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    if not in_fa.exists():
        raise FileNotFoundError(f"input fasta not found: {in_fa}")

    mafft_bin = ensure_tool(args.mafft_bin)
    aligned = out_dir / "aligned.fasta"
    if args.plant_pt_partition:
        if args.partition_dir:
            pdir = Path(args.partition_dir).resolve()
        elif in_fa.is_dir():
            pdir = in_fa / "partitions"
        else:
            pdir = in_fa.parent / "partitions"
        lsc_in = pdir / "LSC_samples.fasta"
        ir_in = pdir / "IR_samples.fasta"
        ssc_in = pdir / "SSC_samples.fasta"
        for pf in (lsc_in, ir_in, ssc_in):
            if not pf.exists():
                raise FileNotFoundError(f"Partition fasta not found: {pf}")
        lsc_aln = out_dir / "LSC.aligned.fasta"
        ir_aln = out_dir / "IR.aligned.fasta"
        ssc_aln = out_dir / "SSC.aligned.fasta"
        run_alignment_with_direction(
            multifasta=lsc_in,
            aligned=lsc_aln,
            mafft_bin=mafft_bin,
            adjust_direction=args.auto_reverse,
            threads=args.threads,
        )
        run_alignment_with_direction(
            multifasta=ir_in,
            aligned=ir_aln,
            mafft_bin=mafft_bin,
            adjust_direction=args.auto_reverse,
            threads=args.threads,
        )
        run_alignment_with_direction(
            multifasta=ssc_in,
            aligned=ssc_aln,
            mafft_bin=mafft_bin,
            adjust_direction=args.auto_reverse,
            threads=args.threads,
        )
        write_partition_sample_stats(
            partition_fastas={"LSC": lsc_aln, "IR": ir_aln, "SSC": ssc_aln},
            out_tsv=out_dir / "partition_sample_stats.aligned.tsv",
        )
        filtered_aln = filter_partition_samples_by_missing(
            partition_fastas={"LSC": lsc_aln, "IR": ir_aln, "SSC": ssc_aln},
            out_dir=out_dir,
            max_missing_frac=args.partition_max_missing_frac,
        )
        concatenate_three_partition_alignments(
            lsc_aln=filtered_aln["LSC"],
            ir_aln=filtered_aln["IR"],
            ssc_aln=filtered_aln["SSC"],
            out_fa=aligned,
        )

        if args.trim:
            trimal_bin = ensure_tool(args.trimal_bin)
            lsc_pre = out_dir / "LSC.trimmed.pre.fasta"
            ir_pre = out_dir / "IR.trimmed.pre.fasta"
            ssc_pre = out_dir / "SSC.trimmed.pre.fasta"
            run_command([trimal_bin, "-automated1", "-in", str(filtered_aln["LSC"]), "-out", str(lsc_pre)])
            run_command([trimal_bin, "-automated1", "-in", str(filtered_aln["IR"]), "-out", str(ir_pre)])
            run_command([trimal_bin, "-automated1", "-in", str(filtered_aln["SSC"]), "-out", str(ssc_pre)])

            trimmed_pre = out_dir / "trimmed.pre.fasta"
            concatenate_three_partition_alignments(
                lsc_aln=lsc_pre,
                ir_aln=ir_pre,
                ssc_aln=ssc_pre,
                out_fa=trimmed_pre,
            )

            lsc_trim = out_dir / "LSC.trimmed.fasta"
            ir_trim = out_dir / "IR.trimmed.fasta"
            ssc_trim = out_dir / "SSC.trimmed.fasta"
            if args.max_missing_frac < 1.0 or args.snp_only:
                k1, t1 = filter_alignment_sites(
                    in_fa=lsc_pre,
                    out_fa=lsc_trim,
                    max_missing_frac=args.max_missing_frac,
                    snp_only=args.snp_only,
                )
                k2, t2 = filter_alignment_sites(
                    in_fa=ir_pre,
                    out_fa=ir_trim,
                    max_missing_frac=args.max_missing_frac,
                    snp_only=args.snp_only,
                )
                k3, t3 = filter_alignment_sites(
                    in_fa=ssc_pre,
                    out_fa=ssc_trim,
                    max_missing_frac=args.max_missing_frac,
                    snp_only=args.snp_only,
                )
                logger.info(
                    "Partition post-trim filters kept LSC %d/%d, IR %d/%d, SSC %d/%d columns (missing<=%.3f, snp_only=%s).",
                    k1,
                    t1,
                    k2,
                    t2,
                    k3,
                    t3,
                    args.max_missing_frac,
                    args.snp_only,
                )
            else:
                shutil.copyfile(lsc_pre, lsc_trim)
                shutil.copyfile(ir_pre, ir_trim)
                shutil.copyfile(ssc_pre, ssc_trim)

            trimmed = out_dir / "trimmed.fasta"
            concatenate_three_partition_alignments(
                lsc_aln=lsc_trim,
                ir_aln=ir_trim,
                ssc_aln=ssc_trim,
                out_fa=trimmed,
            )
            write_partition_sample_stats(
                partition_fastas={"LSC": lsc_trim, "IR": ir_trim, "SSC": ssc_trim},
                out_tsv=out_dir / "partition_sample_stats.tsv",
            )
            logger.info("Align completed (partition mode): %s, %s", aligned, trimmed)
        else:
            write_partition_sample_stats(
                partition_fastas=filtered_aln,
                out_tsv=out_dir / "partition_sample_stats.tsv",
            )
            logger.info("Align completed (partition mode): %s", aligned)
        return 0
    else:
        run_alignment_with_direction(
            multifasta=in_fa,
            aligned=aligned,
            mafft_bin=mafft_bin,
            adjust_direction=args.auto_reverse,
            threads=args.threads,
        )
    if args.trim:
        trimal_bin = ensure_tool(args.trimal_bin)
        pre_trim = out_dir / "trimmed.pre.fasta"
        trimmed = out_dir / "trimmed.fasta"
        run_command([trimal_bin, "-automated1", "-in", str(aligned), "-out", str(pre_trim)])
        if args.max_missing_frac < 1.0 or args.snp_only:
            kept, total = filter_alignment_sites(
                in_fa=pre_trim,
                out_fa=trimmed,
                max_missing_frac=args.max_missing_frac,
                snp_only=args.snp_only,
            )
            logger.info("Align post-trim filter kept %d/%d columns (missing<=%.3f, snp_only=%s).", kept, total, args.max_missing_frac, args.snp_only)
        else:
            shutil.move(str(pre_trim), str(trimmed))
        logger.info("Align completed: %s, %s", aligned, trimmed)
    else:
        logger.info("Align completed: %s", aligned)
    return 0


def cmd_phyview(args: argparse.Namespace) -> int:
    trimmed = Path(args.input).resolve()
    out_dir = Path(args.outdir).resolve()
    mode_count = int(args.run_ml) + int(args.run_nj) + int(args.run_popart) + int(args.run_beast)
    if mode_count != 1:
        raise ValueError("PhyView requires exactly one mode: --run_ml or --run_nj or --run_popart or --run_beast")
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

    if args.run_beast:
        if not args.beast_template:
            raise ValueError("--beast-template is required with --run_beast")
        beast_template = resolve_beast_template_arg(args.beast_template)
        if not beast_template.exists():
            raise FileNotFoundError(f"BEAST template not found: {beast_template}")
        run_beast_pipeline(
            aligned_fasta=trimmed,
            out_dir=out_dir,
            template_xml=beast_template,
            prefix=args.beast_prefix,
            chain_length=args.beast_chain_length,
            store_every=args.beast_store_every,
            beast_bin=args.beast_bin,
            beast_threads=args.beast_threads,
            use_beagle=args.beast_beagle,
            treeannotator_bin=args.treeannotator_bin,
            run_treeannotator=args.run_treeannotator,
            burnin=args.beast_burnin,
            heights=args.beast_heights,
            include_posterior=args.include_posterior,
            old_prefix=args.beast_old_prefix,
        )
        logger.info("PhyView completed (BEAST): %s", out_dir)
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
        help="Run full workflow from VCF+REF or from multifasta to trimmed alignment",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p_run.add_argument("-v", "--vcf", help="Input VCF or VCF.GZ (VCF mode)")
    p_run.add_argument("-r", "--ref", help="Reference fasta (VCF mode)")
    p_run.add_argument("-i", "--multifasta", help="Input multifasta (MSA mode: align+trim only)")
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
    p_run.add_argument(
        "--align-threads",
        default="AUTO",
        help="MAFFT threads used for alignment in run mode (AUTO uses all cores)",
    )
    p_run.add_argument("--trimal-bin", default="trimal", help="Path or name of trimal executable")
    p_run.add_argument(
        "--auto-reverse",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Use MAFFT --adjustdirection during alignment",
    )
    p_run.add_argument(
        "--max-missing-frac",
        type=float,
        default=1.0,
        help="Post-trim site filter: drop columns with missing fraction above this threshold (0-1)",
    )
    p_run.add_argument(
        "--snp-only",
        action="store_true",
        help="Post-trim site filter: keep SNP columns only (useful for highly divergent genomes)",
    )
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
    mode_group.add_argument("--run_beast", "--run-beast", action="store_true", help="Run BEAST tree inference from input alignment")
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
    p_phy.add_argument(
        "--beast-template",
        help="BEAST XML template path (required with --run_beast); use 'default' for built-in template",
    )
    p_phy.add_argument("--beast-prefix", default="organpath_beast", help="Output prefix used by BEAST")
    p_phy.add_argument("--beast-old-prefix", default="ultrametric", help="Template prefix token to replace in BEAST XML")
    p_phy.add_argument("--beast-chain-length", type=int, default=20000000, help="BEAST MCMC chain length")
    p_phy.add_argument("--beast-store-every", type=int, default=5000, help="BEAST storeEvery interval")
    p_phy.add_argument("--beast-bin", default="beast", help="BEAST executable name/path")
    p_phy.add_argument("--beast-threads", type=int, default=8, help="Threads for BEAST")
    p_phy.add_argument("--beast-beagle", action="store_true", help="Enable BEAGLE for BEAST run")
    p_phy.add_argument("--run-treeannotator", action=argparse.BooleanOptionalAction, default=True, help="Run treeannotator on BEAST .trees output")
    p_phy.add_argument("--treeannotator-bin", default="treeannotator", help="treeannotator executable name/path")
    p_phy.add_argument("--beast-burnin", type=int, default=10, help="Burn-in percentage for treeannotator")
    p_phy.add_argument("--beast-heights", default="median", help="treeannotator node heights mode")
    p_phy.add_argument("--include-posterior", action="store_true", help="Keep posterior as internal labels in BEAST Newick output")
    p_phy.set_defaults(func=cmd_phyview)

    p_m2v = subs.add_parser(
        "MSA2VCF",
        help="Align multifasta with MAFFT (auto reverse-complement) and build strict biallelic SNP VCF",
    )
    p_m2v.add_argument("-i", "--input", required=True, help="Input multifasta (e.g. assembled_samples.fasta)")
    p_m2v.add_argument("-o", "--outdir", required=True, help="Output directory")
    p_m2v.add_argument("--mafft-bin", default="mafft", help="MAFFT executable name/path")
    p_m2v.add_argument("--auto-reverse", action=argparse.BooleanOptionalAction, default=True, help="Use MAFFT --adjustdirection")
    p_m2v.add_argument("--trim", action="store_true", help="Trim aligned MSA with trimAl before VCF conversion")
    p_m2v.add_argument("--trimal-bin", default="trimal", help="trimAl executable name/path")
    p_m2v.add_argument("--snp-sites-bin", default="snp-sites", help="snp-sites executable name/path")
    p_m2v.add_argument("--bcftools-bin", default="bcftools", help="bcftools executable name/path")
    p_m2v.add_argument("--norm-ref", help="Reference fasta for bcftools norm -f (default: first MSA sequence)")
    p_m2v.add_argument("--ref-id", help="Reference sequence ID from MSA when --norm-ref is not provided")
    p_m2v.set_defaults(func=cmd_msa2vcf)

    p_align = subs.add_parser(
        "align",
        help="Panel: align multifasta with MAFFT and optional trimAl trimming",
    )
    p_align.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input multifasta; with --plant-pt-partition, can be sortOrgan output directory",
    )
    p_align.add_argument("-o", "--outdir", required=True, help="Output directory")
    p_align.add_argument("--mafft-bin", default="mafft", help="MAFFT executable name/path")
    p_align.add_argument("--threads", default="AUTO", help="MAFFT threads (AUTO uses all cores)")
    p_align.add_argument(
        "--auto-reverse",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Use MAFFT --adjustdirection to auto-handle reverse-complement sequences (default: off)",
    )
    p_align.add_argument(
        "--trim",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Run trimAl after MAFFT",
    )
    p_align.add_argument("--trimal-bin", default="trimal", help="trimAl executable name/path")
    p_align.add_argument(
        "--max-missing-frac",
        type=float,
        default=1.0,
        help="Post-trim site filter: drop columns with missing fraction above this threshold (0-1)",
    )
    p_align.add_argument(
        "--snp-only",
        action="store_true",
        help="Post-trim site filter: keep SNP columns only",
    )
    p_align.add_argument(
        "--plant-pt-partition",
        action="store_true",
        help="Use sortOrgan plant_pt partition FASTAs (LSC/IR/SSC) for separate alignment then concatenate",
    )
    p_align.add_argument(
        "--partition-dir",
        help="Partition directory containing LSC_samples.fasta, IR_samples.fasta, SSC_samples.fasta",
    )
    p_align.add_argument(
        "--partition-max-missing-frac",
        type=float,
        default=0.2,
        help="In --plant-pt-partition mode, drop a sample if any partition missing fraction exceeds this threshold (applied right after MAFFT)",
    )
    p_align.set_defaults(func=cmd_align)

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
        help="Sort and orient assembled organellar contigs by seed; output per-sample fasta (mode defaults built-in)",
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
    p_sort.add_argument(
        "--organelle-mode",
        default="generic",
        choices=["generic", "plant_pt", "plant_mt", "animal_mt"],
        help=(
            "Sorting strategy by organelle type. "
            "Defaults: plant_pt(0.95,1000,100), plant_mt(0.95,3000,100), "
            "animal_mt(0.95,1000,100), generic(0.95,1000,100) for (identity,min_len,gap_n)"
        ),
    )
    p_sort.add_argument(
        "--pt-single-ir",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Plant chloroplast mode: reorder as LSC+singleIR+SSC (default: auto-on in plant_pt)",
    )
    p_sort.add_argument(
        "--pt-keep-ir",
        default="auto",
        choices=["auto", "ira", "irb"],
        help="Which IR copy to keep for --pt-single-ir",
    )
    p_sort.add_argument("--cp-regions", help="Region table from cpstools (contains LSC/SSC/IR coordinates)")
    p_sort.add_argument("--cpstools-bin", default="cpstools", help="cpstools executable name/path")
    p_sort.add_argument(
        "--cpstools-args",
        nargs="*",
        default=[],
        help="Args passed to cpstools. Placeholders: {seed_fasta} {cpstools_out} {cp_regions_tsv}",
    )
    p_sort.add_argument("--pt-fragment-min-len", type=int, default=1000, help="plant_pt fragmented mode: minimum contig hit length")
    p_sort.add_argument("--pt-complete-min-frac", type=float, default=0.85, help="plant_pt complete mode: minimum fraction of expected one-IR length")
    p_sort.add_argument("--min-identity", type=float, default=None, help="Minimum alignment identity (default by --organelle-mode)")
    p_sort.add_argument("--min-len", type=int, default=None, help="Minimum aligned length (default by --organelle-mode)")
    p_sort.add_argument("--gap-n", type=int, default=None, help="Number of Ns between selected contigs (default by --organelle-mode)")
    p_sort.add_argument(
        "--min-non-n-len",
        type=int,
        default=0,
        help="Minimum non-N length required to keep sample in assembled_samples.fasta (default: 0)",
    )
    p_sort.set_defaults(func=cmd_sort_organ)

    p_pan = subs.add_parser(
        "panman",
        help="Run PanGraph -> DIPPER -> TWILIGHT -> panmanUtils from an assembled multifasta",
    )
    p_pan.add_argument("-i", "--input", required=True, help="Input multifasta (e.g. sortOrgan/assembled_samples.fasta)")
    p_pan.add_argument("-o", "--outdir", required=True, help="Output directory")
    p_pan.add_argument("--run-pangraph", action=argparse.BooleanOptionalAction, default=True, help="Run PanGraph stage")
    p_pan.add_argument("--pangraph-bin", default="pangraph", help="PanGraph executable name/path")
    p_pan.add_argument(
        "--pangraph-args",
        nargs="*",
        default=[],
        help="Arguments passed to PanGraph. Supports placeholders: {input_fasta} {pangraph_out} {pangraph_json}",
    )
    p_pan.add_argument("--pangraph-json", help="Existing PanGraph JSON (if already generated)")
    p_pan.add_argument(
        "--run-dipper",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Run DIPPER stage (default: off; provide --guide-tree if skipping TWILIGHT)",
    )
    p_pan.add_argument("--dipper-bin", default="dipper", help="DIPPER executable name/path")
    p_pan.add_argument(
        "--dipper-args",
        nargs="*",
        default=[],
        help="Arguments passed to DIPPER. Supports placeholders: {input_fasta} {dipper_out} {dipper_graph} {aln_fasta}",
    )
    p_pan.add_argument(
        "--run-twilight",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Run TWILIGHT stage (default: off; provide --guide-tree)",
    )
    p_pan.add_argument("--twilight-bin", default="twilight", help="TWILIGHT executable name/path")
    p_pan.add_argument(
        "--twilight-args",
        nargs="*",
        default=[],
        help="Arguments passed to TWILIGHT. Supports placeholders: {dipper_graph} {aln_fasta} {twilight_out} {guide_tree}",
    )
    p_pan.add_argument("--run-panman", action=argparse.BooleanOptionalAction, default=True, help="Run panmanUtils stage")
    p_pan.add_argument("--panman-bin", default="panmanUtils", help="panman executable name/path")
    p_pan.add_argument(
        "--panman-args",
        nargs="*",
        default=[],
        help="Arguments passed to panmanUtils. Supports placeholders: {pangraph_json} {guide_tree} {panman_out} {panman_file}",
    )
    p_pan.add_argument("--aln-file", help="Alignment path for DIPPER/TWILIGHT handoff (optional)")
    p_pan.add_argument(
        "--guide-tree",
        help="Guide tree path for PanMAN handoff (optional). Useful with --no-run-twilight.",
    )
    p_pan.add_argument("--dipper-graph", help="DIPPER graph path (optional)")
    p_pan.add_argument("--blocks-dir", help="Directory intended for block outputs if panman args produce them")
    p_pan.set_defaults(func=cmd_panman)

    p_mt = subs.add_parser(
        "mtBlocks",
        help="Plant mitochondrial route: panman blocks -> block MAFFT+trimAl -> concatenated supermatrix",
    )
    p_mt.add_argument("-i", "--input", required=True, help="Input multifasta (per-sample mt assemblies)")
    p_mt.add_argument("-o", "--outdir", required=True, help="Output directory")
    p_mt.add_argument("--run-pangraph", action="store_true", help="Run PanGraph before TWILIGHT/panman")
    p_mt.add_argument("--pangraph-bin", default="pangraph", help="PanGraph executable name/path")
    p_mt.add_argument(
        "--pangraph-args",
        nargs="*",
        default=[],
        help="Arguments passed to PanGraph. Supports placeholders: {input_fasta} {pangraph_out} {pangraph_json}",
    )
    p_mt.add_argument("--pangraph-json", help="Existing PanGraph JSON (if already generated)")
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
        help="Channel: plant chloroplast (default profile: identity=0.95, min_len=1000, gap_n=100)",
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
    p_ch_pt.add_argument("--min-identity", type=float, default=0.95, help="Minimum identity for contig selection (default: 0.95)")
    p_ch_pt.add_argument("--min-len-pt", type=int, default=1000, help="Minimum length for cp contig selection (default: 1000)")
    p_ch_pt.add_argument("--gap-n", type=int, default=100, help="Ns inserted between selected contigs (default: 100)")
    p_ch_pt.add_argument("--min-non-n-len", type=int, default=0, help="Minimum non-N length to keep sample in merged multifasta")
    p_ch_pt.add_argument(
        "--pt-single-ir",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Reorder chloroplast output to LSC+singleIR+SSC (default: on)",
    )
    p_ch_pt.add_argument("--pt-keep-ir", default="auto", choices=["auto", "ira", "irb"], help="IR copy kept in single-IR mode")
    p_ch_pt.add_argument("--cp-regions", help="Region table from cpstools (LSC/SSC/IR coordinates)")
    p_ch_pt.add_argument("--cpstools-bin", default="cpstools", help="cpstools executable")
    p_ch_pt.add_argument("--cpstools-args", nargs="*", default=[], help="Args passed to cpstools")
    p_ch_pt.add_argument("--pt-fragment-min-len", type=int, default=1000, help="Fragmented cp mode minimum contig length")
    p_ch_pt.add_argument("--pt-complete-min-frac", type=float, default=0.85, help="Complete cp mode minimum fraction of expected one-IR length")
    p_ch_pt.set_defaults(func=cmd_channel_plant_pt)

    p_ch_mt = subs.add_parser(
        "plant_mt",
        help="Channel: plant mitochondria (default profile: identity=0.95, min_len=3000, gap_n=100)",
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
    p_ch_mt.add_argument("--min-identity", type=float, default=0.95, help="Minimum identity for contig selection (default: 0.95)")
    p_ch_mt.add_argument("--min-len-mt", type=int, default=3000, help="Minimum length for mt contig selection (default: 3000)")
    p_ch_mt.add_argument("--gap-n", type=int, default=100, help="Ns inserted between selected contigs (default: 100)")
    p_ch_mt.add_argument("--min-non-n-len", type=int, default=0, help="Minimum non-N length to keep sample in merged multifasta")
    p_ch_mt.add_argument("--run-pangraph", action="store_true", help="Run PanGraph before TWILIGHT/panman")
    p_ch_mt.add_argument("--pangraph-bin", default="pangraph", help="PanGraph executable")
    p_ch_mt.add_argument("--pangraph-args", nargs="*", default=[], help="Arguments passed to PanGraph")
    p_ch_mt.add_argument("--pangraph-json", help="Existing PanGraph JSON for mtBlocks")
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
        help="Channel: animal mitochondria (default profile: identity=0.95, min_len=1000, gap_n=100)",
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
    p_ch_amt.add_argument("--min-identity", type=float, default=0.95, help="Minimum identity for contig selection (default: 0.95)")
    p_ch_amt.add_argument("--min-len-mt", type=int, default=1000, help="Minimum length for mt contig selection (default: 1000)")
    p_ch_amt.add_argument("--gap-n", type=int, default=100, help="Ns inserted between selected contigs (default: 100)")
    p_ch_amt.add_argument("--min-non-n-len", type=int, default=0, help="Minimum non-N length to keep sample in merged multifasta")
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
    mode_group.add_argument("--run_beast", "--run-beast", action="store_true", help="Run BEAST tree inference from input alignment")
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
    parser.add_argument(
        "--beast-template",
        help="BEAST XML template path (required with --run_beast); use 'default' for built-in template",
    )
    parser.add_argument("--beast-prefix", default="organpath_beast", help="Output prefix used by BEAST")
    parser.add_argument("--beast-old-prefix", default="ultrametric", help="Template prefix token to replace in BEAST XML")
    parser.add_argument("--beast-chain-length", type=int, default=20000000, help="BEAST MCMC chain length")
    parser.add_argument("--beast-store-every", type=int, default=5000, help="BEAST storeEvery interval")
    parser.add_argument("--beast-bin", default="beast", help="BEAST executable name/path")
    parser.add_argument("--beast-threads", type=int, default=8, help="Threads for BEAST")
    parser.add_argument("--beast-beagle", action="store_true", help="Enable BEAGLE for BEAST run")
    parser.add_argument("--run-treeannotator", action=argparse.BooleanOptionalAction, default=True, help="Run treeannotator on BEAST .trees output")
    parser.add_argument("--treeannotator-bin", default="treeannotator", help="treeannotator executable name/path")
    parser.add_argument("--beast-burnin", type=int, default=10, help="Burn-in percentage for treeannotator")
    parser.add_argument("--beast-heights", default="median", help="treeannotator node heights mode")
    parser.add_argument("--include-posterior", action="store_true", help="Keep posterior as internal labels in BEAST Newick output")
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
