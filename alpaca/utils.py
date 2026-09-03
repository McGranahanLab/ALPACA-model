import os
import pandas as pd
from pathlib import Path
import json
import importlib
import logging
from datetime import datetime
from typing import Optional
import glob
import sys
import csv
import io
import urllib.request


_GENOME_LENGTH_SOURCES = {
    "hg19": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes",
    "hg38": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes",
}
_CANONICAL_CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
_DEFAULT_GENOME_CACHE = Path.home() / ".cache" / "alpaca" / "genomes"
SUPPORTED_GENOME_BUILDS = tuple(sorted(_GENOME_LENGTH_SOURCES.keys()))


def ensure_chr_table(genome_build: str, cache_dir: Optional[str | Path] = None) -> Path:
    """Download (if needed) and cache chromosome lengths for a genome build.

    Returns the path to a CSV with columns ['chr', 'len'] restricted to
    chr1-22, chrX, chrY. Data are sourced from UCSC chrom.sizes files and
    cached under ~/.cache/alpaca/genomes by default.
    """

    build = (genome_build or "").lower()
    if build not in _GENOME_LENGTH_SOURCES:
        raise ValueError(
            f"Unsupported genome build '{genome_build}'. Choose one of {sorted(_GENOME_LENGTH_SOURCES)}."
        )

    cache_root = Path(cache_dir).expanduser().resolve() if cache_dir else _DEFAULT_GENOME_CACHE
    cache_root.mkdir(parents=True, exist_ok=True)
    cache_file = cache_root / f"{build}_chrom_lengths.csv"
    if cache_file.exists():
        return cache_file

    url = _GENOME_LENGTH_SOURCES[build]
    try:
        with urllib.request.urlopen(url) as response:
            payload = response.read().decode("utf-8")
    except Exception as exc:
        raise RuntimeError(
            f"Failed to download chromosome sizes for {genome_build} from {url}."
        ) from exc

    df = pd.read_csv(
        io.StringIO(payload),
        sep="\t",
        header=None,
        usecols=[0, 1],
        names=["chr", "len"],
        dtype={"chr": str, "len": "int64"},
    )
    df = df[df["chr"].isin(_CANONICAL_CHROMOSOMES)].copy()
    if df.empty:
        raise ValueError(
            f"Downloaded chromosome table for {genome_build} did not contain canonical chromosomes."
        )
    missing = [chrom for chrom in _CANONICAL_CHROMOSOMES if chrom not in df["chr"].values]
    if missing:
        raise ValueError(
            f"Chromosome table for {genome_build} is missing entries: {', '.join(missing)}."
        )

    df["chr"] = pd.Categorical(df["chr"], categories=_CANONICAL_CHROMOSOMES, ordered=True)
    df = df.sort_values("chr")
    df.to_csv(cache_file, index=False)
    return cache_file


def show_version():
    try:
        version = importlib.metadata.version("alpaca")
        print(f"alpaca {version}")
    except importlib.metadata.PackageNotFoundError:
        print("alpaca version unknown (not installed)")


def show_help():
    print(
        "ALPACA = ALlele-specific Phylogenetic Analysis of clone Copy-number Alterations"
    )
    print_logo()
    print("")
    print("Usage:")
    print("  alpaca [command]")
    print("")
    print("Commands:")
    print("  version              Show version")
    print("  help                 Show this help")
    print("  run                  Run ALPACA")
    print("  get-scores           Recalculate D_score and CI_score for a tumour from existing outputs")
    print("  input-conversion     Run input conversion")
    print("  ccd                  Calculate clone copy number diversity")
    print("  wgd                  Calculate parent/child WGD ratios")
    print("  ancestor-delta       Compute copy-number changes to each ancestor")
    print("  plot-tumour          Generate notebooks or PDFs once results exist")
    print("")


def create_logger(name: str, log_dir: Optional[str] = None) -> logging.Logger:
    """
    Create a named logger with both console and file handlers.

    Args:
        name: Name for the logger
        log_dir: Optional directory for log files (defaults to current directory)

    Returns:
        Configured logger instance
    """
    # create logger
    logger = logging.getLogger(name)
    # check for active handlers
    if logger.handlers:
        return logger

    logger.setLevel(logging.INFO)
    formatter = logging.Formatter(
        fmt="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # File handler
    log_time = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = f"{name}_log_{log_time}.log"

    if log_dir:
        os.makedirs(log_dir, exist_ok=True)
        log_path = os.path.join(log_dir, log_filename)
    else:
        log_path = log_filename

    file_handler = logging.FileHandler(log_path)
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    logger.info(f"Logger '{name}' initialized. Log file: {log_path}")

    return logger


def split_to_segments(tumour_dir: str) -> list[str]:
    segments_dir_path = f"{tumour_dir}/segments"
    df_path = f"{tumour_dir}/ALPACA_input_table.csv"
    os.makedirs(segments_dir_path, exist_ok=True)
    df = pd.read_csv(df_path)
    tumour_id = df["tumour_id"].iloc[0]
    assert len(
        df["tumour_id"].unique()
    ), "Found multiple tumour ids. In tummour mode only one tumour_id is allowed per input csv"
    segments = []
    for segment, segment_df in df.groupby("segment"):
        segment_df_path = (
            f"{segments_dir_path}/ALPACA_input_table_{tumour_id}_{segment}.csv"
        )
        segment_df.to_csv(segment_df_path, index=False)
        segments.append(segment_df_path)
    return segments


def concatenate_output(output_dir: str) -> str:
    logger = logging.getLogger("ALPACA")
    # keep only segment files in output files list
    output_files = [
        f
        for f in os.listdir(output_dir)
        if f.endswith(".csv") and (("optimal" in f) or ("all" in f))
    ]
    dfs = [pd.read_csv(f"{output_dir}/{f}") for f in output_files]
    concatenated_df = pd.concat(dfs)
    tumour_id = concatenated_df["tumour_id"].iloc[0]
    output_name = f"{output_dir}/ALPACA_output_{tumour_id}.csv"
    concatenated_df.to_csv(output_name, index=False)
    if os.path.exists(output_name):
        logger.info(f"Combined output saved to {output_name}")
        # remove segment files
        for f in output_files:
            if f != os.path.basename(output_name):
                os.remove(f"{output_dir}/{f}")
    else:
        logger.error(f"Failed to save combined output to {output_name}")
        raise FileNotFoundError(f"Output file not found: {output_name}")
    return output_name


def set_run_mode(config: dict) -> tuple[dict, str]:
    run_mode = config["preprocessing_config"]["mode"]
    if run_mode == "tumour":
        print("Running in tumour mode")
        # create segment files:
        config["preprocessing_config"]["input_files"] = [
            Path(x).name
            for x in split_to_segments(
                config["preprocessing_config"]["input_tumour_directory"]
            )
        ]
        config["preprocessing_config"]["input_data_directory"] = Path(
            config["preprocessing_config"]["input_tumour_directory"]
        ).parent
    return config, run_mode


def read_tree_json(json_path: str) -> list[list[str]]:
    """Load a tree either from a JSON list-of-lists or a Newick file.

    The caller passes the expected JSON path (e.g. tree_paths.json). If that file
    does not exist, a sibling `tree_paths.nwk` file is used instead. Supplying a
    direct path to a `.nwk` file is also supported.
    """
    path = Path(json_path)
    # Accept direct Newick path
    if path.suffix == ".nwk" and path.exists():
        return _read_tree_newick(path)

    if path.exists():
        with open(path, "r") as f:
            tree = json.load(f)
        return tree

    fallback = path.with_suffix(".nwk")
    if fallback.exists():
        return _read_tree_newick(fallback)

    raise FileNotFoundError(
        f"Tree file not found. Expected either '{path}' or '{fallback}' to exist."
    )


def _read_tree_newick(nwk_path: Path) -> list[list[str]]:
    """Parse a Newick file into the internal list-of-paths representation."""
    newick_content = Path(nwk_path).read_text().strip()
    return newick_to_paths(newick_content)


def newick_to_paths(newick: str) -> list[list[str]]:
    """Convert a Newick string into ALPACA's list-of-paths format."""
    newick = newick.strip()
    if not newick:
        raise ValueError("Empty Newick string provided")
    if newick.endswith(";"):
        newick = newick[:-1]

    def consume_whitespace(idx: int) -> int:
        while idx < len(newick) and newick[idx].isspace():
            idx += 1
        return idx

    def parse_label(idx: int) -> tuple[str, int]:
        label_chars = []
        while idx < len(newick) and newick[idx] not in ",()":
            label_chars.append(newick[idx])
            idx += 1
        label = "".join(label_chars).strip()
        if label.startswith("'") and label.endswith("'"):
            label = label[1:-1]
        if label.startswith('"') and label.endswith('"'):
            label = label[1:-1]
        if ":" in label:
            label = label.split(":", 1)[0]
        return label, idx

    def parse_subtree(idx: int) -> tuple[dict, int]:
        idx = consume_whitespace(idx)
        if idx >= len(newick):
            raise ValueError("Unexpected end of Newick string")

        if newick[idx] == "(":
            idx += 1
            children = []
            idx = consume_whitespace(idx)
            while True:
                child, idx = parse_subtree(idx)
                children.append(child)
                idx = consume_whitespace(idx)
                if idx >= len(newick):
                    raise ValueError("Unexpected end of Newick string while parsing children")
                if newick[idx] == ",":
                    idx += 1
                    idx = consume_whitespace(idx)
                    continue
                if newick[idx] == ")":
                    idx += 1
                    break
                raise ValueError(f"Unexpected character '{newick[idx]}' in Newick content")

            idx = consume_whitespace(idx)
            label, idx = parse_label(idx)
            if not label:
                raise ValueError("Newick internal nodes must have names for ALPACA conversion")
            return {"name": label, "children": children}, idx

        label, idx = parse_label(idx)
        if not label:
            raise ValueError("Newick leaves must have names for ALPACA conversion")
        return {"name": label, "children": []}, idx

    def collect_paths(node: dict, prefix: list[str]) -> list[list[str]]:
        current = prefix + [node["name"]]
        if not node["children"]:
            return [current]
        paths: list[list[str]] = []
        for child in node["children"]:
            paths.extend(collect_paths(child, current))
        return paths

    start_idx = consume_whitespace(0)
    root, idx_after = parse_subtree(start_idx)
    idx_after = consume_whitespace(idx_after)
    # Handle chained wrappers like ")parent1)parent2" that wrap the parsed tree
    # into higher-level unary parents. This occurs when the Newick has more
    # closing parentheses than were opened at the start, e.g.
    # "(child)parent1)parent2;". We iteratively wrap the current root inside
    # the newly parsed parent node until the string is fully consumed.
    while idx_after < len(newick):
        if newick[idx_after] != ")":
            raise ValueError("Unexpected trailing characters in Newick string")
        idx_after += 1
        idx_after = consume_whitespace(idx_after)
        parent_label, idx_after = parse_label(idx_after)
        if not parent_label:
            raise ValueError("Newick internal nodes must have names for ALPACA conversion")
        root = {"name": parent_label, "children": [root]}
        idx_after = consume_whitespace(idx_after)

    return collect_paths(root, [])


def find_path_edges(branch, tree_edges):
    branch_edges = []
    for edge in tree_edges:
        if (edge[0] in branch) and (edge[1] in branch):
            branch_edges.append(edge)
    return set(branch_edges)


def get_tree_edges(tree):
    all_edges = list()
    for branch in tree:
        if len(branch) == 2:
            all_edges.append(tuple(branch))
        else:
            for i in range(len(branch) - 1):
                all_edges.append((branch[i], branch[i + 1]))
    unique_edges = set(all_edges)
    return unique_edges


def find_parent(tree, clone_name):
    for branch in tree:
        if branch[0] == clone_name:
            return "diploid"
        if clone_name in branch:
            clone_index = branch.index(clone_name)
            return branch[clone_index - 1]


def flat_list(target_list):
    if isinstance(target_list[0], list):
        return [item for sublist in target_list for item in sublist]
    else:
        return target_list


def get_length_from_name(segment):
    e = int(segment.split("_")[-1])
    s = int(segment.split("_")[-2])
    return e - s


def print_logo():
    print(
        """
     _____ __    _____ _____ _____ _____
    |  _  |  |  |  _  |  _  |     |  _  |
    |     |  |__|   __|     |   --|     |
    |__|__|_____|__|  |__|__|_____|__|__|
    /\\⌒⌒⌒/\\
    (⦿   ⦿)
    ( 'Y' )
     (   )
     (   )
     (   )
     (~ ~~~~~~~~~~)
     ( ~ ~~   ~~  )
     ( ~  ~ ~  ~  )
     (~  ~~~~~   ~)
     │ │      │ │
     │ │      │ │
    """
    )


def save_dataframe_to_csv(df: pd.DataFrame, output_dir: str, output_filename: str):
    """
    Save a DataFrame to a CSV file.

    Args:
        df: DataFrame to save
        output_dir: Path to directory where the CSV file will be saved
        output_filename: Name of the output CSV file
    """
    output_path = os.path.join(output_dir, output_filename)
    os.makedirs(os.path.dirname(output_dir), exist_ok=True)
    df.to_csv(output_path, index=False)


def process_ci_reports(dirpath: str, delete: bool = False, outpath: Optional[str] = None) -> pd.DataFrame:
    pattern = os.path.join(dirpath, "*_ci_report.json")
    files = sorted(glob.glob(pattern))
    if outpath is None:
        outpath = os.path.join(dirpath, "ci_modified_report.csv")

    header = [
        "tumour_id",
        "segment",
        "affected_sample",
        "affected_allele",
        "min_ci",
        "timestamp",
        "source_file",
    ]

    # write header even if no files
    with open(outpath, "w", newline="") as csvf:
        writer = csv.writer(csvf)
        writer.writerow(header)
        for fp in files:
            try:
                with open(fp) as fh:
                    data = json.load(fh)
            except Exception as e:
                print(f"Warning: failed to read {fp}: {e}", file=sys.stderr)
                continue

            tumour = data.get("tumour_id")
            segment = data.get("segment")
            min_ci = data.get("min_ci")
            ts = data.get("timestamp")
            affected_samples = data.get("affected_samples") or []
            affected_alleles = data.get("affected_alleles") or []

            # produce a row per sample x allele combination
            if not affected_samples and not affected_alleles:
                writer.writerow(
                    [tumour, segment, "", "", min_ci, ts, os.path.basename(fp)]
                )
                continue

            for s in affected_samples:
                if affected_alleles:
                    for a in affected_alleles:
                        writer.writerow(
                            [tumour, segment, s, a, min_ci, ts, os.path.basename(fp)]
                        )
                else:
                    writer.writerow(
                        [tumour, segment, s, "", min_ci, ts, os.path.basename(fp)]
                    )

            if delete and bool(files):
                try:
                    os.remove(fp)
                except Exception as e:
                    print(f"Warning: failed to remove {fp}: {e}", file=sys.stderr)


def process_monoclonal_reports(dirpath: str, delete: bool = False, outpath: Optional[str] = None) -> pd.DataFrame:
    pattern = os.path.join(dirpath, "*_monoclonal_samples_report.csv")
    files = sorted(glob.glob(pattern))
    if outpath is None:
        combined_df = os.path.join(dirpath, "monoclonal_samples_report.csv")
    # each file is a csv with header, concatenate them or make an empty dataframe
    if not files:
        combined_df = pd.DataFrame(columns=['tumour_id', 'sample', 'segment', 'cpnA', 'cpnB', 'distance_to_integer_A', 'distance_to_integer_B'])
    else:
        dfs = []
        for fp in files:
            try:
                df = pd.read_csv(fp)
                dfs.append(df)
            except Exception as e:
                print(f"Warning: failed to read {fp}: {e}", file=sys.stderr)
                continue
        if dfs:
            combined_df = pd.concat(dfs, ignore_index=True)
    combined_df.to_csv(outpath, index=False)
    if delete and bool(files):
        try:
            os.remove(fp)
        except Exception as e:
            print(f"Warning: failed to remove {fp}: {e}", file=sys.stderr)


def process_elbow_increase_reports(dirpath: str, delete: bool = False, outpath: Optional[str] = None) -> pd.DataFrame:
    pattern = os.path.join(dirpath, "*_elbow_increase_report.csv")
    files = sorted(glob.glob(pattern))
    if outpath is None:
        combined_df = os.path.join(dirpath, "elbow_increase_report.csv")
    # each file is a csv with header, concatenate them or make an empty dataframe
    if not files:
        combined_df = pd.DataFrame(columns=['complexity', 'D_score', 'CI_score', 'allowed_complexity', 'issue', 'tumour_id', 'segment'])
    else:
        dfs = []
        for fp in files:
            try:
                df = pd.read_csv(fp)
                dfs.append(df)
            except Exception as e:
                print(f"Warning: failed to read {fp}: {e}", file=sys.stderr)
                continue
        if dfs:
            combined_df = pd.concat(dfs, ignore_index=True)
    combined_df.to_csv(outpath, index=False)
    if delete and bool(files):
        try:
            os.remove(fp)
        except Exception as e:
            print(f"Warning: failed to remove {fp}: {e}", file=sys.stderr)


def process_run_summary_reports(dirpath: str, delete: bool = False, outpath: Optional[str] = None) -> pd.DataFrame:
    """Combine segment-level run gap summary reports into a single file.

    Only includes segments with non-zero gap (non-optimal solutions).
    Reports gap value and reason (time_limit, gap_tolerance, or other).
    """
    pattern = os.path.join(dirpath, "*_run_gap_summary.csv")
    files = sorted(glob.glob(pattern))
    if outpath is None:
        outpath = os.path.join(dirpath, "run_gap_summary.csv")
    
    # Concatenate all segment reports
    if not files:
        combined_df = pd.DataFrame(columns=[
            'tumour_id', 'segment', 'max_gap', 'gap_reason', 
            'runtime', 'optimal_complexity', 'strict_gap_enabled'
        ])
    else:
        dfs = []
        for fp in files:
            try:
                df = pd.read_csv(fp)
                dfs.append(df)
            except Exception as e:
                print(f"Warning: failed to read {fp}: {e}", file=sys.stderr)
                continue
        if dfs:
            combined_df = pd.concat(dfs, ignore_index=True)
        else:
            combined_df = pd.DataFrame(columns=[
                'tumour_id', 'segment', 'max_gap', 'gap_reason', 
                'runtime', 'optimal_complexity', 'strict_gap_enabled'
            ])

    if combined_df.empty:
        message = ["All segments reached optimality; no gaps to report"]
        pd.DataFrame([message]).to_csv(outpath, index=False, header=False)
    else:
        combined_df.to_csv(outpath, index=False)
    
    if delete and files:
        for fp in files:
            try:
                os.remove(fp)
            except Exception as e:
                print(f"Warning: failed to remove {fp}: {e}", file=sys.stderr)


class SegmentInfeasibleError(RuntimeError):
    """Raised when a segment's ILP model has no feasible solution.

    Callers should catch this, skip the offending segment, and continue
    processing the remaining segments.
    """

    def __init__(self, tumour_id: str, segment: str, solver_status=None):
        self.tumour_id = tumour_id
        self.segment = segment
        self.solver_status = solver_status
        super().__init__(
            f"Segment infeasible: tumour={tumour_id}, segment={segment}, "
            f"solver_status={solver_status}"
        )


def process_infeasibility_reports(
    dirpath: str, delete: bool = False, outpath: Optional[str] = None
) -> pd.DataFrame:
    """Combine per-segment infeasibility JSON reports into a single CSV.

    Each report is written by the model backend when a segment yields an
    infeasible ILP and is named ``<tumour_id>_<segment>_infeasibility_report.json``.
    The combined CSV has one row per infeasible segment.
    """
    pattern = os.path.join(dirpath, "*_infeasibility_report.json")
    files = sorted(glob.glob(pattern))
    if outpath is None:
        outpath = os.path.join(dirpath, "infeasibility_report.csv")

    header = [
        "tumour_id",
        "segment",
        "solver_status",
        "timestamp",
        "source_file",
    ]

    with open(outpath, "w", newline="") as csvf:
        writer = csv.writer(csvf)
        writer.writerow(header)
        for fp in files:
            try:
                with open(fp) as fh:
                    data = json.load(fh)
            except Exception as e:
                print(f"Warning: failed to read {fp}: {e}", file=sys.stderr)
                continue
            writer.writerow(
                [
                    data.get("tumour_id", ""),
                    data.get("segment", ""),
                    data.get("solver_status", ""),
                    data.get("timestamp", ""),
                    os.path.basename(fp),
                ]
            )

    if delete and files:
        for fp in files:
            try:
                os.remove(fp)
            except Exception as e:
                print(f"Warning: failed to remove {fp}: {e}", file=sys.stderr)

    # return a DataFrame for callers that want programmatic access
    try:
        return pd.read_csv(outpath)
    except Exception:
        return pd.DataFrame(columns=header)
