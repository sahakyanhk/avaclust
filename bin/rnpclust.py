#!/usr/bin/env python3
"""Python entrypoint for the rnpclust pipeline."""

import argparse
import os
import subprocess
import sys
from pathlib import Path
from typing import Optional, Sequence

DEFAULT_OUTDIR = "rnpclust_out"
DEFAULT_CUTOFF = "auto"
DEFAULT_METHOD = "0"
DEFAULT_CHAINS = "A,B"


def strip_trailing_slash(value: str) -> str:
    return value[:-1] if value.endswith("/") else value


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the rnpclust structural clustering pipeline.",
        epilog=(
            "Method notes:\n"
            "  0: hierarchical clustering; omit --cutoff for auto-detection\n"
            "  1: greedy set-cover clustering; requires --cutoff"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        allow_abbrev=False,
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="pdb_dir",
        required=True,
        help="Input directory with .pdb files (required)",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="out_dir",
        default=DEFAULT_OUTDIR,
        help=f"Output directory (default: {DEFAULT_OUTDIR})",
    )
    parser.add_argument(
        "-c",
        "--cutoff",
        default=DEFAULT_CUTOFF,
        help="TM-score cutoff (default: auto)",
    )
    parser.add_argument(
        "-I",
        "--interface-cutoff",
        dest="interface_cutoff",
        default="",
        help="Cutoff in Angstroms for extracting interface residues between chains",
    )
    parser.add_argument(
        "-m",
        "--method",
        default=DEFAULT_METHOD,
        help="Clustering method: 0 = hierarchical, 1 = greedy set-cover (default: 0)",
    )
    parser.add_argument(
        "-a",
        "--chains",
        default=DEFAULT_CHAINS,
        help=f"Chains to consider for alignment (default: {DEFAULT_CHAINS})",
    )
    parser.add_argument(
        "-t",
        "--threads",
        default="",
        help="Number of threads for parallel USalign (default: all cores)",
    )
    return parser


def validate_args(args: argparse.Namespace) -> None:
    args.pdb_dir = strip_trailing_slash(args.pdb_dir)
    args.out_dir = strip_trailing_slash(args.out_dir)

    if not args.pdb_dir:
        print("ERROR: -i <pdb_dir> is required", file=sys.stderr)
        raise SystemExit(1)

    if not Path(args.pdb_dir).is_dir():
        print(f"ERROR: {args.pdb_dir} is not a directory", file=sys.stderr)
        raise SystemExit(1)

    if args.method == "1" and args.cutoff == DEFAULT_CUTOFF:
        print("ERROR: set-cover clustering (-m 1) requires a cutoff (-c <value>)", file=sys.stderr)
        raise SystemExit(1)


def print_run_summary(args: argparse.Namespace) -> None:
    method_label = "hierarchical" if args.method == "0" else "set-cover"
    thread_count = args.threads or str(os.cpu_count() or 1)

    print(f"PDB dir:    {args.pdb_dir}", file=sys.stderr)
    print(f"Output dir: {args.out_dir}", file=sys.stderr)
    print(f"Method:     {method_label}", file=sys.stderr)
    print(f"Cutoff:     {args.cutoff}", file=sys.stderr)
    print(f"ICutoff:    {args.interface_cutoff}", file=sys.stderr)
    print(f"Chains:     {args.chains}", file=sys.stderr)
    print(f"Threads:    {thread_count}", file=sys.stderr)


def run_command(cmd: Sequence[str]) -> None:
    subprocess.run(list(cmd), check=True)


def maybe_extract_interface(
    pdb_dir: str,
    interface_cutoff: str,
    script_dir: Path,
) -> str:
    if not interface_cutoff:
        return pdb_dir

    print(
        f"=== Extracting interface residues with cutoff {interface_cutoff} Angstrom ===",
        file=sys.stderr,
    )
    output_dir = f"{pdb_dir}_i{interface_cutoff}"
    run_command(
        [
            sys.executable,
            str(script_dir / "extract_interface.py"),
            "-i",
            pdb_dir,
            "-c",
            "A",
            "--cut",
            interface_cutoff,
            "-o",
            output_dir,
        ]
    )
    return output_dir


def count_pdb_files(pdb_dir: str) -> int:
    return sum(1 for path in Path(pdb_dir).glob("*.pdb") if path.is_file())


def count_lines(path: Path) -> int:
    with path.open("rb") as handle:
        return sum(1 for _ in handle)


def ensure_all_vs_all(
    pdb_dir: str,
    out_dir: Path,
    chains: str,
    threads: str,
    script_dir: Path,
) -> Path:
    allvsall = out_dir / "all_vs_all.tsv"
    nstructs = count_pdb_files(pdb_dir)
    expected_lines = nstructs * (nstructs - 1) // 2 + 1

    if allvsall.is_file() and count_lines(allvsall) == expected_lines:
        print(f"=== Reusing existing {allvsall} ({expected_lines} lines) ===", file=sys.stderr)
        return allvsall

    print("=== Running all-vs-all alignments ===", file=sys.stderr)
    # Preserve the shell entrypoint behavior: an empty threads argument lets
    # usalign_all_vs_all.bash fall back to all available cores.
    run_command(
        [
            "bash",
            str(script_dir / "usalign_all_vs_all.bash"),
            pdb_dir,
            str(allvsall),
            chains,
            threads,
        ]
    )
    return allvsall


def save_histogram(allvsall: Path, cutoff: str, out_dir: Path) -> None:
    import pandas as pd
    from matplotlib import pyplot as plt

    df = pd.read_csv(allvsall, sep="\t")
    plt.style.use("bmh")
    plt.hist(df.TM, bins=50, alpha=0.8, label=f"cutoff={cutoff}")
    plt.xlabel("TM Score")
    plt.ylabel("Frequency")
    plt.title("Distribution of TM Scores")
    plt.savefig(out_dir / f"all_vs_all_i{cutoff}.png")
    plt.close()


def run_clustering(
    allvsall: Path,
    pdb_dir: str,
    out_dir: Path,
    cutoff: str,
    method: str,
    script_dir: Path,
) -> None:
    clusters_path = out_dir / "clusters.dat"
    aligned_clusters_dir = out_dir / "aligned_clusters"

    if method == "0":
        cmd = [
            sys.executable,
            str(script_dir / "hierarchical_cluster.py"),
            str(allvsall),
            "--col",
            "2",
        ]
        if cutoff != DEFAULT_CUTOFF:
            print(f"=== Hierarchical clustering (cutoff={cutoff}) ===", file=sys.stderr)
            cmd.extend(["--cutoff", cutoff])
        else:
            print("=== Hierarchical clustering (auto k) ===", file=sys.stderr)
        cmd.extend(
            [
                "--pdb-dir",
                pdb_dir,
                "-o",
                str(clusters_path),
                "-sc",
                str(aligned_clusters_dir),
            ]
        )
        run_command(cmd)
        return

    print(f"=== Set-cover clustering (cutoff={cutoff}) ===", file=sys.stderr)
    run_command(
        [
            sys.executable,
            str(script_dir / "setcover_cluster.py"),
            str(allvsall),
            "--col",
            "2",
            "--cutoff",
            cutoff,
            "--centroid",
            "--fmt",
            "wide",
            "--pdb-dir",
            pdb_dir,
            "-o",
            str(clusters_path),
            "-sc",
            str(aligned_clusters_dir),
        ]
    )


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(args)
    print_run_summary(args)

    script_dir = Path(__file__).resolve().parent
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    pdb_dir = maybe_extract_interface(args.pdb_dir, args.interface_cutoff, script_dir)
    allvsall = ensure_all_vs_all(pdb_dir, out_dir, args.chains, args.threads, script_dir)
    save_histogram(allvsall, args.cutoff, out_dir)
    run_clustering(allvsall, pdb_dir, out_dir, args.cutoff, args.method, script_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
