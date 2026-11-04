"""Per-iteration timing report for PWSCF (``relax`` / ``scf``) output files.

Answers the two questions that come up every time a run script or a mixing
parameter is changed:

1. **How long did each SCF iteration take?**  PWSCF prints a running
   ``total cpu time spent up to now`` counter after every iteration, but never
   the increment, so the cost of a single iteration has to be reconstructed by
   differencing consecutive counters.  That increment is the number to compare
   when deciding whether an ``-nk`` (k-point pool) change on the ``mpiexec``
   line actually helped.
2. **Did the change reduce the number of SCF iterations needed to converge?**
   That is the number to compare when judging a ``mixing_beta`` change.

Both are reported per SCF block -- a ``relax`` run restarts the iteration
counter for every BFGS geometry step -- next to a side-by-side comparison
against a baseline for the same system.

A previous run is picked up as the baseline automatically when an output file
has a sibling whose name starts with its own name plus ``.failed_`` (the
convention used when preserving a crashed run before re-running it).  A
baseline can also be named explicitly with ``--baseline``.

Caveat: PWSCF reports *cpu* time per process, not wall time.  For these runs
the two differ by only ~2 % (``1d13h 9m CPU`` against ``1d13h47m WALL``), so the
per-iteration increment is a good proxy for wall time.  The definitive
wall-clock figure for a finished job is the final ``PWSCF : ... WALL`` line,
which is reported here as well.

Usage::

    python scf_timing_check.py                     # scan the Pt4 adsorption tree
    python scf_timing_check.py <file-or-dir> ...   # one or more outputs or folders
    python scf_timing_check.py --iters 10          # cap the per-iteration table
    python scf_timing_check.py --csv timing.csv    # dump every iteration to CSV

While a job is still running::

    watch -n 300 python scf_timing_check.py
"""

import argparse
import csv
import re
from pathlib import Path
from statistics import mean, median

PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_SCAN_DIR = PROJECT_ROOT / "05_metals" / "Pt4" / "Tetra" / "Adsorbed" / "3_layers"

OUTPUT_GLOB = "OUTPUT_*.out*"
BASELINE_MARKER = ".failed_"

# --- header / echo of the input -------------------------------------------
STARTS_RE = re.compile(r"Program PWSCF v\.(\S+)\s+starts on\s+(\S+)\s+at\s+(\S+)")
NPROC_RE = re.compile(r"running on\s+(\d+)\s+processors")
NODES_RE = re.compile(r"MPI processes distributed on\s+(\d+)\s+nodes")
NPOOL_RE = re.compile(r"npool\s*=\s*(\d+)")
MIXBETA_RE = re.compile(r"mixing beta\s*=\s*([\d.]+)")
MIXNDIM_RE = re.compile(r"number of iterations used\s*=\s*(\d+)\s+(\S+)\s+mixing")
NSTEP_RE = re.compile(r"^\s+nstep\s*=\s*(\d+)")
CONVTHR_RE = re.compile(r"scf convergence threshold\s*=\s*(\S+)")
NKS_RE = re.compile(r"number of k points=\s*(\d+)")

# --- the SCF loop itself ---------------------------------------------------
ITER_RE = re.compile(
    r"^\s+iteration #\s*(\d+)\s+ecut=\s*([-\d.]+)\s+Ry\s+beta=\s*([-\d.]+)"
)
TIMING_RE = re.compile(r"total cpu time spent up to now is\s+([\d.]+)\s+secs")
ENERGY_RE = re.compile(r"^\s*!?\s*total energy\s+=\s+(-?\d+\.\d+)\s+Ry")
ACCURACY_RE = re.compile(r"estimated scf accuracy\s+<\s+([\d.eE+-]+)\s+Ry")
SCF_OK_RE = re.compile(r"convergence has been achieved in\s+(\d+)\s+iterations")
SCF_FAIL_RE = re.compile(r"convergence NOT achieved after\s+(\d+)\s+iterations")
# Closes the window in which energy/accuracy lines belong to the last
# iteration: everything after this point describes the finished SCF, and its
# ``!    total energy`` line must not be mistaken for an iteration value.
END_SCF_RE = re.compile(r"End of self-consistent calculation|the Fermi energy is")

# --- ionic loop and final accounting --------------------------------------
BFGS_STEPS_RE = re.compile(r"number of bfgs steps\s+=\s+(\d+)")
GRADIENT_RE = re.compile(r"Gradient error\s+=\s+([\d.eE+-]+)\s+Ry/Bohr")
BFGS_CONVERGED_RE = re.compile(r"bfgs converged")
FINAL_SUMMARY_RE = re.compile(r"^PWSCF\s+:\s+(.+?)\s+CPU\s+(.+?)\s+WALL")
SCF_TIME_RE = re.compile(r"^\s+electrons\s+:\s+([\d.]+)s CPU\s+([\d.]+)s WALL")


def parse_output(path):
    """Parse one PWSCF output file into a report dict.

    Alongside the scalar header fields the result carries

    ``events``
        every ``total cpu time spent up to now`` line in file order, each
        tagged with the iteration it belongs to (``None`` for the init and
        post-SCF lines) and its ``dt`` increment over the previous line;
    ``blocks``
        one entry per SCF block, holding that block's per-iteration events
        and its outcome.
    """
    try:
        text = path.read_text(encoding="utf-8", errors="ignore")
    except OSError as exc:
        return {"path": path, "error": str(exc)}

    info = {
        "path": path,
        "error": None,
        "version": None,
        "started": None,
        "nproc": None,
        "nodes": None,
        "npool": None,
        "nks": None,
        "mixing_beta": None,
        "mixing_ndim": None,
        "mixing_mode": None,
        "nstep": None,
        "conv_thr": None,
        "bfgs_steps": None,
        "gradients": [],
        "bfgs_converged": False,
        "final_cpu": None,
        "final_wall": None,
        "scf_cpu": None,
        "scf_wall": None,
        "events": [],
        "blocks": [],
    }

    block = 0
    iteration = None
    window = None  # the event whose energy/accuracy lines are still expected
    blocks = {}

    for line in text.splitlines():
        match = ITER_RE.match(line)
        if match:
            number = int(match.group(1))
            if number == 1 or block == 0:
                block += 1
                blocks.setdefault(
                    block, {"index": block, "iters": [], "outcome": "in progress"}
                )
            iteration = number
            if info["mixing_beta"] is None:
                info["mixing_beta"] = float(match.group(3))  # "beta=" echo
            window = None
            continue

        match = TIMING_RE.search(line)
        if match:
            event = {
                "block": block,
                "iteration": iteration,
                "cum_cpu": float(match.group(1)),
                "dt": None,
                "energy": None,
                "accuracy": None,
            }
            info["events"].append(event)
            window = event
            iteration = None
            continue

        match = SCF_FAIL_RE.search(line)
        if match:
            if block:
                blocks[block]["outcome"] = f"ABORTED after {int(match.group(1))}"
            window = None
            continue

        match = SCF_OK_RE.search(line)
        if match:
            if block:
                blocks[block]["outcome"] = f"converged in {int(match.group(1))}"
            window = None
            continue

        if END_SCF_RE.search(line):
            window = None
            continue

        # --- metadata, captured once per file -----------------------------
        if info["started"] is None:
            match = STARTS_RE.search(line)
            if match:
                info["version"] = match.group(1)
                info["started"] = f"{match.group(2)} {match.group(3)}"
                continue
        if info["nproc"] is None:
            match = NPROC_RE.search(line)
            if match:
                info["nproc"] = int(match.group(1))
                continue
        if info["nodes"] is None:
            match = NODES_RE.search(line)
            if match:
                info["nodes"] = int(match.group(1))
                continue
        if info["npool"] is None:
            match = NPOOL_RE.search(line)
            if match:
                info["npool"] = int(match.group(1))
                continue
        if info["nks"] is None:
            match = NKS_RE.search(line)
            if match:
                info["nks"] = int(match.group(1))
                continue
        if info["mixing_ndim"] is None:
            match = MIXNDIM_RE.search(line)
            if match:
                info["mixing_ndim"] = int(match.group(1))
                info["mixing_mode"] = match.group(2)
                continue
        if info["mixing_beta"] is None:
            match = MIXBETA_RE.search(line)
            if match:
                info["mixing_beta"] = float(match.group(1))
                continue
        if info["nstep"] is None:
            match = NSTEP_RE.match(line)
            if match:
                info["nstep"] = int(match.group(1))
                continue
        if info["conv_thr"] is None:
            match = CONVTHR_RE.search(line)
            if match:
                info["conv_thr"] = match.group(1)
                continue
        match = BFGS_STEPS_RE.search(line)
        if match:
            info["bfgs_steps"] = int(match.group(1))  # keep the latest count
            continue
        if info["final_wall"] is None:
            match = FINAL_SUMMARY_RE.search(line)
            if match:
                info["final_cpu"] = match.group(1)
                info["final_wall"] = match.group(2)
                continue
        if info["scf_wall"] is None:
            match = SCF_TIME_RE.search(line)
            if match:
                info["scf_cpu"] = float(match.group(1))
                info["scf_wall"] = float(match.group(2))
                continue
        match = GRADIENT_RE.search(line)
        if match:
            info["gradients"].append(float(match.group(1)))
            continue
        if BFGS_CONVERGED_RE.search(line):
            info["bfgs_converged"] = True
            continue

        # Energy / accuracy belong to the iteration whose timing line was seen
        # just above.  The post-SCF "!  total energy" line is deliberately not
        # picked up here, because END_SCF_RE has already closed the window.
        if window is not None and window["iteration"] is not None:
            match = ENERGY_RE.match(line)
            if match and window["energy"] is None:
                window["energy"] = float(match.group(1))
                continue
            match = ACCURACY_RE.search(line)
            if match and window["accuracy"] is None:
                window["accuracy"] = float(match.group(1))
                continue

    previous = None
    for event in info["events"]:
        if previous is not None:
            event["dt"] = round(event["cum_cpu"] - previous, 1)
        previous = event["cum_cpu"]
        if event["dt"] is not None and event["iteration"] is not None:
            entry = blocks.get(event["block"])
            if entry is not None:
                entry["iters"].append(event)

    # PWSCF only prints the "K-points division" line when npool > 1.
    if info["npool"] is None and info["nproc"] is not None:
        info["npool"] = 1

    info["init_cpu"] = info["events"][0]["cum_cpu"] if info["events"] else None
    info["blocks"] = [blocks[key] for key in sorted(blocks)]
    return info


def fmt_duration(seconds):
    """Compact human-readable duration, e.g. ``3m10.5s`` or ``1h13m``."""
    if seconds is None:
        return "-"
    if seconds < 90:
        return f"{seconds:.1f}s"
    minutes, rest = divmod(seconds, 60)
    if minutes < 90:
        return f"{int(minutes)}m{rest:04.1f}s"
    hours, minutes = divmod(minutes, 60)
    return f"{int(hours)}h{int(minutes):02d}m"


def fmt_num(value):
    """Format a number with a sensible number of significant digits."""
    if value is None:
        return "-"
    if value == 0:
        return "0"
    if abs(value) < 1e-3 or abs(value) >= 1e5:
        return f"{value:.2e}"
    return f"{value:.4f}".rstrip("0").rstrip(".")


def iteration_times(entry, steady_only=False):
    """Per-iteration ``dt`` values across every SCF block.

    With ``steady_only`` the first iteration of each block is dropped: it
    carries the wavefunction re-read plus the geometry update of the BFGS step,
    so it is not representative of the cost of a single iteration.
    """
    values = []
    for block in entry["blocks"]:
        times = [item["dt"] for item in block["iters"] if item["dt"] is not None]
        values.extend(times[1:] if steady_only else times)
    return values


def parallel_label(entry):
    """One-line description of how the job was laid out."""
    bits = []
    if entry["nproc"]:
        bits.append(f"{entry['nproc']} procs")
    if entry["nodes"]:
        bits.append(f"{entry['nodes']} node(s)")
    if entry["npool"]:
        bits.append(f"npool={entry['npool']}")
    if entry["nks"]:
        bits.append(f"{entry['nks']} k-pts")
    return ", ".join(bits) if bits else "unknown"


def mixing_label(entry):
    """One-line description of the mixing settings."""
    if entry["mixing_beta"] is None:
        return "unknown"
    text = f"beta={entry['mixing_beta']:g}"
    if entry["mixing_ndim"]:
        text += f", {entry['mixing_ndim']}x {entry['mixing_mode'] or 'local-TF'}"
    return text


def status_label(entry):
    """Overall verdict on the run described by ``entry``."""
    if entry["bfgs_converged"]:
        return "bfgs converged"
    outcomes = [block["outcome"] for block in entry["blocks"]]
    if not outcomes:
        return "no SCF iterations yet"
    if any(item.startswith("ABORTED") for item in outcomes):
        return "ABORTED (SCF did not converge)"
    if outcomes[-1] == "in progress":
        return "in progress"
    return outcomes[-1]


def trim_rows(items, limit):
    """Return ``(item, is_ellipsis)`` pairs: a head, an ellipsis, then the tail."""
    if limit <= 0 or len(items) <= limit:
        return [(item, False) for item in items]
    head = max(limit - 2, 1)
    return (
        [(item, False) for item in items[:head]]
        + [(None, True)]
        + [(item, False) for item in items[-1:]]
    )


def print_report(entry, iter_limit=12, all_blocks=False):
    """Print the header, the per-block summary and the per-iteration table."""
    print("=" * 78)
    print(f"FILE  {entry['path']}")
    print("=" * 78)
    if entry["error"]:
        print(f"  could not read file: {entry['error']}\n")
        return

    version = f" (PWSCF v.{entry['version']})" if entry["version"] else ""
    print(f"  started        : {entry['started'] or '?'}{version}")
    print(f"  parallel       : {parallel_label(entry)}")
    print(f"  mixing         : {mixing_label(entry)}")
    print(f"  scf threshold  : {entry['conv_thr'] or '?'}")
    if entry["init_cpu"] is not None:
        print(f"  init cost      : {fmt_duration(entry['init_cpu'])}")
    if entry["bfgs_steps"] is not None:
        print(f"  bfgs steps     : {entry['bfgs_steps']}")
    if entry["gradients"]:
        print(
            f"  gradient error : {fmt_num(entry['gradients'][0])}"
            f" -> {fmt_num(entry['gradients'][-1])} Ry/Bohr"
        )
    print(f"  status         : {status_label(entry)}")
    if entry["final_wall"]:
        print(f"  total runtime  : {entry['final_cpu']} CPU / {entry['final_wall']} WALL")
    print()

    if not entry["blocks"]:
        print("  (no completed SCF iterations parsed yet)\n")
        return

    print("  SCF block summary  (dt = time for one iteration, cpu time)")
    print(
        f"    {'blk':>3} {'iters':>6}  {'outcome':<22}"
        f" {'mean dt':>10} {'min':>10} {'max':>10} {'sum':>10}"
    )
    for block in entry["blocks"]:
        times = [item["dt"] for item in block["iters"] if item["dt"] is not None]
        cells = (
            (fmt_duration(mean(times)), fmt_duration(min(times)),
             fmt_duration(max(times)), fmt_duration(sum(times)))
            if times else ("-", "-", "-", "-")
        )
        print(
            f"    {block['index']:>3} {len(block['iters']):>6}  {block['outcome']:<22}"
            f" {cells[0]:>10} {cells[1]:>10} {cells[2]:>10} {cells[3]:>10}"
        )
    print()

    chosen = entry["blocks"] if all_blocks else entry["blocks"][-1:]
    for block in chosen:
        if not block["iters"]:
            continue
        print(f"  Per-iteration cost -- block {block['index']} ({block['outcome']})")
        print(
            f"    {'iter':>5} {'dt':>10} {'cum':>11}"
            f" {'energy (Ry)':>20} {'scf acc (Ry)':>14}"
        )
        for item, ellipsis in trim_rows(block["iters"], iter_limit):
            if ellipsis:
                print(f"    {'...':>5}")
                continue
            print(
                f"    {item['iteration']:>5} {fmt_duration(item['dt']):>10}"
                f" {fmt_duration(item['cum_cpu']):>11}"
                f" {fmt_num(item['energy']):>20} {fmt_num(item['accuracy']):>14}"
            )
        print()


def print_comparison(current, baseline):
    """Side-by-side table of the metrics that answer 'did the change help?'."""
    cur_all = iteration_times(current)
    base_all = iteration_times(baseline)
    if not cur_all or not base_all:
        print(f"  comparison with {baseline['path'].name}: not enough data yet\n")
        return

    cur_steady = iteration_times(current, steady_only=True)
    base_steady = iteration_times(baseline, steady_only=True)

    def delta(new, old):
        if not new or not old:
            return "-"
        ratio = mean(new) / mean(old)
        sign = "faster" if ratio < 1 else "slower"
        return f"{ratio:.2f}x {sign}"

    def block_counts(entry):
        return ", ".join(str(len(block["iters"])) for block in entry["blocks"]) or "-"

    rows = [
        ("parallelisation", parallel_label(baseline), parallel_label(current), ""),
        ("mixing", mixing_label(baseline), mixing_label(current), ""),
        ("iters per block", block_counts(baseline), block_counts(current), ""),
        ("mean dt", fmt_duration(mean(base_all)), fmt_duration(mean(cur_all)),
         delta(cur_all, base_all)),
        ("steady mean dt", fmt_duration(mean(base_steady)) if base_steady else "-",
         fmt_duration(mean(cur_steady)) if cur_steady else "-",
         delta(cur_steady, base_steady)),
        ("median dt", fmt_duration(median(base_all)), fmt_duration(median(cur_all)),
         delta([median(cur_all)], [median(base_all)])),
        ("total scf time", fmt_duration(sum(base_all)), fmt_duration(sum(cur_all)), ""),
    ]

    print(f"  COMPARISON  vs  {baseline['path'].name}")
    print(f"    {'metric':<17} {'baseline':<38} {'current':<38} {'change'}")
    print(f"    {'-' * 17} {'-' * 38} {'-' * 38} {'-' * 12}")
    for name, old, new, change in rows:
        print(f"    {name:<17} {str(old):<38} {str(new):<38} {change}")

    # Prefer the steady-state mean for the verdict: the first iteration of each
    # block also pays for the wavefunction restart and the BFGS geometry update,
    # which would swamp the parallelisation signal over a short run.
    if base_steady and len(cur_steady) >= 3:
        primary_cur, primary_base, used = cur_steady, base_steady, "steady mean dt"
    else:
        primary_cur, primary_base, used = cur_all, base_all, "mean dt"

    ratio = mean(primary_cur) / mean(primary_base)
    if ratio < 0.97:
        verdict = f"per-iteration cost DOWN {100 * (1 - ratio):.0f}% -- the change helped"
    elif ratio > 1.03:
        verdict = f"per-iteration cost UP {100 * (ratio - 1):.0f}% -- the change did not help"
    else:
        verdict = "per-iteration cost unchanged -- parallelisation made no difference"

    print(f"    => [{used}] {verdict}")
    if len(cur_steady) < 8:
        print(
            f"       PROVISIONAL: only {len(cur_steady)} steady iteration(s) in the "
            "current run; re-check once more have finished."
        )
    print(
        "    note: 'steady mean dt' drops the first iteration of each block, which\n"
        "          also carries the BFGS step and the wavefunction re-read.\n"
    )


def collect_outputs(paths):
    """Resolve the CLI paths into ``(path, explicitly_named)`` pairs.

    Files handed over on the command line are always reported, even when they
    carry the ``.failed_`` marker; inside a scanned directory that marker means
    "baseline", so those are used only for the comparison.
    """
    roots = list(paths) if paths else [DEFAULT_SCAN_DIR]
    found, seen = [], set()

    def add(path, explicit):
        key = path.resolve()
        if key not in seen:
            seen.add(key)
            found.append((path, explicit))

    for root in roots:
        if root.is_dir():
            for item in sorted(root.glob(f"**/{OUTPUT_GLOB}")):
                add(item, False)
        elif root.exists():
            add(root, True)
        else:
            print(f"warning: no such file or directory: {root}")
    return found


def write_csv(destination, entries):
    """Dump every parsed iteration to a tidy CSV, one row per iteration."""
    with destination.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "system", "file", "block", "iteration", "cum_cpu_s", "dt_s",
            "energy_Ry", "scf_accuracy_Ry", "npool", "nodes", "mixing_beta",
        ])
        for entry in entries:
            if entry.get("error"):
                continue
            for block in entry["blocks"]:
                for item in block["iters"]:
                    writer.writerow([
                        entry["path"].parent.name, entry["path"].name,
                        block["index"], item["iteration"], item["cum_cpu"],
                        item["dt"], item["energy"], item["accuracy"],
                        entry["npool"], entry["nodes"], entry["mixing_beta"],
                    ])
    print(f"wrote {destination}")


def main():
    parser = argparse.ArgumentParser(
        description="Per-iteration timing report for PWSCF relax/scf outputs."
    )
    parser.add_argument(
        "paths", nargs="*", type=Path,
        help="output files or directories to scan "
             "(default: the Pt4/Tetra/Adsorbed/3_layers tree)",
    )
    parser.add_argument(
        "--iters", type=int, default=12,
        help="max per-iteration rows shown per SCF block (default: 12)",
    )
    parser.add_argument(
        "--all-blocks", action="store_true",
        help="show the per-iteration table for every SCF block, not just the last",
    )
    parser.add_argument(
        "--no-baseline", action="store_true",
        help="skip the comparison against preserved .failed_* runs",
    )
    parser.add_argument(
        "--csv", type=Path, metavar="FILE",
        help="also write every parsed iteration to FILE as CSV",
    )
    args = parser.parse_args()

    discovered = collect_outputs(args.paths)
    if not discovered:
        print("no OUTPUT_*.out files found")
        return 1

    baselines = [item for item, _ in discovered if BASELINE_MARKER in item.name]
    currents = [
        item for item, explicit in discovered
        if explicit or BASELINE_MARKER not in item.name
    ]

    cache = {}

    def load(path):
        if path not in cache:
            cache[path] = parse_output(path)
        return cache[path]

    for path in currents:
        entry = load(path)
        print_report(entry, iter_limit=args.iters, all_blocks=args.all_blocks)
        if not args.no_baseline:
            for base_path in baselines:
                if base_path.name.startswith(path.name + BASELINE_MARKER):
                    print_comparison(entry, load(base_path))

    if args.csv:
        everything = [load(item) for item in currents + baselines]
        write_csv(args.csv, everything)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
