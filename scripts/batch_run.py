"""Run main.py over many cases / rules / orders / criteria settings, one process per run.

Each run is a separate `python main.py ...` process with a hard timeout, so a crash or a
runaway run only ends that run. Finished runs are skipped (main.py checks its own output
file; this script also reads it to decide what to skip). If a setting fails to finish on
a case of size n, larger n of the same family are skipped for that setting.

Examples (from the repo root):
    python scripts/batch_run.py --cases "scripts/cases/katsura*.json" "scripts/cases/cyclic*.json"
    python scripts/batch_run.py --cases "scripts/cases/binomial_3-20-4_*" --orders grlex --rules normal sugar
    python scripts/batch_run.py --cases "scripts/cases/*.json" --dry-run
"""
import argparse
import glob
import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent          # repo root (where main.py lives)
MAIN = ROOT / "main.py"
DEFAULT_OUT = ROOT / "output"
ALL_RULES = ["first", "random", "normal", "degree", "sugar"]


def run_id(case, order, rule, seed, crit, rep):
    """Must match the run ID format in main.py."""
    return f"{case.stem}|{order}|{rule}|{seed}|crit{crit}|rep{rep}"


def read_records(out_dir, case):
    """All records already written for this case, keyed by run_id."""
    path = out_dir / f"{case.stem}.jsonl"
    records = {}
    if path.exists():
        with open(path) as fh:
            for line in fh:
                if not line.strip():
                    continue
                try:
                    rec = json.loads(line)
                except json.JSONDecodeError:
                    continue          # half-written line from an interrupted run; it will be redone
                records[rec.get("run_id")] = rec
    return records


def append_record(out_dir, case, record):
    """For runs main.py could not record itself (killed / skipped)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / f"{case.stem}.jsonl", "a") as fh:
        fh.write(json.dumps(record) + "\n")


def load_meta(case):
    with open(case) as fh:
        c = json.load(fh)
    return c["family"], c["n"]


def main():
    parser = argparse.ArgumentParser(description="Batch-run Buchberger experiments.")
    parser.add_argument("--cases", nargs="+", required=True, help="Case files or glob patterns.")
    parser.add_argument("--rules", nargs="+", default=ALL_RULES, choices=ALL_RULES)
    parser.add_argument("--orders", nargs="+", default=["lex", "grlex"])
    parser.add_argument("--criteria", nargs="+", type=int, default=[0, 1], choices=[0, 1])
    parser.add_argument("--random-seeds", type=int, default=5, help="Seeds 0..k-1 for the random rule (default: 5).")
    parser.add_argument("--repeats", type=int, default=3, help="Timing repeats per run (default: 3).")
    parser.add_argument("--max-pairs", type=int, default=5000, help="Pair cap passed to main.py (default: 5000).")
    parser.add_argument("--max-seconds", type=float, default=600, help="Time cap passed to main.py (default: 600).")
    parser.add_argument("--kill-margin", type=float, default=60, help="Extra seconds before the hard kill (default: 60).")
    parser.add_argument("-O", "--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--dry-run", action="store_true", help="List the runs without running them.")
    args = parser.parse_args()

    # Collect cases; sort by (family, n) so smaller cases run first (needed for skip-larger).
    cases = sorted({Path(p).resolve() for pat in args.cases for p in glob.glob(pat)})
    if not cases:
        sys.exit("No case files matched.")
    meta = {c: load_meta(c) for c in cases}
    cases.sort(key=lambda c: (meta[c][0], meta[c][1], c.name))

    failed = {}       # (family, order, rule, seed, crit) -> smallest n that did not finish
    counts = {"ran": 0, "done_before": 0, "unfinished": 0, "killed": 0, "skipped": 0, "error": 0}
    t_start = time.perf_counter()

    for case in cases:
        family, n = meta[case]
        records = read_records(args.output, case)
        for order in args.orders:
            for crit in args.criteria:
                for rule in args.rules:
                    seeds = range(args.random_seeds) if rule == "random" else [0]
                    for seed in seeds:
                        setting = (family, order, rule, seed, crit)
                        for rep in range(args.repeats):
                            rid = run_id(case, order, rule, seed, crit, rep)

                            if rid in records:
                                counts["done_before"] += 1
                                if records[rid].get("finished") is False:
                                    failed.setdefault(setting, n)
                                continue

                            # Skip larger n once a smaller n failed (binomials share n, so only strictly larger).
                            if setting in failed and n > failed[setting]:
                                rec = {"run_id": rid, "case": case.stem, "family": family, "n": n,
                                       "order": order, "rule": rule, "seed": seed, "criteria": bool(crit),
                                       "rep": rep, "finished": False,
                                       "reason": f"skipped: n={failed[setting]} did not finish"}
                                if not args.dry_run:
                                    append_record(args.output, case, rec)
                                counts["skipped"] += 1
                                continue

                            cmd = [sys.executable, str(MAIN), 
                                   "-tc", str(case), 
                                   "-r", rule, 
                                   "-o", order, 
                                   "-s", str(seed), 
                                   "-c", str(crit), 
                                   "-R", str(rep), 
                                   "-O", str(args.output), 
                                   "--paircutoff", str(args.max_pairs),
                                   "--timeout", str(int(args.max_seconds))]
                            if args.dry_run:
                                print("would run:", rid)
                                continue

                            print(f"[{time.perf_counter() - t_start:7.0f}s] {rid} ...", end=" ", flush=True)
                            try:
                                proc = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, timeout=args.max_seconds + args.kill_margin)
                            except subprocess.TimeoutExpired:
                                append_record(args.output, case, {
                                    "run_id": rid, "case": case.stem, "family": family, "n": n,
                                    "order": order, "rule": rule, "seed": seed, "criteria": bool(crit),
                                    "rep": rep, "finished": False, "reason": "killed: hard timeout"})
                                failed.setdefault(setting, n)
                                counts["killed"] += 1
                                print("KILLED")
                                continue

                            if proc.returncode != 0:
                                counts["error"] += 1
                                print("ERROR")
                                print(proc.stderr.strip().splitlines()[-1] if proc.stderr.strip() else "(no stderr)")
                                continue

                            rec = read_records(args.output, case).get(rid, {})
                            counts["ran"] += 1
                            if rec.get("finished") is False:
                                failed.setdefault(setting, n)
                                counts["unfinished"] += 1
                                print(f"did not finish ({rec.get('pairs_processed')} pairs)")
                            else:
                                print(f"{rec.get('poly_addition')} additions, {rec.get('time_total_elapsed', 0):.2f}s, "f"correct={rec.get('correct')}")

    print("\nSummary:", ", ".join(f"{k} {v}" for k, v in counts.items()), f"| {time.perf_counter() - t_start:.0f}s total")
    return 0


if __name__ == "__main__":
    sys.exit(main())