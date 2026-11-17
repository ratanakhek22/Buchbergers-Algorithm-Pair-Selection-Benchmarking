import argparse
import sys
import random
import json
import itertools
import time

from pathlib import Path
import platform
import subprocess

import sympy
from sympy import sympify, Symbol
from sympy.polys.rings import ring
from sympy.polys.domains import GF
from sympy.polys.orderings import lex, grlex
from sympy.polys.monomials import monomial_lcm, monomial_div

# local imports 
from rules import RULES_MAP
from rules.sugar import total_degree, pair_sugar
from progress import Progress

RULES = list(RULES_MAP)
ORDERS = ["lex", "grlex"]
ORDER_MAP = {"lex": lex, "grlex": grlex}
DEFAULT_OUT = Path(__file__).resolve().parent / "scripts" / "output"

def buchberger_algorithm(R, F, rule, state, progress, track_sugar, use_criterion):
    G = list(F)
    B = list(itertools.combinations(range(len(G)), 2))
    Bset = set(B)
    stats = {"pairs_processed": 0, "pairs_eliminated": 0, "zero_reductions": 0, "poly_addition": 0, "time_rule_overhead": 0, "pair_order": []}
    sugars  = [total_degree(f) for f in F] if track_sugar else None
    state["sugar"] = sugars

    pair = None
    while B:
        start_rule_overhead = time.perf_counter()
        k = rule(B, G, state)                           # index of the chosen pair in B
        end_rule_overhead = time.perf_counter()
        stats["time_rule_overhead"] += (end_rule_overhead - start_rule_overhead)

        i, j = B[k]
        pair = (i, j)
        B[k] = B[-1]
        B.pop()
        progress.update(B, G, pair, stats)

        if use_criterion:
            lcm = monomial_lcm(G[i].LM, G[j].LM)
            coprime = lcm == tuple(a + b for a, b in zip(G[i].LM, G[j].LM))
            if coprime or criterion2(i, j, lcm, G, Bset):
                Bset.discard((i, j))
                stats["pairs_eliminated"] += 1
                continue
        Bset.discard((i, j))

        stats["pairs_processed"] += 1
        stats["pair_order"].append((i + 1, j + 1))      # 1-based, to match the paper

        S = sPolynomial(G[i], G[j])
        stats["poly_addition"] += 1
        s_ij = pair_sugar(i, j, G, sugars) if track_sugar else None
        remainder, div_stats, sugar = division(S, G, sugars, s_ij)
        stats["poly_addition"] += div_stats["poly_addition"]

        if remainder:
            G.append(remainder)
            if track_sugar:
                sugars.append(sugar)
            t = len(G) - 1
            B.extend((l, t) for l in range(t))
        else:
            stats["zero_reductions"] += 1

    progress.close(B, G, pair, stats)
    return G, stats

def criterion2(i, j, lcm, G, Bset):
    for k in range(len(G)):
        if k in (i, j):
            continue
        if (min(i, k), max(i, k)) in Bset or (min(j, k), max(j, k)) in Bset:
            continue                                  # (i,k) or (j,k) not handled yet
        if monomial_div(lcm, G[k].LM) is not None:    # LM(g_k) divides the lcm
            return True
    return False

def sPolynomial(g_i, g_j):
    lcm = monomial_lcm(g_i.LM, g_j.LM)
    a = g_i.mul_term((monomial_div(lcm, g_i.LM), 1 / g_i.LC))
    b = g_j.mul_term((monomial_div(lcm, g_j.LM), 1 / g_j.LC))
    return a - b

def division(dividend, divisors, sugars, sugar):
    p = dividend
    if sugars is None:
        sugars = [None] * len(divisors)
    G = [(g, *g.LT, s) for g, s in zip(divisors, sugars)]   # (g, g_lm, g_lc, g_sugar)
    r = dividend.ring.zero
    stats = {"poly_addition": 0}
    while p:
        lm, lc = p.LT                                       # LT of current p
        div_flag = False                                    # Track if G finds something to divide by

        for g, g_lm, g_lc, g_s in G:
            mult = monomial_div(lm, g_lm)
            if mult is not None:
                div_flag = True
                p = p - g.mul_term((mult, lc / g_lc))
                if sugar is not None:
                    sugar = max(sugar, sum(mult) + g_s)     # deg t + sugar(g_k)
                stats["poly_addition"] += 1
                break

        if not div_flag:                        # If p's LT is not divisible, move it to r
            lt = p.ring({lm: lc})
            r = r + lt
            p = p - lt

    return r, stats, sugar

def git_commit():
    try:
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True, cwd=Path(__file__).resolve().parent).stdout.strip() or None
    except OSError:
        return None

def load_case(path, order_name):
    with open(path) as fh:
        case = json.load(fh)

    # Build ring from file
    R, *gens = ring(",".join(case["variables"]), GF(case["modulus"]), ORDER_MAP[order_name])

    # Note: Symbol() ensures no special characters cause errors
    names = {v: Symbol(v) for v in case["variables"]}
    F = [R.from_expr(sympify(p, locals=names)) for p in case["polynomials"]]
    return R, F, case

def main():
    parser = argparse.ArgumentParser(description="Run Buchberger's algorithm with a chosen pair selection strategy.")

    # Required
    parser.add_argument("-tc", "--testcase", type=Path, required=True, help="Path to the test case (input ideal) file.")

    # Optional
    parser.add_argument("-v", "--verbose", action="store_true", help="Show live progress.")
    parser.add_argument("-r", "--rule", type=str, default="normal", choices=RULES, help="Pair selection strategy (default: normal).")
    parser.add_argument("-o", "--order", type=str, default="lex", choices=ORDERS, help="Monomial order (default: lex).")
    parser.add_argument("-s", "--seed", type=int, default=0, help="Random seed, used by the random strategy (default: 0).")
    parser.add_argument("-c", "--criterion", type=int, default=0, choices=(0, 1), help="Whether to use 2 Criterion from CLO to eliminate pairs.")
    parser.add_argument("-O", "--output", type=Path, default=DEFAULT_OUT, help="Output directory (default: scripts/output).")

    args = parser.parse_args()

    print(f"Beginning Test with rule: {args.rule}")

    # Setup
    R, F, case = load_case(args.testcase, args.order)
    rule = RULES_MAP[args.rule]
    track_sugar = (args.rule == "sugar")
    use_criterion = args.criterion == 1
    state = {"ring": R, "rng": random.Random(args.seed)}
    progress = Progress(enabled=args.verbose)

    run_id = f"{args.testcase.stem}|{args.order}|{args.rule}|{args.seed}|crit{args.criterion}"
    out_file = args.output / f"{args.testcase.stem}.jsonl"
    args.output.mkdir(parents=True, exist_ok=True)
    if out_file.exists():
        with open(out_file) as fh:
            done = {json.loads(line).get("run_id") for line in fh if line.strip()}
        if run_id in done:
            print(f"Run {run_id} already in {out_file}, skipping.")
            return 0

    start = time.perf_counter()     # Start Time

    G, stats = buchberger_algorithm(R, F, rule, state, progress, track_sugar, use_criterion)

    end = time.perf_counter()       # End Time
    stats["time_total_elapsed"] = end - start

    # Reduced basis, used both for its size and as a correctness check
    X = [Symbol(v) for v in case["variables"]]
    ours = groebner([g.as_expr() for g in G], *X, order=args.order, modulus=case["modulus"])
    ref = groebner([f.as_expr() for f in F], *X, order=args.order, modulus=case["modulus"])

    record = {
        "run_id": run_id,
        "case": args.testcase.stem,
        "family": case["family"], "n": case["n"], "case_seed": case.get("seed"),
        "order": args.order, "rule": args.rule, "seed": args.seed, "criteria": use_criterion,
        **stats,
        "time_total_elapsed": end - start,
        "G_size": len(G),
        "G_size_reduced": len(ours.exprs),
        "correct": ours.exprs == ref.exprs,
        "reduced_basis": [str(e) for e in ours.exprs],
        "commit": git_commit(),
        "python": platform.python_version(),
        "sympy": sympy.__version__,
    }
    with open(out_file, "a") as fh:                     # append: one JSON object per line
        fh.write(json.dumps(record) + "\n")

    print(f"{run_id}: {stats['pairs_processed']} pairs, {stats['poly_addition']} additions, "
          f"{record['time_total_elapsed']:.3f}s, correct={record['correct']}")
    return 0

if __name__ == "__main__":
    sys.exit(main())