import argparse
import sys
import json
import random

from pathlib import Path
from sympy import sympify, expand
from math import comb

DEFAULT_OUT = Path(__file__).resolve().parent     # ./

def generate_cyclic(n):
    polys = []
    xs = [f"x_{i}" for i in range(0, n)]
    for k in range(1, n):
        terms = ["*".join(xs[(i + j) % n] for j in range(k)) for i in range(n)]
        polys.append(" + ".join(terms))
    polys.append("*".join(xs) + " - 1")
    return xs, polys

def generate_katsura(n):
    def U(l):
        l = abs(l)                          # mirror: u_{-l} = u_l
        return us[l] if l <= n else None    # out of range: drop the term

    polys = []
    us = [f"u_{i}" for i in range(0, n + 1)]
    for m in range(n):
        terms = [f"{U(l)}*{U(m - l)}" for l in range(-n, n + 1) if U(l) and U(m - l)]
        polys.append(" + ".join(terms) + f" - {us[m]}")
    polys.append(" + ".join([us[0]] + [f"2*{u}" for u in us[1:]]) + " - 1")
    return us, polys

def _monomial_of_degree(rng, n, k):
    bars = sorted(rng.sample(range(k + n - 1), n - 1))
    cuts = [-1] + bars + [k + n - 1]
    return tuple(cuts[i + 1] - cuts[i] - 1 for i in range(n))

def _random_monomial(rng, n, d, dist):
    degrees = range(1, d + 1)                       # no constant terms
    if dist == "weighted":                          # degree uniform, then monomial uniform within it
        k = rng.choice(degrees)
    else:                                           # uniform over all monomials of degree 1..d
        k = rng.choices(degrees, weights=[comb(k + n - 1, n - 1) for k in degrees])[0]
    return _monomial_of_degree(rng, n, k)
 
def _term_str(coeff, expv, variables):
    factors = [v if e == 1 else f"{v}**{e}" for v, e in zip(variables, expv) if e > 0]
    return f"{coeff}*" + "*".join(factors)

def generate_binomial(n, d, s, dist, seed, modulus=32003):
    rng = random.Random(seed)
    variables = [f"x_{i}" for i in range(n)]
    polys = []
    for _ in range(s):
        a = _random_monomial(rng, n, d, dist)
        b = _random_monomial(rng, n, d, dist)
        while b == a:                                                       # two distinct terms, or it isn't a binomial
            b = _random_monomial(rng, n, d, dist)
        c1, c2 = rng.randrange(1, modulus), rng.randrange(1, modulus)       # nonzero coefficients
        polys.append(f"{_term_str(c1, a, variables)} + {_term_str(c2, b, variables)}")
    return variables, polys

def write_case(path, family, n, variables, polys, seed=None, modulus=32003, **params):
    case = {
        "family": family,
        "n": n,
        "seed": seed,
        **params,
        "modulus": modulus,
        "variables": variables,
        "polynomials": [str(expand(sympify(p))) for p in polys],
    }
    with open(path, "w") as fh:
        json.dump(case, fh, indent=2)

def main():
    parser = argparse.ArgumentParser(description="Generate test cases (input ideals) as JSON files.")

    # Required
    parser.add_argument("-o", "--output", type=Path, default=DEFAULT_OUT, help="Output directory for /cases (default: scripts/).")
    sub = parser.add_subparsers(dest="family", required=True, help="Input family to generate.")

    # Cyclic-n Options
    p = sub.add_parser("cyclic", help="Cyclic-n.")
    p.add_argument("-n", type=int, nargs="+", required=True, help="One or more n values, e.g. -n 3 4 5.")

    # Katsura-n Options
    p = sub.add_parser("katsura", help="Katsura-n.")
    p.add_argument("-n", type=int, nargs="+", required=True, help="One or more n values, e.g. -n 2 3 4.")

    # Random-Binomial Options
    p = sub.add_parser("binomial", help="Random binomial ideals (Peifer et al.).")
    p.add_argument("-n", type=int, required=True, help="Number of variables.")
    p.add_argument("-d", type=int, required=True, help="Maximum monomial degree.")
    p.add_argument("-s", type=int, required=True, help="Number of generators.")
    p.add_argument("--dist", choices=["weighted", "uniform"], default="weighted", help="Exponent distribution (default: weighted).")
    p.add_argument("--seed", type=int, default=0, help="First seed (default: 0).")
    p.add_argument("--count", type=int, default=1, help="How many ideals: seeds seed..seed+count-1 (default: 1).")

    args = parser.parse_args()

    # generation loop
    out_path = args.output / "cases"
    out_path.mkdir(parents=True, exist_ok=True)      # create <path>/cases if missing

    if args.family == "cyclic":
        for n in args.n:
            variables, polys = generate_cyclic(n)
            write_case(out_path / f"cyclic{n}.json", "cyclic", n, variables, polys)

    elif args.family == "katsura":
        for n in args.n:
            variables, polys = generate_katsura(n)
            write_case(out_path / f"katsura{n}.json", "katsura", n, variables, polys)

    elif args.family == "binomial":
        for seed in range(args.seed, args.seed + args.count):
            variables, polys = generate_binomial(args.n, args.d, args.s, args.dist, seed)
            name = f"binomial_{args.n}-{args.d}-{args.s}_{args.dist}_seed{seed}.json"
            write_case(out_path / name, "binomial", args.n, variables, polys, seed=seed, d=args.d, s=args.s, dist=args.dist)

    return 0

if __name__ == "__main__":
    sys.exit(main())