import argparse
import sys
import random
import json
import itertools

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

def buchberger_algorithm(R, F, rule, state, progress, track_sugar):
    G = list(F)
    B = list(itertools.combinations(range(len(G)), 2))
    stats = {"pairs_processed": 0, "zero_reductions": 0, "poly_addition": 0, "pair_order": []}
    sugars  = [total_degree(f) for f in F] if track_sugar else None
    state["sugar"] = sugars

    pair = None
    while B:
        k = rule(B, G, state)                           # index of the chosen pair in B
        i, j = B[k]
        pair = (i, j)
        B[k] = B[-1]
        B.pop()
        progress.update(B, G, pair, stats)

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
            state["sugar"].append(sugar)
            t = len(G) - 1
            B.extend((l, t) for l in range(t))
        else:
            stats["zero_reductions"] += 1

    progress.close(B, G, pair, stats)
    return G, stats

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
    parser.add_argument("-tc", "--testcase", type=str, required=True, help="Path to the test case (input ideal) file.")

    # Optional
    parser.add_argument("-v", "--verbose", action="store_true", help="Show live progress.")
    parser.add_argument("-r", "--rule", type=str, default="normal", choices=RULES, help="Pair selection strategy (default: normal).")
    parser.add_argument("-o", "--order", type=str, default="lex", choices=ORDERS, help="Monomial order (default: lex).")
    parser.add_argument("-s", "--seed", type=int, default=0, help="Random seed, used by the random strategy (default: 0).")

    args = parser.parse_args()

    print(f"Beginning Test with rule: {args.rule}")

    # Setup
    R, F, case = load_case(args.testcase, args.order)
    rule = RULES_MAP[args.rule]
    track_sugar = (args.rule == "sugar")
    state = {"ring": R, "rng": random.Random(args.seed)}

    progress = Progress(enabled=args.verbose)
    G, stats = buchberger_algorithm(R, F, rule, state, progress, track_sugar)

    # TEMP
    print(stats)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())