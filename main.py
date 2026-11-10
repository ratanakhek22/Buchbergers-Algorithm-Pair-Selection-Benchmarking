import argparse
import sys
import random
import json
import itertools

from sympy import sympify, Symbol
from sympy.polys.rings import ring
from sympy.polys.domains import GF
from sympy.polys.orderings import lex, grlex

# importing all rules
from rules import RULES_MAP

RULES = list(RULES_MAP)
ORDERS = ["lex", "grlex"]
ORDER_MAP = {"lex": lex, "grlex": grlex}

def buchberger_algorithm(R, F, rule, state):
    G = list(F)
    B = list(itertools.combinations(range(len(G)), 2))
    stats = {"pairs_processed": 0, "zero_reductions": 0, "pair_order": []}
    
    while B:
        k = rule(B, G, state)       # index of the chosen pair in B
        i, j = B[k]
        B[k] = B[-1]
        B.pop()

        stats["pairs_processed"] += 1
        stats["pair_order"].append((i + 1, j + 1))      # 1-based, to match the paper

        S = sPolynomial(G[i], G[j])
        remainder = division(S, G)

        if remainder:
            G.append(remainder)
            t = len(G) - 1
            B.extend((l, t) for l in range(t))
        else:
            stats["zero_reductions"] += 1

    return G, stats

def sPolynomial(g_i, g_j):
    return 0        # TODO

def division(dividend, divisor):
    return 0        # TODO

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
    parser.add_argument("-r", "--rule", type=str, default="normal", choices=RULES, help="Pair selection strategy (default: normal).")
    parser.add_argument("-o", "--order", type=str, default="lex", choices=ORDERS, help="Monomial order (default: lex).")
    parser.add_argument("-s", "--seed", type=int, default=0, help="Random seed, used by the random strategy (default: 0).")

    args = parser.parse_args()

    # Setup
    random.seed(args.seed)
    R, F, case = load_case(args.testcase, args.order)
    rule = RULES_MAP[args.rule]

    buchberger_algorithm(R, F, rule)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())