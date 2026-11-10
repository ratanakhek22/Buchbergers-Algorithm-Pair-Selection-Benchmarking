import argparse
import sys

STRATEGIES = ["first", "random", "normal", "degree", "sugar"]
ORDERS = ["lex", "grlex"]

def main():
    parser = argparse.ArgumentParser(description="Run Buchberger's algorithm with a chosen pair selection strategy.")

    # Required
    parser.add_argument("-tc", "--testcase", type=str, required=True, help="Path to the test case (input ideal) file.")

    # Optional
    parser.add_argument("-r", "--strategy", type=str, default="normal", choices=STRATEGIES, help="Pair selection strategy (default: normal).")
    parser.add_argument("-o", "--order", type=str, default="lex", choices=ORDERS, help="Monomial order (default: lex).")
    parser.add_argument("-s", "--seed", type=int, default=0, help="Random seed, used by the random strategy (default: 0).")

    args = parser.parse_args()

    args = parser.parse_args()

    print(f"Testing {args.testcase} as input, using rule: {args.strategy}")
    return 0


if __name__ == "__main__":
    sys.exit(main())