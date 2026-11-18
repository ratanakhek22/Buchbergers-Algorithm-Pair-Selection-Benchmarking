"""Summarize the experiment results in output/*.jsonl.

Reads every result line written by main.py and batch_run.py, cleans the data (every
cleaning step is printed with how many lines it touched), and writes:

    <out>/clean_runs.csv          one row per kept run (after cleaning)
    <out>/cleaning_log.txt        what was dropped or changed, and why
    <out>/fixed_cases.csv         Cyclic / Katsura (and the example): one row per case x setting x rule
    <out>/binomial_summary.csv    binomials: one row per setting x order x criteria x rule
    <out>/binomial_wins.csv       binomials: how often each rule had the fewest additions
    <out>/criteria_effect.csv     additions with criteria / without, per rule
    <out>/fixed_cases.png         additions vs n for Cyclic and Katsura
    <out>/binomial_additions.png  median additions per rule for each binomial setting

and prints the main tables. Default <out> is output/summary/ (ignored by git with output/*).

How runs are combined (see the README section of the cleaning log for details):
  * Counts (additions, pairs, ...) are deterministic for a given case/order/rule/seed/criteria,
    so the timing repeats are collapsed: counts from one repeat, time = median over repeats.
  * Random has several seeds per case: its counts are the median over seeds.
  * A setting only counts as "finished" if every repeat (and every Random seed) finished.

Usage (from the repo root):
    python scripts/summarize.py
    python scripts/summarize.py -i output -O output/summary --include-example
"""
import argparse
import json
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RULES = ["first", "random", "normal", "degree", "sugar"]
COUNT_COLS = ["pairs_processed", "pairs_eliminated", "zero_reductions", "poly_addition", "G_size", "G_size_reduced"]

# Categorical colors (fixed order, one per rule) and ink colors for the plots.
RULE_COLORS = {"first": "#2a78d6", "random": "#eb6834", "normal": "#1baf7a", "degree": "#eda100", "sugar": "#e87ba4"}
INK, INK_2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"

RUN_ID = re.compile(r"^(?P<case>[^|]+)\|(?P<order>[^|]+)\|(?P<rule>[^|]+)\|(?P<seed>-?\d+)\|crit(?P<crit>[01])(\|rep(?P<rep>\d+))?$")
BINOMIAL = re.compile(r"^binomial_(?P<n>\d+)-(?P<d>\d+)-(?P<s>\d+)_(?P<dist>weighted|uniform)_seed(?P<cseed>\d+)$")


class Log:
    """Collects cleaning messages so they are printed and saved together."""
    def __init__(self):
        self.lines = []

    def __call__(self, msg):
        print(msg)
        self.lines.append(msg)


# ----------------------------------------------------------------------------------------------
# 1. Load
# ----------------------------------------------------------------------------------------------
def load(in_dir, log):
    rows, bad = [], 0
    files = sorted(in_dir.glob("*.jsonl"))
    for path in files:
        with open(path) as fh:
            for lineno, line in enumerate(fh, 1):
                if not line.strip():
                    continue
                try:
                    rec = json.loads(line)
                except json.JSONDecodeError:
                    bad += 1
                    continue
                rec["_file"], rec["_line"] = path.name, lineno
                rows.append(rec)
    log(f"Loaded {len(rows)} lines from {len(files)} files in {in_dir}.")
    log(f"  [clean 1] dropped {bad} lines that are not valid JSON (half-written lines from interrupted runs).")
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------------------------
# 2. Clean
# ----------------------------------------------------------------------------------------------
def clean(df, include_example, log):
    n0 = len(df)

    # 2. Lines without a run ID (e.g. the old {"TODO": ...} placeholder).
    no_id = df["run_id"].isna() if "run_id" in df else pd.Series(True, index=df.index)
    df = df[~no_id].copy()
    log(f"  [clean 2] dropped {int(no_id.sum())} lines with no run_id (placeholders, not runs).")

    # 3. Parse the run ID; drop old-format IDs from before the |rep suffix was added.
    parts = df["run_id"].str.extract(RUN_ID)
    unparsable = parts["case"].isna()
    legacy = parts["rep"].isna() & ~unparsable
    log(f"  [clean 3] dropped {int(legacy.sum())} old-format run IDs without '|rep' (from before repeats existed; "
        f"those runs were redone with the new IDs) and {int(unparsable.sum())} IDs that don't match the format.")
    df = df[~legacy & ~unparsable].copy()
    parts = parts[~legacy & ~unparsable]
    for col in ["case", "order", "rule"]:
        df[col] = parts[col]
    df["seed"] = parts["seed"].astype(int)
    df["criteria"] = parts["crit"].astype(int)
    df["rep"] = parts["rep"].astype(int)

    # 4. Duplicate run IDs: keep the last line written (the most recent run).
    dup = df.duplicated("run_id", keep="last")
    log(f"  [clean 4] dropped {int(dup.sum())} duplicate lines for the same run_id (kept the most recent).")
    df = df[~dup].copy()

    # 5. The Section 3 worked example is a sanity check, not part of the experiment.
    ex = df["case"] == "example"
    if not include_example:
        df = df[~ex].copy()
        log(f"  [clean 5] set aside {int(ex.sum())} runs of the worked example (use --include-example to keep them).")
    else:
        log(f"  [clean 5] kept {int(ex.sum())} runs of the worked example (--include-example).")

    # 6. Status of each run.
    reason = df["reason"] if "reason" in df else pd.Series(None, index=df.index)
    finished = df["finished"] if "finished" in df else pd.Series(True, index=df.index)
    df["status"] = "finished"
    df.loc[finished.eq(False), "status"] = "capped"                        # main.py hit the pair/time cap
    df.loc[reason.fillna("").str.startswith("killed"), "status"] = "killed"   # batch_run hard timeout
    df.loc[reason.fillna("").str.startswith("skipped"), "status"] = "skipped" # smaller n already failed
    log("  [clean 6] run status: " + ", ".join(f"{k} {v}" for k, v in df["status"].value_counts().items()) + ".")

    # 7. Finished runs whose basis is NOT a Groebner basis: a bug, never a result.
    wrong = (df["status"] == "finished") & df.get("correct", pd.Series(True, index=df.index)).eq(False)
    if wrong.any():
        log(f"  [clean 7] WARNING: dropped {int(wrong.sum())} finished runs marked correct=False (not a Groebner basis). "
            f"Investigate these: {sorted(df.loc[wrong, 'run_id'])[:5]}")
    else:
        log("  [clean 7] all finished runs are verified Groebner bases (correct=True); nothing dropped.")
    df = df[~wrong].copy()

    # 8. Code versions: reported, not dropped. Mixing is fine unless results-affecting code changed in between.
    commits = df["commit"].fillna("none").value_counts() if "commit" in df else pd.Series(dtype=int)
    log("  [clean 8] commits in the data (not filtered): " + ", ".join(f"{c} ({v} runs)" for c, v in commits.items()) + ".")

    # Family / binomial setting columns.
    b = df["case"].str.extract(BINOMIAL)
    df["d"] = pd.to_numeric(b["d"]); df["s"] = pd.to_numeric(b["s"]); df["dist"] = b["dist"]
    is_b = b["n"].notna()
    df.loc[is_b, "family"] = "binomial"
    df.loc[is_b, "n"] = pd.to_numeric(b.loc[is_b, "n"])
    df["n"] = pd.to_numeric(df["n"]).astype("Int64")
    df["setting"] = df["family"]
    df.loc[is_b, "setting"] = (b.loc[is_b, "n"] + "-" + b.loc[is_b, "d"] + "-" + b.loc[is_b, "s"] + " " + b.loc[is_b, "dist"])

    log(f"Kept {len(df)} of {n0} lines.")
    return df


# ----------------------------------------------------------------------------------------------
# 3. Collapse repeats and Random seeds
# ----------------------------------------------------------------------------------------------
def collapse(df, log):
    """One row per case x order x criteria x rule x seed (repeats collapsed).

    Counts are deterministic, so a run counts as finished if ANY of its repeats finished:
    with the pair cap all repeats agree, but the time cap can stop one repeat (e.g. on a busy
    machine) while another finishes. Time = median over the finished repeats only.
    """
    keys = ["family", "setting", "case", "n", "order", "criteria", "rule", "seed"]
    fin = df[df["status"] == "finished"]
    g = df.groupby(keys, dropna=False)
    out = g.agg(reps=("rep", "nunique"),
                status=("status", lambda s: "finished" if (s == "finished").any() else s.iloc[0])).reset_index()
    mixed = g["status"].agg(lambda s: (s == "finished").any() and not (s == "finished").all())
    log(f"  [clean 9] {int(mixed.sum())} runs finished in some repeats but hit the time cap in others; "
        f"counted as finished (counts are identical across repeats), timed over the finished repeats only.")
    times = fin.groupby(keys, dropna=False).agg(time_median=("time_total_elapsed", "median"),
                                                rule_time_median=("time_rule_overhead", "median")).reset_index()
    first_fin = fin.sort_values("rep").drop_duplicates(keys)
    out = out.merge(times, on=keys, how="left")
    out = out.merge(first_fin[keys + [c for c in COUNT_COLS if c in df]], on=keys, how="left")
    out["all_finished"] = out["status"] == "finished"
    return out


def per_rule(runs):
    """One row per case x order x criteria x rule (Random's seeds combined by median)."""
    keys = ["family", "setting", "case", "n", "order", "criteria", "rule"]
    agg = {c: (c, "median") for c in COUNT_COLS if c in runs}
    g = runs.groupby(keys, dropna=False)
    out = g.agg(seeds=("seed", "nunique"),
                finished=("all_finished", "all"),
                status=("status", lambda s: "finished" if (s == "finished").all() else s[s != "finished"].iloc[0]),
                time_median=("time_median", "median"),
                **agg).reset_index()
    out.loc[~out["finished"], COUNT_COLS[:4]] = float("nan")   # no counts for runs that did not finish
    return out


# ----------------------------------------------------------------------------------------------
# 4. Tables
# ----------------------------------------------------------------------------------------------
def fixed_table(pr):
    fx = pr[pr["family"] != "binomial"].copy()
    if fx.empty:
        return fx, pd.DataFrame()
    fx["n"] = fx["n"].astype(int)
    fx["shown"] = [f"{v:,.0f}" if ok else st for v, ok, st in zip(fx["poly_addition"], fx["finished"], fx["status"])]
    shown = fx.pivot_table(index=["family", "n", "order", "criteria"], columns="rule", values="shown",
                           aggfunc="first").fillna("-")
    return fx, shown[[r for r in RULES if r in shown.columns]]


def binomial_tables(pr):
    bi = pr[pr["family"] == "binomial"].copy()
    if bi.empty:
        return bi, pd.DataFrame(), pd.DataFrame()
    keys = ["setting", "order", "criteria", "rule"]
    summary = bi.groupby(keys).agg(
        ideals=("case", "nunique"),
        finished_pct=("finished", lambda s: 100 * s.mean()),
        additions_median=("poly_addition", "median"),
        additions_mean=("poly_addition", "mean"),
        additions_q1=("poly_addition", lambda s: s.quantile(0.25)),
        additions_q3=("poly_addition", lambda s: s.quantile(0.75)),
        pairs_median=("pairs_processed", "median"),
        time_median_s=("time_median", "median"),
    ).reset_index()

    # Fair comparison: medians over only the ideals that EVERY rule finished (same ideals for all rules).
    n_rules = bi["rule"].nunique()
    done = bi[bi["finished"]].groupby(["setting", "order", "criteria", "case"])["rule"].transform("nunique")
    common = bi[bi["finished"]][done == n_rules]
    cm = common.groupby(keys).agg(ideals_all_finished=("case", "nunique"),
                                  additions_median_common=("poly_addition", "median")).reset_index()
    summary = summary.merge(cm, on=keys, how="left")

    # Wins: for each ideal, which rule(s) had the fewest additions. Ties share the win.
    fin = bi[bi["finished"]]
    best = fin.groupby(["setting", "order", "criteria", "case"])["poly_addition"].transform("min")
    winners = fin[fin["poly_addition"] == best]
    share = 1 / winners.groupby(["setting", "order", "criteria", "case"])["rule"].transform("count")
    wins = (winners.assign(win=share).groupby(["setting", "order", "criteria", "rule"])["win"].sum()
            .unstack("rule").fillna(0).reindex(columns=[r for r in RULES if r in fin["rule"].unique()]))
    return bi, summary, wins


def criteria_effect(pr):
    p = pr[pr["finished"]].pivot_table(index=["family", "setting", "case", "order", "rule"], columns="criteria",
                                       values="poly_addition", aggfunc="first").dropna()
    if p.empty or 0 not in p or 1 not in p:
        return pd.DataFrame()
    p["ratio_with_over_without"] = p[1] / p[0]
    return (p.groupby(["family", "setting", "order", "rule"])["ratio_with_over_without"]
            .median().unstack("rule").reindex(columns=RULES))


# ----------------------------------------------------------------------------------------------
# 5. Plots
# ----------------------------------------------------------------------------------------------
def style(ax, title):
    ax.set_title(title, fontsize=10, color=INK, loc="left")
    ax.set_facecolor(SURFACE)
    for side in ["top", "right"]:
        ax.spines[side].set_visible(False)
    for side in ["left", "bottom"]:
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK_2, labelsize=8)
    ax.grid(axis="y", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def plot_fixed(fx, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fams = [f for f in ["cyclic", "katsura"] if f in fx["family"].unique()]
    if not fams:
        return
    combos = [(o, c) for o in ["lex", "grlex"] for c in [0, 1]]
    fig, axes = plt.subplots(len(fams), len(combos), figsize=(3.2 * len(combos), 2.8 * len(fams)),
                             squeeze=False, facecolor=SURFACE)
    for r, fam in enumerate(fams):
        for c, (order, crit) in enumerate(combos):
            ax = axes[r][c]
            d = fx[(fx["family"] == fam) & (fx["order"] == order) & (fx["criteria"] == crit)]
            for rule in RULES:
                s = d[(d["rule"] == rule) & d["finished"]].sort_values("n")
                if not s.empty:
                    ax.plot(s["n"], s["poly_addition"], color=RULE_COLORS[rule], linewidth=2,
                            marker="o", markersize=5, markeredgecolor=SURFACE, markeredgewidth=1.5, label=rule)
            ax.set_yscale("log")
            ns = sorted(d["n"].dropna().unique())
            if ns:
                ax.set_xticks(ns)
            style(ax, f"{fam.capitalize()}, {order}, criteria {'on' if crit else 'off'}")
            if c == 0:
                ax.set_ylabel("additions (log scale)", fontsize=8, color=INK_2)
            if r == len(fams) - 1:
                ax.set_xlabel("n", fontsize=8, color=INK_2)
    handles, labels = axes[0][0].get_legend_handles_labels()
    for ax in axes.flat:                                   # legend from whichever panel has all rules
        h, l = ax.get_legend_handles_labels()
        if len(l) > len(labels):
            handles, labels = h, l
    fig.legend(handles, labels, loc="upper center", ncol=len(labels), frameon=False, fontsize=8)
    fig.text(0.5, 0.005, "Missing points: the run did not finish within the cap.", ha="center", fontsize=7, color=INK_2)
    fig.tight_layout(rect=(0, 0.03, 1, 0.93))
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)


def plot_binomial(summary, path):
    if summary.empty:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    settings = sorted(summary["setting"].unique())
    combos = [(o, c) for o in ["lex", "grlex"] for c in [0, 1]]
    fig, axes = plt.subplots(len(settings), len(combos), figsize=(3.2 * len(combos), 2.6 * len(settings)),
                             squeeze=False, facecolor=SURFACE)
    for r, setting in enumerate(settings):
        for c, (order, crit) in enumerate(combos):
            ax = axes[r][c]
            d = summary[(summary["setting"] == setting) & (summary["order"] == order) & (summary["criteria"] == crit)]
            d = d.set_index("rule").reindex(RULES)          # fixed slot per rule, even if it never finished
            x = list(range(len(RULES)))
            has = d["additions_median"].notna()
            ax.bar([i for i, h in zip(x, has) if h], d.loc[has, "additions_median"], width=0.6,
                   color=[RULE_COLORS[k] for k in d.index[has]], edgecolor=SURFACE, linewidth=2)
            ax.errorbar([i for i, h in zip(x, has) if h], d.loc[has, "additions_median"],
                        yerr=[d.loc[has, "additions_median"] - d.loc[has, "additions_q1"],
                              d.loc[has, "additions_q3"] - d.loc[has, "additions_median"]],
                        fmt="none", ecolor=INK_2, elinewidth=1, capsize=3)
            if not has.any():
                ax.set_yticks([])
                ax.text(0.5, 0.5, "no rule finished", transform=ax.transAxes, ha="center", fontsize=8, color=INK_2)
            ax.set_xticks(x)
            ax.set_xticklabels([f"{k}\n{'DNF' if pd.isna(v) or v == 0 else f'{v:.0f}%'}"
                                for k, v in zip(d.index, d["finished_pct"])], fontsize=7)
            style(ax, f"{setting}, {order}, criteria {'on' if crit else 'off'}")
            if c == 0:
                ax.set_ylabel("median additions", fontsize=8, color=INK_2)
    fig.text(0.5, 0.005, "Bars: median additions over the ideals each rule finished; whiskers: 25th-75th percentile. Under each rule: % of ideals it finished (DNF = none).",
             ha="center", fontsize=7, color=INK_2)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)


# ----------------------------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Clean and summarize the experiment results.")
    parser.add_argument("-i", "--input", type=Path, default=ROOT / "output", help="Folder with *.jsonl (default: output/).")
    parser.add_argument("-O", "--out", type=Path, default=ROOT / "output" / "summary", help="Where to write the summary.")
    parser.add_argument("--include-example", action="store_true", help="Keep the Section 3 worked example.")
    args = parser.parse_args()

    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 20)
    pd.set_option("display.max_rows", 200)
    log = Log()

    print("=" * 80 + "\nDATA CLEANING\n" + "=" * 80)
    df = load(args.input, log)
    if df.empty:
        sys.exit("No results found.")
    df = clean(df, args.include_example, log)

    runs = collapse(df, log)
    pr = per_rule(runs)
    args.out.mkdir(parents=True, exist_ok=True)
    df.drop(columns=["reduced_basis", "pair_order"], errors="ignore").to_csv(args.out / "clean_runs.csv", index=False)

    print("\n" + "=" * 80 + "\nCYCLIC / KATSURA: polynomial additions\n"
          "(capped/killed/skipped = did not finish; - = no record for that run)\n" + "=" * 80)
    fx, fixed_shown = fixed_table(pr)
    if not fx.empty:
        print(fixed_shown.to_string())
        fx.to_csv(args.out / "fixed_cases.csv", index=False)
        plot_fixed(fx, args.out / "fixed_cases.png")

    print("\n" + "=" * 80 + "\nRANDOM BINOMIALS: per setting, order, criteria, rule\n"
          "additions_median etc. use each rule's finished ideals; additions_median_common uses only the\n"
          "ideals that every rule finished, so it is the fair column for comparing rules.\n" + "=" * 80)
    bi, bsum, wins = binomial_tables(pr)
    if not bsum.empty:
        print(bsum.round(1).to_string(index=False))
        print("\nWins (ideals where the rule had the fewest additions; ties split):")
        print(wins.round(1).to_string())
        bsum.to_csv(args.out / "binomial_summary.csv", index=False)
        wins.to_csv(args.out / "binomial_wins.csv")
        plot_binomial(bsum, args.out / "binomial_additions.png")

    print("\n" + "=" * 80 + "\nCRITERIA EFFECT: median of (additions with criteria / without), finished pairs only\n" + "=" * 80)
    ce = criteria_effect(pr)
    if not ce.empty:
        print(ce.round(3).to_string())
        ce.to_csv(args.out / "criteria_effect.csv")

    (args.out / "cleaning_log.txt").write_text("\n".join(log.lines) + "\n")
    print(f"\nWrote tables, plots and cleaning_log.txt to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())