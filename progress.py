import shutil
import sys
import time


def _preview(items, total, width=None):
    """Format the first few items as [a, b, c, +n more]."""
    shown = [str(x) for x in items]
    if width:
        shown = [s if len(s) <= width else s[:width - 3] + "..." for s in shown]
    extra = total - len(items)
    return "[" + ", ".join(shown) + (f", +{extra} more" if extra > 0 else "") + "]"


class Progress:
    """Live status block that redraws in place (no scrolling)."""

    def __init__(self, enabled=True, k=5, poly_width=30, interval=0.1):
        self.enabled = enabled
        self.k = k                    # how many items of B and G to show
        self.poly_width = poly_width  # max characters per polynomial
        self.interval = interval      # seconds between redraws
        self._last = 0.0
        self._lines = 0               # lines drawn last time

    def update(self, B, G, pair, stats, force=False):
        if not self.enabled:
            return
        now = time.perf_counter()
        if not force and now - self._last < self.interval:
            return
        self._last = now

        B_show = [(i + 1, j + 1) for i, j in B[:self.k]]        # 1-based, like the paper
        rows = [
            f"pair {pair[0] + 1, pair[1] + 1}   processed {stats['pairs_processed']}"
            f"   zero {stats['zero_reductions']}   adds {stats.get('polynomial_additions', 0)}",
            f"B ({len(B)}): {_preview(B_show, len(B))}",
            f"G ({len(G)}): {_preview(G[:self.k], len(G), self.poly_width)}",
        ]
        cols = shutil.get_terminal_size().columns
        rows = [r[:cols - 1] for r in rows]   # never wrap, or the redraw drifts

        out = f"\033[{self._lines}F" if self._lines else ""   # cursor up to first line
        out += "".join("\033[2K" + r + "\n" for r in rows)     # clear line, write row
        sys.stderr.write(out)
        sys.stderr.flush()
        self._lines = len(rows)

    def close(self, B, G, pair, stats):
        """Draw the final state once, regardless of the interval."""
        if pair is not None:
            self.update(B, G, pair, stats, force=True)