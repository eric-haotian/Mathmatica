#!/usr/bin/env python3
"""Figure 1: modulus of the exact response difference (Proposition 3.2) for the
centered binomial pairs with h = 1/8 and r = 1, ..., 4.

The curve is evaluated from the closed form
    |Delta Z(omega)| = (2r)! (omega h / 2)^{2r} / (2^{2r-1} prod_l sqrt(1 + omega^2 t_l^2)),
    t_l = 1 + x_l / 2,  x_l = c + (l - r) h,  l = 0, ..., 2r,
and the unique maximizer omega_* = sqrt(v_*) is located by bisection on the
monotone equation sum_l t_l^2 v / (1 + t_l^2 v) = 2r.  Both quantities are
cross-checked against the archived rational enclosures of the reproducibility
package (code/inherited/fs1/.../results/summary.csv, cases r*_h3) before the
figure is written.  No new experiment is performed; the script only draws the
closed form.

Usage:  python3 make_figure_response_difference.py  [--check path/to/summary.csv]
Output: response_difference.pdf (next to this script).
"""
from __future__ import annotations

import argparse
import csv
import math
from fractions import Fraction
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
DEFAULT_CSV = (
    HERE.parents[2]
    / "reproducibility/Supplementary_Materials/code/inherited/fs1"
    / "EIS_PEM_Math_FS1_Full_Spectrum/results/summary.csv"
)

H = Fraction(1, 8)
C = Fraction(0)
RS = (1, 2, 3, 4)

# Colour-blind safe (Okabe-Ito) hues in fixed order, paired with distinct line
# styles so that the figure also reads in black-and-white print.
STYLE = {
    1: dict(color="#0072B2", ls="-"),
    2: dict(color="#D55E00", ls="--"),
    3: dict(color="#009E73", ls="-."),
    4: dict(color="#CC79A7", ls=":"),
}


def nodes(r: int) -> list[Fraction]:
    return [1 + (C + (l - r) * H) / 2 for l in range(2 * r + 1)]


def modulus(r: int, omega: np.ndarray) -> np.ndarray:
    t = np.array([float(v) for v in nodes(r)])
    num = math.factorial(2 * r) * (omega * float(H) / 2) ** (2 * r)
    den = 2 ** (2 * r - 1) * np.prod(np.sqrt(1 + (omega[:, None] * t[None, :]) ** 2), axis=1)
    return num / den


def uniform_bound(r: int) -> float:
    return math.factorial(2 * r) / 2 ** (2 * r - 1) * float(H) ** (2 * r)


def maximizer(r: int) -> tuple[float, float]:
    """Bisection on the strictly increasing map v -> sum t^2 v/(1+t^2 v) - 2r."""
    t = [float(v) for v in nodes(r)]
    lo, hi = 2 * r / max(t) ** 2, 2 * r / min(t) ** 2  # bracket from (eq:maximizer)

    def g(v: float) -> float:
        return sum(ti * ti * v / (1 + ti * ti * v) for ti in t) - 2 * r

    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if g(mid) < 0:
            lo = mid
        else:
            hi = mid
    v = 0.5 * (lo + hi)
    w = math.sqrt(v)
    return w, float(modulus(r, np.array([w]))[0])


def check_against_archive(path: Path) -> None:
    rows = {row["case"]: row for row in csv.DictReader(path.open())}
    for r in RS:
        row = rows[f"r{r}_h3"]
        assert row["h"] == "1/8" and row["c"] == "0", row
        w, m = maximizer(r)
        lo, hi = float(row["omega_lo"]), float(row["omega_hi"])
        nlo, nhi = float(row["response_norm_lo"]), float(row["response_norm_hi"])
        assert lo - 1e-9 <= w <= hi + 1e-9, (r, w, lo, hi)
        assert nlo * (1 - 1e-8) <= m <= nhi * (1 + 1e-8), (r, m, nlo, nhi)
        print(f"r={r}: omega_* = {w:.10f} in [{lo}, {hi}];  max |Delta Z| = {m:.8e} in [{nlo:.8e}, {nhi:.8e}]")
    print("PASS: figure quantities agree with the archived enclosures")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", type=Path, default=DEFAULT_CSV)
    ap.add_argument("--output", type=Path, default=HERE / "response_difference.pdf")
    args = ap.parse_args()
    if args.check.is_file():
        check_against_archive(args.check)
    else:
        print(f"note: archive summary {args.check} not found; skipping the cross-check")

    plt.rcParams.update(
        {
            "font.family": "STIXGeneral",
            "mathtext.fontset": "stix",
            "axes.formatter.use_mathtext": True,
            "font.size": 9,
            "axes.labelsize": 9,
            "legend.fontsize": 8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "axes.linewidth": 0.6,
            "lines.linewidth": 1.1,
            "pdf.fonttype": 42,
        }
    )
    omega = np.logspace(-2, 3, 1200)
    fig, ax = plt.subplots(figsize=(4.9, 3.25))
    for r in RS:
        y = modulus(r, omega)
        ax.plot(omega, y, label=rf"$r={r}$", **STYLE[r])
        ax.axhline(uniform_bound(r), color=STYLE[r]["color"], lw=0.6, ls=(0, (1, 2)), alpha=0.9)
        w, m = maximizer(r)
        ax.plot([w], [m], marker="o", ms=4, color=STYLE[r]["color"], mec="white", mew=0.6, zorder=5)
        ax.annotate(rf"$r={r}$", xy=(1e3, float(y[-1])), xytext=(3, 0), textcoords="offset points", fontsize=7, color="black", va="center", ha="left", annotation_clip=False)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1e-2, 1e3)
    ax.set_ylim(1e-14, 1)
    ax.set_xlabel(r"$\omega$")
    ax.set_ylabel(r"$|\Delta Z(\omega)|$")
    ax.grid(True, which="major", color="#d0d0d0", lw=0.4)
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=4, handlelength=2.6, columnspacing=1.6)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    fig.tight_layout(pad=0.4)
    fig.savefig(args.output)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
