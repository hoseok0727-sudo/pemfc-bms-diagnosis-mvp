"""Recreate the full-pulse scalar sufficient-bound figure.

Each point is a separate interval-evaluated sufficient bound; connecting lines
are visualization, not a certificate over intervals of parameter values.
"""
from pathlib import Path
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from repair_budget import full_certificate


def main(out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    lambdas = np.linspace(0.001, 1 / 120, 45)
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    for B in (0.1, 0.5, 1.0):
        margins = [full_certificate(str(lam), str(lam / 100000), str(B))[
            'fraction_min_equilibrium_lower'] for lam in lambdas]
        ax.plot(lambdas, margins, label=f'Entry amplitude B = {B}')
    ax.axhline(0, linestyle='--', linewidth=1)
    ax.set_xlabel('Dimensionless pulse parameter λ')
    ax.set_ylabel('Lower bound / minimum equilibrium')
    ax.set_title('Conditional scalar positivity over the full pulse')
    ax.legend()
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path,
        default=Path(__file__).resolve().parents[1] / 'results' / 'full_pulse_margin.png')
    main(parser.parse_args().out)
