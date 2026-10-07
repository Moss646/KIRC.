# fig_Figure6_main_v1.py
# Main-text Figure 6: 30-gene x 7-immune-feature correlation heatmap + core
# association dot plot.
#   Panel A: heatmap of Spearman rho adjusted for ImmuneScore + StromaScore;
#            bold italic = FDR < 0.05 (BH across 210 pairs).
#   Panel B: ten core associations across four analytic strata
#            (full / adjusted / S1 / S2).
# Reads results/metabolic_immune_corr_30genes_v1.csv.
# Output: results/Figure6_metab_immune_corr_v1.{png,svg}

import sys, os, io, atexit
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

import config


def _hold():
    if os.environ.get('NOPAUSE'):
        return
    try:
        if os.name == 'nt':
            import msvcrt
            msvcrt.getch()
        else:
            input('Press Enter to exit...')
    except (EOFError, KeyboardInterrupt):
        pass


atexit.register(_hold)

if sys.platform == 'win32' and sys.stdout is not None:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                                  errors='replace')

CSV = os.path.join(config.RESULTS_DIR, "metabolic_immune_corr_30genes_v1.csv")
OUT_PNG = os.path.join(config.RESULTS_DIR, "Figure6_metab_immune_corr_v1.png")
OUT_SVG = os.path.join(config.RESULTS_DIR, "Figure6_metab_immune_corr_v1.svg")

plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["axes.linewidth"] = 0.8
plt.rcParams["mathtext.default"] = "regular"


def main():
    df = pd.read_csv(CSV)
    assert len(df) == 210, "expected 210 rows, got %d" % len(df)

    GENES_14 = ["CPT1A", "CPT2", "ACOX1", "ACADM", "ACADSB", "ACAA2",
                "SLC27A2", "CD36", "PLIN2", "PPARG", "PPARGC1A", "FASN",
                "SCD", "ITPKA"]
    GENES_16_ENZ = ["LPCAT1", "LPCAT2", "LPCAT3", "LPCAT4",
                    "MBOAT1", "MBOAT2", "LCLAT1", "MBOAT7",
                    "PLA2G4A", "PLA2G6", "PNPLA8", "PLA2G4C"]
    GENES_16_SPLA2 = ["PLA2G2A", "PLA2G1B", "PLA2G2D", "PLA2G4F"]
    GENES = GENES_14 + GENES_16_ENZ + GENES_16_SPLA2
    VARS = ["TIDE", "PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "CD274"]
    assert df.metabolic_gene.nunique() == 30 and df.immune_variable.nunique() == 7

    rho = df.pivot(index="metabolic_gene", columns="immune_variable",
                   values="rho_adjPurity").loc[GENES, VARS]
    fdr = df.pivot(index="metabolic_gene", columns="immune_variable",
                   values="fdr_adjPurity").loc[GENES, VARS]
    assert rho.shape == (30, 7)

    KEY14 = [
        ("ITPKA", "TIDE"), ("PPARGC1A", "PDCD1"), ("PPARGC1A", "LAG3"),
        ("SLC27A2", "TIDE"), ("CPT1A", "CD274"),
    ]
    KEY16 = [
        ("PNPLA8", "CD274"), ("LCLAT1", "HAVCR2"), ("PLA2G6", "CTLA4"),
        ("PLA2G2D", "TIGIT"), ("PLA2G2A", "TIDE"),
    ]
    KEY = KEY14 + KEY16
    N14 = len(KEY14)
    STRATA = [("full", "rho_all", "fdr_all"),
              ("adj", "rho_adjPurity", "fdr_adjPurity"),
              ("S1", "rho_S1", "p_S1"), ("S2", "rho_S2", "p_S2")]
    sig_col = {"full": "fdr_all", "adj": "fdr_adjPurity",
               "S1": "p_S1", "S2": "p_S2"}

    MARK = {"full": ("o", "#000000"), "adj": ("^", "#E69F00"),
            "S1": ("s", "#56B4E9"), "S2": ("D", "#009E73")}
    STR_ORDER = ["full", "adj", "S1", "S2"]

    cmap = matplotlib.colormaps["RdBu_r"]
    vmax = 0.55
    norm = Normalize(vmin=-vmax, vmax=vmax)

    fig = plt.figure(figsize=(3.35, 7.4))
    gs = fig.add_gridspec(2, 1, height_ratios=[2.90, 1.00], hspace=0.28,
                          left=0.15, right=0.92, top=0.965, bottom=0.115)
    axA = fig.add_subplot(gs[0])
    axB = fig.add_subplot(gs[1])

    # Panel A: heatmap
    axA.imshow(rho.values, cmap=cmap, norm=norm, aspect="auto")

    X = np.arange(len(VARS))
    Y = np.arange(len(GENES))
    for i, g in enumerate(GENES):
        for j, v in enumerate(VARS):
            r = rho.loc[g, v]
            f = fdr.loc[g, v]
            rgb = cmap(norm(r))[:3]
            lum = 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2]
            tcol = "black" if lum > 0.62 else "white"
            sig = f < 0.05
            axA.text(j, i, "%.2f" % r, ha="center", va="center",
                     fontsize=5.5, color=tcol,
                     fontweight="bold" if sig else "normal",
                     fontstyle="italic" if sig else "normal")

    axA.axhline(len(GENES_14) - 0.5, color="black", linewidth=1.0)
    axA.axhline(len(GENES_14) + len(GENES_16_ENZ) - 0.5, color="black",
                linewidth=0.6, linestyle="--")

    axA.set_xticks(X)
    axA.set_xticklabels([v if v == "TIDE" else r"$\it{%s}$" % v for v in VARS],
                        fontsize=6)
    axA.set_yticks(Y)
    axA.set_yticklabels([r"$\it{%s}$" % g for g in GENES], fontsize=6)
    axA.tick_params(length=0)
    axA.set_title("Spearman rho adjusted for immune and stromal content\n"
                  "(bold italic: FDR < 0.05, BH across 210 pairs)",
                  fontsize=6.5, pad=4)

    img = axA.images[0]
    cb = fig.colorbar(img, ax=axA, fraction=0.055, pad=0.02)
    cb.ax.tick_params(labelsize=6, pad=1)

    # Panel B: core associations
    pair_pos = np.arange(len(KEY))
    offsets = {"full": -0.21, "adj": -0.07, "S1": 0.07, "S2": 0.21}
    for j, (g, v) in enumerate(KEY):
        row = df[(df.metabolic_gene == g) & (df.immune_variable == v)].iloc[0]
        for st in STR_ORDER:
            rcol = {"full": "rho_all", "adj": "rho_adjPurity",
                    "S1": "rho_S1", "S2": "rho_S2"}[st]
            sc = sig_col[st]
            r = row[rcol]
            p = row[sc]
            mk, col = MARK[st]
            filled = p < 0.05
            axB.scatter(pair_pos[j] + offsets[st], r, marker=mk, s=20,
                        facecolor=col if filled else "white",
                        edgecolor=col, linewidths=0.9, zorder=3)

    axB.axhline(0, color="0.55", linewidth=0.7, zorder=1)
    axB.axvline(N14 - 0.5, color="0.80", linewidth=0.6,
                linestyle="--", zorder=1)
    axB.set_xticks(pair_pos)
    axB.set_xticklabels(
        [r"$\it{%s}$" % g + "\u2013" + r"$\it{%s}$" % v for g, v in KEY],
        fontsize=5.5, rotation=45, ha="right", rotation_mode="anchor")
    axB.set_ylabel("Spearman rho", fontsize=6.5)
    axB.tick_params(axis="y", labelsize=6)
    axB.tick_params(axis="x", length=0)
    axB.set_ylim(-0.9, 0.9)
    axB.margins(x=0.08)

    handles = []
    for st in STR_ORDER:
        mk, col = MARK[st]
        handles.append(plt.Line2D([], [], marker=mk, linestyle="None",
                                  markersize=5, markerfacecolor=col,
                                  markeredgecolor=col,
                                  label={"full": "Full cohort",
                                         "adj": "Adjusted",
                                         "S1": "Within S1",
                                         "S2": "Within S2"}[st]))
    handles.append(plt.Line2D([], [], marker="o", linestyle="None",
                              markersize=5, markerfacecolor="white",
                              markeredgecolor="0.35",
                              label="Open, P \u2265 0.05"))
    axB.legend(handles=handles, fontsize=5.5, ncol=3, frameon=False,
               loc="upper center", bbox_to_anchor=(0.5, 1.30),
               columnspacing=0.9, handletextpad=0.3)

    for ax_, lab in [(axA, "A"), (axB, "B")]:
        ax_.text(-0.06, 1.05, lab, transform=ax_.transAxes, fontsize=10,
                 fontweight="bold", va="top", ha="right")

    fig.savefig(OUT_PNG, dpi=600)
    fig.savefig(OUT_SVG)
    print("SAVED:", OUT_PNG)
    print("SAVED:", OUT_SVG)
    print("sig cells Panel A: %d/210 (FDR<0.05)" % (fdr.values < 0.05).sum())

    for g, v in KEY:
        row = df[(df.metabolic_gene == g) & (df.immune_variable == v)].iloc[0]
        print("KEY %-9s %-8s rho_all=%.3f | adjImmune=%+.3f | "
              "adjPurity=%+.3f (FDR %.1e) | S1 %+.3f | S2 %+.3f"
              % (g, v, row.rho_all, row.rho_adjImmune, row.rho_adjPurity,
                 row.fdr_adjPurity, row.rho_S1, row.rho_S2))


if __name__ == '__main__':
    main()