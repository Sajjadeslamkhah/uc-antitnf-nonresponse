"""Figure for the alternative-discovery analysis (reads the output of reverse_discovery.py).

Usage: python revision/plot_reverse_discovery.py --out OUT_DIR
"""
import argparse
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

LABEL = {'GOL_GSE92415': 'Golimumab GSE92415', 'PBO_GSE92415': 'Placebo GSE92415',
         'IFX_GSE16879': 'Infliximab GSE16879', 'IFX_GSE23597': 'Infliximab GSE23597'}
DISC = {'GOL_GSE92415 (primary)': 'Discovery: golimumab GSE92415 (primary)',
        'IFX_GSE16879': 'Discovery: infliximab GSE16879',
        'IFX_GSE23597': 'Discovery: infliximab GSE23597',
        'IFX_pooled': 'Discovery: infliximab cohorts pooled'}
COLOR = {'test': '#1F4E79', 'negative control (placebo)': '#B03A2E', 'discovery (optimistic)': '#9E9E9E'}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    out = ap.parse_args().out
    r = pd.read_csv(os.path.join(out, 'T_reverse_discovery_AUC.csv'))
    s = pd.read_csv(os.path.join(out, 'T_reverse_discovery_signatures.csv')).set_index('discovery')

    fig, axes = plt.subplots(1, 4, figsize=(13, 3.4), sharex=True)
    for ax, (d, title) in zip(axes, DISC.items()):
        sub = r[r.discovery == d].reset_index(drop=True)
        for i, row in sub.iterrows():
            c = COLOR[row.role]
            ax.plot([row.CI95_low, row.CI95_high], [i, i], color=c, lw=1.6)
            ax.plot(row.AUC_NR, i, 'o', color=c, ms=6, mfc='white' if row.role.startswith('discovery') else c)
            ax.text(1.05, i, f'{row.AUC_NR:.2f}', va='center', fontsize=7)
        ax.axvline(0.5, color='k', lw=0.7, ls=':')
        ax.set_yticks(range(len(sub)))
        ax.set_yticklabels([LABEL[k] for k in sub.scored_cohort] if ax is axes[0] else [], fontsize=7.5)
        ax.invert_yaxis()
        ax.set_xlim(0.2, 1.03)
        n = s.loc[d]
        ax.set_title(f'{title}\n{int(n.signature_genes)} genes (n = {int(n.n_NR)} NR / {int(n.n_R)} R)', fontsize=8)
        ax.set_xlabel('AUC for non-response (95% CI)', fontsize=7.5)
        ax.tick_params(axis='x', labelsize=7)
        ax.spines[['top', 'right']].set_visible(False)
    handles = [plt.Line2D([], [], color=c, marker='o', lw=1.6, label=l,
                          mfc='white' if l.startswith('discovery') else c) for l, c in COLOR.items()]
    fig.legend(handles=handles, loc='lower center', ncol=3, fontsize=7.5, frameon=False, bbox_to_anchor=(0.5, -0.04))
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(out, f'FigureS13_alternative_discovery.{ext}'), dpi=600, bbox_inches='tight')


if __name__ == '__main__':
    main()
