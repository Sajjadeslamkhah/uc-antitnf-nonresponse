"""
Revision analysis (Reviewer 2, comment 7): alternative discovery cohorts.

The primary analysis derives the baseline non-response signature in the Week-0
golimumab arm of GSE92415 and tests it in two infliximab cohorts. This script
re-runs the same pipeline with each infliximab cohort, and with the two
infliximab cohorts pooled, as the discovery set, and scores every resulting
signature in all remaining cohorts, including the GSE92415 placebo arm.

All thresholds and scoring rules are those of the primary analysis
(Supplementary Methods 2.6-2.8) and are not re-tuned:
  * low-expression filter: genes below the 25th percentile of mean expression
    within a contrast are removed before limma
  * limma lmFit + eBayes(trend = TRUE, robust = TRUE)
  * signature: nominal P < 0.01 and |log2 fold change| >= 0.5
  * score: mean within-cohort z (up in NR) - mean within-cohort z (down in NR)
  * AUC with 2,000 bootstrap resamples, Mann-Whitney U test
  * inflammatory burden: mean z of up-regulated core UC genes (Table S4)
    after removal of the signature genes; contribution of the signature
    beyond burden tested by likelihood-ratio test

Usage:
    python revision/reverse_discovery.py --data DATA_DIR --out OUT_DIR

DATA_DIR must contain
    expression/GSE92415_gene_level_all_samples.csv.gz
    expression/GSE16879_gene_level_all_samples.csv.gz
    expression/GSE23597_gene_level_all_samples.csv.gz
    samplesheets/samplesheet_<contrast>.csv  (written by Cell 1.6 of the notebook)
    Supplementary_Tables_S1-S12.xlsx         (sheet S4_core_UC_genes_687)
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import warnings
import datetime

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import chi2, mannwhitneyu
from sklearn.metrics import roc_auc_score

L2_P, L2_LFC = 0.01, 0.5
REPL_P = 0.05
LOW_EXPR_QUANTILE = 0.25
N_BOOT = 2000
SEED = 2026

COHORTS = {
    'GOL_GSE92415': ('GSE92415', 'samplesheet_L2_disc_GOL_GSE92415.csv'),
    'PBO_GSE92415': ('GSE92415', 'samplesheet_L2_placebo_GSE92415.csv'),
    'IFX_GSE16879': ('GSE16879', 'samplesheet_L2_val_IFX_GSE16879.csv'),
    'IFX_GSE23597': ('GSE23597', 'samplesheet_L2_val_IFX_GSE23597.csv'),
}
# discovery set -> cohorts used for testing (the discovery cohort itself is
# also scored, but reported as optimistic)
DISCOVERY = {
    'GOL_GSE92415 (primary)': ['GOL_GSE92415'],
    'IFX_GSE16879': ['IFX_GSE16879'],
    'IFX_GSE23597': ['IFX_GSE23597'],
    'IFX_pooled': ['IFX_GSE16879', 'IFX_GSE23597'],
}

LIMMA_R = r'''
suppressMessages(library(limma))
a    <- commandArgs(TRUE)
e    <- as.matrix(read.csv(a[1], row.names = 1, check.names = FALSE))
p    <- read.csv(a[2], row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
e    <- e[, rownames(p), drop = FALSE]
p$group <- factor(p$group, levels = c(a[4], a[5]))
covs <- if (length(a) >= 6 && nchar(a[6]) > 0) strsplit(a[6], ",")[[1]] else character(0)
d    <- model.matrix(as.formula(paste("~", paste(c("group", covs), collapse = " + "))), data = p)
fit  <- eBayes(lmFit(e, d), trend = TRUE, robust = TRUE)
tt   <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
write.csv(tt, a[3])
cat(as.character(packageVersion("limma")))
'''


def load(data):
    expr, sheets = {}, {}
    for k, (gse, ss) in COHORTS.items():
        if gse not in expr:
            expr[gse] = pd.read_csv(os.path.join(data, 'expression', f'{gse}_gene_level_all_samples.csv.gz'),
                                    index_col=0)
        sheets[k] = pd.read_csv(os.path.join(data, 'samplesheets', ss), index_col=0)
    core = pd.read_excel(os.path.join(data, 'Supplementary_Tables_S1-S12.xlsx'),
                         sheet_name='S4_core_UC_genes_687', index_col=0)
    return expr, sheets, core


def run_limma(g, pheno, covariates, tmp, tag):
    ef, pf, of, rf = (os.path.join(tmp, f'{tag}_{x}') for x in ('expr.csv', 'pheno.csv', 'out.csv', 'limma.R'))
    open(rf, 'w').write(LIMMA_R)
    g.to_csv(ef)
    pheno[['group', *covariates]].to_csv(pf)
    r = subprocess.run(['Rscript', rf, ef, pf, of, 'R', 'NR', ','.join(covariates)],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(r.stderr)
    return pd.read_csv(of, index_col=0), r.stdout.strip()


def low_expr_filter(g):
    m = g.mean(axis=1)
    return g[m > m.quantile(LOW_EXPR_QUANTILE)]


def discovery_table(members, expr, sheets, tmp, tag):
    if len(members) == 1:
        k = members[0]
        s = sheets[k]
        g = low_expr_filter(expr[COHORTS[k][0]][s.index].dropna())
        return run_limma(g, s, (), tmp, tag)
    # pooled: genes present in every member cohort; filter applied within each
    # cohort, cohort entered as a blocking covariate
    parts, ph = [], []
    for k in members:
        s = sheets[k]
        parts.append(low_expr_filter(expr[COHORTS[k][0]][s.index].dropna()))
        ph.append(s[['group']].assign(cohort=k))
    genes = sorted(set.intersection(*[set(p.index) for p in parts]))
    g = pd.concat([p.loc[genes] for p in parts], axis=1)
    return run_limma(g, pd.concat(ph), ('cohort',), tmp, tag)


def zscore(g):
    return g.sub(g.mean(axis=1), axis=0).div(g.std(axis=1).replace(0, np.nan), axis=0)


def score(z, up, dn):
    u = [x for x in up if x in z.index]
    d = [x for x in dn if x in z.index]
    sc = (z.loc[u].mean() if u else 0) - (z.loc[d].mean() if d else 0)
    return sc, len(u) + len(d)


def auc_ci(y, x, rng):
    y, x = np.asarray(y), np.asarray(x)
    a = roc_auc_score(y, x)
    boots = []
    for _ in range(N_BOOT):
        b = rng.choice(len(y), len(y), replace=True)
        if len(np.unique(y[b])) == 2:
            boots.append(roc_auc_score(y[b], x[b]))
    return a, np.percentile(boots, 2.5), np.percentile(boots, 97.5)


def fit_logit(y, X):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return sm.Logit(y, sm.add_constant(X, has_constant='add')).fit(disp=0, maxiter=200)


def lrt(y, Xb, Xf):
    m0, m1 = fit_logit(y, Xb), fit_logit(y, Xf)
    stat = 2 * (m1.llf - m0.llf)
    return m1, stat, chi2.sf(stat, Xf.shape[1] - Xb.shape[1])


def or_text(m, term):
    b, ci = m.params[term], m.conf_int().loc[term]
    if abs(b) > 8 or not np.isfinite(ci).all() or (ci.iloc[1] - ci.iloc[0]) > 16:
        return 'unstable (separation)'
    return f'{np.exp(b):.2f} [{np.exp(ci.iloc[0]):.2f}-{np.exp(ci.iloc[1]):.2f}]'


def std(s):
    return (s - s.mean()) / s.std()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    assert shutil.which('Rscript'), 'Rscript not found'

    expr, sheets, core = load(args.data)
    core_up = list(core.index[core['logFC_disc'] > 0])
    z = {k: zscore(expr[COHORTS[k][0]][s.index]) for k, s in sheets.items()}
    y = {k: (s.group == 'NR').astype(int) for k, s in sheets.items()}

    tmp = tempfile.mkdtemp()
    tables, sigs, limma_version = {}, {}, None
    for name, members in DISCOVERY.items():
        tt, limma_version = discovery_table(members, expr, sheets, tmp, name.split()[0])
        tables[name] = tt
        sig = tt[(tt['P.Value'] < L2_P) & (tt['logFC'].abs() >= L2_LFC)]
        sigs[name] = sig
        tt.to_csv(os.path.join(args.out, f'DEG_discovery_{name.split()[0]}.csv'))
        sig.to_csv(os.path.join(args.out, f'signature_{name.split()[0]}.csv'))

    primary = sigs['GOL_GSE92415 (primary)']
    single_tables = {k: tables[n] for n, k in [('GOL_GSE92415 (primary)', 'GOL_GSE92415'),
                                               ('IFX_GSE16879', 'IFX_GSE16879'),
                                               ('IFX_GSE23597', 'IFX_GSE23597')]}
    # placebo DE table for gene-level replication in the placebo arm
    s = sheets['PBO_GSE92415']
    single_tables['PBO_GSE92415'], _ = run_limma(low_expr_filter(expr['GSE92415'][s.index].dropna()), s, (),
                                                 tmp, 'PBO')

    rng = np.random.default_rng(SEED)
    rows, sig_rows = [], []
    for name, members in DISCOVERY.items():
        sig = sigs[name]
        up, dn = list(sig.index[sig.logFC > 0]), list(sig.index[sig.logFC < 0])
        infl = [g for g in core_up if g not in set(sig.index)]
        sig_rows.append({
            'discovery': name, 'n_R': int(sum((sheets[k].group == 'R').sum() for k in members)),
            'n_NR': int(sum((sheets[k].group == 'NR').sum() for k in members)),
            'genes_tested': len(tables[name]), 'signature_genes': len(sig), 'up_in_NR': len(up), 'down_in_NR': len(dn),
            'n_FDR<0.05': int((tables[name]['adj.P.Val'] < 0.05).sum()),
            'overlap_with_primary_signature': len(set(sig.index) & set(primary.index)),
            'overlap_with_core_UC_genes': len(set(sig.index) & set(core.index)),
            'overlap_genes_with_primary': ';'.join(sorted(set(sig.index) & set(primary.index))),
        })
        for k in COHORTS:
            role = 'discovery (optimistic)' if k in members else (
                'negative control (placebo)' if k == 'PBO_GSE92415' else 'test')
            sc, n_used = score(z[k], up, dn)
            yy = y[k]
            a, lo, hi = auc_ci(yy.values, sc.values, rng)
            p = mannwhitneyu(sc[yy == 1], sc[yy == 0]).pvalue
            # gene-level replication in this cohort (same direction, nominal P < 0.05)
            t = single_tables[k].reindex(sig.index)
            repl = int(((np.sign(t.logFC) == np.sign(sig.logFC)) & (t['P.Value'] < REPL_P)).sum())
            # beyond inflammatory burden
            inf_sc, n_inf = score(z[k], infl, [])
            X = pd.DataFrame({'burden': std(inf_sc), 'signature': std(sc)})
            try:
                m1, stat, p_lrt = lrt(yy, X[['burden']], X)
                or_ = or_text(m1, 'signature')
            except Exception as e:  # noqa: BLE001
                stat, p_lrt, or_ = np.nan, np.nan, f'failed: {e}'
            rows.append({
                'discovery': name, 'scored_cohort': k, 'role': role,
                'n_R': int((yy == 0).sum()), 'n_NR': int(yy.sum()), 'genes_used': n_used,
                'AUC_NR': round(a, 3), 'CI95_low': round(lo, 3), 'CI95_high': round(hi, 3),
                'MannWhitney_P': p, 'genes_replicated_P<0.05_same_direction': repl,
                'AUC_burden_alone': round(roc_auc_score(yy, inf_sc), 3), 'burden_genes_used': n_inf,
                'OR_signature_per_SD_adj_burden': or_, 'LRT_stat_beyond_burden': stat,
                'LRT_P_beyond_burden': p_lrt,
            })

    res = pd.DataFrame(rows)
    sig_tab = pd.DataFrame(sig_rows)
    res.to_csv(os.path.join(args.out, 'T_reverse_discovery_AUC.csv'), index=False)
    sig_tab.to_csv(os.path.join(args.out, 'T_reverse_discovery_signatures.csv'), index=False)
    log = {
        'date': datetime.datetime.now().isoformat(timespec='seconds'),
        'script': 'revision/reverse_discovery.py', 'seed': SEED, 'n_boot': N_BOOT,
        'thresholds': {'L2_P': L2_P, 'L2_LFC': L2_LFC, 'REPL_P': REPL_P, 'low_expr_quantile': LOW_EXPR_QUANTILE},
        'discovery_sets': DISCOVERY, 'core_up_genes': len(core_up),
        'limma': limma_version,
        'R': subprocess.run(['Rscript', '-e', 'cat(R.version.string)'], capture_output=True, text=True).stdout,
        'python': sys.version.split()[0],
        'sample_counts': {k: s.group.value_counts().to_dict() for k, s in sheets.items()},
    }
    json.dump(log, open(os.path.join(args.out, 'run_log_reverse_discovery.json'), 'w'), indent=2, default=str)
    pd.set_option('display.width', 250)
    pd.set_option('display.max_columns', 30)
    print(sig_tab.drop(columns='overlap_genes_with_primary').to_string(index=False))
    print(res.drop(columns=['burden_genes_used']).to_string(index=False))


if __name__ == '__main__':
    main()
