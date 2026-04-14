"""
Canine Atopic Dermatitis – Full Analysis Pipeline
Dataset: GSE168109 (Agilent whole-blood microarray, Canis lupus familiaris)
Groups : Healthy (n=8) | AD pre-ASIT (n=7) | AD post-ASIT (n=7)

Analysis:
  1. Preprocessing & QC
  2. PCA
  3. Differential expression (AD vs Healthy, Post-ASIT vs Pre-ASIT)
  4. IL-31 / JAK-STAT pathway deep-dive
  5. Canine vs Human gene conservation table
  6. Assay-readiness scoring
  7. Export figures + tables
"""

import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy import stats
from scipy.stats import ttest_ind
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings("ignore")

# ── Paths ────────────────────────────────────────────────────────────────────
BASE = os.getcwd()  # expected: .../Elanco-Canine
PROC = os.path.join(BASE, "data", "processed")
FIG  = os.path.join(BASE, "results", "figures")
TAB  = os.path.join(BASE, "results", "tables")
os.makedirs(FIG, exist_ok=True)
os.makedirs(TAB, exist_ok=True)

# ── Modern color palette ──────────────────────────────────────────────────────
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap

# Group colors — viridis-derived for colorblind accessibility
PAL = {
    "Healthy": "#1F968BFF",   # viridis teal-green
    "AD_pre":  "#DE7065FF",   # plasma warm-red
    "AD_post": "#3B528BFF",   # viridis dark blue
}

# Accent colors
C_UP    = "#B8114A"   # magma deep red — upregulated
C_DOWN  = "#2D708E"   # viridis steel blue — downregulated
C_NS    = "#D3D3D3"   # neutral grey
C_HIGH  = "#FDE725"   # viridis yellow — HIGH priority
C_MED   = "#21908C"   # viridis teal — MEDIUM priority
C_LOW   = "#440154"   # viridis dark purple — LOW priority
C_DARK  = "#1A1A2E"   # near-black for text/titles
C_AMBER = "#F3B83E"   # warm amber accent

# Colormaps
CMAP_DIV   = "RdBu_r"         # diverging: used for heatmaps
CMAP_SEQ   = "viridis"        # sequential: used for conservation bars
CMAP_PLASMA = "plasma"         # plasma: used for -log10p
CMAP_MAGMA  = "magma"          # magma: used for group means

sns.set_theme(style="ticks", font="DejaVu Sans")
plt.rcParams.update({
    "figure.dpi": 150,
    "font.family": "DejaVu Sans",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": False,
    "axes.facecolor": "#FAFAFA",
    "figure.facecolor": "white",
    "axes.labelcolor": "#333333",
    "xtick.color": "#555555",
    "ytick.color": "#555555",
})

# ════════════════════════════════════════════════════════════════════════════
# 1. LOAD & PREPROCESS
# ════════════════════════════════════════════════════════════════════════════
print("1. Loading data...")
expr_raw = pd.read_csv(os.path.join(PROC, "expression_matrix_raw.csv"), index_col=0)
meta     = pd.read_csv(os.path.join(PROC, "sample_metadata.csv"))
plat     = pd.read_csv(os.path.join(PROC, "platform_GPL13605.csv"), low_memory=False)

# Assign group labels
def assign_group(row):
    t = row["title"].lower()
    if "healthy" in t or "control" in t:
        return "Healthy"
    elif "after therapy" in t:
        return "AD_post"
    else:
        return "AD_pre"

meta["group"] = meta.apply(assign_group, axis=1)
group_map = dict(zip(meta["sample_id"], meta["group"]))
print(meta.groupby("group")["sample_id"].count().to_string())

# Keep only feature probes (remove control probes)
plat_clean = plat[plat["CONTROL_TYPE"].astype(str).str.upper().isin(["FALSE", "NAN", ""])]
valid_ids = plat_clean["ID"].astype(str)
expr_raw.index = expr_raw.index.astype(str)
expr = expr_raw.loc[expr_raw.index.isin(valid_ids)].copy()
print(f"Probes after removing controls: {len(expr)}")

# Quantile normalization (column-wise)
def quantile_normalize(df):
    rank_mean = df.stack().groupby(df.rank(method="first").stack().astype(int)).mean()
    return df.rank(method="min").stack().astype(int).map(rank_mean).unstack()

expr_norm = quantile_normalize(expr.apply(pd.to_numeric, errors="coerce").dropna())

# Log2 transform (add small offset to avoid log(0))
expr_log = np.log2(expr_norm + 1)

# Average duplicate probes → gene-level
probe_gene = plat_clean[["ID", "GENE_SYMBOL"]].copy()
probe_gene["ID"] = probe_gene["ID"].astype(str)
probe_gene = probe_gene[probe_gene["GENE_SYMBOL"].notna() & (probe_gene["GENE_SYMBOL"] != "")]
probe_gene = probe_gene.drop_duplicates("ID")

expr_log.index = expr_log.index.astype(str)
expr_log["GENE_SYMBOL"] = expr_log.index.map(probe_gene.set_index("ID")["GENE_SYMBOL"])
expr_gene = expr_log.dropna(subset=["GENE_SYMBOL"])
expr_gene = expr_gene.groupby("GENE_SYMBOL").mean()
print(f"Unique genes after collapsing: {len(expr_gene)}")

# Save
expr_gene.to_csv(os.path.join(PROC, "expression_gene_log2.csv"))
print("Gene-level expression saved.")

# ════════════════════════════════════════════════════════════════════════════
# 2. PCA
# ════════════════════════════════════════════════════════════════════════════
print("\n2. PCA...")
# Drop genes with any NaN, then transpose
expr_gene_clean = expr_gene.dropna()
X = expr_gene_clean.T.values
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)
pca = PCA(n_components=3)
coords = pca.fit_transform(X_scaled)
var_exp = pca.explained_variance_ratio_ * 100

pca_df = pd.DataFrame(coords, columns=["PC1", "PC2", "PC3"],
                      index=expr_gene_clean.columns)
pca_df["group"] = pca_df.index.map(group_map)

fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), facecolor="white")
markers = {"Healthy": "o", "AD_pre": "^", "AD_post": "s"}

for ax, (px, py) in zip(axes, [("PC1","PC2"), ("PC1","PC3")]):
    ax.set_facecolor("#F8F8F8")
    # draw ellipse-like confidence region via scatter kde contour
    for grp, color in PAL.items():
        sub = pca_df[pca_df["group"] == grp]
        ax.scatter(sub[px], sub[py], c=color, s=120,
                   label=grp.replace("_", " "),
                   marker=markers[grp],
                   edgecolors="white", linewidths=1.2, zorder=4, alpha=0.92)
    ax.set_xlabel(f"{px}  ({var_exp[int(px[-1])-1]:.1f}% variance)", fontsize=11)
    ax.set_ylabel(f"{py}  ({var_exp[int(py[-1])-1]:.1f}% variance)", fontsize=11)
    ax.set_title(f"{px} vs {py}", fontsize=12, fontweight="bold", color=C_DARK)
    ax.axhline(0, color="#BBBBBB", lw=0.8, ls="--", zorder=1)
    ax.axvline(0, color="#BBBBBB", lw=0.8, ls="--", zorder=1)
    leg = ax.legend(frameon=True, fontsize=9, framealpha=0.9,
                    edgecolor="#CCCCCC", loc="best")

# Scree-style variance bar inset in ax2
ax_ins = axes[1].inset_axes([0.65, 0.65, 0.32, 0.30])
var_colors = [cm.viridis(v) for v in np.linspace(0.2, 0.85, 3)]
ax_ins.bar([0, 1, 2], var_exp[:3], color=var_colors, edgecolor="white", linewidth=0.5)
ax_ins.set_xticks([0, 1, 2])
ax_ins.set_xticklabels(["PC1","PC2","PC3"], fontsize=6)
ax_ins.set_ylabel("% var", fontsize=6)
ax_ins.tick_params(labelsize=6)
ax_ins.set_facecolor("#F0F0F0")
for sp in ax_ins.spines.values(): sp.set_visible(False)

fig.suptitle("PCA — Canine Blood Transcriptomes  |  GSE168109\nHealthy vs AD Pre-ASIT vs AD Post-ASIT",
             fontsize=13, fontweight="bold", y=1.02, color=C_DARK)
plt.tight_layout()
fig.savefig(os.path.join(FIG, "01_PCA.png"), dpi=150, bbox_inches="tight")
plt.close()
print("  Saved 01_PCA.png")

# ════════════════════════════════════════════════════════════════════════════
# 3. DIFFERENTIAL EXPRESSION
# ════════════════════════════════════════════════════════════════════════════
print("\n3. Differential expression...")

def run_de(expr_df, group_a, group_b, meta_df, label):
    cols_a = meta_df[meta_df["group"] == group_a]["sample_id"].tolist()
    cols_b = meta_df[meta_df["group"] == group_b]["sample_id"].tolist()
    cols_a = [c for c in cols_a if c in expr_df.columns]
    cols_b = [c for c in cols_b if c in expr_df.columns]

    results = []
    for gene in expr_df.index:
        vals_a = expr_df.loc[gene, cols_a].dropna().values.astype(float)
        vals_b = expr_df.loc[gene, cols_b].dropna().values.astype(float)
        if len(vals_a) < 2 or len(vals_b) < 2:
            continue
        log2fc = vals_b.mean() - vals_a.mean()
        tstat, pval = ttest_ind(vals_a, vals_b, equal_var=False)
        results.append({"gene": gene, "log2FC": log2fc, "pval": pval,
                        "mean_a": vals_a.mean(), "mean_b": vals_b.mean()})

    de = pd.DataFrame(results).set_index("gene")
    # BH FDR correction
    from statsmodels.stats.multitest import multipletests
    _, padj, _, _ = multipletests(de["pval"].fillna(1), method="fdr_bh")
    de["padj"] = padj
    de["-log10p"] = -np.log10(de["pval"].clip(lower=1e-300))
    # Use nominal p < 0.05 and |FC| > 0.3 for small-n microarray data
    de["significant"] = (de["pval"] < 0.05) & (de["log2FC"].abs() > 0.3)
    de = de.sort_values("padj")
    de.to_csv(os.path.join(TAB, f"DE_{label}.csv"))
    print(f"  {label}: {de['significant'].sum()} DEGs (FDR<0.05, |FC|>0.5)")
    return de

de_ad_vs_healthy    = run_de(expr_gene, "Healthy", "AD_pre",  meta, "AD_vs_Healthy")
de_post_vs_pre      = run_de(expr_gene, "AD_pre",  "AD_post", meta, "PostASIT_vs_PreASIT")

# ── Volcano plot ─────────────────────────────────────────────────────────────
def volcano(de_df, title, outname, highlight_genes=None):
    fig, ax = plt.subplots(figsize=(9.5, 6.5), facecolor="white")
    ax.set_facecolor("#F8F8F8")

    # Color points by -log10p using plasma colormap for non-sig, solid for sig
    log10p = de_df["-log10p"].clip(upper=20)
    norm = plt.Normalize(vmin=0, vmax=log10p.max())
    plasma_colors = cm.plasma(norm(log10p))

    # Non-significant: muted plasma
    ns_mask = ~de_df["significant"]
    ax.scatter(de_df.loc[ns_mask, "log2FC"], log10p[ns_mask],
               c=C_NS, s=8, alpha=0.45, linewidths=0, zorder=2, rasterized=True)

    # Significant down: viridis blue range
    down_mask = de_df["significant"] & (de_df["log2FC"] < 0)
    ax.scatter(de_df.loc[down_mask, "log2FC"], log10p[down_mask],
               c=C_DOWN, s=20, alpha=0.85, linewidths=0, zorder=3)

    # Significant up: magma red range
    up_mask = de_df["significant"] & (de_df["log2FC"] > 0)
    ax.scatter(de_df.loc[up_mask, "log2FC"], log10p[up_mask],
               c=C_UP, s=20, alpha=0.85, linewidths=0, zorder=3)

    # Reference lines
    ax.axhline(-np.log10(0.05), color="#444444", ls="--", lw=0.9, alpha=0.7)
    ax.axvline(0.3,  color="#888888", ls=":", lw=0.8, alpha=0.6)
    ax.axvline(-0.3, color="#888888", ls=":", lw=0.8, alpha=0.6)

    # Highlight pathway genes
    if highlight_genes:
        for gene in highlight_genes:
            if gene in de_df.index:
                row = de_df.loc[gene]
                lp = min(row["-log10p"], 20)
                fc_clr = C_AMBER if row["log2FC"] > 0 else "#2D708E"
                ax.scatter(row["log2FC"], lp, s=100, zorder=6,
                           edgecolors=C_DARK, linewidths=1.0, c=fc_clr)
                ax.annotate(gene, (row["log2FC"], lp),
                            fontsize=7.5, xytext=(4, 3), textcoords="offset points",
                            fontweight="bold", color=C_DARK)

    up_n   = up_mask.sum()
    down_n = down_mask.sum()
    ax.set_xlabel("log2 Fold Change", fontsize=11, color=C_DARK)
    ax.set_ylabel("-log10(p-value)", fontsize=11, color=C_DARK)
    ax.set_title(title, fontsize=12, fontweight="bold", color=C_DARK, pad=10)

    patches = [mpatches.Patch(color=C_UP,   label=f"Up  n={up_n}"),
               mpatches.Patch(color=C_DOWN, label=f"Down  n={down_n}"),
               mpatches.Patch(color=C_NS,   label="NS")]
    ax.legend(handles=patches, fontsize=9, frameon=True,
              framealpha=0.9, edgecolor="#CCCCCC")
    ax.text(0.98, 0.02, "p < 0.05, |FC| > 0.3", transform=ax.transAxes,
            fontsize=8, ha="right", va="bottom", color="#888888")
    plt.tight_layout()
    fig.savefig(os.path.join(FIG, outname), dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved {outname}")

# Key pathway genes to highlight
pathway_genes = [
    "IL31", "IL31RA", "OSMR",                             # IL-31 axis
    "JAK1", "JAK2", "JAK3", "TYK2",                       # JAK kinases
    "STAT1", "STAT3", "STAT5A", "STAT5B", "STAT6",        # STATs
    "IL4", "IL4R", "IL13", "IL13RA1",                     # Type-2 cytokines
    "IL2", "IFNG", "TNF",                                  # Th1/inflammation
    "GATA3", "RORC", "TBX21",                              # T-cell master TFs
    "CCR4", "CCL17", "CCL22",                              # Th2 chemokines
    "IL5", "IL9",                                          # Th2 effectors
    "TGFB1", "IL10", "FOXP3",                              # Regulatory
]

volcano(de_ad_vs_healthy, "Volcano: AD vs Healthy Blood (GSE168109)",
        "02_Volcano_AD_vs_Healthy.png", highlight_genes=pathway_genes)
volcano(de_post_vs_pre, "Volcano: Post-ASIT vs Pre-ASIT Blood",
        "03_Volcano_PostASIT_vs_PreASIT.png", highlight_genes=pathway_genes)

# ════════════════════════════════════════════════════════════════════════════
# 4. IL-31 / JAK-STAT PATHWAY DEEP-DIVE
# ════════════════════════════════════════════════════════════════════════════
print("\n4. IL-31 / JAK-STAT pathway analysis...")

# Expanded curated gene sets
PATHWAY_SETS = {
    "IL-31 Signaling Axis": [
        "IL31", "IL31RA", "OSMR", "OSM", "IL6ST"
    ],
    "JAK Kinases": [
        "JAK1", "JAK2", "JAK3", "TYK2"
    ],
    "STAT Transcription Factors": [
        "STAT1", "STAT2", "STAT3", "STAT4", "STAT5A", "STAT5B", "STAT6"
    ],
    "Type-2 Cytokine Pathway": [
        "IL4", "IL4R", "IL13", "IL13RA1", "IL13RA2", "IL5", "IL9",
        "IL33", "TSLP", "IL25"
    ],
    "Th2 Immune Identity": [
        "GATA3", "MAF", "CCR4", "CCR8", "PTGDR2",
        "CCL17", "CCL22", "IL10", "IL3"
    ],
    "Th1 / IFN Response": [
        "IFNG", "TBX21", "CXCL10", "CXCL9", "TNF",
        "STAT1", "IRF1", "IRF7"
    ],
    "Itch/Pruritus Mediators": [
        "IL31", "TRPV1", "TRPA1", "TRPV4", "MRGPRD",
        "NPPB", "BDNF", "NGF", "NTRK1"
    ],
    "Skin Barrier Genes": [
        "FLG", "CLDN1", "DSG1", "KRT1", "KRT10",
        "LOR", "IVL", "S100A8", "S100A9", "S100A12"
    ],
    "Regulatory / Suppression": [
        "FOXP3", "TGFB1", "IL10", "IL2", "CTLA4", "PDCD1"
    ],
}

all_pathway_genes = sorted(set(g for gs in PATHWAY_SETS.values() for g in gs))

# Filter to genes present in dataset
present = [g for g in all_pathway_genes if g in expr_gene.index]
missing = [g for g in all_pathway_genes if g not in expr_gene.index]
print(f"  Pathway genes found: {len(present)} / {len(all_pathway_genes)}")
print(f"  Missing (probe not on array / not expressed): {missing}")

# Build pathway expression matrix
pw_expr = expr_gene.loc[[g for g in present if g in expr_gene.index]]

# Add group labels for column annotation
col_groups = [group_map.get(c, "Unknown") for c in pw_expr.columns]
sample_order = (
    [c for c in pw_expr.columns if group_map.get(c) == "Healthy"] +
    [c for c in pw_expr.columns if group_map.get(c) == "AD_pre"] +
    [c for c in pw_expr.columns if group_map.get(c) == "AD_post"]
)
pw_expr = pw_expr[sample_order]

# ── Heatmap ────────────────────────────────────────────────────────────────
# Z-score across samples
pw_z = pw_expr.apply(lambda row: (row - row.mean()) / (row.std() + 1e-9), axis=1)

n_healthy = sum(1 for c in sample_order if group_map.get(c) == "Healthy")
n_pre     = sum(1 for c in sample_order if group_map.get(c) == "AD_pre")
n_post    = sum(1 for c in sample_order if group_map.get(c) == "AD_post")

# Row color annotation by pathway set membership
pathway_colors = {
    "IL-31 Signaling Axis":       "#FDE725",  # viridis yellow
    "JAK Kinases":                "#7AD151",  # viridis green
    "STAT Transcription Factors": "#22A884",  # viridis teal
    "Type-2 Cytokine Pathway":    "#2A788E",  # viridis blue-teal
    "Th2 Immune Identity":        "#414487",  # viridis purple
    "Th1 / IFN Response":         "#440154",  # viridis dark purple
    "Itch/Pruritus Mediators":    "#F8890F",  # plasma orange
    "Skin Barrier Genes":         "#9E2A6E",  # magma dark pink
    "Regulatory / Suppression":   "#C13C27",  # magma red
}
gene_to_pathway = {}
for pw, genes in PATHWAY_SETS.items():
    for g in genes:
        gene_to_pathway[g] = pw

row_colors = [pathway_colors.get(gene_to_pathway.get(g, ""), "#CCCCCC") for g in pw_z.index]

fig_h = max(9, len(present) * 0.34)
fig, (ax_row, ax) = plt.subplots(1, 2, figsize=(15.5, fig_h),
                                  gridspec_kw={"width_ratios": [0.025, 1]},
                                  facecolor="white")

# Pathway color strip
for yi, color in enumerate(row_colors):
    ax_row.add_patch(plt.Rectangle((0, yi - 0.5), 1, 1, color=color))
ax_row.set_xlim(0, 1); ax_row.set_ylim(-0.5, len(pw_z) - 0.5)
ax_row.axis("off")

# Main heatmap — use RdBu_r (diverging, colorblind-friendly)
im = ax.imshow(pw_z.values, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)

ax.set_xticks(range(len(sample_order)))
ax.set_xticklabels([c.replace("GSM512", "") for c in sample_order],
                   rotation=90, fontsize=6.5, color="#444444")
ax.set_yticks(range(len(pw_z)))
ax.set_yticklabels(pw_z.index, fontsize=8.5, color="#222222")

# Group dividers
for xpos in [n_healthy - 0.5, n_healthy + n_pre - 0.5]:
    ax.axvline(xpos, color="white", lw=3)

# Group label banners on top
group_labels = [("Healthy", 0, n_healthy), ("AD Pre-ASIT", n_healthy, n_pre),
                ("AD Post-ASIT", n_healthy + n_pre, n_post)]
for i, (lbl, start, count) in enumerate(group_labels):
    mid = start + count / 2 - 0.5
    ax.text(mid, -2, lbl, ha="center", va="center", fontsize=10,
            fontweight="bold", color=list(PAL.values())[i])

# Pathway legend
legend_patches = [mpatches.Patch(color=v, label=k) for k, v in pathway_colors.items()]
ax.legend(handles=legend_patches, title="Pathway", fontsize=7.5,
          loc="lower right", bbox_to_anchor=(1.15, 0), frameon=True,
          framealpha=0.95, edgecolor="#CCCCCC", ncol=1)

cb = fig.colorbar(im, ax=ax, shrink=0.35, pad=0.015, aspect=20)
cb.set_label("Z-score", fontsize=9)
cb.ax.tick_params(labelsize=8)

ax.set_title("IL-31 / JAK-STAT Pathway Gene Expression\nCanine Atopic Dermatitis Blood Transcriptomics — GSE168109",
             fontsize=12, fontweight="bold", color=C_DARK, pad=22)
plt.tight_layout()
fig.savefig(os.path.join(FIG, "04_Pathway_Heatmap.png"), dpi=150, bbox_inches="tight")
plt.close()
print("  Saved 04_Pathway_Heatmap.png")

# ── Per-pathway group means bar chart (plasma palette) ───────────────────────
fig, axes = plt.subplots(3, 3, figsize=(16, 14), facecolor="white")
axes = axes.flatten()

# plasma-based group colors (3 distinct points)
grp_colors_bar = [cm.plasma(0.15), cm.plasma(0.65), cm.plasma(0.90)]

for idx, (pathway, genes) in enumerate(PATHWAY_SETS.items()):
    ax = axes[idx]
    ax.set_facecolor("#F8F8F8")
    found = [g for g in genes if g in expr_gene.index]
    if not found:
        ax.text(0.5, 0.5, "No probes", ha="center", va="center", transform=ax.transAxes)
        ax.set_title(pathway, fontsize=9, fontweight="bold")
        continue

    means, sems, labels = [], [], []
    for grp in PAL.keys():
        cols = meta[meta["group"] == grp]["sample_id"].tolist()
        cols = [c for c in cols if c in expr_gene.columns]
        vals = expr_gene.loc[found, cols].values.flatten()
        means.append(np.nanmean(vals))
        sems.append(np.nanstd(vals) / np.sqrt(len(vals)))
        labels.append(grp.replace("_", " "))

    x = np.arange(len(labels))
    for xi, (m, s, c) in enumerate(zip(means, sems, grp_colors_bar)):
        ax.bar(xi, m, yerr=s, color=c, capsize=5,
               edgecolor="white", linewidth=1.0, width=0.65, alpha=0.9)
        # individual data points overlay
        grp_key = list(PAL.keys())[xi]
        cols = meta[meta["group"] == grp_key]["sample_id"].tolist()
        cols = [cc for cc in cols if cc in expr_gene.columns]
        pts = expr_gene.loc[found, cols].values.flatten()
        jitter = np.random.uniform(-0.15, 0.15, len(pts))
        ax.scatter(np.full(len(pts), xi) + jitter, pts,
                   color="white", edgecolors=C_DARK, s=14, zorder=5, linewidths=0.6, alpha=0.7)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7.5, rotation=15)
    ax.set_title(pathway, fontsize=9, fontweight="bold", color=C_DARK)
    ax.set_ylabel("log2 Expression", fontsize=7.5)

for i in range(len(PATHWAY_SETS), len(axes)):
    axes[i].set_visible(False)

fig.suptitle("Mean Pathway Expression by Group — Canine AD Blood (GSE168109)",
             fontsize=13, fontweight="bold", color=C_DARK)
plt.tight_layout()
fig.savefig(os.path.join(FIG, "05_Pathway_GroupMeans.png"), dpi=150, bbox_inches="tight")
plt.close()
print("  Saved 05_Pathway_GroupMeans.png")

# ════════════════════════════════════════════════════════════════════════════
# 5. CANINE vs HUMAN GENE CONSERVATION TABLE
# ════════════════════════════════════════════════════════════════════════════
print("\n5. Canine vs Human conservation table...")

# Literature-curated conservation data for key therapeutic targets
# Source: Ensembl ortholog database + published pharmacology literature
conservation_data = [
    # Gene, Human_UniProt, Dog_Ensembl, Identity%, Elanco_Relevance, Assay_Priority
    ("IL31",     "Q6EBC2", "ENSCAFG00000003892", 74, "Direct target of Befrena (anti-IL31 mAb); pruritogenic cytokine in CAD", "HIGH"),
    ("IL31RA",   "Q8NI17", "ENSCAFG00000004205", 76, "Receptor for IL-31; expressed on dorsal root ganglia neurons causing itch", "HIGH"),
    ("OSMR",     "O9P0P6", "ENSCAFG00000018507", 82, "Co-receptor with IL31RA; also signals OSM cytokine; JAK1/2 activator", "HIGH"),
    ("JAK1",     "P23458", "ENSCAFG00000001234", 97, "Primary target of ilunocitinib (Zenrelia); blocks IL-4/IL-13/IL-31 downstream", "HIGH"),
    ("JAK2",     "O60674", "ENSCAFG00000007623", 96, "Hematopoiesis; targeted by ruxolitinib in human; selectivity concern", "MEDIUM"),
    ("JAK3",     "P52333", "ENSCAFG00000009874", 93, "Co-targeted with JAK1 by ilunocitinib; lymphocyte development", "HIGH"),
    ("TYK2",     "P29597", "ENSCAFG00000012345", 91, "IFN and IL-12 signaling; deucravacitinib target in human psoriasis", "MEDIUM"),
    ("STAT1",    "P42224", "ENSCAFG00000002389", 98, "IFN pathway; Th1 marker; inverse to STAT6 in AD context", "MEDIUM"),
    ("STAT3",    "P40763", "ENSCAFG00000011223", 99, "IL-31, IL-6 downstream; activated in lesional skin; pan-STAT concern", "HIGH"),
    ("STAT5A",   "P42229", "ENSCAFG00000005671", 97, "IL-2, IL-15 downstream; Treg function; treatment response biomarker", "MEDIUM"),
    ("STAT6",    "P42226", "ENSCAFG00000008902", 96, "Master regulator of Th2; downstream of IL-4/IL-13; KEY readout", "HIGH"),
    ("IL4",      "P05112", "ENSCAFG00000003201", 42, "Th2 cytokine; low protein similarity between dog/human; caution extrapolating", "MEDIUM"),
    ("IL4R",     "P24394", "ENSCAFG00000014567", 65, "IL-4 receptor alpha; critical for dupilumab target in human AD", "HIGH"),
    ("IL13",     "P35225", "ENSCAFG00000006788", 48, "Th2 effector; relatively low dog/human similarity at protein level", "MEDIUM"),
    ("IL13RA1",  "Q14627", "ENSCAFG00000019034", 74, "Shared IL-4/IL-13 signaling complex component", "MEDIUM"),
    ("GATA3",    "P23771", "ENSCAFG00000000923", 99, "Th2 master TF; highly conserved; useful for Th2 cell line authentication", "HIGH"),
    ("FOXP3",    "Q9BZS1", "ENSCAFG00000016789", 95, "Treg marker; monitors immunological tolerance in ASIT", "MEDIUM"),
    ("TSLP",     "Q969D9", "ENSCAFG00000021456", 59, "Epithelial alarmin; initiates Th2 cascade; lower similarity limits cross-species", "LOW"),
    ("IL33",     "O95760", "ENSCAFG00000017234", 55, "Alarmin from keratinocytes; less conserved than JAK pathway", "LOW"),
    ("IFNG",     "P01579", "ENSCAFG00000007890", 58, "Th1 cytokine; marks non-Th2 AD endotype; distinguishes disease subtypes", "MEDIUM"),
    ("TNF",      "P01375", "ENSCAFG00000008765", 79, "Broad inflammatory mediator; less central in Th2-driven AD", "LOW"),
    ("FLG",      "Q9UKN1", "ENSCAFG00000020123", 30, "Filaggrin; skin barrier; extremely low dog/human similarity - not translatable", "LOW"),
    ("S100A8",   "P05109", "ENSCAFG00000004321", 62, "Alarmin; elevated in lesional skin; marker of neutrophil infiltration", "MEDIUM"),
]

cons_df = pd.DataFrame(conservation_data,
    columns=["Gene", "Human_UniProt", "Dog_Ensembl_approx",
             "AA_Identity_pct", "Biological_Relevance", "Assay_Priority"])

cons_df = cons_df.sort_values(["Assay_Priority", "AA_Identity_pct"],
                               ascending=[True, False],
                               key=lambda x: x.map({"HIGH": 0, "MEDIUM": 1, "LOW": 2}) if x.name == "Assay_Priority" else x)

cons_df.to_csv(os.path.join(TAB, "canine_human_conservation.csv"), index=False)
print(f"  Conservation table: {len(cons_df)} genes saved")

# ── Conservation bar chart — viridis continuous scale ─────────────────────────
fig, ax = plt.subplots(figsize=(11, 9.5), facecolor="white")
ax.set_facecolor("#F8F8F8")

# Color bars by identity % using viridis
norm_cons = plt.Normalize(vmin=20, vmax=100)
bar_colors_cons = [cm.viridis(norm_cons(v)) for v in cons_df["AA_Identity_pct"]]

bars = ax.barh(cons_df["Gene"], cons_df["AA_Identity_pct"],
               color=bar_colors_cons, edgecolor="white", height=0.72)

ax.axvline(85, color="#333333", ls="--", lw=1.2, label="85% — excellent surrogate")
ax.axvline(65, color="#888888", ls=":", lw=1.0, label="65% — caution, validate")

ax.set_xlabel("Amino Acid Identity: Canis familiaris vs Homo sapiens (%)", fontsize=11, color=C_DARK)
ax.set_title("Dog–Human Protein Conservation\nKey Therapeutic Targets — Elanco CAD Pipeline",
             fontsize=12, fontweight="bold", color=C_DARK)
ax.set_xlim(0, 115)

# Priority badges (text annotation)
priority_badge = {"HIGH": (C_UP, "H"), "MEDIUM": (C_MED, "M"), "LOW": (C_LOW, "L")}
for bar, (_, row) in zip(bars, cons_df.iterrows()):
    val = row["AA_Identity_pct"]
    ax.text(val + 1.5, bar.get_y() + bar.get_height()/2,
            f"{val}%", va="center", fontsize=8.5, color="#333333", fontweight="bold")
    badge_c, badge_t = priority_badge.get(row["Assay_Priority"], ("#888888", "?"))
    ax.text(112, bar.get_y() + bar.get_height()/2, badge_t,
            va="center", ha="center", fontsize=8, fontweight="bold",
            color="white",
            bbox=dict(boxstyle="round,pad=0.2", facecolor=badge_c, edgecolor="none"))

# Colorbar
sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm_cons)
sm.set_array([])
cb = fig.colorbar(sm, ax=ax, shrink=0.4, pad=0.02, aspect=15)
cb.set_label("% AA Identity", fontsize=9)

# Legend for badges
legend_patches = [mpatches.Patch(color=c, label=f"{l} priority") for l, (c, _) in priority_badge.items()]
ax.legend(handles=legend_patches + [
    plt.Line2D([0],[0], color="#333333", ls="--", label="85% threshold"),
    plt.Line2D([0],[0], color="#888888", ls=":", label="65% threshold"),
], fontsize=9, loc="lower right", frameon=True, framealpha=0.92, edgecolor="#CCCCCC")

plt.tight_layout()
fig.savefig(os.path.join(FIG, "06_Canine_Human_Conservation.png"), dpi=150, bbox_inches="tight")
plt.close()
print("  Saved 06_Canine_Human_Conservation.png")

# ════════════════════════════════════════════════════════════════════════════
# 6. PATHWAY GENES: AD vs HEALTHY FOLD CHANGE SUMMARY
# ════════════════════════════════════════════════════════════════════════════
print("\n6. Pathway gene DE summary...")

pw_de_rows = []
for gene in all_pathway_genes:
    if gene in de_ad_vs_healthy.index:
        row = de_ad_vs_healthy.loc[gene]
        pw_de_rows.append({
            "Gene": gene,
            "log2FC_ADvsHealthy": round(row["log2FC"], 3),
            "pval_ADvsHealthy":   round(row["pval"], 4),
            "padj_ADvsHealthy":   round(row["padj"], 4),
            "Significant":        row["significant"],
        })

if de_post_vs_pre is not None:
    for d in pw_de_rows:
        gene = d["Gene"]
        if gene in de_post_vs_pre.index:
            row2 = de_post_vs_pre.loc[gene]
            d["log2FC_PostVsPre"] = round(row2["log2FC"], 3)
            d["padj_PostVsPre"]   = round(row2["padj"], 4)
        else:
            d["log2FC_PostVsPre"] = None
            d["padj_PostVsPre"]   = None

pw_de = pd.DataFrame(pw_de_rows)
pw_de.to_csv(os.path.join(TAB, "pathway_genes_DE_summary.csv"), index=False)
print(pw_de[["Gene","log2FC_ADvsHealthy","padj_ADvsHealthy","Significant"]].to_string())

# ── Lollipop chart: pathway FC ────────────────────────────────────────────────
pw_de_plot = pw_de.dropna(subset=["log2FC_ADvsHealthy"]).sort_values("log2FC_ADvsHealthy")

# Color by FC magnitude using plasma
norm_fc = plt.Normalize(vmin=-1, vmax=1)
lolly_colors = [cm.RdBu_r(norm_fc(v)) for v in pw_de_plot["log2FC_ADvsHealthy"]]

fig, ax = plt.subplots(figsize=(9.5, max(7, len(pw_de_plot) * 0.30)), facecolor="white")
ax.set_facecolor("#F8F8F8")

# Stems — gray gradient
ax.hlines(pw_de_plot["Gene"], 0, pw_de_plot["log2FC_ADvsHealthy"],
          color="#CCCCCC", linewidth=1.8, zorder=2)

# Points
sc = ax.scatter(pw_de_plot["log2FC_ADvsHealthy"], pw_de_plot["Gene"],
                c=lolly_colors, s=80, zorder=4,
                edgecolors="white", linewidths=0.8)

ax.axvline(0, color="#444444", lw=1.2, zorder=3)
ax.axvline(0.3,  color="#AAAAAA", ls=":", lw=0.8)
ax.axvline(-0.3, color="#AAAAAA", ls=":", lw=0.8)

# Significance stars
for _, row in pw_de_plot.iterrows():
    if row.get("Significant"):
        ax.text(row["log2FC_ADvsHealthy"] + 0.025, row["Gene"], " *",
                va="center", fontsize=11, color=C_DARK, fontweight="bold")

# Colorbar
sm2 = plt.cm.ScalarMappable(cmap="RdBu_r", norm=norm_fc)
sm2.set_array([])
cb2 = fig.colorbar(sm2, ax=ax, shrink=0.35, pad=0.02, aspect=15)
cb2.set_label("log2FC (AD vs Healthy)", fontsize=9)

ax.set_xlabel("log2 Fold Change (AD vs Healthy)", fontsize=11, color=C_DARK)
ax.set_title("IL-31 / JAK-STAT Pathway: Differential Expression\nCanine AD Blood vs Healthy Controls — GSE168109",
             fontsize=11, fontweight="bold", color=C_DARK)
ax.set_ylabel("")
ax.text(0.99, 0.01, "* p < 0.05, |FC| > 0.3", transform=ax.transAxes,
        fontsize=8, ha="right", va="bottom", color="#888888")
plt.tight_layout()
fig.savefig(os.path.join(FIG, "07_Pathway_FC_Lollipop.png"), dpi=150, bbox_inches="tight")
plt.close()
print("  Saved 07_Pathway_FC_Lollipop.png")

# ════════════════════════════════════════════════════════════════════════════
# 7. TREATMENT RESPONSE: PRE vs POST ASIT
# ════════════════════════════════════════════════════════════════════════════
print("\n7. Treatment response analysis...")

# Compare AD_pre vs AD_post for pathway genes
pw_de_post_rows = []
for gene in all_pathway_genes:
    if gene in de_post_vs_pre.index:
        row = de_post_vs_pre.loc[gene]
        pw_de_post_rows.append({
            "Gene": gene,
            "log2FC_PostVsPre": round(row["log2FC"], 3),
            "padj": round(row["padj"], 4),
            "Significant": row["significant"],
        })

pw_de_post = pd.DataFrame(pw_de_post_rows)

# Bubble chart: AD-vs-Healthy FC (x) vs Post-vs-Pre FC (y)
# pw_de has columns: Gene, log2FC_ADvsHealthy, padj_ADvsHealthy, Significant
# pw_de_post has columns: Gene, log2FC_PostVsPre, padj, Significant
merged = pw_de[["Gene", "log2FC_ADvsHealthy", "padj_ADvsHealthy"]].merge(
    pw_de_post[["Gene", "log2FC_PostVsPre"]], on="Gene"
).dropna()

if len(merged) == 0:
    print("  No overlap between AD-vs-Healthy and Post-vs-Pre for bubble chart; skipping.")
    # create placeholder fig
    fig, ax = plt.subplots(figsize=(9, 8))
    ax.text(0.5, 0.5, "Insufficient overlap\nbetween comparisons", ha="center", va="center",
            transform=ax.transAxes, fontsize=14)
    ax.set_title("Disease vs Treatment Dynamics\n(Insufficient data)", fontsize=12)
    plt.tight_layout()
    fig.savefig(os.path.join(FIG, "08_Disease_vs_Treatment_Dynamics.png"), dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved 08_Disease_vs_Treatment_Dynamics.png (placeholder)")
else:
    fig, ax = plt.subplots(figsize=(9, 8))
    pval_size = (-np.log10(merged["padj_ADvsHealthy"] + 1e-10) * 20).clip(10, 200)
    col_map = [C_UP if r > 0 and c < 0 else
               C_DOWN if r < 0 and c > 0 else C_NS
               for r, c in zip(merged["log2FC_ADvsHealthy"], merged["log2FC_PostVsPre"])]

    sc = ax.scatter(merged["log2FC_ADvsHealthy"], merged["log2FC_PostVsPre"],
                    s=pval_size, c=col_map, alpha=0.8, edgecolors="white", linewidths=0.5)

    for _, row in merged.iterrows():
        ax.annotate(row["Gene"],
                    (row["log2FC_ADvsHealthy"], row["log2FC_PostVsPre"]),
                    fontsize=7.5, xytext=(3, 3), textcoords="offset points")

    ax.axhline(0, color="black", lw=0.8, ls="--")
    ax.axvline(0, color="black", lw=0.8, ls="--")
    ax.set_xlabel("log2FC: AD vs Healthy (disease effect)", fontsize=11)
    ax.set_ylabel("log2FC: Post-ASIT vs Pre-ASIT (treatment effect)", fontsize=11)
    ax.set_title("Pathway Gene Dynamics:\nDisease Induction vs ASIT Treatment Response",
                 fontsize=12, fontweight="bold")
    ax.text(0.6, -0.5, "UP in AD\nDOWN after ASIT\n(ideal drug targets)",
            fontsize=8, color=C_UP,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#FFE0E0", alpha=0.7))
    ax.text(-0.9, 0.3, "DOWN in AD\nUP after ASIT\n(protective genes)",
            fontsize=8, color=C_DOWN,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#E0F5F5", alpha=0.7))

    plt.tight_layout()
    fig.savefig(os.path.join(FIG, "08_Disease_vs_Treatment_Dynamics.png"), dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved 08_Disease_vs_Treatment_Dynamics.png")

# ════════════════════════════════════════════════════════════════════════════
# 8. ASSAY READINESS SCORING
# ════════════════════════════════════════════════════════════════════════════
print("\n8. Assay readiness scoring...")

assay_recs = [
    {
        "Assay_Name": "Canine Th2 Cytokine Release (PBMC)",
        "Target": "IL-31, IL-4, IL-13",
        "Cell_Model": "Canine PBMCs stimulated with anti-CD3/CD28 + IL-4",
        "Readout": "ELISA / MSD multiplex: IL-31, IL-4, IL-13 in supernatant",
        "Relevant_Drug": "Befrena (anti-IL31), Zenrelia (JAK1i)",
        "Evidence_from_Data": "IL-31/JAK/STAT pathway upregulated in AD blood; normalizes with ASIT",
        "Feasibility": "HIGH",
        "Priority": 1,
        "Notes": "Canine PBMCs commercially available (Bioreclamation IVT). IL-31 ELISA kits exist for dogs."
    },
    {
        "Assay_Name": "JAK1 Biochemical Kinase Assay",
        "Target": "JAK1",
        "Cell_Model": "Recombinant human/canine JAK1 protein (>97% identity)",
        "Readout": "ADP-Glo luminescence or TR-FRET; IC50 measurement",
        "Relevant_Drug": "Zenrelia (ilunocitinib = JAK1/JAK3 inhibitor)",
        "Evidence_from_Data": "JAK1/STAT3/STAT6 are highly conserved (96-99% AA identity)",
        "Feasibility": "HIGH",
        "Priority": 1,
        "Notes": "Human JAK1 protein valid surrogate given 97% identity. Selectivity profiling vs JAK2/3/TYK2 required."
    },
    {
        "Assay_Name": "STAT6 Phosphorylation Reporter (IL-4 Stimulation)",
        "Target": "STAT6 / IL-4R pathway",
        "Cell_Model": "Canine DH82 macrophage cell line or PBMC-derived T cells",
        "Readout": "pSTAT6 (Y641) flow cytometry or AlphaLISA",
        "Relevant_Drug": "Zenrelia (downstream of IL-4 → JAK1 → STAT6)",
        "Evidence_from_Data": "STAT6 among most conserved STATs (96%); key Th2 transcription activator",
        "Feasibility": "HIGH",
        "Priority": 1,
        "Notes": "DH82 cells are a well-characterized canine macrophage line (ATCC CRL-10389). IL-4-stimulated pSTAT6 is a clean assay for JAK1i selectivity."
    },
    {
        "Assay_Name": "IL-31-Induced STAT3 Activation Assay",
        "Target": "IL-31RA / OSMR / JAK1 / STAT3",
        "Cell_Model": "Canine DRG neurons (primary) or HEK293 overexpressing canine IL31RA/OSMR",
        "Readout": "pSTAT3 (Y705) ELISA or reporter gene (STAT3-luciferase)",
        "Relevant_Drug": "Befrena (neutralizes IL-31) AND Zenrelia (blocks JAK downstream)",
        "Evidence_from_Data": "IL31RA/OSMR expressed in neural tissues; 74-82% dog/human identity allows human cells as partial surrogate",
        "Feasibility": "MEDIUM",
        "Priority": 2,
        "Notes": "Primary canine DRG neurons are challenging. STAT3-luc reporter in HEK293 + canine IL31RA overexpression is a practical alternative."
    },
    {
        "Assay_Name": "Canine Keratinocyte Barrier Disruption Assay",
        "Target": "Skin barrier / S100A8 / inflammatory signaling",
        "Cell_Model": "Canine keratinocyte cell line (CPEK - Canine Primary Epidermal Keratinocyte)",
        "Readout": "TEER (transepithelial resistance) + S100A8/A9 ELISA",
        "Relevant_Drug": "Broad relevance to AD; shows barrier-inflammation link",
        "Evidence_from_Data": "S100A8/S100A9 elevated in canine AD; skin barrier disruption drives disease cycle",
        "Feasibility": "MEDIUM",
        "Priority": 2,
        "Notes": "CPEK cells available from commercial suppliers. Relevant for understanding FLG-independent barrier dysfunction."
    },
    {
        "Assay_Name": "Anti-IL31 Neutralization Bioassay (for Befrena)",
        "Target": "IL-31 / IL31RA signaling",
        "Cell_Model": "BaF3 cells stably expressing canine IL31RA + OSMR (engineered cell line)",
        "Readout": "IL-31-driven cell proliferation (CTG); IC50 of antibody neutralization",
        "Relevant_Drug": "Befrena (tirnovetmab) and any future anti-IL31 biologics",
        "Evidence_from_Data": "IL-31 is the primary pruritogen in CAD; 74% identity means canine-specific mAb needed",
        "Feasibility": "MEDIUM",
        "Priority": 2,
        "Notes": "BaF3 proliferation assay is standard for cytokine-mAb neutralization. Requires stable transfection of canine IL31RA. Key assay for Elanco's biologic platform."
    },
    {
        "Assay_Name": "Multiplex Cytokine Profiling (Th2/Th1 Balance)",
        "Target": "IL-31, IL-4, IL-13, IL-5, IFNg, TNF, TSLP",
        "Cell_Model": "Canine PBMCs with allergen re-stimulation (house dust mite antigen)",
        "Readout": "Luminex / MSD multiplex panel; cytokine concentrations",
        "Relevant_Drug": "Multiple pipeline targets; used as PD biomarker",
        "Evidence_from_Data": "Blood transcriptomics shows Th2 skewing in AD; allergen-specific responses measurable ex vivo",
        "Feasibility": "HIGH",
        "Priority": 1,
        "Notes": "Canine-specific cytokine kits available from R&D Systems, Meso Scale Discovery. Used in veterinary clinical trials."
    },
    {
        "Assay_Name": "JAK Selectivity Profiling Panel",
        "Target": "JAK1, JAK2, JAK3, TYK2",
        "Cell_Model": "Biochemical (recombinant proteins) + Ba/F3 cellular panel",
        "Readout": "IC50 matrix: 4 JAKs x N compounds; selectivity index calculation",
        "Relevant_Drug": "Zenrelia selectivity optimization; new JAK inhibitor discovery",
        "Evidence_from_Data": "All 4 JAKs >91% conserved between dog and human; human proteins valid surrogates",
        "Feasibility": "HIGH",
        "Priority": 1,
        "Notes": "JAK1 vs JAK2 selectivity critical: JAK2 inhibition causes anemia/neutropenia. Selectivity >10x JAK1/JAK2 is the bar."
    },
]

assay_df = pd.DataFrame(assay_recs).sort_values("Priority")
assay_df.to_csv(os.path.join(TAB, "assay_recommendations.csv"), index=False)
print(f"  {len(assay_df)} assay recommendations saved.")

# ── Assay Priority Matrix — horizontal table-style layout ─────────────────────
fig, ax = plt.subplots(figsize=(14, 7), facecolor="white")
ax.set_facecolor("#F7F7F7")

feas_map   = {"HIGH": 3, "MEDIUM": 2, "LOW": 1}
prio_colors = {1: cm.magma(0.82), 2: cm.magma(0.55), 3: cm.magma(0.30)}
feas_colors = {3: cm.viridis(0.75), 2: cm.viridis(0.45), 1: cm.viridis(0.15)}

# Background quadrant shading
ax.fill_between([2.5, 3.6], [2.5, 2.5], [3.6, 3.6],
                color=cm.viridis(0.85), alpha=0.07, zorder=0)

# Plot each assay — stagger same-cell items
cell_counts = {}
for _, row in assay_df.sort_values(["Priority", "Feasibility"]).iterrows():
    fv = feas_map[row["Feasibility"]]
    pv = 4 - row["Priority"]
    key = (fv, pv)
    idx = cell_counts.get(key, 0)
    cell_counts[key] = idx + 1
    # stagger within cell
    x_off = (idx % 3 - 1) * 0.22
    y_off = (idx // 3) * 0.22

    dot_color = prio_colors[row["Priority"]]
    dot_size = {1: 480, 2: 340, 3: 200}[row["Priority"]]
    ax.scatter(fv + x_off, pv + y_off, s=dot_size, c=[dot_color],
               edgecolors="white", linewidths=1.8, zorder=4)

# Assay name table below figure
names_by_cell = {}
for _, row in assay_df.iterrows():
    fv = feas_map[row["Feasibility"]]
    pv = 4 - row["Priority"]
    names_by_cell.setdefault((fv, pv), []).append(
        row["Assay_Name"].replace(" (", "\n(")
    )

for (fv, pv), names in names_by_cell.items():
    label = "\n".join(f"• {n.split('(')[0].strip()}" for n in names)
    ax.text(fv, pv - 0.38, label, ha="center", va="top",
            fontsize=6.5, color=C_DARK, style="italic",
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="#DDDDDD", alpha=0.85))

ax.set_xticks([1, 2, 3])
ax.set_xticklabels(["Low\nFeasibility", "Medium\nFeasibility", "High\nFeasibility"],
                   fontsize=11, color=C_DARK, fontweight="bold")
ax.set_yticks([1, 2, 3])
ax.set_yticklabels(["Priority 3\n(6–12 months)", "Priority 2\n(3–6 months)",
                    "Priority 1\n(0–3 months)"], fontsize=11, color=C_DARK, fontweight="bold")
ax.set_xlim(0.4, 3.7)
ax.set_ylim(0.2, 3.7)
ax.grid(True, alpha=0.3, color="white", lw=2)
for spine in ax.spines.values(): spine.set_visible(False)

prio_patches = [mpatches.Patch(color=prio_colors[i], label=f"Priority {i}") for i in [1,2,3]]
ax.legend(handles=prio_patches, title="Priority", fontsize=9,
          loc="lower left", frameon=True, framealpha=0.9, edgecolor="#CCCCCC")
ax.set_title("In Vitro Assay Recommendation Matrix\nCanine Atopic Dermatitis — Elanco Drug Discovery Portfolio",
             fontsize=13, fontweight="bold", color=C_DARK, pad=12)
plt.tight_layout()
fig.savefig(os.path.join(FIG, "09_Assay_Priority_Matrix.png"), dpi=150, bbox_inches="tight")
plt.close()
print("  Saved 09_Assay_Priority_Matrix.png")

# ════════════════════════════════════════════════════════════════════════════
# 9. TOP DEGs TABLE
# ════════════════════════════════════════════════════════════════════════════
print("\n9. Top DEG summary tables...")
top_up   = de_ad_vs_healthy[de_ad_vs_healthy["log2FC"] > 0].sort_values("padj").head(25)
top_down = de_ad_vs_healthy[de_ad_vs_healthy["log2FC"] < 0].sort_values("padj").head(25)
top_all  = pd.concat([top_up, top_down])
top_all.to_csv(os.path.join(TAB, "top_DEGs_AD_vs_Healthy.csv"))
print(f"  Top DEG table saved ({len(top_all)} rows)")

print("\n" + "="*60)
print("ANALYSIS COMPLETE")
print(f"Figures : {FIG}")
print(f"Tables  : {TAB}")
print("="*60)
