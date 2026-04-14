"""
Download GEO datasets for canine atopic dermatitis analysis.
Datasets:
  - GSE168109: Canine AD microarray (healthy vs AD pre-ASIT vs AD post-ASIT)
  - Human pathway gene lists are embedded (literature-curated)
"""

import GEOparse
import pandas as pd
import os

RAW_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "raw")
PROC_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "processed")
os.makedirs(RAW_DIR, exist_ok=True)
os.makedirs(PROC_DIR, exist_ok=True)

# ── 1. Download GSE168109 ────────────────────────────────────────────────────
print("Downloading GSE168109 (Canine AD microarray)...")
gse = GEOparse.get_GEO(geo="GSE168109", destdir=RAW_DIR, silent=False)

print(f"\nDataset title : {gse.metadata['title'][0]}")
print(f"Organism      : {gse.metadata.get('organism', ['unknown'])[0]}")
print(f"Samples (n)   : {len(gse.gsms)}")
print(f"Platforms     : {list(gse.gpls.keys())}")

# ── 2. Inspect sample metadata ───────────────────────────────────────────────
records = []
for gsm_name, gsm in gse.gsms.items():
    title       = gsm.metadata["title"][0]
    source      = gsm.metadata.get("source_name_ch1", ["NA"])[0]
    description = gsm.metadata.get("description", ["NA"])[0]
    char_raw    = gsm.metadata.get("characteristics_ch1", [])
    characteristics = " | ".join(char_raw)
    records.append({
        "sample_id": gsm_name,
        "title": title,
        "source": source,
        "description": description,
        "characteristics": characteristics
    })

meta_df = pd.DataFrame(records)
meta_path = os.path.join(PROC_DIR, "sample_metadata.csv")
meta_df.to_csv(meta_path, index=False)
print(f"\nSample metadata saved -> {meta_path}")
print(meta_df[["sample_id", "title", "source"]].to_string())

# ── 3. Extract expression matrix ─────────────────────────────────────────────
print("\nExtracting expression matrix...")
pivot_dfs = []
for gsm_name, gsm in gse.gsms.items():
    tbl = gsm.table.copy()
    if tbl.empty:
        print(f"  WARNING: {gsm_name} table is empty — skipping")
        continue
    # The value column is typically 'VALUE'
    val_col = "VALUE" if "VALUE" in tbl.columns else tbl.columns[-1]
    tbl = tbl[["ID_REF", val_col]].rename(columns={val_col: gsm_name})
    tbl["ID_REF"] = tbl["ID_REF"].astype(str)
    pivot_dfs.append(tbl.set_index("ID_REF"))

if pivot_dfs:
    expr_df = pd.concat(pivot_dfs, axis=1)
    expr_df = expr_df.apply(pd.to_numeric, errors="coerce")
    expr_path = os.path.join(PROC_DIR, "expression_matrix_raw.csv")
    expr_df.to_csv(expr_path)
    print(f"Expression matrix shape: {expr_df.shape}")
    print(f"Saved -> {expr_path}")
else:
    print("ERROR: No expression data extracted!")

# ── 4. Extract probe-to-gene mapping from platform ───────────────────────────
print("\nExtracting probe-gene mapping...")
for gpl_name, gpl in gse.gpls.items():
    print(f"  Platform: {gpl_name}")
    print(f"  Columns : {list(gpl.table.columns)}")
    gpl_path = os.path.join(PROC_DIR, f"platform_{gpl_name}.csv")
    gpl.table.to_csv(gpl_path, index=False)
    print(f"  Saved -> {gpl_path}")

print("\nDownload complete.")
