"""Shared helpers for the ANCA-GN Xenium kidney course.

Students receive only the raw Xenium output folder for slide 0011695 plus this
repo (which ships ``data/sample_bounding_boxes.csv`` and ``data/marker_genes.csv``).
These helpers read the raw slide, split it into the 8 tissue samples by their
bounding boxes, and load the Lake KPMP single-cell reference.
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "data"

# Data locations are resolved without a server-specific path: the XENIUM_RAW_DIR /
# LAKE_REF environment variables first, then a shared course location, then the
# user's home. Set the environment variables if your paths differ.
RAW_FOLDER = "output-XETG00088__0011695__Region_1__20240202__104242"


def xenium_raw_dir():
    """Folder with the raw Xenium output for slide 0011695."""
    cands = []
    if os.environ.get("XENIUM_RAW_DIR"):
        cands.append(Path(os.environ["XENIUM_RAW_DIR"]))
    cands += [Path("/home/shared") / RAW_FOLDER,
              Path("/home/shared/xenium_anca_kidney") / RAW_FOLDER,
              Path.home() / RAW_FOLDER]
    for c in cands:
        if (c / "cell_feature_matrix.h5").exists():
            return c
    raise FileNotFoundError(
        f"Xenium raw folder not found. Set XENIUM_RAW_DIR to the '{RAW_FOLDER}' "
        "folder for slide 0011695, or place it under /home/shared.")


def lake_reference_path():
    """Path to the Lake/KPMP single-cell reference h5ad."""
    cands = []
    if os.environ.get("LAKE_REF"):
        cands.append(Path(os.environ["LAKE_REF"]))
    cands += [Path("/home/shared/lake_reduced.h5ad"), Path.home() / "lake_reduced.h5ad"]
    for c in cands:
        if c.exists():
            return c
    raise FileNotFoundError(
        "Lake reference not found. Set LAKE_REF to lake_reduced.h5ad, or place it "
        "at /home/shared/lake_reduced.h5ad.")


def load_sample_boxes():
    """The per-sample bounding boxes (microns) that split the slide."""
    return pd.read_csv(DATA / "sample_bounding_boxes.csv")


def load_markers():
    """Kidney cell-type marker lists, as a dict {celltype: [genes]}."""
    df = pd.read_csv(DATA / "marker_genes.csv")
    return {r.celltype: r.markers.split(";") for r in df.itertuples()}


def read_xenium_raw(genes_only=True):
    """Read the raw slide into an AnnData (cell x gene), with the cell table in
    ``obs`` and control features separated out as QC columns."""
    d = xenium_raw_dir()
    adata = sc.read_10x_h5(d / "cell_feature_matrix.h5")
    adata.var_names_make_unique()
    cells = pd.read_parquet(d / "cells.parquet").set_index("cell_id")
    cells.index = cells.index.astype(str)
    adata.obs_names = adata.obs_names.astype(str)
    adata.obs = adata.obs.join(cells)
    adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].to_numpy()

    is_gene = adata.var.get("feature_types", "Gene Expression") == "Gene Expression"
    adata.obs["total_counts_all"] = np.asarray(adata.X.sum(1)).ravel()
    # Controls: this run keeps only Gene Expression in the matrix, but the cell
    # table carries the control counts, so build the QC fraction from there.
    ctrl_cols = [c for c in ["control_probe_counts", "control_codeword_counts",
                             "unassigned_codeword_counts", "deprecated_codeword_counts"]
                 if c in adata.obs]
    if ctrl_cols:
        adata.obs["control_counts"] = adata.obs[ctrl_cols].sum(1)
        denom = (adata.obs["total_counts_all"] + adata.obs["control_counts"]).clip(lower=1)
        adata.obs["control_frac"] = adata.obs["control_counts"] / denom
    elif (~is_gene).any():
        adata.obs["control_counts"] = np.asarray(adata[:, (~is_gene).values].X.sum(1)).ravel()
        adata.obs["control_frac"] = adata.obs["control_counts"] / adata.obs["total_counts_all"].clip(lower=1)
    if genes_only:
        adata = adata[:, is_gene.values].copy()
    return adata


def assign_samples(adata, boxes=None, drop_unassigned=True):
    """Add ``sample``, ``patient_sample_id`` and ``disease`` to ``adata.obs`` by
    locating each cell's centroid inside a sample bounding box. Cells in the gaps
    between tissue cores are dropped by default."""
    if boxes is None:
        boxes = load_sample_boxes()
    x = adata.obs["x_centroid"].to_numpy()
    y = adata.obs["y_centroid"].to_numpy()
    samp = np.full(adata.n_obs, "", dtype=object)
    psid = np.full(adata.n_obs, "", dtype=object)
    dis = np.full(adata.n_obs, "", dtype=object)
    for r in boxes.itertuples():
        inb = (x >= r.x_min) & (x <= r.x_max) & (y >= r.y_min) & (y <= r.y_max)
        samp[inb], psid[inb], dis[inb] = r.sample, r.patient_sample_id, r.disease
    adata.obs["sample"] = pd.Categorical(samp)
    adata.obs["patient_sample_id"] = pd.Categorical(psid)
    adata.obs["disease"] = pd.Categorical(dis)
    if drop_unassigned:
        adata = adata[adata.obs["sample"] != ""].copy()
        for c in ["sample", "patient_sample_id", "disease"]:
            adata.obs[c] = adata.obs[c].cat.remove_unused_categories()
    return adata


def load_or_build_qc(path="results/slide_0011695_qc.h5ad",
                     min_counts=10, min_genes=5, max_control=0.05):
    """Return the QC'd slide. Loads ``path`` if it exists (written by notebook 0),
    otherwise reads the raw slide, splits samples, filters, and caches it. This lets
    each notebook run standalone."""
    p = Path(path)
    if p.exists():
        return sc.read_h5ad(p)
    adata = assign_samples(read_xenium_raw())
    sc.pp.calculate_qc_metrics(adata, percent_top=None, inplace=True)
    if "control_frac" not in adata.obs:
        adata.obs["control_frac"] = 0.0
    keep = ((adata.obs["total_counts"] >= min_counts) &
            (adata.obs["n_genes_by_counts"] >= min_genes) &
            (adata.obs["control_frac"] <= max_control))
    adata = adata[keep].copy()
    p.parent.mkdir(parents=True, exist_ok=True)
    adata.write(p)
    return adata


def load_annotated(path="results/slide_0011695_annotated.h5ad"):
    """Return the slide with a working cell-type annotation in ``obs['celltype']``
    and log-normalised ``X`` (raw counts in ``layers['counts']``). Loads the cached
    object written by notebook 2 if present; otherwise builds a quick marker-based
    annotation so notebooks 3-5 run standalone."""
    p = Path(path)
    if p.exists():
        return sc.read_h5ad(p)
    adata = load_or_build_qc()
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    mk = {k: [g for g in v if g in adata.var_names] for k, v in load_markers().items()}
    mk = {k: v for k, v in mk.items() if v}
    for ct, gs in mk.items():
        sc.tl.score_genes(adata, gs, score_name=f"s::{ct}")
    S = adata.obs[[f"s::{ct}" for ct in mk]].to_numpy()
    cats = list(mk)
    adata.obs["celltype"] = pd.Categorical([cats[i] for i in S.argmax(1)])
    p.parent.mkdir(parents=True, exist_ok=True)
    adata.write(p)
    return adata


def load_lake_reference(shared_with=None, subclass="subclass.l1"):
    """Load the Lake KPMP reference. If ``shared_with`` (a gene list) is given,
    subset to genes shared with it (e.g. the Xenium panel)."""
    ref = sc.read_h5ad(lake_reference_path())
    if shared_with is not None:
        shared = [g for g in ref.var_names if g in set(shared_with)]
        ref = ref[:, shared].copy()
    return ref
