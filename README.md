# Spatial transcriptomics of kidney: a Xenium glomerulonephritis course

A hands-on course on analysing imaging-based single-cell spatial transcriptomics
(10x **Xenium**) data, worked end to end on a human kidney biopsy slide with several
glomerulonephritis samples. One day is face to face; the notebooks are an end-to-end
reference you can keep working through afterwards.

Robin Khatri, Institute of Medical Systems Bioinformatics, UKE.

## The data

One Xenium slide (ID `0011695`, 480-gene multi-tissue panel, FFPE) carrying **8 tissue
samples**:

| Disease | Samples |
| --- | --- |
| ANCA-GN | X33, X34, X35, X36, X37 |
| anti-GBM (GBM) | X38 |
| Lupus nephritis (SLE) | X39, X40 |

You are given the **raw Xenium output folder** for this slide only. The slide is split
into the 8 samples by the bounding boxes in `data/sample_bounding_boxes.csv` (coordinates
in microns; see `scripts/course_utils.py::assign_samples`). The single-cell **reference**
for classification and integration is the Lake / KPMP kidney atlas (panel-aligned).

## Layout

```
data/sample_bounding_boxes.csv   per-sample bounding boxes (x/y microns) + disease labels
data/marker_genes.csv            kidney cell-type marker lists (panel genes only)
scripts/course_utils.py          read the slide, split samples, load the reference
notebooks/                       worked notebooks with blanks for the exercises
notebooks_solutions/             the same notebooks with exercises filled in (revealed later)
slides/                          lecture slides (PDF + sources + figures)
```

## Notebooks

0. **Setup and samples** read the raw slide, split it into the 8 samples, QC and controls.
1. **Transcriptomic analysis** normalise, cluster, and annotate one sample with markers.
2. **Cell-type classification** marker-based vs reference-based (Lake) label transfer.
3. **Integration** combine the samples, correct batch, joint clustering and composition.
4. **Differential expression and enrichment** ANCA-GN vs the other glomerulonephritis types.
5. **Spatial analysis** neighbourhood enrichment, niches, and ligand-receptor interactions.

Each notebook carries method notes that frame the choices, and interpretation of the
results, so it works as a standalone reference.

## Setup

```bash
conda env create -f environment.yml
conda activate env-xenium-course
```

Point the helpers at your data if your paths differ:

```bash
export XENIUM_RAW_DIR=/path/to/output-XETG00088__0011695__...
export LAKE_REF=/path/to/lake_reduced.h5ad
```
