# Gene Expression Visualizer

A Flask web app for visualizing gene expression data. Enter gene names,
log2 fold changes, and p-values to generate an interactive bar chart,
volcano plot, and fold-change heatmap.

## Features

- Bar chart of fold change per gene, colored by regulation direction
- Volcano plot (log2 fold change vs. -log10 p-value) with configurable
  significance and fold-change thresholds
- Fold-change heatmap with a diverging color scale
- Sortable summary table of all parsed genes
- Input validation with clear, line-specific error messages

## Input format

One gene per line: `gene, fold change, p-value`, separated by commas,
tabs, or spaces. An optional header row is detected automatically.

```
Gene, Log2FC, PValue
TP53, 2.8, 0.0009
BRCA1, -3.1, 0.0002
EGFR, 1.9, 0.012
```

## Setup

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
python app.py
```

Then open http://127.0.0.1:5000 in your browser.
