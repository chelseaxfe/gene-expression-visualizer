const dataInput = document.getElementById("dataInput");
const rowCount = document.getElementById("rowCount");
const fcThresholdInput = document.getElementById("fcThreshold");
const pThresholdInput = document.getElementById("pThreshold");
const visualizeBtn = document.getElementById("visualizeBtn");
const exampleBtn = document.getElementById("exampleBtn");
const clearBtn = document.getElementById("clearBtn");
const errorMessage = document.getElementById("errorMessage");
const results = document.getElementById("results");

const EXAMPLE_DATA = `Gene, Log2FC, PValue
TP53, 2.8, 0.0009
BRCA1, -3.1, 0.0002
EGFR, 1.9, 0.012
MYC, 3.4, 0.00005
PTEN, -2.2, 0.003
KRAS, 0.4, 0.45
CDKN2A, -1.8, 0.021
VEGFA, 2.1, 0.008
AKT1, 0.6, 0.31
STAT3, -0.9, 0.19
CCND1, 1.2, 0.041
BCL2, -2.6, 0.0015
NOTCH1, 0.3, 0.62
JUN, 1.5, 0.033
FOS, -1.1, 0.09`;

const REGULATION_COLORS = {
  up: "#34d399",
  down: "#f87171",
  not_significant: "#64748b",
};

const REGULATION_LABELS = {
  up: "Up",
  down: "Down",
  not_significant: "Not Significant",
};

const PLOTLY_LAYOUT_BASE = {
  paper_bgcolor: "transparent",
  plot_bgcolor: "transparent",
  font: { color: "#e8ecf5", family: "Inter, sans-serif" },
  margin: { t: 20, r: 20, b: 60, l: 60 },
  xaxis: { gridcolor: "#253049", zerolinecolor: "#37476b" },
  yaxis: { gridcolor: "#253049", zerolinecolor: "#37476b" },
};

const PLOTLY_CONFIG = { displaylogo: false, responsive: false };

const MIN_PX_PER_BAR = 40;
const MIN_PX_PER_ROW = 28;

// Sizes a plot's inner div in pixels against its (fixed-size, scrollable)
// viewport, growing beyond the viewport when there are many data points so
// the viewport scrolls instead of the plot overflowing onto other sections.
function sizePlotElement(plotId, { width, height }) {
  const el = document.getElementById(plotId);
  const viewport = el.parentElement;
  const targetWidth = Math.max(viewport.clientWidth, width || 0);
  const targetHeight = Math.max(viewport.clientHeight, height || 0);
  el.style.width = `${targetWidth}px`;
  el.style.height = `${targetHeight}px`;
  return { width: targetWidth, height: targetHeight };
}

function countGenes() {
  const lines = dataInput.value
    .split("\n")
    .map((l) => l.trim())
    .filter(Boolean);
  return lines.length;
}

function updateRowCount() {
  const n = countGenes();
  rowCount.textContent = `${n} line${n === 1 ? "" : "s"}`;
}

function showError(message) {
  errorMessage.textContent = message;
  errorMessage.hidden = false;
  results.hidden = true;
}

function hideError() {
  errorMessage.hidden = true;
}

function renderBarChart(genes) {
  const trace = {
    type: "bar",
    x: genes.map((g) => g.gene),
    y: genes.map((g) => g.log2fc),
    marker: { color: genes.map((g) => REGULATION_COLORS[g.regulation]) },
    hovertemplate: "<b>%{x}</b><br>Log2FC: %{y}<extra></extra>",
  };

  const { width, height } = sizePlotElement("barChart", { width: genes.length * MIN_PX_PER_BAR });

  const layout = {
    ...PLOTLY_LAYOUT_BASE,
    width,
    height,
    autosize: false,
    xaxis: { ...PLOTLY_LAYOUT_BASE.xaxis, title: "Gene", tickangle: -45 },
    yaxis: { ...PLOTLY_LAYOUT_BASE.yaxis, title: "Log2 Fold Change" },
  };

  Plotly.newPlot("barChart", [trace], layout, PLOTLY_CONFIG);
}

function renderVolcanoPlot(genes, fcThreshold, pThreshold) {
  const negLogPThreshold = -Math.log10(pThreshold);
  const groups = ["up", "down", "not_significant"];

  const traces = groups.map((group) => {
    const groupGenes = genes.filter((g) => g.regulation === group);
    return {
      type: "scatter",
      mode: "markers",
      name: REGULATION_LABELS[group],
      x: groupGenes.map((g) => g.log2fc),
      y: groupGenes.map((g) => g.neg_log10_p),
      text: groupGenes.map((g) => g.gene),
      marker: { color: REGULATION_COLORS[group], size: 10, opacity: 0.85 },
      hovertemplate: "<b>%{text}</b><br>Log2FC: %{x}<br>-Log10(p): %{y}<extra></extra>",
    };
  });

  const maxAbsFc = Math.max(1, ...genes.map((g) => Math.abs(g.log2fc))) * 1.15;
  const maxNegLogP = Math.max(negLogPThreshold, ...genes.map((g) => g.neg_log10_p)) * 1.15;
  const { width, height } = sizePlotElement("volcanoChart", {});

  const layout = {
    ...PLOTLY_LAYOUT_BASE,
    width,
    height,
    autosize: false,
    xaxis: { ...PLOTLY_LAYOUT_BASE.xaxis, title: "Log2 Fold Change", range: [-maxAbsFc, maxAbsFc] },
    yaxis: { ...PLOTLY_LAYOUT_BASE.yaxis, title: "-Log10(p-value)", range: [0, maxNegLogP] },
    legend: { orientation: "h", y: -0.25 },
    shapes: [
      {
        type: "line",
        x0: fcThreshold,
        x1: fcThreshold,
        y0: 0,
        y1: maxNegLogP,
        line: { color: "#37476b", width: 1, dash: "dash" },
      },
      {
        type: "line",
        x0: -fcThreshold,
        x1: -fcThreshold,
        y0: 0,
        y1: maxNegLogP,
        line: { color: "#37476b", width: 1, dash: "dash" },
      },
      {
        type: "line",
        x0: -maxAbsFc,
        x1: maxAbsFc,
        y0: negLogPThreshold,
        y1: negLogPThreshold,
        line: { color: "#37476b", width: 1, dash: "dash" },
      },
    ],
  };

  Plotly.newPlot("volcanoChart", traces, layout, PLOTLY_CONFIG);
}

function renderHeatmap(genes) {
  const sorted = [...genes].sort((a, b) => b.log2fc - a.log2fc);
  const maxAbs = Math.max(1, ...sorted.map((g) => Math.abs(g.log2fc)));

  const trace = {
    type: "heatmap",
    x: ["Log2FC"],
    y: sorted.map((g) => g.gene),
    z: sorted.map((g) => [g.log2fc]),
    zmin: -maxAbs,
    zmax: maxAbs,
    colorscale: [
      [0, "#38bdf8"],
      [0.5, "#141b2d"],
      [1, "#f87171"],
    ],
    colorbar: { title: "Log2FC", titlefont: { color: "#e8ecf5" }, tickfont: { color: "#8b96ad" } },
    hovertemplate: "<b>%{y}</b><br>Log2FC: %{z}<extra></extra>",
  };

  const { width, height } = sizePlotElement("heatmapChart", { height: sorted.length * MIN_PX_PER_ROW });

  const layout = {
    ...PLOTLY_LAYOUT_BASE,
    width,
    height,
    autosize: false,
    margin: { t: 20, r: 20, b: 60, l: 100 },
    xaxis: { ...PLOTLY_LAYOUT_BASE.xaxis, title: "" },
    yaxis: { ...PLOTLY_LAYOUT_BASE.yaxis, title: "", autorange: "reversed" },
  };

  Plotly.newPlot("heatmapChart", [trace], layout, PLOTLY_CONFIG);
}

function renderTable(genes) {
  const tbody = document.querySelector("#geneTable tbody");
  tbody.innerHTML = "";

  genes.forEach((g) => {
    const row = document.createElement("tr");
    row.innerHTML = `
      <td>${g.gene}</td>
      <td>${g.log2fc}</td>
      <td>${g.pvalue}</td>
      <td>${g.neg_log10_p}</td>
      <td><span class="badge ${g.regulation}">${REGULATION_LABELS[g.regulation]}</span></td>
    `;
    tbody.appendChild(row);
  });
}

function renderResults(data) {
  document.getElementById("statTotal").textContent = data.summary.total;
  document.getElementById("statUp").textContent = data.summary.upregulated;
  document.getElementById("statDown").textContent = data.summary.downregulated;
  document.getElementById("statNs").textContent = data.summary.not_significant;

  // Charts are sized from their viewport's on-screen dimensions, so the
  // section must be visible (and laid out) before Plotly measures it.
  results.hidden = false;

  renderBarChart(data.genes);
  renderVolcanoPlot(data.genes, data.thresholds.fc_threshold, data.thresholds.p_threshold);
  renderHeatmap(data.genes);
  renderTable(data.genes);
}

async function visualize() {
  visualizeBtn.disabled = true;
  visualizeBtn.textContent = "Visualizing...";

  try {
    const response = await fetch("/api/analyze", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        data: dataInput.value,
        fc_threshold: fcThresholdInput.value,
        p_threshold: pThresholdInput.value,
      }),
    });
    const data = await response.json();

    if (!response.ok) {
      showError(data.error || "Something went wrong. Please try again.");
      return;
    }

    hideError();
    renderResults(data);
  } catch (err) {
    showError("Could not reach the server. Please try again.");
  } finally {
    visualizeBtn.disabled = false;
    visualizeBtn.textContent = "Visualize";
  }
}

dataInput.addEventListener("input", updateRowCount);

visualizeBtn.addEventListener("click", visualize);

exampleBtn.addEventListener("click", () => {
  dataInput.value = EXAMPLE_DATA;
  updateRowCount();
  visualize();
});

clearBtn.addEventListener("click", () => {
  dataInput.value = "";
  updateRowCount();
  hideError();
  results.hidden = true;
});

updateRowCount();
