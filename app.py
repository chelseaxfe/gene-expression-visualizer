import math
import re

from flask import Flask, jsonify, render_template, request

app = Flask(__name__)

SPLIT_PATTERN = re.compile(r"[,\t]+")


def is_number(value: str) -> bool:
    try:
        float(value)
        return True
    except ValueError:
        return False


def split_row(line: str) -> list[str]:
    parts = [p.strip() for p in SPLIT_PATTERN.split(line) if p.strip()]
    if len(parts) < 3:
        parts = line.split()
    return parts


def parse_input(raw_text: str):
    rows = []
    errors = []

    for i, raw_line in enumerate(raw_text.strip().splitlines(), start=1):
        line = raw_line.strip()
        if not line:
            continue

        parts = split_row(line)

        if i == 1 and len(parts) >= 3 and not is_number(parts[1]) and not is_number(parts[2]):
            continue  # header row

        if len(parts) < 3:
            errors.append(f"Line {i}: expected \"gene, fold change, p-value\", got {len(parts)} value(s).")
            continue

        gene, fc_raw, p_raw = parts[0], parts[1], parts[2]

        if not gene:
            errors.append(f"Line {i}: missing gene name.")
            continue

        try:
            log2fc = float(fc_raw)
        except ValueError:
            errors.append(f"Line {i}: invalid fold change value \"{fc_raw}\".")
            continue

        try:
            pvalue = float(p_raw)
        except ValueError:
            errors.append(f"Line {i}: invalid p-value \"{p_raw}\".")
            continue

        if not (0 <= pvalue <= 1):
            errors.append(f"Line {i}: p-value must be between 0 and 1, got {pvalue}.")
            continue

        rows.append({"gene": gene, "log2fc": log2fc, "pvalue": pvalue})

    return rows, errors


def classify(log2fc: float, pvalue: float, fc_threshold: float, p_threshold: float) -> str:
    if pvalue < p_threshold and log2fc >= fc_threshold:
        return "up"
    if pvalue < p_threshold and log2fc <= -fc_threshold:
        return "down"
    return "not_significant"


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/api/analyze", methods=["POST"])
def analyze():
    payload = request.get_json(silent=True) or {}
    raw_data = payload.get("data", "")

    try:
        fc_threshold = float(payload.get("fc_threshold", 1.0))
        p_threshold = float(payload.get("p_threshold", 0.05))
    except (TypeError, ValueError):
        return jsonify({"error": "Thresholds must be numeric."}), 400

    if fc_threshold < 0 or not (0 < p_threshold <= 1):
        return jsonify({"error": "Fold-change threshold must be ≥ 0 and p-value threshold must be between 0 and 1."}), 400

    rows, errors = parse_input(raw_data)

    if errors:
        preview = errors[:6]
        message = "\n".join(preview)
        if len(errors) > len(preview):
            message += f"\n...and {len(errors) - len(preview)} more error(s)."
        return jsonify({"error": message}), 400

    if not rows:
        return jsonify({"error": "Please enter at least one gene (name, fold change, p-value)."}), 400

    genes = []
    counts = {"up": 0, "down": 0, "not_significant": 0}

    for row in rows:
        neg_log10_p = -math.log10(max(row["pvalue"], 1e-300))
        regulation = classify(row["log2fc"], row["pvalue"], fc_threshold, p_threshold)
        counts[regulation] += 1
        genes.append({
            "gene": row["gene"],
            "log2fc": round(row["log2fc"], 4),
            "pvalue": row["pvalue"],
            "neg_log10_p": round(neg_log10_p, 4),
            "regulation": regulation,
        })

    genes.sort(key=lambda g: g["log2fc"], reverse=True)

    return jsonify({
        "genes": genes,
        "summary": {
            "total": len(genes),
            "upregulated": counts["up"],
            "downregulated": counts["down"],
            "not_significant": counts["not_significant"],
        },
        "thresholds": {"fc_threshold": fc_threshold, "p_threshold": p_threshold},
    })


if __name__ == "__main__":
    app.run(debug=True)
