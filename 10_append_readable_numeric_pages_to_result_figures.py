#!/usr/bin/env python3
from pathlib import Path
import csv
from tempfile import NamedTemporaryFile

from pypdf import PdfReader, PdfWriter
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4


PROJECT_DIR = Path(__file__).resolve().parents[1]
RESULT_DIR = PROJECT_DIR / "result_file"
FIG_DIR = RESULT_DIR / "figures"
REPORT_DIR = PROJECT_DIR / "reports"
REPORT_DIR.mkdir(parents=True, exist_ok=True)


def read_csv_rel(rel_path):
    path = RESULT_DIR / rel_path
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def fnum(x, digits=3):
    if x is None or x == "":
        return "NA"
    try:
        return f"{float(x):.{digits}f}"
    except Exception:
        return str(x)


def one_row(rows, key, value):
    for row in rows:
        if row.get(key) == value:
            return row
    return {}


def lines_perf(rows, title):
    lines = [title]
    for row in rows:
        lines.extend([
            f"{row.get('dataset','')} | AUC={fnum(row.get('auc'))} (95% CI {fnum(row.get('ci_low'))}-{fnum(row.get('ci_high'))}); threshold={fnum(row.get('threshold'))}",
            f"  TP={row.get('tp','')}, FP={row.get('fp','')}, FN={row.get('fn','')}, TN={row.get('tn','')}; Acc={fnum(row.get('accuracy'))}; Sen={fnum(row.get('sensitivity'))}; Spe={fnum(row.get('specificity'))}",
        ])
    return lines


def lines_candidate(rows, title, dataset_filter=None, limit=None):
    rows_use = [r for r in rows if dataset_filter is None or r.get("dataset") == dataset_filter]
    if limit:
        rows_use = rows_use[:limit]
    lines = [title]
    for row in rows_use:
        th = row.get("threshold", "")
        lines.append(
            f"{row.get('method','')} | {row.get('dataset','')} | AUC={fnum(row.get('auc'))} "
            f"(95% CI {fnum(row.get('ci_low'))}-{fnum(row.get('ci_high'))}); threshold={fnum(th)}"
        )
    return lines


def lines_feature_table(rows, title, gene_col="genesymbol", limit=None):
    vals = [r.get(gene_col, "") for r in rows if r.get(gene_col, "")]
    if limit:
        vals = vals[:limit]
    return [title, ", ".join(vals)]


def table_subset_lines(rows, title, columns, limit=10):
    out = [title]
    for row in rows[:limit]:
        parts = [f"{col}={row.get(col, '')}" for col in columns]
        out.append(" | ".join(parts))
    if len(rows) > limit:
        out.append(f"... total rows: {len(rows)}")
    return out


def build_sections():
    qc = read_csv_rel("tables/cohort_and_qc_summary.csv")
    de = read_csv_rel("tables/differential_analysis_summary.csv")
    ma = read_csv_rel("modelA/modelA_final_performance.csv")
    ma_cand = read_csv_rel("modelA/modelA_candidate_model_performance.csv")
    ma_cv = read_csv_rel("modelA/modelA_rf_feature_set_cross_validation.csv")
    ma_feat = read_csv_rel("modelA/modelA_feature_table.csv")
    mb = read_csv_rel("modelB/modelB_final_performance.csv")
    mb_cand = read_csv_rel("modelB/modelB_candidate_model_performance.csv")
    mb_rep = read_csv_rel("modelB/modelB_repeated_split_summary.csv")
    mb_feat = read_csv_rel("modelB/modelB_feature_table.csv")
    mc = read_csv_rel("modelC/modelC_final_performance.csv")
    mc_oof = read_csv_rel("modelC/modelC_discovery_5fold_oof_auc.csv")
    mc_feat = read_csv_rel("modelC/modelC_feature_table.csv")

    qcd = {r["item"]: r["value"] for r in qc}
    ded = {r["comparison"]: r for r in de}

    fig1 = [
        "Cohort and QC key values",
        f"Samples: total={fnum(qcd.get('total_samples'),0)}, Healthy={fnum(qcd.get('healthy'),0)}, Lung cancer={fnum(qcd.get('lung_cancer'),0)}",
        f"NSCLC={fnum(qcd.get('NSCLC'),0)}, SCLC={fnum(qcd.get('SCLC'),0)}, proteins={fnum(qcd.get('proteins'),0)}",
        f"NSCLC N0/N+={fnum(qcd.get('NSCLC_N0'),0)}/{fnum(qcd.get('NSCLC_Nplus'),0)}; SCLC N0/N+={fnum(qcd.get('SCLC_N0'),0)}/{fnum(qcd.get('SCLC_Nplus'),0)}",
        f"Protein detection median={fnum(qcd.get('protein_detection_fraction_median'))}; PCA PC1={fnum(qcd.get('PCA_PC1_percent'))}%, PC2={fnum(qcd.get('PCA_PC2_percent'))}%",
    ]

    fig2 = []
    fig2 += lines_perf(ma, "Model A final five-protein score")
    fig2 += [""]
    fig2 += lines_candidate(ma_cand, "Model A discovery candidate model screen", dataset_filter="Discovery")
    fig2 += [""]
    fig2 += table_subset_lines(ma_cv, "Model A feature-set CV AUC", ["feature_set", "mean_cv_auc", "low_2.5", "high_97.5"], limit=10)
    fig2 += [""]
    fig2 += lines_feature_table(ma_feat, "Model A features")

    fig3 = [
        "NSCLC progression differential analysis",
        f"NSCLC N0 vs Healthy: Up={ded['NSCLC_N0_vs_Healthy']['nominal_positive_up']}, Down={ded['NSCLC_N0_vs_Healthy']['nominal_negative_up']}, FDR significant={ded['NSCLC_N0_vs_Healthy']['fdr_significant']}",
        f"NSCLC N+ vs N0: Up={ded['NSCLC_Nplus_vs_N0']['nominal_positive_up']}, Down={ded['NSCLC_Nplus_vs_N0']['nominal_negative_up']}, FDR significant={ded['NSCLC_Nplus_vs_N0']['fdr_significant']}",
        "Figure3D/F protein-group scores are recalculated from expression matrix using curated display protein groups.",
    ]

    fig4 = [
        "Subtype nodal differential analysis",
        f"LUAD N+ vs N0: Up={ded['LUAD_Nplus_vs_N0']['nominal_positive_up']}, Down={ded['LUAD_Nplus_vs_N0']['nominal_negative_up']}, FDR significant={ded['LUAD_Nplus_vs_N0']['fdr_significant']}",
        f"LUSC N+ vs N0: Up={ded['LUSC_Nplus_vs_N0']['nominal_positive_up']}, Down={ded['LUSC_Nplus_vs_N0']['nominal_negative_up']}, FDR significant={ded['LUSC_Nplus_vs_N0']['fdr_significant']}",
        f"NSCLC N+ vs N0 reference: Up={ded['NSCLC_Nplus_vs_N0']['nominal_positive_up']}, Down={ded['NSCLC_Nplus_vs_N0']['nominal_negative_up']}, FDR significant={ded['NSCLC_Nplus_vs_N0']['fdr_significant']}",
    ]

    fig5 = []
    fig5 += lines_perf(mb, "Model B final RF top18")
    fig5 += [""]
    fig5 += lines_candidate(mb_cand, "Model B discovery candidate model screen", dataset_filter="Discovery")
    rep = mb_rep[0]
    fig5 += ["", f"Repeated split robustness: mean AUC={fnum(rep['auc_mean'])}; 2.5%-97.5%={fnum(rep['auc_low_2.5'])}-{fnum(rep['auc_high_97.5'])}"]
    fig5 += [""] + lines_feature_table(mb_feat, "Model B top18 features")

    fig6 = []
    fig6 += lines_perf(mc, "Model C exploratory equal-direction score")
    fig6 += [""]
    fig6 += lines_candidate(mc_oof, "Model C discovery 5-fold OOF candidate screen")
    fig6 += [""] + lines_feature_table(mc_feat, "Model C top9 features")

    supp1 = [
        "Supplementary Figure 1 QC key values",
        f"Total samples={fnum(qcd.get('total_samples'),0)}; proteins={fnum(qcd.get('proteins'),0)}",
        f"Protein detection median={fnum(qcd.get('protein_detection_fraction_median'))}",
        "Sample-level detection/abundance plotted data are in result_file/figures/plotted_data/Supplementary_Figure_1*.csv",
    ]
    supp2 = [
        "Supplementary Figure 2 lung cancer vs healthy",
        f"LC vs Healthy: Up={ded['LC_vs_Healthy']['nominal_positive_up']}, Down={ded['LC_vs_Healthy']['nominal_negative_up']}, FDR significant={ded['LC_vs_Healthy']['fdr_significant']}",
        "Enrichment universe: all 3058 quantified urine proteins mapped to ENTREZID.",
    ]
    supp3 = [
        "Supplementary Figure 3 external validation and marker consistency",
        f"Zhang validation Model A: AUC={fnum(one_row(ma, 'dataset', 'Zhang validation').get('auc'))} (95% CI {fnum(one_row(ma, 'dataset', 'Zhang validation').get('ci_low'))}-{fnum(one_row(ma, 'dataset', 'Zhang validation').get('ci_high'))})",
        "Marker direction consistency table: result_file/modelA/modelA_marker_direction_consistency_original_vs_Zhang.csv",
    ]
    supp4 = [
        "Supplementary Figure 4 Model B diagnostics",
        f"Model B validation AUC={fnum(one_row(mb, 'dataset', 'Validation holdout').get('auc'))}; threshold={fnum(one_row(mb, 'dataset', 'Validation holdout').get('threshold'))}",
        f"Repeated split mean AUC={fnum(rep['auc_mean'])}; 2.5%-97.5%={fnum(rep['auc_low_2.5'])}-{fnum(rep['auc_high_97.5'])}",
    ]
    supp5 = [
        "Supplementary Figure 5 SCLC exploratory biology",
        f"SCLC N+ vs N0 all exploratory: Up={ded['SCLC_Nplus_vs_N0_all_exploratory']['nominal_positive_up']}, Down={ded['SCLC_Nplus_vs_N0_all_exploratory']['nominal_negative_up']}, FDR significant={ded['SCLC_Nplus_vs_N0_all_exploratory']['fdr_significant']}",
        f"Model C validation AUC={fnum(one_row(mc, 'dataset', 'Validation holdout').get('auc'))}; threshold={fnum(one_row(mc, 'dataset', 'Validation holdout').get('threshold'))}",
    ]

    return {
        "Figure1.pdf": fig1,
        "Figure2.pdf": fig2,
        "Figure3.pdf": fig3,
        "Figure4.pdf": fig4,
        "Figure5.pdf": fig5,
        "Figure6.pdf": fig6,
        "Supplementary_Figure_1.pdf": supp1,
        "Supplementary_Figure_2.pdf": supp2,
        "Supplementary_Figure_3.pdf": supp3,
        "Supplementary_Figure_4.pdf": supp4,
        "Supplementary_Figure_5.pdf": supp5,
    }


def wrap(line, width=106):
    words = str(line).split()
    if not words:
        return [""]
    out, cur = [], words[0]
    for word in words[1:]:
        if len(cur) + 1 + len(word) <= width:
            cur += " " + word
        else:
            out.append(cur)
            cur = word
    out.append(cur)
    return out


def make_pages(pdf_name, lines):
    pages = []
    current = [f"Readable numeric key: {pdf_name}", ""]
    for line in lines:
        for wrapped in wrap(line):
            if len(current) >= 39:
                pages.append(current)
                current = [f"Readable numeric key: {pdf_name} (continued)", ""]
            current.append(wrapped)
    pages.append(current)
    return pages


def create_numeric_page(lines):
    tmp = NamedTemporaryFile(suffix=".pdf", delete=False)
    tmp.close()
    c = canvas.Canvas(tmp.name, pagesize=A4)
    width, height = A4
    y = height - 42
    for i, line in enumerate(lines):
        if i == 0:
            c.setFont("Helvetica-Bold", 14)
        else:
            c.setFont("Courier", 8.8)
        c.drawString(36, y, line)
        y -= 19 if i == 0 else 17
    c.save()
    return Path(tmp.name)


def append_numeric_pages(pdf_name, lines):
    target = FIG_DIR / pdf_name
    if not target.exists():
        return None
    reader = PdfReader(str(target))
    writer = PdfWriter()
    for page in reader.pages:
        writer.add_page(page)
    temp_pages = []
    for page_lines in make_pages(pdf_name, lines):
        numeric_pdf = create_numeric_page(page_lines)
        temp_pages.append(numeric_pdf)
        writer.add_page(PdfReader(str(numeric_pdf)).pages[0])
    with target.open("wb") as fh:
        writer.write(fh)
    for page in temp_pages:
        page.unlink(missing_ok=True)
    return {
        "figure": pdf_name,
        "source_pdf": str(target.relative_to(PROJECT_DIR)),
        "output_pdf": str(target.relative_to(PROJECT_DIR)),
        "numeric_pages_added": len(make_pages(pdf_name, lines)),
    }


def main():
    sections = build_sections()
    rows = []
    for pdf_name, lines in sections.items():
        row = append_numeric_pages(pdf_name, lines)
        if row:
            rows.append(row)
    report = REPORT_DIR / "result_figures_readable_numeric_pages.csv"
    with report.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["figure", "source_pdf", "output_pdf", "numeric_pages_added"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"Readable numeric pages appended to {len(rows)} result figure PDFs.")
    print(f"Report: {report}")


if __name__ == "__main__":
    main()
