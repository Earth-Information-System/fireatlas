"""
Cross-branch FEDS output comparison script.

Loads perimeter snapshot FGB files from S3 for three branches (prod, staging, dev),
computes WKT-based geometric matches, prints a summary, and saves folium HTML maps
and a text report to S3 under E2E_test_outputs/{timestamp}/.

Usage (called from compare_branches.sh):
    python3 compare_outputs.py \
        --prod-nrt-regnm=CONUS_TEST_... \
        --staging-nrt-regnm=CONUS_TEST_... \
        --dev-nrt-regnm=CONUS_TEST_... \
        --prod-archive-regnm=CONUS_TEST_... \
        --staging-archive-regnm=CONUS_TEST_... \
        --dev-archive-regnm=CONUS_TEST_... \
        --nrt-date-string=20260417AM \
        --archive-date-string=20240417AM \
        --timestamp=20260501120000
"""

import argparse
import os
import subprocess
import tempfile
import warnings

warnings.filterwarnings("ignore", category=FutureWarning, module="geopandas")

import folium
import geopandas as gpd
import s3fs

BRANCH_BASES = {
    "prod":    "s3://maap-ops-workspace/shared/gsfc_landslides/FEDSoutput-v3",
    "staging": "s3://maap-ops-workspace/shared/zbecker/FEDSstaging/FEDSoutput-v3",
    "dev":     "s3://maap-ops-workspace/shared/gsfc_landslides/FEDS-staging/FEDSoutput-v3",
}

OUTPUT_BASE = "s3://maap-ops-workspace/shared/gsfc_landslides/FEDSoutput-v3/E2E_test_outputs"


def find_snapshot_file(s3: s3fs.S3FileSystem, base: str, regnm: str, date_string: str):
    year = date_string[:4]
    pattern = f"{base.replace('s3://', '')}/{regnm}/{year}/Snapshot/*/*perimeter.fgb"
    candidates = s3.glob(pattern)
    matches = [f for f in candidates if date_string in f]
    if not matches:
        return None
    return "s3://" + matches[0]


def load_perimeters(filepath: str) -> gpd.GeoDataFrame:
    gdf = gpd.read_file(filepath)
    gdf = gdf.sort_values("farea", ascending=False).reset_index(drop=True)
    gdf["wkt"] = gdf.geometry.to_wkt()
    return gdf


def save_to_s3(local_path: str, s3_path: str) -> bool:
    subprocess.run(["aws", "s3", "cp", local_path, s3_path], check=True)
    s3 = s3fs.S3FileSystem()
    return s3.exists(s3_path.replace("s3://", ""))


def compare_mode(
    s3: s3fs.S3FileSystem,
    mode: str,
    prod_regnm: str,
    staging_regnm: str,
    dev_regnm: str,
    date_string: str,
    output_dir: str,
    report_lines: list,
) -> None:
    header = f"\n[{mode.upper()} mode — {date_string}]\n  perimeter:"
    print(header)
    report_lines.append(header)

    def log(msg):
        print(msg)
        report_lines.append(msg)

    prod_file = find_snapshot_file(s3, BRANCH_BASES["prod"], prod_regnm, date_string)
    staging_file = find_snapshot_file(s3, BRANCH_BASES["staging"], staging_regnm, date_string)
    dev_file = find_snapshot_file(s3, BRANCH_BASES["dev"], dev_regnm, date_string)

    regnms = {"prod": prod_regnm, "staging": staging_regnm, "dev": dev_regnm}
    for label, path in [("prod", prod_file), ("staging", staging_file), ("dev", dev_file)]:
        if path is None:
            log(f"    ERROR: No snapshot file found for {label} ({regnms[label]})")
            return

    prod_df = load_perimeters(prod_file)
    staging_df = load_perimeters(staging_file)
    dev_df = load_perimeters(dev_file)

    prod_wkts = set(prod_df["wkt"])
    staging_wkts = set(staging_df["wkt"])
    dev_wkts = set(dev_df["wkt"])

    three_way = prod_wkts & staging_wkts & dev_wkts
    prod_staging = prod_wkts & staging_wkts
    prod_dev = prod_wkts & dev_wkts

    n_prod = len(prod_df)
    n_staging = len(staging_df)
    n_dev = len(dev_df)
    n_three = len(three_way)

    pct_staging = 100.0 * len(prod_staging) / n_prod if n_prod else 0.0
    pct_dev = 100.0 * len(prod_dev) / n_prod if n_prod else 0.0

    prod_unmatched = n_prod - n_three
    staging_unmatched = n_staging - n_three
    dev_unmatched = n_dev - n_three

    log(f"    prod: {n_prod} features | staging: {n_staging} | dev: {n_dev}")
    log(f"    3-way matches: {n_three}")
    log(f"    prod→staging match: {pct_staging:.1f}% ({len(prod_staging)}/{n_prod})")
    log(f"    prod→dev match:     {pct_dev:.1f}% ({len(prod_dev)}/{n_prod})")
    log(f"    prod unmatched: {prod_unmatched} | staging unmatched: {staging_unmatched} | dev unmatched: {dev_unmatched}")

    prod_nm = prod_df[~prod_df["wkt"].isin(three_way)]
    staging_nm = staging_df[~staging_df["wkt"].isin(three_way)]
    dev_nm = dev_df[~dev_df["wkt"].isin(three_way)]

    map_cols = [c for c in ["geometry", "t", "farea", "meanFRP", "n_pixels"] if c in prod_df.columns]

    def build_map(layers):
        m = None
        for df, name, color in layers:
            if df.empty:
                continue
            if m is None:
                m = df[map_cols].explore(name=name, color=color)
            else:
                df[map_cols].explore(m=m, name=name, color=color)
        if m is not None:
            folium.LayerControl().add_to(m)
        return m

    def save_map(m, filename):
        s3_path = f"{output_dir}/{filename}"
        with tempfile.NamedTemporaryFile(suffix=".html", delete=False, mode="w", encoding="utf-8") as tmp:
            tmp_path = tmp.name
            m.save(tmp_path)
        try:
            exists = save_to_s3(tmp_path, s3_path)
            status = "saved" if exists else "UPLOAD FAILED"
            log(f"    → {status}: {s3_path}")
        finally:
            os.unlink(tmp_path)

    # Matches map
    prod_match = prod_df[prod_df["wkt"].isin(three_way)]
    staging_match = staging_df[staging_df["wkt"].isin(three_way)]
    dev_match = dev_df[dev_df["wkt"].isin(three_way)]

    matches_map = build_map([
        (prod_match, "Prod (conus-dps)", "blue"),
        (staging_match, "Staging", "red"),
        (dev_match, "Dev", "purple"),
    ])
    if matches_map is not None:
        save_map(matches_map, f"{mode}_perimeter_matches.html")

    # Non-matches map
    if prod_nm.empty and staging_nm.empty and dev_nm.empty:
        log("    All features match — no non-matches map generated.")
        return

    nm_map = build_map([
        (prod_nm, "Prod (conus-dps)", "blue"),
        (staging_nm, "Staging", "red"),
        (dev_nm, "Dev", "purple"),
    ])
    if nm_map is not None:
        save_map(nm_map, f"{mode}_perimeter_nonmatches.html")


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare FEDS perimeter outputs across branches.")
    parser.add_argument("--prod-nrt-regnm", required=True)
    parser.add_argument("--staging-nrt-regnm", required=True)
    parser.add_argument("--dev-nrt-regnm", required=True)
    parser.add_argument("--prod-archive-regnm", required=True)
    parser.add_argument("--staging-archive-regnm", required=True)
    parser.add_argument("--dev-archive-regnm", required=True)
    parser.add_argument("--nrt-date-string", required=True, help="e.g. 20260417AM")
    parser.add_argument("--archive-date-string", required=True, help="e.g. 20240417AM")
    parser.add_argument("--timestamp", required=True, help="Run timestamp (YYYYMMDDHHMMSS)")
    parser.add_argument("--output-dir", default=OUTPUT_BASE)
    args = parser.parse_args()

    s3 = s3fs.S3FileSystem()
    output_dir = f"{args.output_dir.rstrip('/')}/{args.timestamp}"
    report_lines = ["=== Cross-Branch Comparison ===", f"Timestamp: {args.timestamp}"]

    print("\n=== Cross-Branch Comparison ===")

    compare_mode(
        s3=s3,
        mode="nrt",
        prod_regnm=args.prod_nrt_regnm,
        staging_regnm=args.staging_nrt_regnm,
        dev_regnm=args.dev_nrt_regnm,
        date_string=args.nrt_date_string,
        output_dir=output_dir,
        report_lines=report_lines,
    )

    compare_mode(
        s3=s3,
        mode="archive",
        prod_regnm=args.prod_archive_regnm,
        staging_regnm=args.staging_archive_regnm,
        dev_regnm=args.dev_archive_regnm,
        date_string=args.archive_date_string,
        output_dir=output_dir,
        report_lines=report_lines,
    )

    # Save report.txt to S3
    report_lines.append("\n=== Generated HTML Maps ===")
    for mode in ("nrt", "archive"):
        for suffix in ("matches", "nonmatches"):
            path = f"{output_dir}/{mode}_perimeter_{suffix}.html"
            exists = s3.exists(path.replace("s3://", ""))
            status = "EXISTS" if exists else "NOT FOUND"
            line = f"  [{status}] {path}"
            print(line)
            report_lines.append(line)

    report_s3_path = f"{output_dir}/report.txt"
    with tempfile.NamedTemporaryFile(suffix=".txt", delete=False, mode="w", encoding="utf-8") as tmp:
        tmp_path = tmp.name
        tmp.write("\n".join(report_lines) + "\n")
    try:
        save_to_s3(tmp_path, report_s3_path)
        print(f"\n  Report saved: {report_s3_path}")
    finally:
        os.unlink(tmp_path)


if __name__ == "__main__":
    main()
