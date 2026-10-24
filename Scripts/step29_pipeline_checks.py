"""step29_pipeline_checks.py — consistency checks run at the end of the
pipeline. Any FAIL stops `run_pipeline.sh` with exit code 1.

Everything downstream rests on the units of Steps 01 – 05 (populations,
demes, BLs). These checks verify that every output still agrees with
them, so a stale table or a legacy helper cannot slip through:

 1. Units      — the event → deme table rebuilds identically from the
                 crosswalk and matches Step 01's deme table.
 2. BL order   — BL names follow Step 05's rule (BL{rank}, area DESC →
                 connectivity DESC) and `srk_bl_constants.BL_ORDER` is
                 that order.
 3. BL labels  — every table with populationID + BL (or locationCode +
                 BL) agrees with the current BLs.
 4. IDs        — every populationID / (population, year, deme) in the
                 outputs exists in the units.
 5. Targets    — Step 17's deme targets equal `srk_detection`'s rule
                 (the prediction model), for both years.
 6. Year       — the single-year sampling tables hold 2025 (run last).
 7. Code       — no script uses the legacy EO-level BL lookup except
                 Step 01's placeholder; no hard-coded BL order lists.
 8. Tables     — every file in tables/ is written or read by a script
                 (or is a listed external input); tracked tables carry
                 no coordinate columns.
"""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pandas as pd

from srk_bl_constants import BL_ORDER, locationCode_to_bl
from srk_detection import mothers_for_coverage
from srk_units import load_event_units, sampling_year

TABLES = Path("tables")
SCRIPTS = Path("Scripts")
# Copied from SRK_bioinformatics at the repo split; kept for reference,
# not produced here. Only step26i_L1_carrier_inventory is an input.
EXTERNAL_PREFIXES = ("step26", "step27_Sdomain", "SRK_allele_")
COORD_COL = re.compile(r"(^|_)(lat|lon|latitude|longitude)($|_)|Latitude|Longitude",
                       re.IGNORECASE)

results: list[tuple[str, str, str]] = []


def check(name: str, ok: bool, detail: str = "", warn: bool = False) -> None:
    results.append((("PASS" if ok else ("WARN" if warn else "FAIL")), name, detail))


def read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", encoding="utf-8-sig", low_memory=False)


def main() -> None:
    # 1. Units
    units = load_event_units(validate=True)  # raises on mismatch with Step 01
    saved = read(TABLES / "step29a_event_to_deme.tsv")
    cols = ["event_year", "eventID", "populationID", "BL", "demeID"]
    same = (units[cols].sort_values(cols[:2]).reset_index(drop=True)
            .equals(saved[cols].sort_values(cols[:2]).reset_index(drop=True)))
    check("units: event → deme table matches a fresh rebuild", same)
    pop_bl = units.drop_duplicates("populationID").set_index("populationID")["BL"]
    demes = {(int(y), int(p), int(d)) for y, p, d in
             units[["event_year", "populationID", "demeID"]].drop_duplicates()
             .itertuples(index=False)}

    # 2. BL order
    bl_def = read(TABLES / "step30h_bl_definition.tsv").sort_values("BL_rank")
    names_ok = (bl_def["BL"] == [f"BL{r}" for r in bl_def["BL_rank"]]).all()
    rule_ok = bl_def["convex_hull_area_km2"].is_monotonic_decreasing
    check("BL order: names are BL{rank} (Step 05)", bool(names_ok))
    check("BL order: rank follows area DESC", bool(rule_ok),
          "ties broken by connectivity" if not rule_ok else "")
    check("BL order: BL_ORDER == Step 05 order",
          BL_ORDER == bl_def["BL"].tolist(), str(BL_ORDER))
    check("BL labels: every population's BL is a Step 05 BL",
          set(pop_bl) <= set(bl_def["BL"]))

    # 3 + 4. Labels and IDs in every table
    for f in sorted(TABLES.glob("*.tsv")):
        if f.name.startswith(EXTERNAL_PREFIXES):
            continue
        try:
            df = read(f)
        except Exception as e:  # noqa: BLE001
            check(f"read {f.name}", False, str(e))
            continue
        if "populationID" in df.columns and pd.api.types.is_numeric_dtype(df["populationID"]):
            ids = set(df["populationID"].dropna().astype(int))
            unknown = ids - set(pop_bl.index)
            check(f"IDs: {f.name} populations exist", not unknown,
                  f"unknown {sorted(unknown)[:8]}")
            if "BL" in df.columns:
                m = df.dropna(subset=["populationID", "BL"])
                m = m[m["populationID"].astype(int).isin(pop_bl.index)]
                bad = m[m["BL"] != m["populationID"].astype(int).map(pop_bl)]
                check(f"BL labels: {f.name}", bad.empty,
                      f"{len(bad)} rows, e.g. P{bad['populationID'].iloc[0]} "
                      f"{bad['BL'].iloc[0]}" if len(bad) else "")
            ycol = next((c for c in ("year", "event_year") if c in df.columns), None)
            if "demeID" in df.columns and ycol:
                keys = {(int(y), int(p), int(d)) for y, p, d in
                        df[[ycol, "populationID", "demeID"]].dropna()
                        .drop_duplicates().itertuples(index=False)}
                check(f"IDs: {f.name} demes exist", keys <= demes,
                      f"{len(keys - demes)} unknown")
        elif {"locationCode", "BL"} <= set(df.columns):
            cur = locationCode_to_bl(df["locationCode"]).to_numpy()
            bad = (df["BL"].to_numpy() != cur) & pd.notna(cur)
            check(f"BL labels: {f.name} (by location)", not bad.any(),
                  f"{int(bad.sum())} rows")

    # 5. Deme targets = srk_detection rule
    for yr in (2025, 2026):
        ev = read(TABLES / f"step29c_sampling_frag_aware_per_event_{yr}.tsv")
        d = ev.drop_duplicates(["populationID", "demeID"])
        exp = d["deme_N_fertile"].map(lambda n: mothers_for_coverage(int(n)) if n > 0 else 1)
        bad = d[d["M_paternal_only"].to_numpy() != exp.to_numpy()]
        check(f"targets {yr}: deme M equals srk_detection rule", bad.empty,
              f"{len(bad)} demes")

    # 6. Sampling year
    check("year: single-year sampling tables hold 2025", sampling_year() == 2025,
          str(sampling_year()))

    # 7. Code
    src = {p.name: p.read_text() for p in SCRIPTS.glob("*.py")}
    legacy = [n for n, t in src.items()
              if "legacy_locationCode_to_bl(" in t
              and n not in ("srk_bl_constants.py", "step01_populations.py",
                            "step29_pipeline_checks.py")]
    check("code: legacy EO-level BL lookup used only by Step 01", not legacy, str(legacy))
    hard = [n for n, t in src.items()
            if re.search(r'\[\s*"BL\d"\s*,\s*"BL\d"\s*,\s*"BL\d"', t)
            and n not in ("srk_bl_constants.py",)]
    check("code: no hard-coded BL order lists", not hard, str(hard))

    # 8. Tables
    text = "\n".join(src.values())
    orphans = []
    for f in sorted(TABLES.glob("*.*")):
        stem = re.sub(r"_(2025|2026)$", "", f.stem)
        # Step 22 --demo writes step30_B{tag}_<name> with tag = "_DEMO"
        demo_part = stem.split("_DEMO_", 1)[1] if "_DEMO_" in stem else None
        if (f.name.startswith(EXTERNAL_PREFIXES) or stem in text
                or f.stem in text or (demo_part and f"{{tag}}_{demo_part}" in text)):
            continue
        orphans.append(f.name)
    check("tables: every table is written or read by a script", not orphans,
          ", ".join(orphans))
    tracked = subprocess.run(["git", "ls-files", "tables"], capture_output=True,
                             text=True).stdout.split()
    with_coords = []
    for t in tracked:
        p = Path(t)
        if p.suffix not in (".tsv", ".csv") or not p.exists():
            continue
        sep = "\t" if p.suffix == ".tsv" else ","
        cols = pd.read_csv(p, sep=sep, nrows=0, encoding="utf-8-sig").columns
        if any(COORD_COL.search(c) for c in cols):
            with_coords.append(p.name)
    check("tables: no coordinate columns in git-tracked tables", not with_coords,
          ", ".join(with_coords))

    # Documents: every relative link points to a file in the repository, so
    # nothing renders broken on GitHub (local-only figures must not be linked)
    all_tracked = set(subprocess.run(["git", "ls-files"], capture_output=True,
                                     text=True).stdout.split())
    broken = []
    for md in [p for p in all_tracked if p.endswith(".md")]:
        base = Path(md).parent
        for m in re.finditer(r"\]\(([^)#\s]+)\)", Path(md).read_text(encoding="utf-8")):
            link = m.group(1)
            if link.startswith(("http", "mailto:")):
                continue
            q = str(Path(base, link)).replace("\\", "/")
            q = str(Path(q).resolve().relative_to(Path.cwd().resolve())) if Path(q).exists() else q
            if q not in all_tracked and not any(x.startswith(q.rstrip("/") + "/") for x in all_tracked):
                broken.append(f"{md} → {link}")
    live = [b for b in broken if not b.startswith("archive/")]
    old = [b for b in broken if b.startswith("archive/")]
    check("docs: every link points to a file in the repository", not live, "; ".join(live))
    check("docs: archived documents (not maintained) — dead links", not old,
          f"{len(old)} dead links in archive/", warn=True)

    n_fail = sum(r[0] == "FAIL" for r in results)
    for status, name, detail in results:
        if status != "PASS" or name.startswith(("units", "BL order", "targets", "year", "code", "tables")):
            print(f"[step29] {status:4}  {name}" + (f" — {detail}" if detail and status != "PASS" else ""))
    print(f"[step29] {len(results)} checks: "
          f"{sum(r[0] == 'PASS' for r in results)} pass, "
          f"{sum(r[0] == 'WARN' for r in results)} warn, {n_fail} fail")
    if n_fail:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
