"""print_figures.py — every figure embedded in Seed_bulking.md, in document
order, ready to print: print_figures/FigNN_<name>.pdf (vector PDF where one
exists, else the PNG/JPG) plus print_figures/Seed_bulking_figures.pdf, all of
them merged, one figure per page. The figures are read from the markdown
itself, so the folder follows the document. Rebuilt from scratch on every run;
gitignored (local copies of figures_seed_bulking/). Run from the repository root:
    python Scripts/print_figures.py
"""
from __future__ import annotations

import re
import shutil
from pathlib import Path

from PIL import Image
from pypdf import PdfWriter

DOC = Path("Seed_bulking.md")
OUT = Path("print_figures")


def main() -> None:
    figs = re.findall(r"!\[Figure (\d+)\]\(([^)]+)\)", DOC.read_text(encoding="utf-8"))
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir()
    merged = PdfWriter()
    for n, src in figs:
        src = Path(src)
        pdf = src.with_suffix(".pdf")
        dst = OUT / f"Fig{int(n):02d}_{src.stem}.pdf"
        if pdf.exists():
            shutil.copy(pdf, dst)
        else:                                            # raster only: one page at 200 dpi
            Image.open(src).convert("RGB").save(dst, resolution=200)
        merged.append(str(dst))
    merged.write(OUT / "Seed_bulking_figures.pdf")
    print(f"[print_figures] {len(figs)} figures → {OUT}/ (+ Seed_bulking_figures.pdf)")


if __name__ == "__main__":
    main()
