"""talk_figures.py — copy every figure of the seed-bulking talk into
talk_figures/, numbered in presentation order (zero-padded prefix), so the
slides can be built by sorting the folder.

Order follows Seed_bulking.md (Figures 1 – 14, plus a field-board photo before
Figure 8); where a figure has slide or
build-up versions, those come instead of (or before) the document version:
Figure 1 is split into three slides (the challenge alone, panels 1 – 2,
panels 3 – 4, each pair built panel by panel); Figure 3 builds up step by step;
Figure 2 opens on the three types (then what they mean for restoration)
and is split into two landscape slides (objectives 1 – 2, then 3a – 3b); Figure 8 into two 16:9 slides, Figures 9 and 10
start from an empty frame, and Figure 11 is followed by the BL1 close-up
(kept local, like every per-BL map). File names: NN_FigXX[-part]_name.png.

White borders are trimmed to a slim edge (EDGE pixels), for a colored slide
background; frames of one build-up share the same crop, so that nothing jumps
between clicks. The folder is rebuilt from scratch on every run and is gitignored
(cropped copies of figures_seed_bulking/, which stay as they are). Run from the
repository root:
    python Scripts/talk_figures.py
"""
from __future__ import annotations

import re
import shutil
from pathlib import Path

import numpy as np
from PIL import Image

SRC = Path("figures_seed_bulking")
OUT = Path("talk_figures")
EDGE = 30                                                # white edge kept, pixels (at 200 dpi)

# (figure number, part label, source file stem) in presentation order
TALK = [
    (1, "slide1", "breeding_system_cartoon_slide1"),
    (1, "slide2a", "breeding_system_cartoon_slide2_step1"),
    (1, "slide2b", "breeding_system_cartoon_slide2"),
    (1, "slide3a", "breeding_system_cartoon_slide3_step1"),
    (1, "slide3b", "breeding_system_cartoon_slide3"),
    (2, "typesa", "restoration_program_cartoon_types"),
    (2, "typesb", "restoration_program_cartoon_types_use"),
    (2, "slide1", "restoration_program_cartoon_slide1"),
    (2, "slide2a", "restoration_program_cartoon_slide2_step1"),
    (2, "slide2b", "restoration_program_cartoon_slide2"),
    (3, "a", "source_screening_cartoon_step1"),
    (3, "b", "source_screening_cartoon_step2"),
    (3, "c", "source_screening_cartoon"),
    (4, "", "EO_sampling_map_black"),
    (5, "", "BL_fragmentation_network_map"),
    (6, "", "units_cartoon_landscape"),
    (7, "", "population_EO_map_BL"),
    (7, "BL1BL3", "BL1_BL3_populations_EO_map"),   # EO27 (BL1), EO24 (BL3); local only
    (8, "photo", "plant_photos/size_medium_OC2592_crown13cm"),   # field board, whole photo
    (8, "slide1", "plant_size_examples_slide1"),
    (8, "slide2", "plant_size_examples_slide2"),
    (9, "a", "seed_by_size_points"),
    (9, "b", "seed_by_size_plants"),
    (10, "a", "population_seed_yield_frame"),
    (10, "b", "population_seed_yield"),
    (11, "", "seed_bulking_status_map"),
    (11, "BL1", "BL1_seed_bulking_status_map"),
    (12, "", "genetic_example_P42_P34"),
    (13, "", "bulking_filter_funnel"),
    (14, "", "seed_bulking_model_map"),
]


def ink_box(im: Image.Image) -> tuple[int, int, int, int]:
    """Bounding box (left, top, right, bottom) of everything that is not white."""
    flat = Image.new("RGBA", im.size, "white")
    flat.alpha_composite(im.convert("RGBA"))
    ink = np.asarray(flat.convert("RGB")).min(axis=2) < 250
    rows, cols = np.where(ink.any(axis=1))[0], np.where(ink.any(axis=0))[0]
    return cols[0], rows[0], cols[-1] + 1, rows[-1] + 1


def main() -> None:
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir()
    missing, groups = [], {}
    for n, (fig, part, stem) in enumerate(TALK, start=1):
        src = next((f for f in (SRC / f"{stem}.png", SRC / f"{stem}.jpg") if f.exists()), None)
        if src is None:
            missing.append(stem)
            continue
        tag = f"Fig{fig:02d}" + (f"-{part}" if part else "")
        key = (fig, re.sub(r"[a-e]$", "", part))         # the frames of one build-up
        groups.setdefault(key, []).append((src, OUT / f"{n:02d}_{tag}_{src.name}"))
    for frames in groups.values():
        images = [Image.open(src) for src, _ in frames]
        boxes = np.array([ink_box(im) for im in images])
        w, h = images[0].size
        box = (max(boxes[:, 0].min() - EDGE, 0), max(boxes[:, 1].min() - EDGE, 0),
               min(boxes[:, 2].max() + EDGE, w), min(boxes[:, 3].max() + EDGE, h))
        for im, (_, dst) in zip(images, frames):
            im.crop(tuple(int(v) for v in box)).save(dst, dpi=(200, 200), quality=95)
    print(f"[talk_figures] {len(TALK) - len(missing)} figures → {OUT}/")
    if missing:
        print("  missing (run the pipeline part that makes them):", ", ".join(missing))


if __name__ == "__main__":
    main()
