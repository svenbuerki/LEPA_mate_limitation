"""talk_figures.py — copy every figure of the seed-bulking talk into
talk_figures/, numbered in presentation order (zero-padded prefix), so the
slides can be built by sorting the folder.

Order follows Seed_bulking.md (Figures 1 – 16, plus a field-board photo before
Figure 9); where a figure has slide or build-up versions, those come instead
of (or before) the document version: Figure 1 is the traffic light (Recovery
Criterion 2); Figure 2 is split into three slides (the challenge alone, panels
1 – 2, panels 3 – 4, each pair built panel by panel); Figure 3 opens on the
three types (then what they mean for restoration) and is split into two
landscape slides (objectives 1 – 2, then 3a – 3b); Figure 4 builds up step by
step; Figure 9 is split into two 16:9 slides; Figures 10 and 11 start from an
empty frame; Figure 12 is the traffic light filled with our populations;
Figure 13 (talk version, no EO27 zoom) is followed by the BL1 / BL3 close-ups
(kept local, like every per-BL map); Figure 16 (validation with genetics)
builds up panel 1, then adds panel 2. File names: NN_FigXX[-part]_name.png.

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
    (1, "", "traffic_light_cartoon"),                    # recovery plan (Criterion 2) → traffic light
    (2, "slide1", "breeding_system_cartoon_slide1"),
    (2, "slide2a", "breeding_system_cartoon_slide2_step1"),
    (2, "slide2b", "breeding_system_cartoon_slide2"),
    (2, "slide3a", "breeding_system_cartoon_slide3_step1"),
    (2, "slide3b", "breeding_system_cartoon_slide3"),
    (3, "typesa", "restoration_program_cartoon_types"),
    (3, "typesb", "restoration_program_cartoon_types_use"),
    (3, "slide1", "restoration_program_cartoon_slide1"),
    (3, "slide2a", "restoration_program_cartoon_slide2_step1"),
    (3, "slide2b", "restoration_program_cartoon_slide2"),
    (4, "a", "source_screening_cartoon_step1"),
    (4, "b", "source_screening_cartoon_step2"),
    (4, "c", "source_screening_cartoon_step3"),
    (4, "d", "source_screening_cartoon"),
    (5, "", "EO_sampling_map_BL_outline"),
    (6, "", "BL_fragmentation_network_map"),
    (7, "", "units_cartoon_landscape"),
    (8, "", "population_EO_map_BL"),
    (8, "BL1BL3", "BL1_BL3_populations_EO_map"),   # EO27 (BL1), EO24 (BL3); local only
    (9, "photo", "plant_photos/size_medium_OC2592_crown13cm"),   # field board, whole photo
    (9, "slide1", "plant_size_examples_slide1"),
    (9, "slide2", "plant_size_examples_slide2"),
    (10, "a", "seed_by_size_points"),
    (10, "b", "seed_by_size_plants"),
    (11, "a", "population_seed_yield_frame"),
    (11, "b", "population_seed_yield"),
    (12, "", "traffic_light_results"),                    # our populations under the traffic light
    (13, "", "seed_bulking_status_map_talk"),                  # no EO27 zoom (talk only)
    (13, "BL1BL3", "BL1_BL3_seed_bulking_status_map"),           # close-ups, local only
    (14, "", "bulking_filter_funnel"),
    (15, "", "seed_bulking_model_map"),
    (16, "a", "genetic_example_P42_P34_step1"),               # validation: S-types and plants only
    (16, "b", "genetic_example_P42_P34_step2"),               # + genotypes, 3 commonest sets
    (16, "c", "genetic_example_P42_P34"),
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
