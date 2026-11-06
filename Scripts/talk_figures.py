"""talk_figures.py — copy every figure of the seed-bulking talk into
talk_figures/, numbered in presentation order (zero-padded prefix), so the
slides can be built by sorting the folder.

Order follows Seed_bulking.md (Figures 1 – 14); where a figure has slide or
build-up versions, those come instead of (or before) the document version:
Figure 1 opens on the challenge alone and builds up step by step, as does
Figure 3; Figure 2 is split into two landscape slides (objectives 1 – 2,
then 3a – 3b) with build-up frames; Figure 8 into two 16:9 slides, Figures 9 and 10
start from an empty frame, and Figure 11 is followed by the BL1 close-up
(kept local, like every per-BL map). File names: NN_FigXX[-part]_name.png.

The folder is rebuilt from scratch on every run and is gitignored (copies of
figures_seed_bulking/). Run from the repository root:
    python Scripts/talk_figures.py
"""
from __future__ import annotations

import shutil
from pathlib import Path

SRC = Path("figures_seed_bulking")
OUT = Path("talk_figures")

# (figure number, part label, source file stem) in presentation order
TALK = [
    (1, "a", "breeding_system_cartoon_step0"),
    (1, "b", "breeding_system_cartoon_step1"),
    (1, "c", "breeding_system_cartoon_step2"),
    (1, "d", "breeding_system_cartoon_step3"),
    (1, "e", "breeding_system_cartoon"),
    (2, "slide1a", "restoration_program_cartoon_slide1_step1"),
    (2, "slide1b", "restoration_program_cartoon_slide1"),
    (2, "slide2a", "restoration_program_cartoon_slide2_step1"),
    (2, "slide2b", "restoration_program_cartoon_slide2"),
    (3, "a", "source_screening_cartoon_step1"),
    (3, "b", "source_screening_cartoon_step2"),
    (3, "c", "source_screening_cartoon"),
    (4, "", "EO_sampling_map_black"),
    (5, "", "BL_fragmentation_network_map"),
    (6, "", "units_cartoon_landscape"),
    (7, "", "population_EO_map_BL"),
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


def main() -> None:
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir()
    missing = []
    for n, (fig, part, stem) in enumerate(TALK, start=1):
        src = SRC / f"{stem}.png"
        if not src.exists():
            missing.append(src.name)
            continue
        tag = f"Fig{fig:02d}" + (f"-{part}" if part else "")
        shutil.copy2(src, OUT / f"{n:02d}_{tag}_{stem}.png")
    print(f"[talk_figures] {len(TALK) - len(missing)} figures → {OUT}/")
    if missing:
        print("  missing (run the pipeline part that makes them):", ", ".join(missing))


if __name__ == "__main__":
    main()
