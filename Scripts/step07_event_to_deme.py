"""Step 07 — event → population → deme table.

Writes `tables/step29a_event_to_deme.tsv`: one row per (event_year,
eventID) with its populationID, BL and demeID. Every later step that
works on demes reads this table (via `srk_units.load_event_units`), so
all analyses share the units defined in Steps 01 – 05.

Run after Step 05 (populationIDs and BLs are final) and Step 06.
The rebuilt demes are checked against
`tables/step29a_demes_per_population_year.tsv`; a mismatch stops the run.
"""
from srk_units import main

if __name__ == "__main__":
    main()
