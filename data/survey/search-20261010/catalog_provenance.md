# Search target provenance

- TIC data: MAST TIC via astroquery `Catalogs.query_criteria`, queried 2026-10-10 in two subsets (M dwarfs: Teff 3000–4000 K, d<40 pc, R<0.6 R☉, Tmag<12; K dwarfs: Teff 4000–5200 K, d<40 pc, R<0.9 R☉, Tmag<12), `objType=STAR`. Strict target cuts applied after retrieval: 2700–5200 K, 0<d<100 pc, 0<R<0.9 R☉, 0<Tmag<14.
- Known transiting planet periods: NASA Exoplanet Archive TAP `ps` table, `tic_id, pl_orbper`, fetched 2026-10-10. TOI periods were queried from the NASA Exoplanet Archive `toi` table and merged with confirmed-planet `ps` periods; candidate preflight also checks current TOIs.
- TIC query returned 7219 rows before deduplication; 7219 unique rows passed strict cuts.
- Sector/cadence coverage is queried for the ranked shortlist before survey launch. Empty coverage fields mean unknown.
- Ranking is a model-based selection heuristic, not a per-star planet probability.
