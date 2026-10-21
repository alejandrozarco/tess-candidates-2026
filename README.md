# TESS transit candidates and other TESS results, 2026

Status: computational results from public data, not peer reviewed. Date of this record: 2026-10-21 (UTC); first version 2026-10-18. Author: A. Zarco.

> [!IMPORTANT]
> **Produced by AI models** under the direction of the repository owner; not peer reviewed, and not checked by a
> professional astronomer. See [`AI_DISCLOSURE.md`](AI_DISCLOSURE.md). The five candidates are leads for independent
> follow-up, not discoveries; none is confirmed. Credit for confirming, refining or refuting them belongs to whoever
> does that work. To refer to this record, please cite the archived repository: Zenodo, [doi:10.5281/zenodo.23453967](https://doi.org/10.5281/zenodo.23453967) (all versions; v1.1: [10.5281/zenodo.23463030](https://doi.org/10.5281/zenodo.23463030)).
> Questions, checks and corrections: [GitHub issues](https://github.com/alejandrozarco/tess-candidates-2026/issues).

This record reports periodic transit-shaped signals on five stars observed in TESS full-frame images, a candidate
second transit signal on one of them, and further TESS results: eclipsing binaries whose eclipses change, nine
eclipsing binaries without a previous eclipse report, dimming events of a young M dwarf, and changes in the transits
of three known planets. Each object has a detail page in [`docs/`](docs/); the section numbers of the earlier single-page
version are kept. The full summary paragraphs are in [docs/summary.md](docs/summary.md).

## Summary of all objects

| object | type | key numbers | status / caveat | earlier public records | details |
|---|---|---|---|---|---|
| TIC 402710280 | planet candidate | $`P = 36.34005`$ d, about 1680 ppm, $`R_\mathrm{p} \approx 3.7\,R_\oplus`$ | not confirmed; mid-transit times deviate from a linear ephemeris by up to 29 min (four epochs); FPP about $`10^{-3}`$ | R. Vora, Zenodo 2026-10-15; duotransit table samgill844/bruce (GitHub), 2025-10-23 | [section 1](docs/1_TIC402710280.md) |
| TIC 324429952 | planet candidate | $`P = 42.37574`$ d, about 1500 ppm, $`R_\mathrm{p} \approx 3.9\,R_\oplus`$ | not confirmed; FPP about 0.03 (above the validation threshold) | one sector-13 transit marked by a Planet Hunters TESS volunteer, 2019-08-30 | [section 2](docs/2_TIC324429952.md) |
| TIC 89247735 (SCR J1820-4533) | planet candidate | $`P = 18.323173`$ d, 1.13-1.22 % (dilution corrected), $`R_\mathrm{p} = 8.7 \pm 0.6\,R_\oplus`$ | not confirmed; a brown dwarf below about 30 $`M_\mathrm{J}`$ is not excluded | F. Molina, Zenodo 2026-10-17 | [section 3](docs/3_TIC89247735.md) |
| TIC 455017950 (HD 134147) | planet candidate | $`P = 40.60566`$ d, about 1310 ppm, $`R_\mathrm{p} = 3.16 \pm 0.12\,R_\oplus`$ | not confirmed; timing variation not established; FPP $`(6 \pm 4) \times 10^{-4}`$ | none listed | [section 4](docs/4_TIC455017950.md) |
| TIC 311179599 b | planet candidate | $`P = 24.82684`$ d, about 1435 ppm, $`R_\mathrm{p} = 2.47 \pm 0.10\,R_\oplus`$ | not confirmed; sector-39 transit 2.2 h early, timing variation not established; FPP 0.025-0.031 | none listed | [section 5](docs/5_TIC311179599.md) |
| TIC 311179599 c | candidate second transit signal | $`P = 73.8672`$ d or $`36.9336`$ d, about 880 ppm, $`R_\mathrm{p} \approx 2.0\,R_\oplus`$ | three events; the transit of 2026-11-02 13:16 UTC occurs only for the shorter period | none listed | [section 5.10](docs/5_TIC311179599.md#510-a-second-transit-signal-candidate-c) |
| TIC 11613676 (CD-30 448) | eclipsing binary, changing eclipses | $`P = 12.1158`$ d; primary absent in 2018, 0.69 % (2020), 5.4 % (2025), 4.0-4.1 % (2026) | inclination change inferred from models; third body not detected | Planet Hunters TESS volunteers, 2020-11/12 (sector 30), PHT eclipsing-binary candidate list; SPOC threshold-crossing events (2021) | [section 6.1](docs/6.1_TIC11613676.md) |
| TOI-540 b | known planet, transits not detected in one sector | $`P = 1.239`$ d; about 1.7 ppt in sectors 4-32, $`83 \pm 268`$ ppm in sector 98 | one sector; a data problem is not excluded; one ground-based transit decides | planet: Ment et al. (2021) | [section 6.2](docs/6.2_TOI-540b.md) |
| TIC 13957664 | eclipsing binary | $`P = 245.66`$ d (Gaia orbit), one 2.0 % eclipse | single eclipse at the conjunction predicted by the Gaia DR3 orbit | none listed | [section 6.3](docs/6.3_eclipsing_binaries.md) |
| TIC 304936165 | eclipsing binary | $`P = 61.7000`$ d, 2.5-2.8 % | two eclipses; period selected by the Gaia SB1 period | none listed | [section 6.3](docs/6.3_eclipsing_binaries.md) |
| TIC 120522352 | eclipsing binary (Gaia SB1) | $`P = 8.517`$ d, 6.3 % / 0.46 % | no spectroscopy | Planet Hunters TESS volunteers, 2020-11; PHT eclipsing-binary candidate list | [section 6.3](docs/6.3_eclipsing_binaries.md) |
| TIC 382523691 | eclipsing binary (Gaia SB1) | $`P = 4.915`$ d, 2.6 % / 1.17 % | no spectroscopy | Planet Hunters TESS volunteers, 2020-11; PHT eclipsing-binary candidate list | [section 6.3](docs/6.3_eclipsing_binaries.md) |
| TIC 261702897 | eclipsing binary (Gaia SB1) | $`P = 6.934`$ d, 7.5 % / 0.55 % | no spectroscopy | none listed | [section 6.3](docs/6.3_eclipsing_binaries.md) |
| TIC 265436769 | eclipsing binary (Gaia SB1) | $`P = 15.475`$ d, 0.52 % | grazing; no secondary eclipse | none listed | [section 6.3](docs/6.3_eclipsing_binaries.md) |
| TIC 433642561 (CPD-62 2890) | eclipsing binary in NGC 4609 | $`P = 356.182/n`$ d, 17.5-27 % | period not determined | none listed | [section 6.3](docs/6.3_eclipsing_binaries.md) |
| TIC 245068282 | eclipsing binary in a hierarchical triple | $`P = 1.7075`$ d, about 50 % / 18 % | association with a G dwarf 13.2" away from Gaia parallax and proper motion | Planet Hunters TESS volunteers, 2026-04/06 | [section 6.3](docs/6.3_eclipsing_binaries.md) |
| TIC 342194997 (TYC 9251-1375-1) | eclipsing binary (known Gaia binary) | $`P = 550.75`$ d (Gaia orbit), one 16.9 % eclipse | single eclipse; predicted occultation not observed | T. Chopard-Vilhem, Zenodo 2026-10-13 | [section 6.3](docs/6.3_eclipsing_binaries.md) |
| TIC 124751008 | eclipsing binary, changing eclipses | $`P = 12.61`$ d; primary 2.25 % (2018), 0.28 % (2024-2025), 0.77 % (2026) | inclination change inferred from models; third body not detected | none listed | [section 6.4](docs/6.4_TIC124751008.md) |
| HIP 23309 (TIC 220473309) | dimming events | three further events of 1.9-5.7 %, about 1 d long | no interpretation established | two earlier events: Gollotti et al. (2024) | [section 6.5](docs/6.5_HIP23309.md) |
| L 98-59 d | known planet, transit-duration change candidate | $`\mathrm{d}T_{14}/\mathrm{d}t = -1.17 \pm 0.12`$ min/yr (about 52 to 44 min, 2018.7-2025.96) | change of the orbit not established; error conditional on a fixed transit shape | earlier duration search, not significant: Bessho et al. (2025) | [section 6.6](docs/6.6_L98-59d.md) |
| TOI-620 b | known planet, transit-duration change candidate (weaker) | $`-2.9 \pm 0.4`$ min/yr (about 80 to 62 min, 2019-2026) | slope depends on the aperture; five epochs | TTV report: Naponiello (2026) | [section 6.7](docs/6.7_TOI-620b.md) |

Earlier public records. Two of the five signals were reported publicly by others before this record: TIC 402710280 by
R. Vora (Zenodo, 2026-10-15) and TIC 89247735 by F. Molina (Zenodo, 2026-10-17); see sections
[1.8](docs/1_TIC402710280.md#18-catalogue-and-literature-checks) and [3.8](docs/3_TIC89247735.md#38-catalogue-and-literature-checks).
TIC 402710280 is also listed in a public duotransit table (2025-10-23), and a single sector-13 transit of TIC 324429952
was marked by a Planet Hunters TESS volunteer in 2019 (section [2.7](docs/2_TIC324429952.md#27-catalogue-and-literature-checks)).
Earlier records of the eclipses of TIC 11613676 and of four of the binaries of section 6.3 are listed in sections
[6.1.7](docs/6.1_TIC11613676.md#617-catalogue-and-literature-checks) and [6.3](docs/6.3_eclipsing_binaries.md).
HIP 23309 has two published dimming events (Gollotti et al. 2024); the three events of section
[6.5](docs/6.5_HIP23309.md) are further ones.

## Data and reproduction

The tables in [`data/`](data/) are recomputed from public MAST, Gaia and ESO products by the scripts in
[`scripts/`](scripts/); the exceptions (values from the vetting record without a script here), the commands and the
software versions are in section [8](docs/8_data_and_reproduction.md). The figures are regenerated with
`python figures/make_figures.py` from the CSV files alone; `shasum -a 256 -c MANIFEST.sha256` checks every file.

## Caveats

- Not peer reviewed; nothing in this record has been checked by anyone outside the production of the record.
- No high-resolution imaging and no radial velocities were obtained for this work; none of the candidates is
  confirmed, and companions close to the targets are constrained only partly (by Gaia DR3 and, for TIC 455017950, a
  published SOAR speckle observation).
- The deviations of the mid-transit times from linear ephemerides rest on few epochs; for TIC 455017950 and
  TIC 311179599 a transit-timing variation is not established.
- The third bodies of TIC 11613676 and TIC 124751008 are inferred from the eclipse histories, not detected; the
  nodal and outer periods are orders of magnitude. The TOI-540 b non-detection rests on one TESS sector.
- The duration changes of L 98-59 d and TOI-620 b rest on few observing campaigns of one instrument (with modest JWST
  support for L 98-59 d), and their significances are not corrected for the search of 802 planets that selected them.
- QLP light curves carry a barycentric timing error in sectors 99-104; no QLP time stamp enters the transit tables.
  Full list: section [7](docs/7_caveats.md).

## How to cite

A. Zarco (2026). TESS transit candidates and other TESS results, 2026. Zenodo. DOI: 10.5281/zenodo.23453967 (all versions).
Citation metadata: `CITATION.cff`.

## Licence

Data, figures and text: CC BY 4.0. Scripts: MIT. See `LICENSE`, which also lists the third-party material not covered.

## Contents

- [`docs/`](docs/): sections 1-8 (one file per section; section 6 split by subsection), the
  [summary paragraphs](docs/summary.md), the [references](docs/references.md) and the
  [per-file list](docs/files.md)
- [`data/`](data/): CSV tables; lines starting with `#` identify the file and its columns
- [`figures/`](figures/): light and dark SVG figures, made by `figures/make_figures.py` from `data/`
- [`scripts/`](scripts/): download, photometry, fit, vetting and table scripts
- `AI_DISCLOSURE.md`, `CITATION.cff`, `.zenodo.json`, `LICENSE`, `requirements.txt`, `MANIFEST.sha256`: disclosure,
  citation metadata, licence, environment, checksums
