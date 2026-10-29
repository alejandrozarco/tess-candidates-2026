# TESS transit candidates and other TESS results, 2026

Status: computational results from public data, not peer reviewed. Date of this record: 2026-10-29 (UTC); first version 2026-10-18. Author: A. Zarco.

> [!IMPORTANT]
> **Produced by AI models** under the direction of the repository owner; not peer reviewed, and not checked by a
> professional astronomer. See [`AI_DISCLOSURE.md`](AI_DISCLOSURE.md). The five candidates are leads for independent
> follow-up, not discoveries; none is confirmed. Credit for confirming, refining or refuting them belongs to whoever
> does that work. To refer to this record, please cite the archived repository: Zenodo, [doi:10.5281/zenodo.23453967](https://doi.org/10.5281/zenodo.23453967) (all versions; v1.1: [10.5281/zenodo.23463030](https://doi.org/10.5281/zenodo.23463030)).
> Questions, checks and corrections: [GitHub issues](https://github.com/alejandrozarco/tess-candidates-2026/issues).

This record reports periodic transit-shaped signals on five stars observed in TESS full-frame images, a candidate
second transit signal on one of them, and further TESS results: eclipsing binaries whose eclipses change, nine
eclipsing binaries without a previous eclipse report, dimming events of a young M dwarf, and changes in the transits
of four known planets and of one TESS planet candidate (TOI-2638.01), and a transit-like signal at the period and phase
of TOI-1117 c, a planet known from radial velocities. Each object has a detail page in [`docs/`](docs/); the section numbers of the earlier single-page
version are kept. The full summary paragraphs are in [docs/summary.md](docs/summary.md).

## Summary of all objects

| object | type | key numbers | status / caveat | earlier public records | details |
|---|---|---|---|---|---|
| TIC 402710280 | planet candidate | $`P = 36.34005`$ d, about 1680 ppm, $`R_\mathrm{p} \approx 3.7\,R_\oplus`$ | not confirmed; mid-transit times deviate from a linear ephemeris by up to 29 min (four epochs); FPP about $`10^{-3}`$ | R. Vora, Zenodo 2026-10-15; duotransit table samgill844/bruce (GitHub), 2025-10-23 | [section 1](docs/1_TIC402710280.md) |
| TIC 324429952 | planet candidate | $`P = 42.37574`$ d, about 1500 ppm, $`R_\mathrm{p} \approx 3.9\,R_\oplus`$ | not confirmed; FPP about 0.03 (above the validation threshold) | one sector-13 transit marked by a Planet Hunters TESS volunteer, 2019-08-30 | [section 2](docs/2_TIC324429952.md) |
| TIC 89247735 (SCR J1820-4533) | planet candidate | $`P = 18.323173`$ d, 1.13-1.22 % (dilution corrected), $`R_\mathrm{p} = 8.7 \pm 0.6\,R_\oplus`$ | not confirmed; a brown dwarf below about 30 $`M_\mathrm{J}`$ is not excluded | F. Molina, Zenodo 2026-10-17 | [section 3](docs/3_TIC89247735.md) |
| TIC 455017950 (HD 134147) | planet candidate | $`P = 40.60566`$ d, about 1310 ppm, $`R_\mathrm{p} = 3.16 \pm 0.12\,R_\oplus`$ | not confirmed; timing variation not established; FPP $`(6 \pm 4) \times 10^{-4}`$ | none listed | [section 4](docs/4_TIC455017950.md) |
| TIC 311179599 b | planet candidate | $`P = 24.82686`$ d, about 1430 ppm, $`R_\mathrm{p} = 2.47 \pm 0.10\,R_\oplus`$ | not confirmed; sector-39 transit 2.2 h early, timing variation not established; FPP 0.025-0.031 | none listed | [section 5](docs/5_TIC311179599.md) |
| TIC 311179599 c | candidate second transit signal | $`P = 73.8672`$ d or $`36.9336`$ d, about 880 ppm, $`R_\mathrm{p} \approx 2.0\,R_\oplus`$ | three events; the transit of 2026-11-02 13:16 UTC occurs only for the shorter period | none listed | [section 5.10](docs/5_TIC311179599.md#510-a-second-transit-signal-candidate-c) |
| TIC 11613676 (CD-30 448) | eclipsing binary, changing eclipses | $`P = 12.1158`$ d; primary absent in 2018, 0.69 % (2020), 5.4 % (2025), 4.0-4.1 % (2026) | inclination change inferred from models; third body not detected | Planet Hunters TESS volunteers, 2020-11/12 (sector 30), PHT eclipsing-binary candidate list; SPOC threshold-crossing events (2021) | [section 6.1](docs/6.1_TIC11613676.md) |
| TOI-540 b | known planet, transits not detected in one sector | $`P = 1.239`$ d; about 1.7 ppt in sectors 4-32, $`81 \pm 268`$ ppm in sector 98 | one sector; a data problem is not excluded; one ground-based transit decides | planet: Ment et al. (2021) | [section 6.2](docs/6.2_TOI-540b.md) |
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
| L 98-59 d | known planet, evidence for a decreasing transit duration | $`\mathrm{d}T_{14}/\mathrm{d}t = -1.185 \pm 0.112`$ min/yr (about 51 to 44 min, 2019-2025; nominal local $`z`$ = 6.7-7.4; an illustrative Gaussian trials adjustment gives $`z`$ = 5.4-6.4, subject to the assumptions in section 6.6.10) | rests on four TESS campaigns (1.8 $`\sigma`$ within campaigns); its two JWST NIRSpec detectors differ by about 0.8 min per visit, with opposite signs; with the post-hoc detector averaging and empirical systematic floor, the NIRSpec pair differs from no change by 1.3 $`\sigma`$, and adding NIRISS gives a slope 2.2 $`\sigma`$ from zero; both are consistent with the TESS trend, and JWST does not independently confirm it; change of the orbit not established | earlier duration search, not significant: Bessho et al. (2025) | [section 6.6](docs/6.6_L98-59d.md) |
| TOI-620 b | known planet, transit-duration change candidate (weaker) | $`-2.9 \pm 0.4`$ min/yr (about 80 to 62 min, 2019-2026) | slope depends on the aperture; five epochs | TTV report: Naponiello (2026) | [section 6.7](docs/6.7_TOI-620b.md) |
| TOI-6482.01 (TIC 189300053) | known TOI, ephemeris update | $`P = 20.047937 \pm 0.000060`$ d (six TESS transits, 2019-2026; TOI table: 1483.5175 d), about 9.4 ppt, $`R_\mathrm{p} \approx 17\,R_\oplus`$ (TIC stellar radius) | planet candidate (TOI); nature not established | TOI (TESS project); TFOP period 20.04753 d in the ExoFOP notes; earlier version of this ephemeris: ExoFOP upload 2026-10-20 | [section 6.8](docs/6.8_TOI-6482.md) |
| TOI-2638.01 (TIC 458342769) | known TOI, transit depth increasing | $`P = 3.92995`$ d; depth 1.8 ppt (2021) to 4.3 ppt (2026) in SPOC PDCSAP, $`+384 \pm 56`$ ppm/yr; grazing; A-type host, $`v \sin i = 110`$ km/s | change seen in every TESS product and on the target in pixel-level fits; the shape change that separates a geometric change from dilution is about 2.6-3.1 $`\sigma`$; cause and companion nature not established; TESS data only | TOI (TESS project); no earlier report of the change found | [section 6.9](docs/6.9_TOI-2638.md) |
| TOI-1117 c (TIC 295541511) | known RV planet, transit-like signal | $`P = 4.5786510`$ d, about 490 ppm, $`T_{14} = 1.85`$ h, $`R_\mathrm{p} = 2.73^{+0.16}_{-0.15}\,R_\oplus`$ if on the target; times within 1.0 $`\sigma`$ of the RV-predicted conjunctions of c | not confirmed; a $`G = 19.6`$ star at 34" disfavoured only at $`\Delta\chi^2 = 7.5`$; impact parameter 0.85 with the stellar-density prior, boxier shape without it ($`\chi^2`$ lower by 3.5, not calibrated as a significance); TRICERATOPS FPP 0.0146-0.0149 with an approximate contrast curve, NFPP 0.013 (above the 0.001 threshold); TESS data and published RVs only | planet c: Lockley et al. (2025), RVs, described as non-transiting; no earlier transit report found (ExoFOP unreachable) | [section 6.10](docs/6.10_TOI-1117c.md) |
| LHS 475 b (TIC 369327947, GJ 4102, TOI-910.01) | known planet (JWST transmission spectrum), candidate decrease of the transit duration | $`P = 2.0291`$ d; $`T_{14}`$ about 42.6 min (2019-2021) to about 35.4 min (2025-2026); $`\mathrm{d}b/\mathrm{d}t = +0.0156 \pm 0.0043`$ yr$`^{-1}`$ with injection-calibrated errors at $`a/R_\star = 15.77`$ (local 3.7 $`\sigma`$) | selected from a screen of 271 known planets, about 2-3 $`\sigma`$ after that selection; TESS data only (ten sectors, per-transit S/N about 4); the rate in $`b`$ depends on the adopted $`a/R_\star`$; own refit of the 2022 JWST transits ($`T_{14} = 39.85`$ min) agrees with the trend but does not confirm it | planet: Ment et al. (2024), Lustig-Yaeger et al. (2023); no earlier report of a change found | [section 6.11](docs/6.11_LHS475b.md) |

Earlier public records. Two of the five signals were reported publicly by others before this record: TIC 402710280 by
R. Vora (Zenodo, 2026-10-15) and TIC 89247735 by F. Molina (Zenodo, 2026-10-17, doi:10.5281/zenodo.23443245); see sections
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
- The duration-change evidence for L 98-59 d and TOI-620 b rests on few observing campaigns of one instrument. For
  L 98-59 d, the within-campaign diagnostic is only 1.8 $`\sigma`$ and JWST does not independently confirm the trend. Its quoted
  trials-adjusted significance assumes Gaussian tails and 801 or 3204 independent trials; TOI-620 b has no
  corresponding adjustment.
- The depth change of TOI-2638.01 rests on seven TESS sectors (one of them QLP only); the shape change that distinguishes
  a change of the transit geometry from a change of dilution is about 2.6-3.1 $`\sigma`$. Only TESS data are used; the
  ground-based follow-up listed on ExoFOP is not used.
- The transit-like signal of TOI-1117 c (section 6.10) rests on TESS data and the published HARPS velocities. It is not
  validated: the pixel data disfavour a $`G = 19.6`$ star at 34" only at $`\Delta\chi^2 = 7.5`$, the transit needs an impact
  parameter of 0.85 with the stellar-density prior, the sector-13 transits are 13-22 min early, and ExoFOP could not be
  checked for earlier reports.
- The duration change of LHS 475 b (section 6.11) rests on ten TESS sectors of one instrument at a per-transit S/N of
  about 4. Its local significance (3.7 $`\sigma`$) uses injection-calibrated errors; the planet was selected from a
  screen of 271 planets, which leaves about 2-3 $`\sigma`$. The rate in $`b`$ depends on the adopted $`a/R_\star`$; the
  JWST epoch of 2022 lies mid-baseline and does not confirm the change.
- QLP light curves carry a barycentric timing notice for sectors 99-104. QLP-derived mid-times enter the ephemerides of
  TOI-6482.01 (section 6.8) and TOI-2638.01 (sector 10 only; section 6.9). The TOI-2638 product comparisons also use QLP
  timestamps, including sector 99. See section [7](docs/7_caveats.md) for the timing qualifications.

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
