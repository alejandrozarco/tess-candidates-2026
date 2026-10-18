# TESS transit candidates 2026: TIC 402710280 and TIC 324429952

Status: computational results from public data, not peer reviewed. Date of this record: 2026-10-18 (UTC). Author: A. Zarco.

> [!IMPORTANT]
> **Produced by AI models** under the direction of the repository owner; not peer reviewed, and not checked by a
> professional astronomer. See [`AI_DISCLOSURE.md`](AI_DISCLOSURE.md). The two candidates are leads for independent
> follow-up, not discoveries; neither is confirmed. Credit for confirming, refining or refuting them belongs to whoever
> does that work. To refer to this record, please cite the archived repository (Zenodo DOI: to be added).
> Questions, checks and corrections: [GitHub issues](https://github.com/alejandrozarco/tess-candidates-2026/issues).

Two stars observed in TESS full-frame images show periodic transit-shaped dips that are not listed as TOI or CTOI, in
the NASA Exoplanet Archive, VSX, SIMBAD, the Gaia DR3 variability and non-single-star tables, VizieR or ADS (checks
dated in sections 1.8 and 2.7). TIC 402710280 (K dwarf, 67 pc, $`T = 9.24`$): $`P = 36.34005`$ d, depth about 1680 ppm,
duration 5.3 h, flat-bottomed, $`R_\mathrm{p} \approx 3.7\,R_\oplus`$, four transits in sectors 39, 93, 100 and 102
whose mid-transit times deviate from a linear ephemeris by up to 29 min; TRICERATOPS false-positive probability about
$`10^{-3}`$; two public HARPS epochs give a radial-velocity change of $`-4.1 \pm 1.8`$ m/s between orbital phases 0.046
and 0.074. TIC 324429952 (G dwarf, 108 pc, $`T = 9.64`$): $`P = 42.37574`$ d, depth about 1500 ppm, duration 3.8-4.3 h,
$`R_\mathrm{p} \approx 3.9\,R_\oplus`$, eight transits in sectors 13 to 103 on a linear ephemeris; TRICERATOPS
false-positive probability about 0.03. No high-resolution imaging and no radial velocities were obtained for this work.
The tables in `data/` are recomputed from public MAST, Gaia and ESO products by the scripts in `scripts/`.

## Contents

| path | content |
|---|---|
| `README.md` | this document |
| `AI_DISCLOSURE.md` | how the content was produced |
| `data/transit_times_<tic>.csv` | mid-transit times, errors, O-C, per-epoch depths |
| `data/folded_<tic>.csv` | detrended light curve of every transit epoch (TTV-corrected) with the model |
| `data/model_<tic>.csv` | transit model curve used in the figures |
| `data/sectors_<tic>.csv` | TESS sectors used, cadence, spans, predicted transits and their coverage |
| `data/predicted_transits_<tic>.csv` | predicted mid-transit times for the 12 months from 2026-10-18 |
| `data/tess_sectors_<tic>.csv` | TESS sectors (observed and future) in which the target is on silicon |
| `data/diffimage_<tic>.csv`, `data/diffimage_sources_<tic>.csv` | per-sector difference images and source pixel positions |
| `data/diffimage_centroids_324429952.csv` | difference-image centroid offsets per sector |
| `data/harps_public_rvs.csv` | the two public HARPS epochs of TIC 402710280 |
| `data/triceratops_input_<tic>.csv` | the binned light curve given to TRICERATOPS |
| `data/triceratops/` | TRICERATOPS scenario tables (`probs_*.csv`, verbatim), run index (`runs.csv`), per-configuration summary (`summary.csv`), Gaia background samples, the assumed contrast curve of the reference runs |
| `figures/` | `fold_`, `epochs_`, `oc_`, `diffimage_` figures, light and dark SVG, made by `figures/make_figures.py` from `data/` |
| `scripts/` | download, photometry, fit, difference images, alias/secondary search, prediction, TRICERATOPS and table scripts |
| `requirements.txt`, `LICENSE`, `CITATION.cff`, `.zenodo.json`, `MANIFEST.sha256` | environment, licence, citation metadata, checksums |

Lines starting with `#` in the CSV files are comments that identify the file and its columns.
Times: BTJD = BJD$`_\mathrm{TDB}`$ $`-`$ 2457000. All mid-transit times come from TESScut cutouts or SPOC light curves
(section 3 on the QLP timing notice).

## 1. TIC 402710280

### 1.1 Host star

| quantity | value | source |
|---|---|---|
| identifiers | TIC 402710280, Gaia DR3 5828441289812582400, TYC 9041-3677-1, 2MASS J16340568-6332588 | SIMBAD, ExoFOP |
| position (ICRS, epoch 2016.0) | RA 248.5242026, Dec $`-`$63.5496191 | Gaia DR3 |
| proper motion (mas/yr) | 59.909, 13.685 | Gaia DR3 |
| parallax (mas) | $`14.9636 \pm 0.0137`$ | Gaia DR3 |
| distance (pc) | $`66.63 \pm 0.14`$ | TIC v8.2 |
| $`G`$, $`BP-RP`$ | 9.781, 1.023 | Gaia DR3 |
| $`T`$, $`V`$, $`B`$, $`J`$, $`H`$, $`K`$ | 9.243, 9.982, 10.828, 8.528, 8.120, 8.026 | TIC v8.2 (ExoFOP) |
| $`T_\mathrm{eff}`$ (K) | $`5182 \pm 121`$ (TIC); 5089 (Gaia GSP-Phot) | TIC v8.2; Gaia DR3 |
| $`R_\star`$ ($`R_\odot`$) | $`0.817 \pm 0.047`$ | TIC v8.2 |
| $`M_\star`$ ($`M_\odot`$) | $`0.88 \pm 0.11`$ | TIC v8.2 |
| $`\log g`$ | $`4.558 \pm 0.083`$ | TIC v8.2 |
| $`\rho_\star`$ (g cm$`^{-3}`$) | $`2.27 \pm 0.53`$ | TIC v8.2 |
| RUWE; `ipd_frac_multi_peak`; `non_single_star` | 1.001; 0; 0 | Gaia DR3 |
| radial velocity (km/s) | $`0.63 \pm 0.17`$ (26 transits); `rv_amplitude_robust` 1.93 | Gaia DR3 |
| TIC contamination ratio | 0.060 | TIC v8.2 |
| nearest Gaia DR3 source | 5828441289812582656, $`G = 16.15`$, separation 2.3" at epoch 2016.0 (up to about 3" at the TESS epochs, the target moving at 61 mas/yr); its parallax places it at about 1.45 kpc | Gaia DR3 |
| X-ray | an eRASS1 source is matched to the target (VizieR J/A+A/682/A34) | eROSITA |

### 1.2 Transits and ephemeris

Common-shape fit of the four transits in the TESScut $`r = 2.0`$ px circular-aperture photometry (`scripts/fit_transits.py`;
batman model, quadratic limb darkening $`u = (0.46, 0.19)`$, one quadratic baseline per epoch solved jointly, quality-0
cadences within $`\pm 1`$ d of each transit):

| quantity | value |
|---|---|
| $`P`$ (d) | $`36.34005 \pm 0.00006`$ |
| $`T_0`$ (BJD$`_\mathrm{TDB}`$) | $`2461128.92695 \pm 0.00140`$ (BTJD $`4128.92695`$) |
| $`k = R_\mathrm{p}/R_\star`$ | $`0.041 \pm 0.004`$ |
| depth $`k^2`$ | about 1680 ppm (per epoch 1619-1759 ppm, table 1.3) |
| $`T_{14}`$ (h) | 5.27 |
| $`a/R_\star`$ | 54.8 |
| $`b`$ | 0 (at the bound of the fit; $`b \approx 0.3`$ for the TIC density) |
| $`\rho_\star`$ from $`a/R_\star`$, circular orbit (g cm$`^{-3}`$) | 2.36 (TIC: $`2.27 \pm 0.53`$) |
| $`R_\mathrm{p} = k R_\star`$ ($`R_\oplus`$) | $`3.7 \pm 0.4`$ (TRICERATOPS best-fit planet radius 3.76) |
| $`a = (a/R_\star) R_\star`$ (au) | 0.21 |
| detection-stage trapezoid (3 $`\times`$ 3-pixel box, five events incl. S66) | depth $`1835 \pm 36`$ ppm, $`T_{14} = 5.03 \pm 0.06`$ h, $`T_\mathrm{flat}/T_{14} = 0.96`$ |

The $`\pm 0.004`$ on $`k`$ covers the spread between apertures ($`r`$ = 1.5, 2.0, 3.0 px) and baseline choices; the formal
error of the single fit is smaller. The QLP detrended light curves show the transits at about 550 ppm (the detrending
suppresses them); the QLP SAP flux and the TESS-SPOC, TARS and TESScut photometry show them at the depths given here.

### 1.3 Transit epochs and O-C

`data/transit_times_402710280.csv`. Errors are the larger of (formal error $`\times \beta`$) and the prayer-bead error,
where $`\beta`$ is the time-averaging factor of the residuals over 20-90 min bins.

| epoch | sector | $`t_c`$ (BTJD) | error (min) | O $`-`$ C (min) | depth (ppm) | note |
|---|---|---|---|---|---|---|
| $`-48`$ | 39 | 2384.609329 | 3.4 | $`+6.9`$ | $`1664 \pm 84`$ | |
| $`-28`$ | 66 | 3111.397293 | 25.6 | $`-11.9`$ | $`1962 \pm 162`$ | scattered light; not used |
| $`-8`$ | 93 | 3838.186504 | 2.8 | $`-28.9`$ | $`1619 \pm 72`$ | |
| $`-1`$ | 100 | 4092.585400 | 4.9 | $`-2.2`$ | $`1643 \pm 95`$ | |
| 0 | 102 | 4128.943317 | 2.7 | $`+23.6`$ | $`1759 \pm 67`$ | |

Linear ephemeris through the four used epochs: $`\chi^2 = 190`$ for 2 degrees of freedom. The pattern is the same for
$`r`$ = 1.5, 2.0 and 3.0 px apertures, for linear, quadratic and cubic baselines, and in the QLP SAP light curves (the QLP
time stamps themselves were not used; section 3). The sector-103 transit (BTJD 4165.27) falls in a data gap; sectors 12
and 101 contain no predicted transit.

<picture><source srcset="figures/fold_402710280_dark.svg" media="(prefers-color-scheme: dark)"><img src="figures/fold_402710280_light.svg" alt="TIC 402710280 folded transit"></picture>

<picture><source srcset="figures/epochs_402710280_dark.svg" media="(prefers-color-scheme: dark)"><img src="figures/epochs_402710280_light.svg" alt="TIC 402710280 individual transits"></picture>

<picture><source srcset="figures/oc_402710280_dark.svg" media="(prefers-color-scheme: dark)"><img src="figures/oc_402710280_light.svg" alt="TIC 402710280 O-C"></picture>

### 1.4 Checks against false positives

- The four transits are present in every public light curve of the star (QLP SAP, TESS-SPOC SAP and PDCSAP, TARS) and in
  the independent TESScut photometry; the sector-66 event (BTJD 3111.397) lies in scattered light and is kept out of the
  ephemeris.
- Odd and even transits agree in depth (difference 0.4 $`\sigma`$).
- Secondary eclipse: $`< 170`$ ppm (3 $`\sigma`$) at phase 0.5 (`scripts/alias_secondary_search.py`).
- Period aliases $`P/k`$, $`k = 2`$-12: the predicted epochs that fall in the data show no transit.
- Transit shape: grazing geometries are disfavoured ($`\Delta\chi^2 = 95`$ over the $`b`$ profile), so a host diluted by
  more than 0.4 mag is excluded; the 2.3" neighbour ($`G = 16.15`$, 6.4 mag fainter, background) is excluded as the host
  on this ground.
- Difference images of sectors 39, 93, 100 and 102 peak on the target (section 1.5); in-transit centroid shifts exclude
  the stars at 39-60" at more than 10 $`\sigma`$.
- Gaia DR3: RUWE 1.001, no multi-peak detections, no non-single-star solution; the robust radial-velocity amplitude of
  1.93 km/s over 26 transits limits a companion at this period to below about 20 $`M_\mathrm{J}`$.
- The archival HARPS epochs (section 1.7) exclude stellar and brown-dwarf companions at $`P = 36.34`$ d.

### 1.5 Difference-image localisation

`data/diffimage_402710280.csv`, `data/diffimage_sources_402710280.csv`, `scripts/difference_images.py`. The difference
image (out-of-transit model minus in-transit mean of every cutout pixel) is fitted with a free 2-D Gaussian and compared
with the out-of-transit image; the dimming is centred on the target in all four sectors, and the stars at 39", 54" and
60" are excluded at more than 10 $`\sigma`$. The 2.3" neighbour is unresolved at the 21" pixel scale and is constrained
only by the transit shape (section 1.4) and by the TRICERATOPS run of section 1.6.

<picture><source srcset="figures/diffimage_402710280_dark.svg" media="(prefers-color-scheme: dark)"><img src="figures/diffimage_402710280_light.svg" alt="TIC 402710280 difference images"></picture>

### 1.6 TRICERATOPS

TRICERATOPS 1.1.0 (Giacalone et al. 2021). Input: the TTV-corrected fold of the four TESScut epochs, each divided by its
jointly fitted quadratic baseline ($`\pm 1.0`$ d), 10-min inverse-variance-weighted bins within $`\pm 0.5`$ d (144 bins;
out-of-transit bin scatter 220 ppm, used as the single flux error); raw aperture flux with the per-sector $`r < 2.0`$ px
apertures given to `calc_depths` (`dilution_corrected=False`); background population from Gaia DR3 (the TRILEGAL service
was not reachable); $`N = 10^6`$ draws, `nsamples=20`, 10-min exposure time; 10 random seeds. No contrast curve.

| configuration | seeds | FPP mean $`\pm`$ std | FPP median (range) | NFPP |
|---|---|---|---|---|
| nominal | 10 | $`0.0011 \pm 0.0021`$ | $`2.6 \times 10^{-4}`$ ($`2.8 \times 10^{-6}`$ to $`6.5 \times 10^{-3}`$) | $`\sim 10^{-25}`$ |
| flux errors $`\times 1.5`$ | 5 | $`0.0049 \pm 0.0046`$ | $`3.4 \times 10^{-3}`$ ($`1.8 \times 10^{-3}`$ to $`1.3 \times 10^{-2}`$) | $`\sim 10^{-11}`$ |
| 2.3" neighbour set to a K dwarf (0.85 $`M_\odot`$, 0.80 $`R_\odot`$, 5200 K; it has no TIC parameters) | 5 | $`0.0021 \pm 0.0027`$ | $`7.1 \times 10^{-4}`$ ($`1.8 \times 10^{-4}`$ to $`6.5 \times 10^{-3}`$) | $`\sim 10^{-25}`$ |

Mean scenario probabilities of the nominal runs: TP 0.84, PTP 0.125, DTP 0.034, STP 0.0011 (the whole FPP); every
eclipsing-binary scenario on the target or on neighbours is below $`10^{-6}`$. All nominal seeds and four of the five
runs with inflated errors are below the thresholds FPP $`< 0.015`$ and NFPP $`< 0.001`$ of Giacalone et al. (2021). Runs
with an assumed speckle contrast curve are included in `data/triceratops/` for reference only (`cc_hypothetical = True`
in `runs.csv`); they do not represent an observation.

### 1.7 Archival radial velocities

The ESO Science Archive holds five HARPS science exposures of the star: two public (2025-08-25 and 2025-08-26, programme
112.261U.004) and three proprietary (2026-02-23, 2026-02-28 and 2026-03-29, programme 116.29B6.001, not accessed). The
two public epochs were retrieved anonymously (`data/harps_public_rvs.csv`, with the archive identifiers and checksums)
and the DRS 3.8 pipeline values (K5 mask cross-correlation) are used as delivered:

| date (UTC) | BJD$`_\mathrm{TDB}`$ | RV (m/s) | FWHM (km/s) | contrast (%) | bisector span (m/s) | orbital phase |
|---|---|---|---|---|---|---|
| 2025-08-25 01:18:46 | 2460912.562865 | $`800.94 \pm 1.07`$ | 6.329 | 42.7 | 44.9 | 0.046 |
| 2025-08-26 01:32:12 | 2460913.572387 | $`796.85 \pm 1.42`$ | 6.323 | 42.8 | 35.1 | 0.074 |

$`\Delta\mathrm{RV} = -4.1 \pm 1.8`$ m/s. On the transit ephemeris with a circular orbit the two epochs sample
$`\Delta\mathrm{RV} = -0.162\,K`$, i.e. 16 % of the semi-amplitude, so $`K = 25 \pm 11`$ m/s (photon noise only),
$`\lvert K\rvert < 58`$ m/s at 3 $`\sigma`$ and $`M_2 \sin i < 0.87\,M_\mathrm{J}`$ for $`M_\star = 0.88\,M_\odot`$. A
stellar or brown-dwarf companion at $`P = 36.34`$ d would have produced $`\Delta\mathrm{RV}`$ between $`-0.14`$ km/s
(13 $`M_\mathrm{J}`$, $`K = 0.86`$ km/s) and $`-2.8`$ km/s (0.3 $`M_\odot`$, $`K = 17.3`$ km/s) over the same baseline. A
sub-Neptune of 10-20 $`M_\oplus`$ ($`K`$ = 2-4 m/s) is not detectable with two epochs. The DRS cross-correlation function
is a single Gaussian on both nights (identical FWHM and contrast), and a $`\pm 400`$ km/s self-mask cross-correlation
of the spectra shows no second component above about 0.5 % of the continuum; the star is not double-lined. The Ca II K
line shows a central emission reversal on both nights, i.e. the K dwarf is chromospherically active and radial-velocity
jitter of several m/s is to be expected.

### 1.8 Catalogue and literature checks

All queries were run with a positive control (a known object returned by the same query).

| date | source | result |
|---|---|---|
| 2026-10-15 | VSX (15"), SIMBAD (10"; the star is TYC 9041-3677-1, type high-proper-motion star, 3 references) | no variability or planet entry |
| 2026-10-15 | Gaia DR3 `vari_*` tables (10) and `nss_*` tables | no entry |
| 2026-10-15 | ExoFOP TOI list (8151 rows), CTOI list (5138 rows), target page | no TOI, no CTOI |
| 2026-10-15 | NASA Exoplanet Archive (`toi`, `ps`, `pscomppars`, `td`) | no entry |
| 2026-10-15 | ADS full text, 17 identifier strings (TIC, Gaia DR2/DR3, TYC, 2MASS, UCAC4, WISE) | no publication on the star |
| 2026-10-15 | source files of 14 arXiv papers on TESS duo- and single-transit candidates | not listed |
| 2026-10-15 | VizieR, all tables within 5" (116 tables) | catalogue entries only; no planet or eclipse table |
| 2026-10-15 | TESS-SPOC FFI target lists; MAST data validation products | in the TESS-SPOC target list for two sectors; no SPOC threshold-crossing event |
| 2026-10-18 | ESO Science Archive (5" cone, programme pages) | five HARPS exposures (section 1.7); no other instrument |
| 2026-10-18 | ADS (identifier strings, repeat) | no publication |

Not checked: KELT and ASAS-SN light curves; catalogue tables that exist only outside the arXiv source files searched.

### 1.9 Predicted transits (2026-10-18 to 2027-10-18)

`data/predicted_transits_402710280.csv`. Errors are the propagated linear-ephemeris errors; the observed O-C of up to
29 min means a timing window of at least $`\pm 30`$ min should be added. $`T_{14} = 5.27`$ h.

| epoch | mid-transit (UTC) | BJD$`_\mathrm{TDB}`$ | ephemeris error (min) | TESS |
|---|---|---|---|---|
| 6 | 2026-11-02 11:12 | 2461346.9673 | 2.1 | |
| 7 | 2026-12-08 19:22 | 2461383.3073 | 2.1 | |
| 8 | 2027-01-14 03:32 | 2461419.6474 | 2.1 | |
| 9 | 2027-02-19 11:41 | 2461455.9874 | 2.2 | sector 114 |
| 10 | 2027-03-27 19:51 | 2461492.3275 | 2.2 | |
| 11 | 2027-05-03 04:01 | 2461528.6675 | 2.2 | |
| 12 | 2027-06-08 12:10 | 2461565.0076 | 2.3 | |
| 13 | 2027-07-14 20:20 | 2461601.3476 | 2.3 | |
| 14 | 2027-08-20 04:30 | 2461637.6877 | 2.4 | |
| 15 | 2027-09-25 12:39 | 2461674.0277 | 2.4 | |

TESS (tess-point 0.9.5): the star was on silicon in sectors 12, 39, 66, 93, 100, 101, 102 and 103 and will be in sectors
113 and 114; the nominal sector-113 window contains no transit, the sector-114 window contains epoch 9
(`data/tess_sectors_402710280.csv`).

## 2. TIC 324429952

### 2.1 Host star

| quantity | value | source |
|---|---|---|
| identifiers | TIC 324429952, Gaia DR3 6416596855786566528, CD-73 1398, TYC 9300-1959-1, 2MASS J19035007-7306494, RAVE J190350.1-730649 | SIMBAD, ExoFOP |
| position (ICRS, epoch 2016.0) | RA 285.9590554, Dec $`-`$73.1140172 | Gaia DR3 |
| proper motion (mas/yr) | 25.209, $`-`$66.783 | Gaia DR3 |
| parallax (mas) | $`9.2400 \pm 0.0139`$ | Gaia DR3 |
| distance (pc) | $`107.92 \pm 0.32`$ | TIC v8.2 |
| $`G`$, $`BP-RP`$ | 10.093, 0.855 | Gaia DR3 |
| $`T`$, $`V`$, $`B`$, $`J`$, $`H`$, $`K`$ | 9.637, 10.235, 10.988, 9.050, 8.760, 8.652 | TIC v8.2 (ExoFOP) |
| $`T_\mathrm{eff}`$ (K) | $`5737 \pm 132`$ (TIC); 5778 (Gaia GSP-Phot) | TIC v8.2; Gaia DR3 |
| $`R_\star`$ ($`R_\odot`$) | $`0.933 \pm 0.047`$ | TIC v8.2 |
| $`M_\star`$ ($`M_\odot`$) | $`1.02 \pm 0.13`$ | TIC v8.2 |
| $`\log g`$ | $`4.507 \pm 0.078`$ (TIC); 4.52 (Gaia GSP-Phot) | TIC v8.2; Gaia DR3 |
| $`\rho_\star`$ (g cm$`^{-3}`$) | $`1.77 \pm 0.39`$ | TIC v8.2 |
| RUWE; `ipd_frac_multi_peak`; `non_single_star` | 1.090; 0; 0 | Gaia DR3 |
| radial velocity (km/s) | $`-18.94`$ (24 transits); `rv_amplitude_robust` 3.38; `rv_chisq_pvalue` 0.047 | Gaia DR3 |
| variability | Gaia `phot_variable_flag` VARIABLE, classifier SOLAR_LIKE; VSX type ROT without period; rotation period 7.5 d in TESS rotation catalogues (VizieR J/ApJ/935/90: 7.48 d; J/ApJS/276/57: 7.52 d; J/A+A/709/A180: 7.52 d; J/AJ/167/189: 7.65 d; J/ApJS/281/52: 7.73 d) | Gaia DR3, VSX, VizieR |
| TIC contamination ratio | 0.011 | TIC v8.2 |
| common-proper-motion companion | Gaia DR3 6416597615998643200 at 43.6" (epoch 2016.0): $`G = 15.17`$, $`BP-RP = 2.93`$, parallax $`9.30 \pm 0.03`$ mas, proper motion (24.34, $`-`$66.26) mas/yr, $`T_\mathrm{eff}`$ 3289 K (GSP-Phot); shows flares in the cutouts | Gaia DR3 |

### 2.2 Transits and ephemeris

Common-shape fit of the eight transits (`scripts/fit_transits.py`): TESScut $`r = 2.0`$ px photometry for sectors 27 to
103 and the SPOC 2-min SAP light curve for sector 13; batman model with $`u = (0.40, 0.21)`$; quadratic baselines per
epoch over $`\pm 1`$ d. With all three shape parameters free the fit gives $`k = 0.0394 \pm 0.0017`$, $`b = 0.84 \pm 0.06`$,
$`a/R_\star = 50.1 \pm 9.0`$ ($`\rho_\star = 1.32`$ g cm$`^{-3}`$). The $`\chi^2`$ profile over $`b`$ is flat
($`\Delta\chi^2 \le 1.3`$) for $`b \le 0.85`$ and rises to 20 at $`b = 0.95`$ and 31 at $`b = 1.0`$; with $`a/R_\star`$
fixed to the TIC density, $`b = 0.80`$ is within $`\Delta\chi^2 = 1.5`$ of the free minimum while $`b \le 0.70`$ is
worse by $`\Delta\chi^2 \ge 59`$. The adopted shape is the density-constrained solution:

| quantity | value |
|---|---|
| $`P`$ (d) | $`42.37574 \pm 0.00008`$ |
| $`T_0`$ (BJD$`_\mathrm{TDB}`$) | $`2461174.88770 \pm 0.00290`$ (BTJD $`4174.88770`$) |
| $`k = R_\mathrm{p}/R_\star`$ | $`0.0386 \pm 0.0018`$ |
| depth $`k^2`$ | about 1490 ppm (per epoch 634-2215 ppm, table 2.3) |
| $`b`$ | $`0.80 \pm 0.09`$ |
| $`a/R_\star`$ | $`55.9 \pm 11.3`$ ($`\rho_\star = 1.84`$ g cm$`^{-3}`$, circular; TIC $`1.77 \pm 0.39`$) |
| $`T_{14}`$ (h) | 3.84 (model); detection-stage trapezoid $`4.32 \pm 0.15`$ h with $`T_{23} = 2.3`$ h |
| $`R_\mathrm{p} = k R_\star`$ ($`R_\oplus`$) | $`3.9 \pm 0.3`$ (3.6-4.4 across detrending choices; TRICERATOPS best fit 3.8) |
| $`a`$ (au) | 0.24 |
| detection-stage values (3 $`\times`$ 3-pixel box) | depth $`1433 \pm 51`$ ppm; $`P = 42.375818 \pm 0.000071`$ d; QLP detrended flux shows about 800 ppm |

### 2.3 Transit epochs and O-C

`data/transit_times_324429952.csv`; errors as in section 1.3.

| epoch | sector | product | $`t_c`$ (BTJD) | error (min) | O $`-`$ C (min) | depth (ppm) |
|---|---|---|---|---|---|---|
| $`-59`$ | 13 | SPOC 2-min SAP | 1674.722669 | 8.3 | $`+5.2`$ | $`1432 \pm 129`$ |
| $`-50`$ | 27 | TESScut | 2056.100882 | 4.1 | $`+0.3`$ | $`1553 \pm 111`$ |
| $`-25`$ | 66 | TESScut | 3115.491792 | 4.8 | $`-3.5`$ | $`1496 \pm 110`$ |
| $`-8`$ | 93 | TESScut | 3835.883771 | 18.8 | $`+2.9`$ | $`2215 \pm 146`$ |
| $`-7`$ | 94 | TESScut | 3878.255474 | 5.4 | $`-2.9`$ | $`1349 \pm 113`$ |
| $`-2`$ | 100 | TESScut | 4090.139084 | 9.7 | $`+4.1`$ | $`634 \pm 142`$ |
| $`-1`$ | 102 | TESScut | 4132.517034 | 20.5 | $`+7.3`$ | $`1649 \pm 148`$ |
| 0 | 103 | TESScut | 4174.897401 | 11.1 | $`+14.0`$ | $`1620 \pm 138`$ |

Linear ephemeris: $`\chi^2 = 3.1`$ for 6 degrees of freedom; no transit-timing variation is detected. The per-epoch
depths of sectors 93 and 100 differ from the others; both sectors have the largest red-noise factors ($`\beta`$ = 2.4
and 2.0). The predicted transits of sectors 39, 67 and 101 fall outside the data (`data/sectors_324429952.csv`).

<picture><source srcset="figures/fold_324429952_dark.svg" media="(prefers-color-scheme: dark)"><img src="figures/fold_324429952_light.svg" alt="TIC 324429952 folded transit"></picture>

<picture><source srcset="figures/epochs_324429952_dark.svg" media="(prefers-color-scheme: dark)"><img src="figures/epochs_324429952_light.svg" alt="TIC 324429952 individual transits"></picture>

<picture><source srcset="figures/oc_324429952_dark.svg" media="(prefers-color-scheme: dark)"><img src="figures/oc_324429952_light.svg" alt="TIC 324429952 O-C"></picture>

### 2.4 Checks against false positives

- The sector-13 transit is in the SPOC 2-min light curve (SAP and PDCSAP, 1429-1516 ppm, 13 $`\sigma`$) with no
  in-transit shift of the SPOC flux-weighted centroid; all eight epochs are present in independent photometry (QLP SAP
  for every sector; TESS-SPOC SAP and PDCSAP for sectors 13, 27 and 66; eleanor-lite and TARS where available; the
  TESScut circular apertures of this record). No SPOC threshold-crossing event exists.
- Odd and even transits agree (detection stage: $`1421 \pm 63`$ vs $`1445 \pm 63`$ ppm; 0.4 $`\sigma`$).
- Secondary eclipse: $`< 500`$ ppm at phase 0.5.
- Period aliases $`P/k`$, $`k = 2`$-12, are excluded by covered epochs without transit (for $`P/2`$: sector 94 at BTJD
  3857.08, $`-403 \pm 136`$ ppm; sector 67 at 3136.71, $`287 \pm 127`$ ppm; sector 101 at 4111.33, flat).
- Localisation: joint PSF fits to the difference images of six sectors put the dimming on the target at 7-10 $`\sigma`$
  and not on the common-proper-motion M dwarf at 43.6", which would need eclipses of about 35 % and shows none in its
  own pixels (section 2.5).
- Gaia DR3: RUWE 1.09, no multi-peak detections, no non-single-star solution; the robust radial-velocity amplitude of
  3.4 km/s over 24 transits limits a companion at this period to below about 30 $`M_\mathrm{J}`$.
- The transit depth does not depend on the phase of the 7.5-d rotational modulation; all eight transits fall near
  minima of the modulation (2.5 $`\sigma`$ against a random distribution), which is noted and not explained.

### 2.5 Difference-image localisation

`data/diffimage_324429952.csv`, `data/diffimage_sources_324429952.csv`, `data/diffimage_centroids_324429952.csv`. The
difference-image centroid offsets from the out-of-transit centroid, weighted over the eight sectors, are
$`\Delta\mathrm{E} = -2.7 \pm 1.7`$" and $`\Delta\mathrm{N} = +5.5 \pm 3.0`$" ($`\chi^2 = 22.5`$ for 16 degrees of freedom
against zero); the M dwarf lies 43.6" away. The joint PSF fit in the six sectors with 200-s or 10-min cadence assigns
the dimming to the target at 7-10 $`\sigma`$.

<picture><source srcset="figures/diffimage_324429952_dark.svg" media="(prefers-color-scheme: dark)"><img src="figures/diffimage_324429952_light.svg" alt="TIC 324429952 difference images"></picture>

### 2.6 TRICERATOPS

Same set-up as section 1.6 with all eight epochs (sector 13: SPOC aperture and 2-min SAP flux), baselines over
$`\pm 0.6`$ d (nominal), 144 bins, out-of-transit bin scatter 178 ppm, lag-1 autocorrelation of the out-of-transit bins
0.21 (0.72 with $`\pm 1.0`$ d baselines).

| configuration | seeds | FPP mean $`\pm`$ std | FPP median (range) | NFPP |
|---|---|---|---|---|
| nominal ($`\pm 0.6`$ d baselines) | 10 | $`0.033 \pm 0.015`$ | 0.036 (0.011 to 0.061) | $`3 \times 10^{-7}`$ |
| $`\pm 1.0`$ d baselines | 5 | $`0.095 \pm 0.025`$ | 0.105 (0.060 to 0.12) | $`8 \times 10^{-5}`$ |
| flux errors $`\times 1.5`$ | 5 | $`0.073 \pm 0.018`$ | 0.073 (0.049 to 0.099) | $`2 \times 10^{-5}`$ |
| sectors 93 and 100 left out | 5 | $`0.052 \pm 0.019`$ | 0.055 (0.027 to 0.079) | $`3 \times 10^{-6}`$ |

Mean scenario probabilities of the nominal runs: TP 0.805, PTP 0.156, STP 0.0135, SEB 0.0193, DTP 0.0067; scenarios on
the M dwarf (NTP, NEB) are at or below $`10^{-4}`$. The FPP is in the range Giacalone et al. (2021) call "likely planet"
(FPP $`< 0.5`$, NFPP $`< 0.001`$) and above their threshold for statistical validation (0.015). The companion scenarios
STP and SEB carry the FPP, which the V-shaped fold ($`b \approx 0.8`$) leaves room for. Reference runs with an assumed
contrast curve (FPP $`0.0089 \pm 0.004`$) are in `data/triceratops/`, flagged as hypothetical.

### 2.7 Catalogue and literature checks

| date | source | result |
|---|---|---|
| 2026-10-15 | VSX (60") | the star, type ROT, no period |
| 2026-10-15 | SIMBAD (30") | CD-73 1398, high-proper-motion star, 4 references; no planet or eclipse |
| 2026-10-15 | Gaia DR3 `vari_*` and `nss_*` tables | `vari_classifier_result` SOLAR_LIKE; no eclipse, transit or NSS entry |
| 2026-10-14/15 | ExoFOP TOI (8151 rows) and CTOI (5138 rows) lists, target page | no TOI, no CTOI |
| 2026-10-15 | NASA Exoplanet Archive (4 tables) | no entry |
| 2026-10-15 | ADS full text, 27 name strings | no publication on the star beyond catalogue papers |
| 2026-10-15 | VizieR, all tables within 5" (151 tables); Prsa et al. 2022 TESS eclipsing binaries (J/ApJS/258/16); Kostov et al. 2025 (J/ApJS/279/50) | rotation-period entries only; no eclipse or planet entry |
| 2026-10-15 | source files of 20 arXiv transit-candidate papers; RAVEN and T16 candidate tables; MAST data-validation products | not listed; no SPOC threshold-crossing event |
| 2026-10-18 | ESO Science Archive (5" cone) | no data |

### 2.8 Predicted transits (2026-10-18 to 2027-10-18)

`data/predicted_transits_324429952.csv`; $`T_{14}`$ = 3.8-4.3 h.

| epoch | mid-transit (UTC) | BJD$`_\mathrm{TDB}`$ | ephemeris error (min) | TESS |
|---|---|---|---|---|
| 4 | 2026-10-30 21:22 | 2461344.3907 | 4.2 | |
| 5 | 2026-12-12 06:23 | 2461386.7664 | 4.2 | |
| 6 | 2027-01-23 15:24 | 2461429.1421 | 4.2 | |
| 7 | 2027-03-07 00:25 | 2461471.5179 | 4.3 | sector 114 |
| 8 | 2027-04-18 09:26 | 2461513.8936 | 4.3 | |
| 9 | 2027-05-30 18:27 | 2461556.2694 | 4.3 | |
| 10 | 2027-07-12 03:28 | 2461598.6451 | 4.3 | |
| 11 | 2027-08-23 12:30 | 2461641.0208 | 4.4 | |
| 12 | 2027-10-04 21:31 | 2461683.3966 | 4.4 | |

TESS: on silicon in sectors 13, 27, 39, 66, 67, 93, 94, 100, 101, 102 and 103, and in sector 114 (epoch 7).

## 3. Caveats

- No high-resolution imaging exists for either star. For TIC 402710280 the 2.3" neighbour is unresolved by TESS and
  excluded as host only by the transit shape; a companion closer than that is unconstrained. For TIC 324429952 any
  companion within about 1" is unconstrained, and the TRICERATOPS FPP is dominated by bound-companion scenarios.
- No radial velocities were obtained for this work. For TIC 402710280 the two public HARPS epochs sample 16 % of the
  orbital velocity amplitude and exclude only stellar and brown-dwarf companions; three further epochs are proprietary.
  For TIC 324429952 there are no archival spectra.
- QLP timing notice: MAST reports a barycentric timing error in the QLP light curves from the second half of sector 99
  to sector 104 (orbits 210-220), not yet corrected at the date of this record. No QLP time stamp enters the tables; all
  mid-transit times come from TESScut cutouts (SPOC-calibrated full-frame images) or SPOC light curves. For
  TIC 402710280 the QLP-minus-TESScut time-stamp differences measured in sectors 100 and 102 are within $`\pm 1.7`$ min
  (medians $`+0.3`$ and $`-0.3`$ min), far below the 24-29 min deviations of section 1.3.
- The 2.3" neighbour of TIC 402710280 (Gaia DR3 5828441289812582656; TIC 1208600633, without TIC stellar parameters) is
  6.4 mag fainter and lies in the background. TRICERATOPS treats it with solar parameters by default; the run with
  K-dwarf parameters gives the same FPP within the seed-to-seed scatter.
- Stellar parameters are TIC v8.2 values and were not re-derived; the planet radii scale with the TIC stellar radii.
- The per-epoch depths of TIC 324429952 scatter by about 20 % between detrending choices and sectors; the red-noise
  factors reach 2.4. TRICERATOPS uses one flux-error value per light curve and a Gaia DR3 rather than TRILEGAL background.
- The transit-timing variations of TIC 402710280 rest on four epochs with errors of 3-5 min; two further epochs
  (2027-02-19 in TESS sector 114, and any ground-based transit) would test them.
- Nothing in this record has been checked by anyone outside the production of the record.

## 4. Data and reproduction

Public products used (all retrieved 2026-10-15 to 2026-10-18):

- TESScut 21 $`\times`$ 21-pixel cutouts of the SPOC-calibrated full-frame images (MAST), sectors listed in
  `data/sectors_<tic>.csv`, at the TIC v8.2 coordinates.
- MAST light curves for the product-by-product checks: QLP (`hlsp_qlp_tess_ffi_s<sector>-<tic>_tess_v01_llc`), TESS-SPOC
  (`hlsp_tess-spoc_tess_phot_<tic>-s<sector>_tess_v1`), TARS (`hlsp_tars_tess_ffi_s<sector>-<tic>_tess`), GSFC-eleanor-lite
  (sector 13 of TIC 324429952), and the SPOC 2-min light curve of TIC 324429952 in sector 13 (`tess*-s0013-0000000324429952-*_lc.fits`).
- Gaia DR3 (`gaia_source`, `vari_*`, `nss_*`), TIC v8.2 (via ExoFOP and MAST), ExoFOP TOI/CTOI lists, NASA Exoplanet
  Archive, VSX, SIMBAD, VizieR, ADS.
- ESO Science Archive: HARPS ADP products ADP.2025-08-26T01:00:42.037/.038 and ADP.2025-08-27T01:04:15.751/.752
  (checksums in `data/harps_public_rvs.csv`).

Software: Python 3.12.11 with the versions in `requirements.txt`; TRICERATOPS 1.1.0 in its own environment (versions in
`requirements.txt`). Commands, run from the repository root (work files go to `./work/`):

```
python scripts/tesscut_download.py 402710280          # cutouts (MAST)
python scripts/aperture_photometry.py 402710280       # r = 1.0-3.0 px light curves
python scripts/fit_transits.py 402710280              # common-shape fit, timing errors, ephemeris -> work/402710280/fit_r2.0_39_93_100_102.json
python scripts/fit_transits.py 402710280 --epochs 39,66,93,100,102   # the sector-66 event
python scripts/difference_images.py 402710280 --fit fit_r2.0_39_93_100_102.json
python scripts/alias_secondary_search.py 402710280
python scripts/make_tables.py 402710280 --fit fit_r2.0_39_93_100_102.json --extra-fit fit_r2.0_39_66_93_100_102.json --extra-sector 66
python scripts/predict_transits.py 402710280 --start 2026-10-18
python scripts/triceratops_prep.py 402710280 --fit fit_r2.0_39_93_100_102.json --window 1.0 --tag _J1.0
python scripts/triceratops_target.py 402710280 --sectors 39,93,100,102 --depth 0.00168     # TRICERATOPS environment
for s in 1 2 3 4 5 6 7 8 9 10; do python scripts/triceratops_run.py 402710280 $s --tag _J1.0; done
python scripts/triceratops_summary.py
```

For TIC 324429952 use `--fit fit_r2.0_13_27_66_93_94_100_102_103_fix.json` (the shape is fixed through `targets.py`),
`--window 0.6 --tag _J0.6` for the TRICERATOPS fold, `--sectors 13,27,66,93,94,100,102,103 --depth 0.00149`, and the
SPOC 2-min light curve of sector 13 downloaded into `work/324429952/dl/` (MAST `Observations.download_products`). The
figures are regenerated with `python figures/make_figures.py` from the CSV files alone. `MANIFEST.sha256` lists the
SHA-256 checksum of every file (`shasum -a 256 -c MANIFEST.sha256`).

## 5. How to cite

A. Zarco (2026). TESS transit candidates 2026: TIC 402710280 and TIC 324429952. Zenodo. DOI: 10.5281/zenodo.XXXXXXX
(to be added after the deposit). Citation metadata: `CITATION.cff`.

## 6. Licence

Data, figures and text: CC BY 4.0. Scripts: MIT. See `LICENSE`, which also lists the third-party material not covered.

## References

- Giacalone, S., Dressing, C. D., et al. 2021, AJ 161, 24 (TRICERATOPS).
- Kreidberg, L. 2015, PASP 127, 1161 (batman).
- Hippke, M., et al. 2019, AJ 158, 143 (wotan).
- Huang, C. X., et al. 2020, RNAAS 4, 204 (QLP); Caldwell, D. A., et al. 2020, RNAAS 4, 201 (TESS-SPOC); Brasseur, C. E.,
  et al. 2019, ascl:1905.007 (TESScut); Burke, C. J., et al. 2020, ascl:2005.003 (tess-point).
- Stassun, K. G., et al. 2019, AJ 158, 138 (TIC v8); Gaia Collaboration 2023, A&A 674, A1 (Gaia DR3).
- Prsa, A., et al. 2022, ApJS 258, 16; Kostov, V. B., et al. 2025, ApJS 279, 50.
