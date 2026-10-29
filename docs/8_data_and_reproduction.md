[README](../README.md) / 8. Data and reproduction

## 8. Data and reproduction

Public products used (all retrieved 2026-10-15 to 2026-10-19):

- TESScut 21 $`\times`$ 21-pixel cutouts of the SPOC-calibrated full-frame images (MAST), sectors listed in
  `data/sectors_<tic>.csv`, at the TIC v8.2 coordinates (TIC 89247735, TIC 455017950 and TIC 311179599: at the Gaia DR3
  position).
- MAST light curves for the product-by-product checks: QLP (`hlsp_qlp_tess_ffi_s<sector>-<tic>_tess_v01_llc`), TESS-SPOC
  (`hlsp_tess-spoc_tess_phot_<tic>-s<sector>_tess_v1`), TARS (`hlsp_tars_tess_ffi_s<sector>-<tic>_tess`), GSFC-eleanor-lite
  (sector 13 of TIC 324429952 and TIC 89247735; sector 12 of TIC 455017950), the SPOC 2-min light curve of TIC 324429952 in sector 13 (`tess*-s0013-0000000324429952-*_lc.fits`)
  and those of TIC 455017950 in sectors 12 and 65; for TIC 311179599 the QLP light curves of all nine sectors, the
  TESS-SPOC light curves of sectors 12, 13, 39 and 66 and the GSFC-eleanor-lite light curves of sectors 12 and 13.
- Gaia DR3 (`gaia_source`, `vari_*`, `nss_*`), TIC v8.2 (via ExoFOP and MAST), ExoFOP TOI/CTOI lists, NASA Exoplanet
  Archive, VSX, SIMBAD, VizieR, ADS.
- SVO Theoretical Spectra service (Castelli & Kurucz 2003 ODFNEW spectra) and SVO Filter Profile Service (transmission
  curves, zero points), for the SED fits of TIC 89247735 and TIC 455017950.
- Published SOAR speckle contrast limits of TIC 455017950 (Tokovinin et al. 2020, VizieR J/AJ/160/7, table 3), as
  [`data/triceratops/cc_HIP74344_SOAR2019_I.txt`](../data/triceratops/cc_HIP74344_SOAR2019_I.txt).
- ESO Science Archive: HARPS ADP products ADP.2025-08-26T01:00:42.037/.038 and ADP.2025-08-27T01:04:15.751/.752
  (checksums in [`data/harps_public_rvs.csv`](../data/harps_public_rvs.csv)).

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

For TIC 89247735 (epochs are labelled 66, 93a, 93b and 104 in the work files):

```
python scripts/tesscut_download.py 89247735
python scripts/aperture_photometry.py 89247735
python scripts/fit_transits.py 89247735 --window 0.8            # -> work/89247735/fit_r2.0_66_93a_93b_104_w0.8.json
python scripts/fit_transits.py 89247735 --window 0.8 --aperture 1.5   # and --aperture 3.0
python scripts/fit_transits.py 89247735 --window 0.8 --bfix 0.7 --sig-from fit_r2.0_66_93a_93b_104_w0.8.json   # b profile; --arfix 29.4
python scripts/fit_transits.py 89247735 --window 0.6 --order 1
python scripts/difference_images.py 89247735 --fit fit_r2.0_66_93a_93b_104_w0.8.json   # also queries Gaia DR3 (60")
python scripts/dilution_correction.py 89247735 --fits fit_r1.5_66_93a_93b_104_w0.8,fit_r2.0_66_93a_93b_104_w0.8,fit_r3.0_66_93a_93b_104_w0.8
python scripts/gaia_companion_limits.py 89247735 --shape fit_r2.0_66_93a_93b_104_w0.8_shape.json
python scripts/alias_secondary_search.py 89247735 --harmonics
python scripts/event_checks.py 89247735 --events 13:1682.134,66:3111.729,93:3836.585,93:3843.075,93:3849.745,104:4185.066
python scripts/make_tables.py 89247735 --fit fit_r2.0_66_93a_93b_104_w0.8.json --window 0.8
python scripts/predict_transits.py 89247735 --start 2026-10-18
python scripts/triceratops_prep.py 89247735 --fit fit_r2.0_66_93a_93b_104_w0.8.json --window 0.8 --tag _J0.8
python scripts/triceratops_target.py 89247735 --sectors 66,93,104 --depth 0.00826            # TRICERATOPS environment
for s in 1 2 3; do python scripts/triceratops_run.py 89247735 $s --tag _J0.8; done
python scripts/sed_fit.py                                       # SVO spectra and filters -> work/89247735/sed/
python scripts/fit_free_dilution.py 89247735 --fit fit_r2.0_66_93a_93b_104_w0.8.json   # emcee
```

`event_checks.py` uses the QLP light curves of the star in `work/89247735/dl/` (MAST `Observations.download_products`)
when present. Exceptions for TIC 89247735: the cutouts of sectors 66, 93 and 104 used for the tables were retrieved
from TESScut on 2026-10-18 with the same request as `tesscut_download.py` (Gaia DR3 position, 21 $`\times`$ 21
pixels); the re-download of these three sectors timed out at MAST on 2026-10-18 (sector 13 was re-downloaded and is
identical). [`data/triceratops_input_89247735.csv`](../data/triceratops_input_89247735.csv) is the binned light curve that was given to TRICERATOPS; it was built
from the common-shape fit itself, and `triceratops_prep.py`, which refits the shape with the mid-transit times fixed,
reproduces it to within 9 ppm per bin (bin scatter 542.0 instead of 541.7 ppm). The TRICERATOPS tables
(`data/triceratops/probs_89247735_*.csv`, `gaia_background_89247735.csv`) are the outputs of those runs, copied
verbatim. The product-by-product box depths of section [3.4](3_TIC89247735.md#34-checks-against-false-positives) (QLP, TARS) and the catalogue checks of section [3.8](3_TIC89247735.md#38-catalogue-and-literature-checks) have no
script in this repository. [`scripts/sed_fit.py`](../scripts/sed_fit.py) holds the catalogue magnitudes it fits (copied from VizieR: Gaia DR3,
SkyMapper DR4, APASS DR9, 2MASS, AllWISE); it does not query them. The other $`T_\mathrm{eff}`$ and [Fe/H] values,
the Mann et al. (2015, 2019) and Kervella et al. (2004) relations, the $`M_G`$-based luminosity and the adopted stellar
parameters of section [3.1](3_TIC89247735.md#31-host-star), the search for additional transiting signals and its injection-recovery test (section [3.10](3_TIC89247735.md#310-search-for-additional-transiting-signals))
and the archival photometry of section [3.11](3_TIC89247735.md#311-archival-ground-based-photometry) have no script in this repository. [`scripts/fit_free_dilution.py`](../scripts/fit_free_dilution.py) was run
on the same cutouts of sectors 66, 93 and 104 as the tables, after `aperture_photometry.py` and `fit_transits.py
--window 0.8`, in a Python 3.11 environment with emcee (versions in `requirements.txt`).

For TIC 455017950:

```
python scripts/tesscut_download.py 455017950
python scripts/aperture_photometry.py 455017950
python scripts/fit_transits.py 455017950 --window 0.8          # -> work/455017950/fit_r2.0_38_65_99_100_102_w0.8.json
python scripts/fit_transits.py 455017950 --window 0.8 --aperture 1.5   # and --aperture 3.0
python scripts/fit_transits.py 455017950 --window 0.8 --bfix 0.5 --sig-from fit_r2.0_38_65_99_100_102_w0.8.json   # b profile; --arfix 60.6
python scripts/fit_transits.py 455017950 --window 0.6 --order 3  # and --window 1.0 --order 3; --epochs 38,65,100,102
python scripts/difference_images.py 455017950 --fit fit_r2.0_38_65_99_100_102_w0.8.json   # also queries Gaia DR3 (60")
python scripts/dilution_correction.py 455017950 --fits fit_r1.5_38_65_99_100_102_w0.8,fit_r2.0_38_65_99_100_102_w0.8,fit_r3.0_38_65_99_100_102_w0.8
python scripts/gaia_companion_limits.py 455017950 --shape fit_r2.0_38_65_99_100_102_w0.8_shape.json
python scripts/alias_secondary_search.py 455017950
python scripts/event_checks.py 455017950 --events 38:2346.001,65:3089.995,99:4073.054,101:4119.794,102:4133.134,102:4138.374
python scripts/make_tables.py 455017950 --fit fit_r2.0_38_65_99_100_102_w0.8.json --window 0.8
python scripts/predict_transits.py 455017950 --start 2026-10-18
python scripts/triceratops_prep.py 455017950 --fit fit_r2.0_38_65_99_100_102_w0.8.json --window 0.8 --tag _J0.8
python scripts/triceratops_target.py 455017950 --sectors 38,65,99,100,102 --depth 0.001207     # TRICERATOPS environment
for s in 1 2 3 4 5; do python scripts/triceratops_run.py 455017950 $s --tag _J0.8 --label SOARcc --cc data/triceratops/cc_HIP74344_SOAR2019_I.txt --cc-measured; done
for s in 1 2; do python scripts/triceratops_run.py 455017950 $s --tag _J0.8 --label nocc; done
python scripts/sed_fit.py --tic 455017950                       # SVO spectra and filters -> work/455017950/sed/
python scripts/fit_free_dilution.py 455017950 --fit fit_r2.0_38_65_99_100_102_w0.8.json --walkers 44 --steps 6000   # emcee
python scripts/fit_free_dilution.py 455017950 --fit fit_r2.0_38_65_99_100_102_w0.8.json --walkers 44 --steps 6000 --rho-prior
```

`event_checks.py` uses the QLP light curves of the star in `work/455017950/dl/` (MAST `Observations.download_products`)
when present. Exceptions for TIC 455017950: the eight cutouts were re-downloaded from TESScut on 2026-10-18 with
`tesscut_download.py`; their pixel data, times, quality flags and WCS are identical to those of the cutouts first used
(only the request coordinates in the header differ in the 13th digit), and the tables regenerated from them are
identical. [`data/triceratops_input_455017950.csv`](../data/triceratops_input_455017950.csv) is the binned light curve that was given to TRICERATOPS; it was built
from a common-shape fit of the same photometry ($`b = 0.16`$, $`a/R_\star = 67.3`$, $`k^2 = 1207`$ ppm), and
`triceratops_prep.py` reproduces it to within 6 ppm per bin (bin scatter 180.7 instead of 180.2 ppm). The TRICERATOPS
tables (`data/triceratops/probs_455017950_*.csv`, `gaia_background_455017950.csv`) are the outputs of those runs, copied
verbatim; their period (40.605653 d) is that of the first fit, 0.1 $`\sigma`$ from the adopted one. Without script in
this repository: the product-by-product box depths and the SPOC 2-min fit (section [4.4](4_TIC455017950.md#44-checks-against-false-positives)), the diluted-eclipse fits
(section [4.2](4_TIC455017950.md#42-transits-and-ephemeris)), the robust-spline detrending, the red-noise and phase-shift tests (section [4.4](4_TIC455017950.md#44-checks-against-false-positives)), the transit-timing
analysis (section [4.3](4_TIC455017950.md#43-transit-epochs-and-o-c)), the simulated Gaia radial-velocity and astrometric limits and the companion grid (section [4.6](4_TIC455017950.md#46-companion-limits)),
the rotation periods, gyrochronology and kinematics and the Monte Carlo SED fit with free extinction (section [4.1](4_TIC455017950.md#41-host-star)), the
catalogue checks (section [4.8](4_TIC455017950.md#48-catalogue-and-literature-checks)), the ground-based observability (section [4.9](4_TIC455017950.md#49-predicted-transits-2026-10-18-to-2027-10-18)), the search for additional transiting
signals and its injection-recovery test (section [4.10](4_TIC455017950.md#410-search-for-additional-transiting-signals)) and the archival photometry (section [4.11](4_TIC455017950.md#411-archival-photometry)). [`scripts/sed_fit.py`](../scripts/sed_fit.py)
holds the catalogue magnitudes it fits (copied from VizieR: Gaia DR3, Tycho-2, SkyMapper DR4, 2MASS, AllWISE).
[`scripts/fit_free_dilution.py`](../scripts/fit_free_dilution.py) was run in the Python 3.11 environment with emcee (versions in `requirements.txt`).

For TIC 311179599:

```
python scripts/tesscut_download.py 311179599
python scripts/aperture_photometry.py 311179599
python scripts/fit_transits.py 311179599 --window 0.6          # -> work/311179599/fit_r2.0_12_13_39_66_93_100_102_103_w0.6.json
python scripts/fit_transits.py 311179599 --window 0.6 --aperture 1.5   # and --aperture 3.0
python scripts/fit_transits.py 311179599 --window 0.4 --order 1        # and --window 1.0 --order 3
python scripts/fit_transits.py 311179599 --window 0.6 --bfix 0.5 --sig-from fit_r2.0_12_13_39_66_93_100_102_103_w0.6.json   # b profile; --arfix 50.1
python scripts/difference_images.py 311179599 --fit fit_r2.0_12_13_39_66_93_100_102_103_w0.6.json   # also queries Gaia DR3 (60")
python scripts/dilution_correction.py 311179599 --fits fit_r1.5_12_13_39_66_93_100_102_103_w0.6,fit_r2.0_12_13_39_66_93_100_102_103_w0.6,fit_r3.0_12_13_39_66_93_100_102_103_w0.6
python scripts/gaia_companion_limits.py 311179599 --shape fit_r2.0_12_13_39_66_93_100_102_103_w0.6_shape.json
python scripts/alias_secondary_search.py 311179599
python scripts/event_checks.py 311179599 --events 13:1654.542,66:3104.798,66:3111.758,66:3122.768,93:3841.907,93:3855.227,100:4086.307,100:4088.457,100:4092.877,100:4100.147,101:4119.823,101:4126.103,102:4138.398,102:4151.308,103:4152.951,103:4158.671,103:4162.411,103:4164.801,103:4166.501,103:4170.941,103:4177.741
python scripts/make_tables.py 311179599 --fit fit_r2.0_12_13_39_66_93_100_102_103_w0.6.json --window 0.6
python scripts/predict_transits.py 311179599 --start 2026-10-19
python scripts/triceratops_prep.py 311179599 --fit fit_r2.0_12_13_39_66_93_100_102_103_w0.6.json --window 0.6 --drop 12,13 --tag _J0.6
python scripts/triceratops_target.py 311179599 --sectors 39,66,93,100,102,103 --depth 0.001279     # TRICERATOPS environment; see below
for s in 1 2 3 4 5; do python scripts/triceratops_run.py 311179599 $s --tag _J0.6; done
python scripts/fit_free_dilution.py 311179599 --fit fit_r2.0_12_13_39_66_93_100_102_103_w0.6.json --window 0.6 --walkers 44 --steps 6000 --rho-prior   # emcee
python scripts/fit_free_dilution.py 311179599 --fit fit_r2.0_12_13_39_66_93_100_102_103_w0.6.json --window 0.6 --walkers 44 --steps 6000
```

Tables written: `make_tables.py` writes [`data/transit_times_311179599.csv`](../data/transit_times_311179599.csv), `folded_311179599.csv`,
`model_311179599.csv`, `sectors_311179599.csv`, `diffimage_311179599.csv`, `diffimage_sources_311179599.csv` and
`diffimage_centroids_311179599.csv`; `predict_transits.py` writes [`data/predicted_transits_311179599.csv`](../data/predicted_transits_311179599.csv) and
`tess_sectors_311179599.csv`; `fit_free_dilution.py` writes [`data/fit_free_dilution_311179599.csv`](../data/fit_free_dilution_311179599.csv) and
`fit_free_dilution_rho_311179599.csv`; `triceratops_run.py` writes [`data/triceratops/probs_311179599_nominal_seed1.csv`](../data/triceratops/probs_311179599_nominal_seed1.csv)
to `_seed5.csv` (summarised in [`data/triceratops/runs.csv`](../data/triceratops/runs.csv) and `summary.csv`); [`data/triceratops_input_311179599.csv`](../data/triceratops_input_311179599.csv)
and [`data/triceratops/gaia_background_311179599.csv`](../data/triceratops/gaia_background_311179599.csv) are the light curve and the background population of the
TRICERATOPS runs. The cutouts of all nine sectors were downloaded with `tesscut_download.py` on 2026-10-19. The sector-39 transit is
fitted with the common shape but left out of the adopted ephemeris (`ephem_exclude` in [`scripts/targets.py`](../scripts/targets.py));
`fit_transits.py` then also writes the ephemeris through all epochs, and `make_tables.py` adds the column
`oc_all_epochs_min`. `event_checks.py` uses the QLP light curves of the star in `work/311179599/dl/mastDownload/HLSP/`
(MAST HLSP files `hlsp_qlp_tess_ffi_s<sector>-0000000311179599_tess_v01_llc.fits`, retrieved 2026-10-19). [`data/triceratops_input_311179599.csv`](../data/triceratops_input_311179599.csv) is the light curve written by `triceratops_prep.py` before the correction of 2026-10-22 (end of this section). The
TRICERATOPS target object was not built by `triceratops_target.py`: on 2026-10-19 the TESScut queries that TRICERATOPS
makes through lightkurve failed for two or three of the six sectors in each of two attempts (MAST "Request failed"),
and the script stops when an image is missing. The object used (`work/311179599/triceratops_target.pkl`) was built with
the same TRICERATOPS 1.1.0 calls (field-star table, Gaia DR3 background population, $`r < 2.0`$ px apertures of 11-13
pixels) from cutouts of the public TESS full-frame-image cubes, which are pixel-identical to TESScut; its
`calc_depths` was recomputed with the depth of the command above, and `triceratops_run.py` was then run unchanged. The
TRICERATOPS tables (`data/triceratops/probs_311179599_*.csv`, `gaia_background_311179599.csv`) are the outputs of
these runs, copied verbatim.
Without script in this repository: the Mann et al. (2015, 2019) relations and the rotation search (section [5.1](5_TIC311179599.md#51-host-star)), the PSF
photometry, the timing of the MAST light curves, the injection-recovery timing tests, the red-noise statistics of
sector 39, the sector-101 constraint, the sinusoidal timing fit and the perturber scaling (section [5.3](5_TIC311179599.md#53-transit-epochs-and-o-c)), the independent
Metropolis fit (section [5.2](5_TIC311179599.md#52-transits-and-ephemeris)), the presence of the transits in the MAST light curves (section [5.4](5_TIC311179599.md#54-checks-against-false-positives)), the TESS-band flux
ratios of the faint Gaia neighbours (section [5.5](5_TIC311179599.md#55-difference-image-localisation)), the derived insolation and equilibrium temperature (section [5.2](5_TIC311179599.md#52-transits-and-ephemeris)) and
the catalogue checks (section [5.8](5_TIC311179599.md#58-catalogue-and-literature-checks)). [`scripts/fit_free_dilution.py`](../scripts/fit_free_dilution.py) was run in the Python 3.11 environment with emcee
(versions in `requirements.txt`).

For the candidate second signal of TIC 311179599 (section [5.10](5_TIC311179599.md#510-a-second-transit-signal-candidate-c)):

```
python scripts/tesscut_download.py 311179599 --sectors 12,100,103   # the three sectors with events
python scripts/aperture_photometry.py 311179599
python scripts/fit_second_signal.py                    # joint fit, r = 2.0 px (about 8 min) -> data/*_311179599c.csv
python scripts/fit_second_signal.py --aperture 1.5     # -> work/311179599/fit_c_r1.5.json
python scripts/predict_transits.py 311179599 --start 2026-10-19 --signal c
python3 figures/make_figures.py c311179599
```

Tables written: `fit_second_signal.py` writes [`data/transit_times_311179599c.csv`](../data/transit_times_311179599c.csv), `fit_311179599c.csv`,
`folded_311179599c.csv` and `model_311179599c.csv`; `predict_transits.py --signal c` writes
[`data/predicted_transits_311179599c.csv`](../data/predicted_transits_311179599c.csv) (ephemerides in `second` of [`scripts/targets.py`](../scripts/targets.py)). The cutouts were downloaded on
2026-10-19. The TESScut time stamps are used as delivered, as for the other tables; the vetting record applied a
target-specific barycentric correction (-31 s, +44 s and +15 s in sectors 12, 100 and 103) and a 60000-step chain, and
gave $`P = 73.86712 \pm 0.00045`$ d, O $`-`$ C of $`-18`$ and $`+23`$ min for A and B, depth 874 ppm and $`T_{14} = 5.25`$ h
($`r = 2.0`$ px), within the errors of the values here. Exported from the vetting record of 2026-10-19 (no script in
this repository): [`data/period_tests_311179599c.csv`](../data/period_tests_311179599c.csv), `epoch_tests_311179599c.csv`, `completeness_311179599c.csv` and
`triceratops_input_311179599c.csv`, the TRICERATOPS tables `data/triceratops/probs_311179599c_n*_seed*.csv` (copied
verbatim; runs with `calc_depths` redone on the target object of section [5.7](5_TIC311179599.md#57-triceratops)) and their rows in `runs.csv` and
`summary.csv`, the PSF-photometry fit, the light-curve, background, centroid, difference-image, solar-system-object,
comparison-star and event-statistics checks (section [5.10.3](5_TIC311179599.md#5103-checks)), the timing comparison (section [5.10.5](5_TIC311179599.md#5105-transit-timing-of-candidate-b-and-the-second-signal)) and the searches
for further signals (section [5.10.6](5_TIC311179599.md#5106-search-for-further-signals-and-completeness)); the check of the 52.26-d peak has no script.

For the other TESS results (section [6](6_other_results.md)), with the MAST products listed below and the work files in `./work/`:

```
python scripts/other_download.py 11613676                # TESScut cutouts of sectors 3-106, SPOC 2-min light curves of sectors 3 and 30
python scripts/eb_eclipses.py                           # photometry, eclipse fits, ephemeris, difference images, predictions
python scripts/eb_model.py                              # geometric models M0, MA, N4 (about 20 min)
python scripts/other_download.py 124751008             # TESScut cutouts of sectors 6-109 (about 170 MB)
python scripts/eb_124751008.py --start 2026-10-19       # bleed-masked photometry, eclipse fits, ephemeris, secondary limits, predictions
python scripts/other_download.py 200322593 --controls   # SPOC 2-min TPFs and light curves, 20-s light curves; control stars
python scripts/toi540_analysis.py                       # photometry, ephemeris, sector-98 tests, injections, b(t), windows
python scripts/toi540_controls.py                       # control planets in sector 98
python scripts/toi540_timing.py                         # sector-98 time stamps; also queries JPL Horizons
python3 figures/make_figures.py other                   # the figures of sections 6.1 and 6.2
python3 figures/make_figures.py eb124751008             # the figure of section 6.4
```

Tables written: `eb_eclipses.py` writes [`data/eclipse_times_11613676.csv`](../data/eclipse_times_11613676.csv), `depths_11613676.csv`, `sectors_11613676.csv`,
`diffimage_11613676.csv`, `checks_11613676.csv`, `profiles_11613676.csv` and `predicted_eclipses_11613676.csv`;
`eb_model.py` writes [`data/eclipse_model_11613676.csv`](../data/eclipse_model_11613676.csv) and `model_history_11613676.csv`; `toi540_analysis.py` writes
[`data/ephemeris_200322593.csv`](../data/ephemeris_200322593.csv), `transit_depths_200322593.csv`, `s98_tests_200322593.csv`, `injection_200322593.csv`,
`bt_fits_200322593.csv`, `perturber_200322593.csv`, `folded_200322593.csv` and `ground_windows_200322593.csv` (run with
`--start 2026-10-20`); `toi540_controls.py` writes [`data/control_planets_s98.csv`](../data/control_planets_s98.csv); `toi540_timing.py` writes
[`data/timing_audit_200322593.csv`](../data/timing_audit_200322593.csv); `eb_124751008.py` writes [`data/sectors_124751008.csv`](../data/sectors_124751008.csv), `eclipse_times_124751008.csv`,
`depths_124751008.csv` and `predicted_eclipses_124751008.csv`. The formats are described in the comment headers of the
files.

MAST products: TESScut 21 $`\times`$ 21-pixel cutouts of TIC 11613676 (Gaia DR3 position) and its SPOC 2-min light
curves `tess2018263035959-s0003-0000000011613676-0123-s_lc.fits` and `tess2020266004630-s0030-0000000011613676-0195-s_lc.fits`;
the SPOC 2-min target pixel files and light curves of TOI-540 (`tess*-s00{04,05,06,31,32,98}-0000000200322593-*-s_{tp,lc}.fits`)
and its 20-s light curves (`*-a_fast-lc.fits`, sectors 31, 32, 98); the SPOC 2-min light curves of the control stars
([`data/control_planets_s98.csv`](../data/control_planets_s98.csv)). The SPOC data-validation summaries of TIC 11613676
(`tess2018206190142-s0001-s0036-0000000011613676-01-00471_dvs.pdf` and `-02-00471_dvs.pdf`) are quoted in section [6.1.7](6.1_TIC11613676.md#617-catalogue-and-literature-checks)
and not used by the scripts. The control-planet ephemerides are the ExoFOP TOI-table values of 2026-10-18, copied into
[`scripts/targets.py`](../scripts/targets.py); `toi540_timing.py` queries the JPL Horizons vectors of TESS (-95). The scripts were run on
2026-10-19 with the versions of `requirements.txt` (MAST downloads of 2026-10-19).

Without script in this repository, for TIC 11613676: the models N1, N2, MB, MC and N3 and the second implementation of
N4 (section [6.1.6](6.1_TIC11613676.md#616-checks-and-geometric-models)), the out-of-eclipse periodograms, the nodal-precession estimate, the Gaia DR3, RAVE and catalogue
values of section [6.1.1](6.1_TIC11613676.md#611-target) and the catalogue checks of section [6.1.7](6.1_TIC11613676.md#617-catalogue-and-literature-checks). For TOI-540 b: the per-pixel S/N maps, the QLP
sector-98 light curve, the MEarth light curve and the rotation periodogram (section [6.2.5](6.2_TOI-540b.md#625-checks-for-instrumental-and-processing-causes)), and the catalogue checks of
section [6.2.8](6.2_TOI-540b.md#628-catalogue-and-literature-checks). The detection-stage values of TIC 11613676 (section [6.1.4](6.1_TIC11613676.md#614-eclipse-depths-and-durations)) come from 3 $`\times`$ 3-pixel box photometry
that is not part of this repository.

For TIC 124751008 (section [6.4](6.4_TIC124751008.md)), without script in this repository (vetting record of 2026-10-19): the Gaia DR3,
catalogue and X-ray values and the luminosity and photocentre estimates of section [6.4.1](6.4_TIC124751008.md#641-target), the difference images, the
depths in other light-curve sources, the phase scans for the secondary eclipse, the out-of-eclipse variability, the
O $`-`$ C periodogram, the geometric models, the per-eclipse inclinations and the $`i(t)`$ fits
([`data/inclinations_124751008.csv`](../data/inclinations_124751008.csv), `model_history_124751008.csv`; the third-light-0.70 values were recomputed on
2026-10-19 with the per-eclipse inclination grid extended from 88.2 to 89.0 deg, because the 2019 values had been
truncated at 88.2 deg), the third-body estimates and the catalogue checks of section [6.4.7](6.4_TIC124751008.md#647-catalogue-and-literature-checks). The tables of
`eb_124751008.py` reproduce the per-eclipse depths and times of the vetting record to within the errors (the same
cutouts; the linear ephemeris agrees to 0.0001 d in $`T_0`$).

Section [6.3](6.3_eclipsing_binaries.md) and [`data/notable_ebs_2026.csv`](../data/notable_ebs_2026.csv) have no script in this repository: the values are exported from the vetting
record of 2026-10-18 (TIC 342194997: 2026-10-19). Method summary: QLP SAP light curves of all available sectors (MAST, `hlsp_qlp_tess_ffi_*`) and
TESScut 21 $`\times`$ 21-pixel cutouts with circular-aperture photometry; eclipse mid-times measured in each eclipse
(half-depth crossings for TIC 304936165, TIC 13957664 and TIC 261702897) and a linear ephemeris fitted to them; secondary-eclipse and alias searches at the covered epochs; stacked difference images fitted with a Gaussian PSF
at every Gaia DR3 source within 60" and at a free position; Gaia DR3 NSS elements from VizieR (I/357), with the
conjunction times of TIC 13957664, TIC 304936165 and TIC 265436769 propagated by Monte Carlo over the NSS errors;
SkyBoT for TIC 13957664; for TIC 342194997 a Monte Carlo over the Gaia DR3 Thiele-Innes elements, eccentricity,
period and periastron time with their correlation matrix (inclination, conjunction times, transverse velocities,
projected separations; the predicted conjunction recomputed on 2026-10-19 with the parameter order of the Gaia DR3
documentation for the correlation vector), a two-sphere eclipse model with quadratic limb darkening fitted to the TESScut photometry of
sector 100, and an occultation-template search of sector 99; and the catalogue queries listed in section [6.3](6.3_eclipsing_binaries.md)
(2026-10-18 and 2026-10-19).

For L 98-59 d and TOI-620 b (sections [6.6](6.6_L98-59d.md) and [6.7](6.7_TOI-620b.md)) the reproduction steps are given in the sections themselves:
`scripts/transit_duration/run_l9859.sh` and `run_toi620.sh` (sections 6.6.3-6.6.9, 6.7),
`scripts/transit_duration/joint/run_joint.sh` (section 6.6.10) and `scripts/transit_duration/jwst_chromatic/run_chromatic.sh`
(section [6.6.11](6.6_L98-59d.md#6611-jwst-durations-per-detector-and-wavelength-2026-10-21), after steps 8-10 of `run_joint.sh`; tables by `make_tables_chromatic.py`, figure by
`python3 figures/make_figures.py jwstchrom`). The per-column JWST light curves of section 6.6.11 (about 90 MB) are not
part of this repository; `run_chromatic.sh` rebuilds them from the public MAST `calints` subsets retrieved by
`run_joint.sh`.

For TOI-6482.01 (section [6.8](6.8_TOI-6482.md)), run from the repository root with the versions of `requirements.txt`:

```
python scripts/toi6482/fit_toi6482.py            # QLP light curves (MAST HLSP) -> data/*_189300053.csv
python scripts/toi6482/timestamp_check_s99.py    # sector-99 QLP time stamps against TESScut (MAST)
python figures/make_figures.py toi6482
```

The QLP files (`hlsp_qlp_tess_ffi_s00{08,35,62,89,99}-0000000189300053_tess_v01_llc.fits`) were retrieved on 2026-10-20 and
again on 2026-10-21 (identical checksums); the run of 2026-10-21 reproduces the mid-times of the ExoFOP upload of
2026-10-20 within 0.11 min. With numpy 2.4.4, scipy 1.17.1 and astropy 7.2.0 the transit-time, phase-window and depth
tables are identical, the ephemeris table differs only beyond the precision quoted in section 6.8, and the light-curve
table differs in the last written digit. The TRAPPIST-South light curve (reduction v2 provided by the observer) is not
in this repository; without `--trappist FILE` the script keeps the TRAPPIST row of the existing table, if present. Without
script in this repository: the Gaia DR3, AllWISE and CatWISE2020 values, the speckle separation limits (interpolated
in the contrast-curve tables on ExoFOP), the tess-point sector check and the catalogue checks of section [6.8.6](6.8_TOI-6482.md#686-catalogue-and-literature-checks).

For TOI-2638.01 (section [6.9](6.9_TOI-2638.md)), run from the repository root with the versions of `requirements.txt` (pyvo 1.9.1 for the
Gaia DR3 query):

```
python scripts/toi2638/download_toi2638.py --tpf    # MAST: SPOC 2-min light curves and target pixel files, TESS-SPOC, QLP, TGLC; QLP of the 17" neighbour
python scripts/toi2638/fit_toi2638.py               # -> data/{ephemeris,depths,pertransit,shape,jointfit,trends,predictions}_458342769.csv
python scripts/toi2638/checks_toi2638.py            # -> data/{sectors,neighbours,checks,precession}_458342769.csv (queries Gaia DR3)
python scripts/toi2638/pixels_toi2638.py            # -> data/pixels_458342769.csv
python scripts/toi2638/sims_toi2638.py              # -> data/simulations_458342769.csv
python figures/make_figures.py toi2638
```

The MAST products (the download script writes their MAST URIs to `work/458342769/dl/products.txt`) were retrieved on 2026-10-22; the light
curves are identical to copies retrieved on 2026-10-21. The scripts use TESS data, Gaia DR3 and tess-point only.
Without script in this repository: the search that selected the planet, the TOI-table and ExoFOP values and the
Gaia DR3 astrophysical parameters of section [6.9.1](6.9_TOI-2638.md#691-target), and the catalogue and literature checks of section
[6.9.7](6.9_TOI-2638.md#697-catalogue-and-literature-checks).

For TOI-1117 c (section [6.10](6.10_TOI-1117c.md)), run from the repository root; `download_toi1117.py`, `fit_toi1117.py`,
`checks_toi1117.py`, `pixels_toi1117.py` and the figure use the versions of `requirements.txt`, and `mcmc_toi1117.py`,
`ephemeris_toi1117.py`, `rv_toi1117.py` and `triceratops_toi1117.py` the emcee/TRICERATOPS environment listed there:

```
python scripts/toi1117/download_toi1117.py      # MAST: SPOC 2-min and 20-s, TESS-SPOC, QLP, eleanor-lite and TARS light curves, target pixel files, DV reports; Gaia DR3 and TIC cones; arXiv source 2506.05521
python scripts/toi1117/fit_toi1117.py           # -> data/{depths,pertransit,depth_tests}_295541511.csv and work/295541511/lc_fit_{pdc,sap}.npz
python scripts/toi1117/mcmc_toi1117.py          # -> data/fit_295541511.csv and work/295541511/chain_*.npz (about 30 min)
python scripts/toi1117/ephemeris_toi1117.py     # -> data/{ephemeris,transit_times,predicted_transits,folded,model}_295541511.csv
python scripts/toi1117/rv_toi1117.py            # -> data/{rv_harps,rvfit,rv_ephemeris}_295541511.csv
python scripts/toi1117/checks_toi1117.py        # -> data/{falsealarm,epochs,checks,timing}_295541511.csv (about 25 min)
python scripts/toi1117/pixels_toi1117.py        # -> data/{pixels,sources,apertures}_295541511.csv
python scripts/toi1117/triceratops_toi1117.py   # -> data/triceratops_295541511.csv, data/triceratops_input_295541511.csv, data/triceratops/probs_295541511_*.csv
python figures/make_figures.py toi1117
```

The MAST products (listed with their URIs in `work/295541511/dl/products.txt`) were retrieved on 2026-10-22 and listed
again on 2026-10-23 (no new or changed product names); the arXiv source of Lockley et al. (2025) has sha256
8076ff3adfbc1699395091e0fe8dfe7e75b64288de82f41dabf539472f3c5e7e. The HARPS velocities in
`data/rv_harps_295541511.csv` are transcribed from that source. The transit models sample each exposure at the centres of
equal parts (batman's `supersample_factor` is not used). The MCMC, false-alarm, injection and TRICERATOPS runs use explicit seeds;
TRICERATOPS also depends on the TIC and Gaia DR3 queries made at run time, and a repeat of the four runs with the same seeds gave the same FPP and NFPP to four decimals. Without script in this repository: the search
that selected the signal, the planet-b ephemeris used for masking, the Gaia DR3 values of section
[6.10.1](6.10_TOI-1117c.md#6101-target) other than those in the Gaia DR3 cone, the values quoted from Lockley et al. (2025),
the approximate contrast curve, and the catalogue and literature checks of section
[6.10.7](6.10_TOI-1117c.md#6107-catalogue-and-literature-checks).

For LHS 475 b (section [6.11](6.11_LHS475b.md)), run from the repository root with the versions of `requirements.txt`:

```
python scripts/lhs475/trend_lhs475.py      # -> data/{trend,forecast}_369327947.csv from data/{sectors,jwst}_369327947.csv
python figures/make_figures.py lhs475
```

Without script in this repository (vetting record of 2026-10-28 and 2026-10-29): the per-sector fits of
`data/sectors_369327947.csv` (SPOC 2-min PDCSAP and SAP light curves of sectors 12-102, QLP light curves and own TESScut
photometry, from MAST, as used in the vetting of 2026-10-28; batman models with `supersample_factor`, i.e. the sampling that the
correction of 2026-10-22 below replaced elsewhere, applied identically to every 120-s sector), the injection-recovery of
`data/injection_369327947.csv`, the tests of `data/checks_369327947.csv` (including the chi-square profiles, the template
variants, the shape and fold tests, the noise and timing checks and the stellar density from the Gaia DR3 parallax and the
2MASS $`K_s`$ magnitude), the fits of the JWST GO 1981 `x1dints` products (CAL_VER 3.0.0, MAST, retrieved 2026-10-29) in
`data/jwst_369327947.csv`, the fits of the ExoFOP TFOP light curves and the MEarth DR11 survey photometry in
`data/ground_369327947.csv`, the companion estimates of `data/perturber_369327947.csv`, the screen of 271 planets that
selected the planet, the TESS sector and ground-based observing windows of section [6.11.8](6.11_LHS475b.md#6118-forecast-and-decisive-tests), and the catalogue and
literature checks of section [6.11.9](6.11_LHS475b.md#6119-catalogue-and-literature-checks).

Correction of 2026-10-22 (exposure sampling; no conclusion changed). batman's `supersample_factor` places the
sub-exposures from one end of the exposure to the other, which smooths the model more than the exposure does (by
factors of 1.07-1.73 in rms for 15 to 2 sub-exposures). The transit models of `fit_transits.py`, `make_tables.py`,
`triceratops_prep.py`, `fit_free_dilution.py`, `fit_second_signal.py`, `toi540_analysis.py`, `toi540_controls.py` and of
`transit_duration/lib.py`, `v8_inject.py`, `joint/t3_windows.py`, `joint/t11_tpf.py`, `joint/j2_fit.py` and
`jwst_chromatic/c3_binfit.py`, `c7_white_diag.py`, `c9_sysvariants.py` now sample each exposure at the centres of its
n equal parts (`exp_time` passed as $`T (n - 1)/n`$); `scripts/toi2638/` uses the same sampling. Regenerated with the
corrected scripts on 2026-10-22, with the commands above: `transit_times_`, `folded_`, `model_` and `sectors_` of TIC 402710280,
324429952, 89247735, 455017950 and 311179599; `predicted_transits_311179599.csv` and `tess_sectors_311179599.csv`;
`fit_free_dilution[_rho]_` of TIC 89247735, 455017950 and 311179599; the four tables of candidate c (and its $`r = 1.5`$ px fit quoted in section 5.10.1); the eight tables of
`toi540_analysis.py` and `control_planets_s98.csv` (with the control light curves of the original selection). The
values quoted in sections 1-6.2 were updated where they changed at the quoted precision. The figure models of `make_tables.py`
are fixed-mid-time refits of the common shape (TIC 324429952: shape fixed in `targets.py`) and can differ slightly from
the simultaneous fits quoted in sections 1-5. The ephemeris of TIC 311179599 b
in `targets.py` was updated to the corrected fit ($`P = 24.826859 \pm 0.000044`$ d, $`T_0`$ = BTJD $`4082.33414 \pm 0.00154`$;
before: $`24.826839 \pm 0.000047`$ d, BTJD $`4082.33390 \pm 0.00167`$). The adopted ephemerides of TIC 402710280,
324429952, 89247735 and 455017950 and of candidate c, which set the O $`-`$ C values and the predicted transits, were
kept; the corrected fits differ from them by at most 0.15 $`\sigma`$ (TIC 402710280: $`T_0`$ = BTJD 4128.92682, $`P = 36.340046`$ d;
TIC 324429952: 4174.88768, 42.375735 d; TIC 89247735: 4196.80900, 18.323174 d; TIC 455017950: 4055.48919, 40.605654 d;
candidate c: 4088.51676, 73.86722 d). Not regenerated: the difference-image tables (`diffimage_*`), which use the
mid-transit times of the earlier fits (shifts of at most 1.7 min); the TRICERATOPS inputs and scenario tables, whose
inputs change by at most 9 ppm per bin (TIC 311179599: the sector-100 mid-time moved across 10-min bin edges, so its
input is re-binned, with bin scatter 279 instead of 285 ppm, depth 1252 instead of 1254 ppm); the FPPs were not re-run;
and the tables of sections 6.6 and 6.7, whose slopes and durations change by at most 0.21 $`\sigma`$ in a check of
the same scripts with both samplings (TOI-620 b PDCSAP slope $`-2.86`$ to $`-2.94`$ min/yr; the L 98-59 d core of
section 6.6.10 does not use batman). The 2-min PDCSAP per-sector timing files used in sections 6.6.3-6.6.9 predate a
later edit of `scripts/transit_duration/lib.py`; a rerun of the scripts gives slopes that differ from the tabulated
ones by up to 0.021 min/yr (0.18 $`\sigma`$).
