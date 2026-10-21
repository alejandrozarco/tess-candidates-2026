#!/bin/sh
# L 98-59 d (TIC 307210830): JWST transit durations per detector and per wavelength bin, systematics and white-light
# diagnostics, empirical detector floor, JWST inventory and predictions (section 6.6.11 of docs/6.6_L98-59d.md).
#
# Run from the repository root, after steps 8-10 of scripts/transit_duration/joint/run_joint.sh (HTTP-range subsets of
# the public calints products and the NIRISS x1dints in $TDJ_WORK/dl/jwst; own and pipeline white light in
# $TDJ_WORK/out/jwst; joint white-light fits ns_main, ns_main_v2, ns_nrs1, ns_nrs2, pipe_nrs2, ctrl_c_nrs1, ni_1201 in
# $TDJ_WORK/out/jwst/fits):
#   sh scripts/transit_duration/jwst_chromatic/run_chromatic.sh
# Work files go to $TDC_WORK (default work/transit_duration_chromatic); $TDJ_WORK defaults to
# work/transit_duration_joint. The tables in data/ are then written by
#   python scripts/transit_duration/jwst_chromatic/make_tables_chromatic.py
# and the figure by python3 figures/make_figures.py jwstchrom.
# Python: numpy, scipy, astropy, astroquery, batman-package, emcee, matplotlib (versions in requirements.txt).
# Downloads: the Kurucz ATLAS9 intensity grid im05k2new.pck (32 MB); MAST queries only (steps 6-7).
# Memory: the column extractions need about 1.4 GB each and run one at a time.
# emcee starting positions are seeded; the emcee moves are not, so per-bin posteriors reproduce to within their Monte
# Carlo error. The deterministic steps (2, 4, 5) reproduce the published outputs exactly from the same inputs.
set -e
PY=${PY:-python}
TDJ=${TDJ_WORK:-work/transit_duration_joint}; export TDJ_WORK=$TDJ
W=${TDC_WORK:-work/transit_duration_chromatic}; export TDC_WORK=$W
D=scripts/transit_duration/jwst_chromatic
J=$TDJ/dl/jwst
mkdir -p $W/out $W/cfg $W/fits $W/dl $W/logs

# 1. per-column extractions (same settings as ../joint/j1_extract.py, box half-width 3 px) and NIRISS SOSS bins
for v in 1224001001 4098002001; do
  s=${v%??????}
  for det in nrs1 nrs2; do
    lo=2.87; hi=3.71; [ $det = nrs2 ] && lo=3.83 && hi=5.15
    $PY $D/c1_extract_cols.py $W/out/cols_${s}_$det.npz $lo $hi 3 $J/jw0${v}_04102_00001-seg001_${det}_calints_sub.fits $J/jw0${v}_04102_00001-seg002_${det}_calints_sub.fits
  done
done
$PY $D/c1_extract_cols.py $W/out/cols_c2512a_nrs1.npz 2.87 3.71 3 $J/jw02512015001_04102_00001-seg00[123]_nrs1_calints_sub.fits
$PY $D/c1_extract_cols.py $W/out/cols_c2512b_nrs1.npz 2.87 3.71 3 $J/jw02512016001_04102_00001-seg00[123]_nrs1_calints_sub.fits
$PY $D/c4_soss_bins.py $W/out/soss_1201_bins.npz $J/jw01201311001_04101_00001-seg00[123]_nis_x1dints.fits

# 2. limb darkening per bin (Castelli & Kurucz 2003 grid, [M/H] = -0.5) and the fit configurations
[ -f $W/dl/im05k2new.pck ] || curl -sS -o $W/dl/im05k2new.pck http://kurucz.harvard.edu/grids/gridm05/im05k2new.pck
$PY $D/c2_ld.py $W/dl/im05k2new.pck $W/out/ld_kurucz.json
$PY $D/c5_make_cfgs.py

# 3. per-bin fits (resumable: c3_binfit.py skips bins already in its output)
for c in d1224_nrs1_prior d1224_nrs2_prior d4098_nrs1_prior d4098_nrs2_prior d1224_nrs1_free d1224_nrs2_free d4098_nrs1_free d4098_nrs2_free \
         c2512a_nrs1_prior c2512b_nrs1_prior d1201_soss1_prior d1201_soss2_prior c2512a_nrs1_free c2512b_nrs1_free d1201_soss1_free d1201_soss2_free \
         d1224_nrs1_ramp d1224_nrs2_ramp; do
  $PY $D/c3_binfit.py $W/cfg/$c.json > $W/logs/fit_$c.log 2>&1
done

# 4. chromatic analysis, systematics-model variants, white-light diagnostics, floor and detector averages
$PY $D/c6_analyze.py $W/fits $W/out/chromatic.json > $W/logs/c6.log
$PY $D/c9_sysvariants.py $W/out/sysvariants.json > $W/logs/c9.log
$PY $D/c7_white_diag.py $W/out/white_diag.json $W/out/white_diag.png > $W/logs/c7.log
$PY $D/c11_summary.py $W/out/jwst_summary.json > $W/logs/c11.log

# 5. predictions for GO 10498
$PY $D/m5_predictions.py $W/out/predictions.json > $W/logs/m5.log

# 6-7. JWST inventory in MAST (results change with time; the published table is the query of 2026-10-21)
$PY $D/m1_mast_inventory.py $W/out/mast_inventory.csv > $W/logs/m1.log
$PY $D/m2_products.py jw10498-o004_t001_nirspec_f290lp-g395h-s1600a1-sub2048 jw01201-o302_t019_niriss_clear-gr700xd-substrip96 \
  jw01201301001_02101_00001-seg001_nis jw03730008001_03101_00001-seg001_mirimage jw03942-o001_t001_nirspec_f290lp-g395h-s1600a1-sub2048 \
  jw03942-o002_t001_nirspec_f290lp-g395h-s1600a1-sub2048 jw03942-o003_t001_nirspec_f290lp-g395h-s1600a1-sub2048 \
  jw03942-o004_t001_nirspec_f290lp-g395h-s1600a1-sub2048 jw01201-o311_t019_niriss_clear-gr700xd-substrip96 \
  jw04098-o001_t001_niriss_clear-gr700xd-substrip96 jw04098-o002_t011_nirspec_f290lp-g395h-s1600a1-sub2048 \
  jw02512-o015_t003_nirspec_f290lp-g395h-s1600a1-sub2048 jw02512-o016_t003_nirspec_f290lp-g395h-s1600a1-sub2048 \
  jw10498001001_xx102_00001_nirspec > $W/logs/m2.log
echo done
