#!/bin/sh
# L 98-59 d (TIC 307210830): joint inference of the transit geometry over all TESS transits, significance diagnostics
# (block bootstraps, null injections into the light curves and into the target pixels, look-elsewhere factor),
# telemetry and activity audit, and own extractions and joint fits of the JWST NIRSpec and NIRISS white-light curves
# (section 6.6.10 of docs/6.6_L98-59d.md).
#
# Run from the repository root, after steps 1 and 7 of scripts/transit_duration/run_l9859.sh (SPOC light curves and
# the product list; JWST pipeline white-light curves):
#   sh scripts/transit_duration/joint/run_joint.sh [survey_summary.csv]
# Work files go to $TDJ_WORK (default work/transit_duration_joint); the inputs of run_l9859.sh are read from $TD_WORK
# (default work/transit_duration). Results: $TDJ_WORK/out/summary.json, out/significance.json, out/jwst_slopes.json;
# the tables in data/ are then written by
#   python scripts/transit_duration/joint/make_tables_joint.py
# The optional argument is the per-planet table of the 802-planet search (not part of this repository); without it
# t22_significance.py skips the survey rank and uses the fixed trial numbers 801 and 3204.
# Python: numpy, scipy, astropy, astroquery, batman-package, emcee, h5py, arviz, celerite (versions in
# requirements.txt). Downloads: 27 SPOC target pixel files (about 1.3 GB), HTTP-range subsets of 20 JWST calints files
# and 6 NIRISS x1dints files (about 5 GB). The steps follow the order of the original runs, which were made in
# separate queues. Bootstraps, injection phases and emcee starting positions are seeded; the emcee moves are not, so
# posterior values reproduce to within their Monte Carlo error.
set -e
PY=${PY:-python}
TD=${TD_WORK:-work/transit_duration}
W=${TDJ_WORK:-work/transit_duration_joint}; export TDJ_WORK=$W
D=scripts/transit_duration/joint
O=$W/out
mkdir -p $O/mcmc $O/jwst/fits $W/dl/tpf $W/dl/jwst
cp $D/claret_ld.json $O/

# 1. SPOC 2-min and 20-s light curves with telemetry; ephemeris and per-sector depth of the 1.0492-d eclipses;
#    per-transit windows of d (b, c, eclipses and flares masked)
$PY $D/t1_extract.py $TD/dl/l9859 $O/lc
$PY $D/t2_eb.py
$PY $D/t3_windows.py 2m pdc mask $O/win_2m_pdc_ebmask.npy

# 2. clipping, per-transit white noise, per-sector red-noise factor, starting point, per-transit mid-times
$PY $D/t5_prep.py $O/win_2m_pdc_ebmask.npy $O/prep_2m_pdc_ebmask.npz

# 3. posteriors: main pair (b linear in time against constant b) and variants, each with its constant partner
run() { $PY $D/t6_mcmc.py $O/prep_2m_pdc_ebmask.npz $2 $3 $4 $5 $O/mcmc/$1; }
run main_lin_b lin_b 4 32 2500
run main_const const 4 32 2500
for v in lin_b_kcamp const_kcamp lin_b_ldfree const_ldfree lin_b_ecc const_ecc lin_i camp camp_slope lin_b_fixtc const_fixtc lin_b_rhofree; do
  run v_$v $v 2 24 1200
done

# 4. per-transit impact parameters and durations; telemetry and activity audit
$PY $D/t8_pertransit.py $O/prep_2m_pdc_ebmask.npz $O/mcmc/main_lin_b.json $O/pertransit_2m_pdc.json
$PY $D/t13_audit.py $O/pertransit_2m_pdc.json $O/prep_2m_pdc_ebmask.npz $O/mcmc/main_lin_b.json $O/win_2m_pdc_ebmask_flares.npy $O/audit_2m_pdc.json

# 5. alternative data sets (least squares, same chain): SAP, eclipses modelled, regressors, 20 s
alt() { $PY $D/t12_ownphot.py $1 $2 $3 $4 $5 $6 $O/alt_$1.json; }
alt pdc2m 2m pdc mask - 0
alt pdc2m_regs 2m pdc mask - 1
alt sap2m 2m sap mask - 0
alt pdc2m_ebmodel 2m pdc model - 0
alt pdc20s 20s pdc mask - 0
alt sap20s 20s sap mask - 0
# own aperture photometry from the target pixel files (no SPOC photometry)
$PY scripts/transit_duration/dl_get.py $TD/dl/l9859_products.json $W/dl/tpf tp
$PY $D/t11_tpf.py compact $W/dl/tpf
for ap in spoc r15 r25; do
  $PY $D/t11_tpf.py phot $ap $O/lc_own_$ap
  alt own_$ap 2m pdc mask $O/lc_own_$ap 0
  alt own_${ap}_regs 2m pdc mask $O/lc_own_$ap 1
done

# 6. block bootstraps and jackknives
$PY $D/t9_boot.py $O/win_2m_pdc_ebmask.npy jack 0 $O/boot_jack.npz
$PY $D/t9_boot.py $O/win_2m_pdc_ebmask.npy transit 500 $O/boot_transit.npz
$PY $D/t9_boot.py $O/win_2m_pdc_ebmask.npy sector 300 $O/boot_sector.npz

# 7. null injections of a constant-geometry copy of d at random phases: into the PDCSAP flux, and into the target pixels
$PY $D/t10_inject.py lc 200 123 $O/inj_lc.jsonl $O/mcmc/main_lin_b.json
$PY $D/t10_inject.py pix 40 456 $O/inj_pix2.jsonl $O/mcmc/main_lin_b.json

# 8. JWST: HTTP-range subsets of the public calints products; stage-2 NIRISS x1dints
J=$W/dl/jwst
$PY $D/jw_hdu_get.py $J SCI,WAVELENGTH,INT_TIMES \
  jw01224001001_04102_00001-seg001_nrs1_calints.fits jw01224001001_04102_00001-seg002_nrs1_calints.fits \
  jw04098002001_04102_00001-seg001_nrs1_calints.fits jw04098002001_04102_00001-seg002_nrs1_calints.fits \
  jw01224001001_04102_00001-seg001_nrs2_calints.fits jw01224001001_04102_00001-seg002_nrs2_calints.fits \
  jw04098002001_04102_00001-seg001_nrs2_calints.fits jw04098002001_04102_00001-seg002_nrs2_calints.fits \
  jw02512015001_04102_00001-seg001_nrs1_calints.fits jw02512015001_04102_00001-seg002_nrs1_calints.fits \
  jw02512015001_04102_00001-seg003_nrs1_calints.fits jw02512016001_04102_00001-seg001_nrs1_calints.fits \
  jw02512016001_04102_00001-seg002_nrs1_calints.fits jw02512016001_04102_00001-seg003_nrs1_calints.fits
$PY $D/jw_hdu_get.py $J SCI,INT_TIMES \
  jw01201311001_04101_00001-seg001_nis_calints.fits jw01201311001_04101_00001-seg002_nis_calints.fits \
  jw01201311001_04101_00001-seg003_nis_calints.fits jw04098001001_04101_00001-seg001_nis_calints.fits \
  jw04098001001_04101_00001-seg002_nis_calints.fits jw04098001001_04101_00001-seg003_nis_calints.fits
$PY $D/jw_get.py $J jw01201311001_04101_00001-seg001_nis_x1dints.fits jw01201311001_04101_00001-seg002_nis_x1dints.fits \
  jw01201311001_04101_00001-seg003_nis_x1dints.fits jw04098001001_04101_00001-seg001_nis_x1dints.fits \
  jw04098001001_04101_00001-seg002_nis_x1dints.fits jw04098001001_04101_00001-seg003_nis_x1dints.fits

# 9. own white-light extractions (NRS1 2.87-3.71 um, NRS2 3.83-5.15 um, box half-width 3 px; variants: 5 px,
#    NRS2 halves); NIRISS order 1 from x1dints and from calints; pipeline NRS2 white light (run_l9859.sh step 7)
for v in 1224001001 4098002001; do
  s=${v%??????}
  for det in nrs1 nrs2; do
    lo=2.87; hi=3.71; [ $det = nrs2 ] && lo=3.83 && hi=5.15
    $PY $D/j1_extract.py $O/jwst/wl_${s}_$det.npz $lo $hi 3 $J/jw0${v}_04102_00001-seg001_${det}_calints_sub.fits $J/jw0${v}_04102_00001-seg002_${det}_calints_sub.fits
    $PY $D/j1_extract.py $O/jwst/wlhw5_${s}_$det.npz $lo $hi 5 $J/jw0${v}_04102_00001-seg001_${det}_calints_sub.fits $J/jw0${v}_04102_00001-seg002_${det}_calints_sub.fits
  done
  $PY $D/j1_extract.py $O/jwst/wlblue_${s}_nrs2.npz 3.83 4.45 3 $J/jw0${v}_04102_00001-seg001_nrs2_calints_sub.fits $J/jw0${v}_04102_00001-seg002_nrs2_calints_sub.fits
  $PY $D/j1_extract.py $O/jwst/wlred_${s}_nrs2.npz 4.45 5.15 3 $J/jw0${v}_04102_00001-seg001_nrs2_calints_sub.fits $J/jw0${v}_04102_00001-seg002_nrs2_calints_sub.fits
done
$PY $D/j1_extract.py $O/jwst/wl_c2512a_nrs1.npz 2.87 3.71 3 $J/jw02512015001_04102_00001-seg00[123]_nrs1_calints_sub.fits
$PY $D/j1_extract.py $O/jwst/wl_c2512b_nrs1.npz 2.87 3.71 3 $J/jw02512016001_04102_00001-seg00[123]_nrs1_calints_sub.fits
$PY $D/j1b_soss.py $O/jwst/wl_1201_soss.npz $J/jw01201311001_04101_00001-seg00[123]_nis_x1dints.fits
$PY $D/j1b_soss.py $O/jwst/wl_4098o1_soss.npz $J/jw04098001001_04101_00001-seg00[123]_nis_x1dints.fits
X=$J/jw01201311001_04101_00001-seg001_nis_x1dints.fits
$PY $D/j1c_soss2d.py $O/jwst/wl2d_4098o1_soss.npz 13 $X $J/jw04098001001_04101_00001-seg00[123]_nis_calints_sub.fits
$PY $D/j1c_soss2d.py $O/jwst/wl2d_1201_soss.npz 13 $X $J/jw01201311001_04101_00001-seg00[123]_nis_calints_sub.fits
cp $TD/jw_d_1224.npz $O/jwst/wlpipe_1224_nrs2.npz
cp $TD/jw_d_4098b.npz $O/jwst/wlpipe_4098_nrs2.npz

# 10. JWST fits (configurations in jwst_cfg/; light-curve paths there are relative to $TDJ_WORK)
for c in ns_main ns_main_v2 ns_const ns_aRfree ns_aR33.6 ns_aR37.1 ns_ldGressier ns_sqrtlaw ns_power2 ns_noregs ns_ramp \
         ns_poly1 ns_poly3 ns_trim30 ns_nrs1 ns_nrs2 ns_gp ns_gp_nrs1 ns_gp_nrs2 ns_matchcov ns_matchcov_nrs1 ns_matchcov_nrs2 \
         pipe_nrs2 pipe_nrs2_ramp ex_hw5_nrs1 ex_hw5_nrs2 ex_blue_nrs2 ex_red_nrs2 ctrl_c_nrs1 \
         ni_1201 ni_4098o1 ni_pair ni2d_pair all_visits all_visits_gp all_visits_2d all_visits_2d_nrs1 all_visits_2d_nrs2; do
  $PY $D/j2_fit.py $D/jwst_cfg/$c.json $O/jwst/fits/$c
done

# 11. summaries
$PY $D/t20_summary.py
$PY $D/t21_jwst_slope.py ns_main_v2 ns_gp ns_nrs1 ns_nrs2 ns_gp_nrs1 ns_gp_nrs2 ns_matchcov ns_matchcov_nrs1 ns_matchcov_nrs2 \
  pipe_nrs2 pipe_nrs2_ramp ex_hw5_nrs1 ex_hw5_nrs2 ex_blue_nrs2 ex_red_nrs2 ns_aRfree ns_aR33.6 ns_aR37.1 ns_ldGressier \
  ns_sqrtlaw ns_power2 ns_noregs ns_ramp ns_poly1 ns_poly3 all_visits all_visits_gp all_visits_2d all_visits_2d_nrs1 \
  all_visits_2d_nrs2 ni2d_pair
$PY $D/t22_significance.py "$@"
echo done
