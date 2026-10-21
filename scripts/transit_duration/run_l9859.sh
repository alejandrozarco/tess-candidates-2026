#!/bin/sh
# L 98-59 (TIC 307210830): per-sector transit durations of planet d and of the controls b and c from the SPOC light
# curves, group fits, injections, and fits of the public JWST white-light curves (README section 6.6).
# Run from the repository root:   sh scripts/transit_duration/run_l9859.sh
# Work files go to $TD_WORK (default work/transit_duration); the tables in data/ are written afterwards by
#   python scripts/transit_duration/make_tables_td.py
# Python: numpy, scipy, astropy, astroquery, batman-package, emcee (versions in requirements.txt). Run time about 2 h
# on one core, most of it in v6_persector.py and j2_fit.py.
set -e
PY=${PY:-python}
W=${TD_WORK:-work/transit_duration}; export TD_WORK=$W
D=scripts/transit_duration
P=7.450729
T="0.0453 34.80 0.9257 0.1895 0.4604"   # template for d: constant-geometry common fit C of v5_groups.py (k, a/R*, b, u1, u2)
mkdir -p $W/dl/l9859 $W/dl/jwst

# 1. MAST downloads: SPOC 2-min and 20-s light curves; JWST stage-3 white-light curves (*_whtlt.ecsv)
$PY $D/dl_list.py 307210830 $W/dl/l9859_products.json
$PY $D/dl_get.py $W/dl/l9859_products.json $W/dl/l9859 lc,fastlc
$PY $D/jwst_get_whtlt.py $W/dl/jwst \
  jw01201-o311_t019_niriss_clear-gr700xd-substrip96 jw01224-o001_t009_nirspec_f290lp-g395h-s1600a1-sub2048 \
  jw04098-o002_t011_nirspec_f290lp-g395h-s1600a1-sub2048 jw01201-o302_t019_niriss_clear-gr700xd-substrip96 \
  jw02512-o015_t003_nirspec_f290lp-g395h-s1600a1-sub2048 jw02512-o016_t003_nirspec_f290lp-g395h-s1600a1-sub2048 \
  jw03942-o001_t001_nirspec_f290lp-g395h-s1600a1-sub2048 jw03942-o002_t001_nirspec_f290lp-g395h-s1600a1-sub2048 \
  jw03942-o003_t001_nirspec_f290lp-g395h-s1600a1-sub2048

# 2. light curves (quality 0) and per-transit windows (+-0.16 d; transits of the other planets masked)
$PY $D/v1_extract.py $W/dl/l9859 $W/lc
$PY $D/v2_windows.py L9859 d 2m pdc $W/win_L9859_d_2m_pdc.npy
$PY $D/v2_windows.py L9859 d 2m sap $W/win_L9859_d_2m_sap.npy
$PY $D/v2_windows.py L9859 d 20s pdc $W/win_L9859_d_20s_pdc.npy
$PY $D/v2_windows.py L9859 b 2m pdc $W/win_L9859_b_2m_pdc.npy
$PY $D/v2_windows.py L9859 c 2m pdc $W/win_L9859_c_2m_pdc.npy

# 3. per-transit mid-times, flare clipping and common shape (three iterations)
$PY $D/v4_timing.py $W/win_L9859_d_2m_pdc.npy $P 120 $W/tim_L9859_d_2m_pdc.npz 0.0453 34.80 0.9257 0.19 0.46
$PY $D/v4_timing.py $W/win_L9859_d_2m_sap.npy $P 120 $W/tim_L9859_d_2m_sap.npz 0.0453 34.80 0.9257 0.19 0.46
$PY $D/v4_timing.py $W/win_L9859_d_20s_pdc.npy $P 20 $W/tim_L9859_d_20s_pdc.npz 0.0453 34.80 0.9257 0.19 0.46
$PY $D/v4_timing.py $W/win_L9859_b_2m_pdc.npy 2.2531140 120 $W/tim_L9859_b_2m_pdc.npz 0.0243 15.15 0.51 0.19 0.46
$PY $D/v4_timing.py $W/win_L9859_c_2m_pdc.npy 3.6906764 120 $W/tim_L9859_c_2m_pdc.npz 0.0385 21.05 0.41 0.23 0.33

# 4. group fits of d (four epoch groups; free a/R*, and a/R* = 33.6 from the stellar density of Cadieux et al. 2025)
$PY $D/v5_groups.py $W/win_L9859_d_2m_pdc.npy $W/tim_L9859_d_2m_pdc.npz $P 120 ttv prior year $W/grp_d_2m_pdc_ttv_prior_year.json 0.17 0.44
$PY $D/v5_groups.py $W/win_L9859_d_2m_pdc.npy $W/tim_L9859_d_2m_pdc.npz $P 120 ttv prior year $W/grp_d_aR33.6.json 0.17 0.44 33.6 0.02

# 5. per-sector T14, b and width w with the fixed template; variants of the time stacking, flux and cadence
for m in ttv pt lin sec; do
  $PY $D/v6_persector.py $W/win_L9859_d_2m_pdc.npy $W/tim_L9859_d_2m_pdc.npz $P 120 $m sector $W/ps_d_2m_pdc_$m.json $T
done
$PY $D/v6_persector.py $W/win_L9859_d_2m_sap.npy $W/tim_L9859_d_2m_sap.npz $P 120 ttv sector $W/ps_d_2m_sap.json $T
$PY $D/v6_persector.py $W/win_L9859_d_20s_pdc.npy $W/tim_L9859_d_20s_pdc.npz $P 20 ttv sector $W/ps_d_20s_pdc.json $T
$PY $D/v7_filter.py $W/win_L9859_d_2m_pdc.npy $W/tim_L9859_d_2m_pdc.npz $W/win_L9859_d_2m_pdc_clean.npy $W/tim_L9859_d_2m_pdc_clean.npz
$PY $D/v6_persector.py $W/win_L9859_d_2m_pdc_clean.npy $W/tim_L9859_d_2m_pdc_clean.npz $P 120 ttv sector $W/ps_d_clean.json $T
# controls b and c (templates: k and b of Cadieux et al. 2025 or of step 3, a/R* from the stellar density; b grid from 0)
BMIN=0.0 $PY $D/v6_persector.py $W/win_L9859_c_2m_pdc.npy $W/tim_L9859_c_2m_pdc.npz 3.6906764 120 ttv sector $W/ps_c_2m_pdc.json 0.0393 21.05 0.41 0.232 0.334
BMIN=0.0 $PY $D/v6_persector.py $W/win_L9859_b_2m_pdc.npy $W/tim_L9859_b_2m_pdc.npz 2.2531140 120 lin sector $W/ps_b_2m_pdc.json 0.0243 15.15 0.51 0.19 0.46

# 6. constant-geometry injections (b = 0.915) at phases 0.30 and 0.62 of d, full chain re-run
for n in 30 62; do
  $PY $D/v8_inject.py 0.$n 0.915 $W/lc_inj$n
  LCDIR=$W/lc_inj$n $PY $D/v2_windows.py L9859inj$n x 2m pdc $W/win_inj${n}_2m_pdc.npy
  $PY $D/v4_timing.py $W/win_inj${n}_2m_pdc.npy $P 120 $W/tim_inj${n}_2m_pdc.npz 0.0453 34.80 0.915 0.19 0.46
  $PY $D/v6_persector.py $W/win_inj${n}_2m_pdc.npy $W/tim_inj${n}_2m_pdc.npz $P 120 ttv sector $W/ps_inj$n.json $T
done

# 7. JWST white-light curves (NIRSpec G395H: NRS2; NIRISS SOSS: order 1), first 3 min removed
J=$W/dl/jwst
$PY $D/j1_prep.py $J/jw01201-o311_t019_niriss_clear-gr700xd-substrip96_whtlt.ecsv whitelight_flux_order_1 $W/jw_d_1201.npz 3
$PY $D/j1_prep.py $J/jw01224-o001_t009_nirspec_f290lp-g395h-s1600a1-sub2048_whtlt.ecsv whitelight_flux_NRS2 $W/jw_d_1224.npz 3
$PY $D/j1_prep.py $J/jw04098-o002_t011_nirspec_f290lp-g395h-s1600a1-sub2048_whtlt.ecsv whitelight_flux_NRS2 $W/jw_d_4098b.npz 3
$PY $D/j1_prep.py $J/jw01224-o001_t009_nirspec_f290lp-g395h-s1600a1-sub2048_whtlt.ecsv whitelight_flux_NRS2 $W/jw_d_1224t10.npz 10
$PY $D/j1_prep.py $J/jw04098-o002_t011_nirspec_f290lp-g395h-s1600a1-sub2048_whtlt.ecsv whitelight_flux_NRS2 $W/jw_d_4098bt10.npz 10
$PY $D/j1_prep.py $J/jw01201-o302_t019_niriss_clear-gr700xd-substrip96_whtlt.ecsv whitelight_flux_order_1 $W/jw_c_1201.npz 3
$PY $D/j1_prep.py $J/jw02512-o015_t003_nirspec_f290lp-g395h-s1600a1-sub2048_whtlt.ecsv whitelight_flux_NRS2 $W/jw_c_2512a.npz 3
$PY $D/j1_prep.py $J/jw02512-o016_t003_nirspec_f290lp-g395h-s1600a1-sub2048_whtlt.ecsv whitelight_flux_NRS2 $W/jw_c_2512b.npz 3
for o in 1:a 2:b 3:c; do
  $PY $D/j1_prep.py $J/jw03942-o00${o%:*}_t001_nirspec_f290lp-g395h-s1600a1-sub2048_whtlt.ecsv whitelight_flux_NRS2 $W/jw_b_3942${o#*:}.npz 3
done
# GO 4098 visit cut to the pre-ingress coverage of the GTO 1224 visit (38.67 min before mid-transit)
$PY $D/j3_cut.py $W/jw_d_4098b.npz $W/jw_d_4098bmim.npz 3381.86015

# 8. JWST fits (emcee, 40 walkers, 3000 steps; errors inflated by the red-noise factor beta)
export NWALK=40 NSTEP=3000
for n in d_1201 d_1224 d_4098b; do
  $PY $D/j2_fit.py $W/jfit_$n.json $P single $n
  AR_FIX=33.6 $PY $D/j2_fit.py $W/jfit_${n}_aR33.6.json $P single $n
done
$PY $D/j2_fit.py $W/jfit_d_nrs2_joint.json $P joint d_1224 d_4098b
$PY $D/j2_fit.py $W/jfit_d_nrs2_joint_t10.json $P joint d_1224t10 d_4098bt10
$PY $D/j2_fit.py $W/jfit_d_nrs2_joint_mim.json $P joint d_1224 d_4098bmim
AR_FIX=33.6 $PY $D/j2_fit.py $W/jfit_d_nrs2_joint_aR33.6.json $P joint d_1224 d_4098b
$PY $D/j2_fit.py $W/jfit_c_2512_joint.json 3.6906764 joint c_2512a c_2512b
$PY $D/j2_fit.py $W/jfit_c_1201.json 3.6906764 single c_1201
$PY $D/j2_fit.py $W/jfit_b_3942_joint.json 2.2531140 joint b_3942a b_3942b b_3942c
echo done
