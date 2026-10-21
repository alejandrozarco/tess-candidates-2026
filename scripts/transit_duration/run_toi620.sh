#!/bin/sh
# TOI-620 b (TIC 296739893): per-sector and per-transit durations from the SPOC light curves and from own photometry of
# the SPOC target pixel files, and the two-star decomposition of the difference images (README section 6.7).
# Run from the repository root:   sh scripts/transit_duration/run_toi620.sh   (after or without run_l9859.sh)
# Tables: python scripts/transit_duration/make_tables_td.py
set -e
PY=${PY:-python}
W=${TD_WORK:-work/transit_duration}; export TD_WORK=$W
D=scripts/transit_duration
P=5.098784; T0=1518.004514           # ExoFOP TOI-table ephemeris (window centres and difference images)
mkdir -p $W/dl/toi620 $W/tpf620
pars() { $PY -c "import numpy as np; print(' '.join('%.5f' % x for x in np.load('$1', allow_pickle=True)['pars']))"; }

# 1. MAST downloads: SPOC 2-min light curves and target pixel files (sectors 8, 35, 62, 89, 99)
$PY $D/dl_list.py 296739893 $W/dl/toi620_products.json
$PY $D/dl_get.py $W/dl/toi620_products.json $W/dl/toi620 lc,tp
$PY $D/v1_extract.py $W/dl/toi620 $W/lc

# 2. SPOC PDCSAP and SAP: windows of +-0.2 d, timing, per-sector and per-transit measurements (b grid from 0.6)
for fl in pdc sap; do
  $PY $D/v2_windows.py TOI620 b 2m $fl $W/win_T620_2m_$fl.npy 0.2
  $PY $D/v4_timing.py $W/win_T620_2m_$fl.npy $P 120 $W/tim_T620_2m_$fl.npz 0.0627 18.87 0.887 0.30 0.33
  BMIN=0.6 $PY $D/v6_persector.py $W/win_T620_2m_$fl.npy $W/tim_T620_2m_$fl.npz $P 120 pt sector $W/ps_T620_2m_${fl}_pt.json $(pars $W/tim_T620_2m_$fl.npz)
done
TP=$(pars $W/tim_T620_2m_pdc.npz)
BMIN=0.6 $PY $D/v6_persector.py $W/win_T620_2m_pdc.npy $W/tim_T620_2m_pdc.npz $P 120 sec sector $W/ps_T620_2m_pdc_sec.json $TP
BMIN=0.6 $PY $D/v6_persector.py $W/win_T620_2m_pdc.npy $W/tim_T620_2m_pdc.npz $P 120 sec transit $W/ps_T620_2m_pdc_sec_tr.json $TP

# 3. own photometry of the target pixel files (circular apertures, local background) and difference images
for f in $W/dl/toi620/*_tp.fits; do
  s=$(basename $f | sed 's/.*-s\(00[0-9][0-9]\)-.*/s\1/')
  $PY $D/p1_tpf.py $f $W/tpf620/$s $P $T0
done
$PY $D/p2_psf.py s0008 s0035 s0062 s0089 s0099

# 4. the same measurements in the own r = 1.5 px and 3.0 px apertures
for ap in 1.5 3.0; do
  tag=r$(echo $ap | tr -d .)
  $PY $D/p3_aplc.py 296739893 r$ap $W/lc_t620_$tag $W/tpf620/s0008 $W/tpf620/s0035 $W/tpf620/s0062 $W/tpf620/s0089 $W/tpf620/s0099
  LCDIR=$W/lc_t620_$tag $PY $D/v2_windows.py TOI620 b 2m pdc $W/win_T620_$tag.npy 0.2
  $PY $D/v4_timing.py $W/win_T620_$tag.npy $P 120 $W/tim_T620_$tag.npz 0.0602 16.64 0.9326 0.294 0.324
  BMIN=0.6 $PY $D/v6_persector.py $W/win_T620_$tag.npy $W/tim_T620_$tag.npz $P 120 sec sector $W/ps_T620_${tag}_sec.json 0.0602 16.64 0.9326 0.294 0.324
done
echo done
