#!/bin/bash
# TOI-421 b (TIC 94986319): JWST white-light fits of section 6.12 (program GO 1935, public). Run from the repository root:
#   bash scripts/toi421/jwst/run_jwst.sh        (PY=<python> to choose the interpreter; about 2 GB of downloads)
# Inputs: Zenodo record 10.5281/zenodo.15176438 (Davenport et al. 2025; file TOI421b.zip, only the NIRSpec Stage-4
# white-light h5 is used) and the MAST stage-2 x1dints of visits jw01935002001 (NIRSpec G395M, NRS1) and jw01935001001
# (NIRISS SOSS). Downloaded files are read as data only (nothing in them is executed).
set -e
PY=${PY:-python}; S=scripts/toi421/jwst; W=work/toi421/jwst
mkdir -p $W/data $W/lc $W/fits $W/plots $W/data/exotic_ld_data
cd $W/data
[ -f TOI421b.zip ] || curl -L -o TOI421b.zip "https://zenodo.org/records/15176438/files/TOI421b.zip?download=1"
H5=$(unzip -Z1 TOI421b.zip | grep 'Stage4/White_Light/S4_toi421b_ap3_bg8_LCData.h5' | head -1)
[ -f zenodo/$H5 ] || unzip -o -q TOI421b.zip "$H5" -d zenodo
for f in jw01935002001_04102_00001-seg00{1,2,3}_nrs1_x1dints.fits jw01935001001_04101_00001-seg00{1,2,3,4}_nis_x1dints.fits; do
  [ -f $f ] || curl -L -o $f "https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:JWST/product/$f"
done
cd - >/dev/null
$PY -I $S/prep_eureka_nrs1.py $W/data/zenodo/$H5 $W/lc/V1_NIRSpecG395M_NRS1_eureka.npz
$PY -I $S/prep_x1dints.py $W/lc/V1_NIRSpecG395M_NRS1_mast.npz V1_NIRSpecG395M_NRS1_mast NIRSpec_G395M_NRS1_2.87-5.17um -1 2.87 5.17 \
  "planet=TOI-421 b (TIC 94986319); visit=2023-09-10 JWST NIRSpec G395M (program 1935, obs 2); detector=NRS1; source=MAST stage-2 x1dints, own white-light sum (cross-check of Eureka LC)" \
  $W/data/jw01935002001_04102_00001-seg00?_nrs1_x1dints.fits
$PY -I $S/prep_x1dints.py $W/lc/V2_NIRISSSOSS_O1_mast.npz V2_NIRISSSOSS_O1_mast NIRISS_SOSS_order1_0.85-2.83um 1 0.85 2.83 \
  "planet=TOI-421 b (TIC 94986319); visit=2023-11-06 JWST NIRISS SOSS SUBSTRIP96 (program 1935, obs 1); detector=NIS order 1; source=MAST stage-2 x1dints, own white-light sum" \
  $W/data/jw01935001001_04101_00001-seg00?_nis_x1dints.fits
$PY -I $S/ld_theory.py $W/data/exotic_ld_data $W/ld_theory.json
E=$W/lc/V1_NIRSpecG395M_NRS1_eureka.npz; M=$W/lc/V1_NIRSpecG395M_NRS1_mast.npz; O=$W/lc/V2_NIRISSSOSS_O1_mast.npz
for j in "$E a 24000" "$E b 24000" "$E c 24000" "$E bs 12000" "$E bl 12000" "$M a 12000" "$M b 12000" \
         "$O a 24000" "$O b 24000" "$O c 24000" "$O bs 12000" "$O bl 12000"; do
  set -- $j
  $PY -I $S/fit_wlc.py $1 $2 $W/ld_theory.json $W/fits $3 > $W/fits/log_$(basename $1 .npz)_$2.txt 2>&1
done
$PY -I $S/summarize.py $W/fits $W/fit_table.md $W/plots
