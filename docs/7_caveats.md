[README](../README.md) / 7. Caveats

## 7. Caveats

- No high-resolution imaging was obtained for this work; the only existing observation is the published SOAR speckle
  observation of TIC 455017950. For TIC 402710280 the 2.3" neighbour is unresolved by TESS and
  excluded as host only by the transit shape; a companion closer than that is unconstrained. For TIC 324429952 any
  companion within about 1" is unconstrained, and the TRICERATOPS FPP is dominated by bound-companion scenarios. For
  TIC 89247735 companions closer than the Gaia DR3 sources at 2-5" are unconstrained, and the TRICERATOPS FPP consists
  of the bound-companion (STP) scenario. For TIC 455017950 companions inside 0.15" or fainter than $`\Delta m = 2.74`$
  (0.15-1") are not excluded by the imaging. For TIC 311179599 no imaging exists; companions within about 9" are
  constrained only by the Gaia DR3 astrometric and radial-velocity data, and the TRICERATOPS FPP consists of the
  bound-companion (STP) scenario.
- No radial velocities were obtained for this work. For TIC 402710280 the two public HARPS epochs sample 16 % of the
  orbital velocity amplitude and exclude only stellar and brown-dwarf companions; thirteen further epochs are proprietary until 2027.
  For TIC 324429952, TIC 89247735, TIC 455017950 and TIC 311179599 there are no archival spectra (for TIC 455017950 there are
  catalogue radial velocities of km/s precision only, section [4.6](4_TIC455017950.md#46-companion-limits)).
- QLP timing notice: MAST reports a barycentric timing error in the QLP light curves from the second half of sector 99
  to sector 104 (orbits 210-220), not yet corrected at the date of this record. The mid-transit times of
  sections 1-5 come from TESScut cutouts or SPOC light curves. TOI-2638.01 also uses a QLP sector-10 mid-time in its
  ephemeris and QLP timestamps in its product comparisons, including sector 99. For
  TIC 402710280 the QLP-minus-TESScut time-stamp differences measured in sectors 100 and 102 are within $`\pm 1.7`$ min
  (medians $`+0.3`$ and $`-0.3`$ min), far below the 24-29 min deviations of section [1.3](1_TIC402710280.md#13-transit-epochs-and-o-c). Exception: the mid-transit times of
  TOI-6482.01 (section [6.8](6.8_TOI-6482.md)) come from QLP light curves; in sector 99 the QLP time stamps differ from
  SPOC full-frame-image time stamps referred to the target position by $`+1.71`$ to $`+2.04`$ s (orbits 209 and 210). For TOI-1117 c (section [6.10](6.10_TOI-1117c.md)) the ephemeris comes from SPOC 2-min data; QLP sectors 100-103 enter only the depth and timing cross-checks.
- The 2.3" neighbour of TIC 402710280 (Gaia DR3 5828441289812582656; TIC 1208600633, without TIC stellar parameters) is
  6.4 mag fainter and lies in the background. TRICERATOPS treats it with solar parameters by default; the run with
  K-dwarf parameters gives the same FPP within the seed-to-seed scatter.
- Stellar parameters are TIC v8.2 values and were not re-derived, except for TIC 89247735 and TIC 455017950 (SED fits,
  sections [3.1](3_TIC89247735.md#31-host-star) and [4.1](4_TIC455017950.md#41-host-star)) and TIC 311179599 (Mann et al. relations, section [5.1](5_TIC311179599.md#51-host-star), which reproduce the TIC values); the
  planet radii scale with the adopted stellar radii.
- The per-epoch depths of TIC 324429952 scatter by about 20 % between detrending choices and sectors; the red-noise
  factors reach 2.4. TRICERATOPS uses one flux-error value per light curve and a Gaia DR3 rather than TRILEGAL background.
- The transit-timing variations of TIC 402710280 rest on four epochs with errors of 3-5 min; two further epochs
  (2027-02-19 in TESS sector 114, and any ground-based transit) would test them.
- TIC 89247735: the four transits come from three sectors (the predicted sector-13 transit falls in a data gap, and
  sector 104 starts just after a predicted transit); the period rests on two consecutive transits in sector 93 and the
  long baseline to sectors 66 and 104. One transit (epoch $`-20`$) lies 0.49 d after the start of a sector.
- TIC 89247735: a brown dwarf below about 30 $`M_\mathrm{J}`$ is not excluded; the Gaia DR3 radial-velocity scatter is
  the only constraint on the companion mass, the TESS data give no usable ellipsoidal or beaming limit, and
  TRICERATOPS does not model a brown-dwarf or low-mass-star companion on the target.
- TIC 89247735: the TIC contamination ratio is 0.274, and the raw depth depends on the aperture (0.70-0.92 %). The
  dilution correction rests on a Gaussian PSF model of the Gaia DR3 sources within 60" (sources beyond 60" are not
  included); in the joint fit the target flux fraction is set by its prior from that model, not by the light curve.
  The adopted radius ($`8.7 \pm 0.6\,R_\oplus`$) uses the SED-based stellar radius, whose error (5 %) is dominated by
  the choice of $`T_\mathrm{eff}`$ (the SED fit is formally poor and does not constrain the extinction); with the TIC
  radius it is 9.1-9.4 $`R_\oplus`$. The orbital eccentricity is not measured ($`e < 0.75`$ at 95 %). The TRICERATOPS
  FPP of three seeds spans $`4 \times 10^{-15}`$ to $`2.7 \times 10^{-3}`$.
- TIC 89247735: the search for further transiting planets (section [3.10](3_TIC89247735.md#310-search-for-additional-transiting-signals)) excludes only planets larger than about
  3-4 $`R_\oplus`$; no archival ground-based light curve constrains the transit, and the KELT-South light curves, which
  cover the star, could not be retrieved (section [3.11](3_TIC89247735.md#311-archival-ground-based-photometry)).
- TIC 455017950: the transit of epoch 0 (sector 99) lies at the end of an orbit, 3.6 h before the data end, and its
  depth depends on the detrending; it agrees with the ephemeris but is not an independent confirmation. The period rests
  on the other four transits, two of them (epochs 1 and 2) one cycle apart.
- TIC 455017950: the departure of the mid-transit times from a linear ephemeris rests on one epoch (epoch 2, +8.0 to
  +17 min); its significance depends on the photometry and the baseline model ($`p`$ = 0.001-0.10), and without that
  epoch all timing sets are linear. A transit-timing variation is not established; further transits (TESS sector 114 on
  2027-02-25, or from the ground) would test it.
- TIC 455017950: the star is chromospherically and coronally active ($`P_\mathrm{rot} = 8.3`$ d, spot modulation
  1.5-2.7 %, $`\log L_X/L_\mathrm{bol} \approx -4.5`$, age 0.3-0.5 Gyr). The expected radial-velocity semi-amplitude
  is about 2.1 m/s (10.3 $`M_\oplus`$ from the Chen & Kipping 2017 mass-radius relation), while the activity jitter is
  expected at 3-20 m/s (FF' method, Aigrain et al. 2012, rotational term only) up to 40-65 m/s (photometric amplitude
  times equatorial velocity); a 3 $`\sigma`$ (5 $`\sigma`$) mass of 10 $`M_\oplus`$ would need of order 18-40 (50-110)
  high-precision epochs if the activity model leaves 2-3 m/s. These are estimates, not exposure-time calculations.
- TIC 455017950: ingress and egress are not resolved, so the impact parameter and the shape are not constrained by the
  light curve alone; the adopted radius uses the stellar-density prior (section [4.2](4_TIC455017950.md#42-transits-and-ephemeris)), and the eccentricity is not
  measured. A blended eclipsing binary is excluded statistically (difference images, stellar density, TRICERATOPS), not
  by the transit shape.
- TIC 455017950: bound companions of $`\lesssim 0.4\,M_\odot`$ beyond about 30 au (0.5-1") are not excluded
  (section [4.6](4_TIC455017950.md#46-companion-limits)), and a brown dwarf below about 10-17 $`M_\mathrm{J}`$ transiting the target is constrained only by the
  Gaia radial velocities. TRICERATOPS was run with the SOAR $`I`$-band contrast curve applied in the TESS band and was not
  re-run with the limits of section [4.6](4_TIC455017950.md#46-companion-limits). The origin of the Hipparcos "suspected non-single" flag is not documented.
- TIC 455017950: the Gaia DR3 rotation period (3.39 d) does not appear in the TESS data; the search for further
  transiting planets (section [4.10](4_TIC455017950.md#410-search-for-additional-transiting-signals)) is limited by the spot modulation and excludes only planets larger than about
  2-3 $`R_\oplus`$ inside 30 d.
- TIC 311179599: the mid-transit times do not follow a linear ephemeris. The sector-39 transit, 2.2 h earlier than the
  ephemeris of the other seven epochs, is a single 10-min-cadence epoch in a sector with strong scattered light; a noise
  dip combined with a transit hidden by noise at the predicted time has a probability of order $`10^{-3}`$ per epoch
  from the red-noise and injection statistics of that sector, which rest on 300 injections and so sample tail
  probabilities below about 1 % poorly. The other offsets (up to 49 min) lie at or near the 95th percentiles of the
  injection distributions, whose tails are not Gaussian. A transit-timing variation is not established; the adopted
  ephemeris excludes sector 39, and predicted times need windows of at least $`\pm 1`$ h. The two transits of TESS
  sectors 113 and 114 (2027-01-24, 2027-02-18) would test it.
- TIC 311179599: the TRICERATOPS runs use six of the eight epochs (the 30-min sectors 12 and 13 are left out) and no
  contrast curve; a brown dwarf below about 20 $`M_\mathrm{J}`$ on the transit period is constrained only by the Gaia
  radial velocities. The TARS light curves were not examined.
- TIC 311179599, candidate second signal c: three events, one of them at 30-min cadence; $`P = 36.93`$ d rests on the
  sector-101 epoch not being counted (scattered light) and is tested by the transit of 2026-11-02. The impact parameter
  is not constrained, the joint-fit posterior has a grazing tail ($`b`$ up to 0.94 at 97.5 %), and the event depths are
  5-6 $`\sigma`$ each. The chance-alignment estimate (about 1 %) rests on zero or one comparable event in 223 d of data and
  is an order of magnitude. The TRICERATOPS NFPP (0.0014-0.0086) exceeds the limit of Giacalone et al. (2021) for a
  likely planet, and the multiplicity boost of Guerrero et al. (2021) was derived for SPOC 2-min TOIs, not for
  candidates found in full-frame images; the boosted values of section [5.10.4](5_TIC311179599.md#5104-triceratops-and-the-multiplicity-of-candidates) are indicative only. The comparison with
  the timing of candidate b uses eight epochs with non-Gaussian errors. No difference image or solar-system-object check
  was made for the sector-12 event. The checks of section [5.10.3](5_TIC311179599.md#5103-checks) and the period tests come from the vetting record.
- TIC 311179599: no rotation period below 15 d is found, so the radial-velocity jitter of the star is not constrained;
  the expected semi-amplitude is about 2.1 m/s (6.7 $`M_\oplus`$ from the Chen & Kipping 2017 mass-radius relation; an
  estimate). The eccentricity is not measured ($`e < 0.78`$ at 95 %), and the adopted radius uses the stellar-density
  prior (section [5.2](5_TIC311179599.md#52-transits-and-ephemeris)).
- TIC 11613676: the depth history rests on five TESS sectors, with no data between 2020.8 and 2025.7. No model of
  section [6.1.6](6.1_TIC11613676.md#616-checks-and-geometric-models) fits the light curves formally ($`\chi^2/N \ge 2.7`$); the adopted model (N4) does not reproduce the
  2020 eclipses, the geometric parameters are degenerate at this grazing geometry, the third light was fixed
  at zero, and the nodal period (about 23 yr) and the outer-period estimates (90-170 d) are orders of magnitude. A third
  body is inferred from the eclipse history, not detected; the Gaia DR3 RUWE and radial-velocity amplitude are consistent
  with one but are not measurements of it. The shift of the secondary eclipse within sector 97 is not explained. The Gaia
  DR3 parallax is not reliable (RUWE 48); the distance of about 100 pc is photometric. The eclipse times do not follow a
  linear ephemeris ($`\chi^2 = 6247`$ for 7 degrees of freedom), so predicted times need windows of at least $`\pm 30`$ min.
- TIC 124751008: the depth history rests on seven TESS sectors with no data between 2021.0 and 2025.0. The source
  position of the shallowest eclipses (sector 87) is only weakly located by the difference images; in sector 108 the
  bleed column of $`\kappa`$ Col passes 1.6 px from the target, four aperture pixels are masked and the depths of the
  second sector-108 eclipse differ by 0.6-1.2 % between photometry methods; the sector-7 eclipse of cycle $`-198`$ lies
  on an orbit-start ramp; sector 109 covers 6.9 d. The geometry is degenerate at this grazing configuration ($`k`$ = 0.57
  and 0.78 fit equally; $`r_1`$, $`J`$ and $`e\sin\omega`$ are correlated), the global smooth $`i(t)`$ models did not
  converge, the per-eclipse inclination errors ignore the shared geometry, and the nodal period (about 20 yr), the
  outer-period estimates (110-185 d) and the depth predictions are orders of magnitude for stated assumptions; the
  outer orbit was not fitted. A third body is inferred, not detected. The third light (0 or 0.70) is assumed. No radial
  velocity of the star exists, and no eclipse depth before 2018 is known.
- TOI-540 b: the non-detection rests on one TESS sector, and a data or processing problem specific to this star and
  sector cannot be excluded by the checks of section [6.2.5](6.2_TOI-540b.md#625-checks-for-instrumental-and-processing-causes); one ground-based transit would test it (section [6.2.9](6.2_TOI-540b.md#629-decisive-test)).
  The change of the impact parameter between 2018 and 2020 is 1.2-2.9 $`\sigma`$ depending on the photometry, and the
  companion masses of section [6.2.7](6.2_TOI-540b.md#627-interpretation) are order-of-magnitude estimates for stated assumptions; no companion is detected.
  Ground-based transit observations made after 2019 but not posted to ExoFOP or published cannot be checked.
- The duration-change evidence for L 98-59 d and TOI-620 b (sections [6.6](6.6_L98-59d.md) and [6.7](6.7_TOI-620b.md)) rests on few observing campaigns of one instrument. For L 98-59 d, the within-campaign diagnostic is only 1.8 $`\sigma`$ and JWST does not independently confirm the trend. Its quoted trials-adjusted significance assumes Gaussian tails and 801 or 3204 independent trials; TOI-620 b has no corresponding adjustment.
- TOI-2638.01 (section [6.9](6.9_TOI-2638.md)): the depth change rests on seven TESS sectors of one instrument (sector 10 from QLP only, which is
  not deblended here), and the sector-99 depth differs between products (3.2-4.4 ppt). The shape change that separates a
  change of the transit geometry from a change of dilution is about 2.6-3.1 $`\sigma`$ when calibrated by simulation; the
  bootstrap error of the duration slope underestimates its scatter. At this grazing geometry $`k`$, $`b`$ and $`a/R_\star`$ are
  degenerate, so the impact parameters, the end of the grazing phase and the predicted depths depend on the model; no
  gravity-darkened model was fitted. The ellipsoidal limit uses PDCSAP fluxes. The precession estimate is an order of
  magnitude for assumed apsidal constant and moment of inertia, and the stellar obliquity is not measured. The
  ground-based photometry, imaging and spectroscopy listed on ExoFOP were not used.
- TOI-1117 c (section [6.10](6.10_TOI-1117c.md)): the signal rests on TESS data (six sectors, 27 transits, S/N 12-14) and the
  published HARPS velocities; it is not validated. In the pixel data a $`G = 19.59`$ star at 34" is disfavoured only at
  $`\Delta\chi^2 = 7.5`$ (summed over six pixel sets; it fits better than the target in two of them), with a Gaussian-mixture
  PSF rather than the TESS PRF. The adopted transit needs $`b = 0.85`$ with the stellar density of Lockley et al. (2025);
  the fit without that prior, also using uniform limb-darkening parameters, favours a boxier shape, with $`\chi^2`$ lower
  by 3.5 (not calibrated as a statistical significance); eccentricity was fixed at zero. Odd and even depths
  differ by 2.1 $`\sigma`$, and the sector-13 transits are 13-22 min early (1.2-2.0 $`\sigma`$ with empirical timing errors).
  The RV conjunction constrains the phase only to $`\pm 0.11`$ d. The TRICERATOPS runs use an approximate contrast curve
  written from the text of Lockley et al. (2025). ExoFOP was unreachable during all checks (2026-10-22 23:22 UTC to
  2026-10-23 00:48 UTC), so ExoFOP notes, files and live TOI/community-TOI lists are unchecked; the arXiv API check returned
  no records even for an alias present in a known abstract. The transit and RV fits use the stellar parameters of Lockley
  et al. (2025); TRICERATOPS uses TIC stellar parameters.
- Eclipsing binaries of section [6.3](6.3_eclipsing_binaries.md): no spectroscopy; the periods come from TESS full-frame images, with QLP eclipse
  times in sectors 99-104 (timing notice above); TIC 13957664 and TIC 342194997 rest on one eclipse each and their
  periods on the Gaia orbits, TIC 304936165 on two eclipses, and the period of CPD-62 2890 is not determined; the
  companion masses are estimates for $`i \approx 90`$ deg and the TIC primary masses (TIC 342194997: 0.68-0.76
  $`M_\odot`$); for TIC 342194997 the companion radius is about 20-25 % larger than expected for its mass, for
  $`R_1`$ = 0.69-0.72 $`R_\odot`$ the eclipse duration implies a transverse velocity 1.6-2.8 $`\sigma`$ above the orbital
  value, and the predicted occultation was not observed (sector-99 data gap). The values come from the vetting record
  and are not recomputed by a script in this repository.
- Nothing in this record has been checked by anyone outside the production of the record.
