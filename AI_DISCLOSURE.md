# AI disclosure

AI models produced the content of this repository: the analysis scripts, the tables, the figures and the text. The
content covers the five transit candidates, a candidate second transit signal on TIC 311179599 (README section 5.10)
and, in README section 6, the eclipsing binaries TIC 11613676 and TIC 124751008 (section 6.4), the sector-98
observation of TOI-540 b and nine new eclipsing binaries (section 6.3, values from the vetting record). The repository owner chose the targets to pursue, directed the work and decided on
scope and release. The owner did not check the code or the numbers line by line, and no professional astronomer has
reviewed them.

**Models.** Claude models (Anthropic), run through Claude Code, for the detection, the vetting, the TRICERATOPS runs,
the archival-spectrum analysis and this package. Internal cross-checks of the numbers were also made by AI models; such
checks are not peer review.

**Corrections made during the work.** The detection-stage transit depth of TIC 402710280 (1835 ± 36 ppm, 3 × 3-pixel box
photometry, five events) was superseded by the circular-aperture common-shape fit on four events (about 1680 ppm); the
sector-66 event was moved out of the ephemeris because it lies in scattered light. Formal mid-transit-time errors were
replaced by time-averaging-scaled and prayer-bead errors. The QLP light curves were dropped as a source of transit times
after the MAST notice of a timing error in sectors 99-104. For TIC 89247735 the detection-stage period (18.323343 ±
0.00012 d, box search on 3 × 3-pixel photometry) was superseded by the common-shape fit of the four transits (18.323173
± 0.000031 d); the detection-stage radius (8.5 R⊕, 3 × 3-pixel box, no correction for the light of nearby stars) was
superseded by the dilution-corrected value (9.1-9.4 R⊕), and the radius 9.1-9.4 R⊕ (TIC stellar radius) was superseded
by 8.7 ± 0.6 R⊕ after an SED fit of the host star; the TRICERATOPS runs were repeated on the circular-aperture
photometry (the detection-stage runs on 3 × 3-pixel photometry gave FPP of 10⁻⁷ to 10⁻⁵). The UTC columns of the predicted-transit tables of all three
targets first showed BJD_TDB printed as a calendar date; they now give geocentric UTC (3-4 min earlier). For TIC 455017950 the
detection-stage rotation period (3.39 d, the Gaia DR3 `vari_rotation_modulation` value) was superseded by 8.3 d, measured
in all eight TESS sectors; the detection-stage description of the transit as flat-bottomed is not supported by the data
(ingress and egress are not resolved, and diluted eclipses of larger bodies fit equally well); the detection-stage FPP
(0.006-0.013 without and 0.0008-0.0012 with the SOAR contrast curve, 3 seeds each) and depth (about 1460 ppm, the
limb-darkened central depth of the raw aperture flux) were superseded by the TRICERATOPS runs of README section 4.7 and
the dilution-corrected density-prior fit of section 4.2, and the detection-stage period (40.605645 ± 0.000052 d) by the
common-shape fit of the five transits (40.605660 ± 0.000074 d). For TIC 311179599 the detection-stage period
(24.826807 ± 0.00016 d, 3 × 3-pixel box photometry, seven epochs) was superseded by the common-shape fit
(24.826839 ± 0.000047 d, the same seven epochs), the detection-stage radius (about 2.3 R⊕, TIC radius, no dilution
correction) by 2.47 ± 0.10 R⊕ (density-prior fit with free dilution, README section 5.2), the detection-stage
odd-even depth difference (1.8 σ, 3 × 3-pixel box) by 0.8 σ, and the detection-stage TRICERATOPS FPP (0.036-0.049, NFPP
about 6 × 10⁻⁶; r < 1.5 px apertures of three sectors) by 0.025-0.031 (NFPP about 7 × 10⁻⁷; r < 2.0 px apertures of
six sectors, README section 5.7; 0.020-0.028 with the photometry of the vetting re-check). For TIC 342194997 (README
section 6.3) the detection-stage description of the eclipse as flat-bottomed was superseded (a model without limb darkening fits worse; the eclipse is limb-darkened), the detection-stage
statement that no occultation was seen at the other conjunction "consistent with the geometry" was superseded (at
that conjunction an occultation of the companion is geometrically required; the orbit conditioned on the sector-99
non-detection places it in the sector-99 data gap), and the detection-stage depth (15.4 %, box) by the model central
depth (16.9 %); its predicted conjunction (observed − predicted +2.0 ± 3.8 d) was recomputed with the parameter order
of the Gaia DR3 documentation for the correlation vector (predicted BTJD 4089.06 ± 3.46, z = +0.57). For TIC 124751008
the detection-stage description of the eclipses as fading was superseded by "decreased and then increased", and the
third-light inclination fits were recomputed after the per-eclipse inclination grid was found to truncate the 2019
values (nodal period 20.2 ± 4.5 yr instead of 12.0 ± 1.8 yr); the candidate-c fit of TIC 311179599 was repeated with a
longer sampler chain on the TESScut time stamps (P = 73.86721 instead of 73.86712 d).

**What the software checks and what is trusted.** The scripts recompute the tables from public products (TESScut, SPOC,
QLP, TESS-SPOC, TARS, Gaia DR3, TIC, ESO archive); exceptions are listed in README.md (Data and reproduction). Stellar
parameters are TIC v8.2 values and are not re-derived, except for TIC 89247735 and TIC 455017950 (SED fits; README
sections 3.1 and 4.1) and TIC 311179599 (Mann et al. relations; README section 5.1). The TRICERATOPS false-positive probabilities depend on that
code's priors and on a single flux-error value per light curve; no high-resolution imaging constraint was available,
except the published SOAR speckle limits for TIC 455017950. The values of the nine eclipsing binaries of README
section 6.3, and the checks, period tests and models of sections 5.10 and 6.4 listed in README section 8, come from the
vetting record and are not recomputed by a script in this repository.

**What this means for use.** The candidates are leads to be checked independently; nothing here is confirmed. Credit
for confirming, refining or refuting any of them belongs to whoever does that work. Questions, checks and corrections:
[GitHub issues](https://github.com/alejandrozarco/tess-candidates-2026/issues).
