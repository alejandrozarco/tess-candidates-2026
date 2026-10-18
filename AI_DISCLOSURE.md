# AI disclosure

AI models produced the content of this repository: the analysis scripts, the tables, the figures and the text. The
repository owner chose the targets to pursue, directed the work and decided on scope and release. The owner did not check
the code or the numbers line by line, and no professional astronomer has reviewed them.

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
targets first showed BJD_TDB printed as a calendar date; they now give geocentric UTC (3-4 min earlier).

**What the software checks and what is trusted.** The scripts recompute the tables from public products (TESScut, SPOC,
QLP, TESS-SPOC, TARS, Gaia DR3, TIC, ESO archive); exceptions are listed in README.md (Data and reproduction). Stellar
parameters are TIC v8.2 values and are not re-derived, except for TIC 89247735 (SED fit; README section 3.1). The TRICERATOPS false-positive probabilities depend on that
code's priors and on a single flux-error value per light curve; no high-resolution imaging constraint was available.

**What this means for use.** The candidates are leads to be checked independently; nothing here is confirmed. Credit
for confirming, refining or refuting any of them belongs to whoever does that work. Questions, checks and corrections:
[GitHub issues](https://github.com/alejandrozarco/tess-candidates-2026/issues).
