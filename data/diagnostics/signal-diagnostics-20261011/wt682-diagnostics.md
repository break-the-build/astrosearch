# Signal diagnostics: TIC 366989133 (WT 682), SPOC, Profile J five-pass

Report-only. Candidate rows are unchanged; verdicts name competing explanations, not classifications.

## Event pattern 1 (11 detections)

| Detection | Period (d) | Events | Reported depth (ppm) | Estimator | Box depth | Local depth |
|---|---|---|---|---|---|---|
| astrosearch:21 | 0.18772020 | 470 | 1.9 | local (transit-masked detrending) | 2837.1 | 1.9 |
| astrosearch:22 | 0.18769245 | 349 | 2694.9 | box (local not positive) | 2694.9 | -115.7 |
| astrosearch:25 | 0.18775143 | 449 | 84.2 | local (transit-masked detrending) | 2089.5 | 84.2 |
| astrosearch:26 | 0.37524452 | 189 | 2947.1 | box (local not positive) | 2947.1 | -587.5 |
| astrosearch:27 | 0.18762111 | 340 | 2232.7 | box (local not positive) | 2232.7 | -136.9 |
| astrosearch:28 | 0.18767808 | 394 | 2090.5 | box (local not positive) | 2090.5 | -1.1 |
| astrosearch:29 | 0.18774796 | 384 | 53.4 | local (transit-masked detrending) | 2055.1 | 53.4 |
| astrosearch:30 | 0.18771971 | 597 | 365.1 | local (transit-masked detrending) | 2277.3 | 365.1 |
| astrosearch:31 | 0.18776184 | 379 | 25.2 | local (transit-masked detrending) | 2044.0 | 25.2 |
| astrosearch:32 | 0.18762210 | 413 | 1947.8 | box (local not positive) | 1947.8 | -186.5 |
| astrosearch:33 | 0.37564298 | 204 | 769.0 | local (transit-masked detrending) | 3286.6 | 769.0 |

- astrosearch:21 and astrosearch:22: 277 shared events (79% of the smaller set), same-period
- astrosearch:21 and astrosearch:25: 357 shared events (80% of the smaller set), same-period
- astrosearch:21 and astrosearch:26: 148 shared events (78% of the smaller set), 2:1 harmonic
- astrosearch:21 and astrosearch:27: 269 shared events (79% of the smaller set), same-period
- astrosearch:21 and astrosearch:28: 256 shared events (65% of the smaller set), same-period
- astrosearch:21 and astrosearch:29: 270 shared events (70% of the smaller set), same-period
- astrosearch:21 and astrosearch:31: 264 shared events (70% of the smaller set), same-period
- astrosearch:21 and astrosearch:32: 313 shared events (76% of the smaller set), same-period
- astrosearch:21 and astrosearch:33: 153 shared events (75% of the smaller set), 2:1 harmonic
- astrosearch:22 and astrosearch:25: 223 shared events (64% of the smaller set), same-period
- astrosearch:22 and astrosearch:26: 104 shared events (55% of the smaller set), 2:1 harmonic
- astrosearch:22 and astrosearch:28: 256 shared events (73% of the smaller set), same-period
- astrosearch:22 and astrosearch:29: 260 shared events (74% of the smaller set), same-period
- astrosearch:22 and astrosearch:31: 257 shared events (74% of the smaller set), same-period
- astrosearch:22 and astrosearch:32: 180 shared events (52% of the smaller set), same-period
- astrosearch:22 and astrosearch:33: 154 shared events (76% of the smaller set), 2:1 harmonic
- astrosearch:25 and astrosearch:27: 239 shared events (70% of the smaller set), same-period
- astrosearch:25 and astrosearch:28: 204 shared events (52% of the smaller set), same-period
- astrosearch:25 and astrosearch:29: 225 shared events (59% of the smaller set), same-period
- astrosearch:25 and astrosearch:31: 215 shared events (57% of the smaller set), same-period
- astrosearch:25 and astrosearch:32: 285 shared events (69% of the smaller set), same-period
- astrosearch:25 and astrosearch:33: 155 shared events (76% of the smaller set), 2:1 harmonic
- astrosearch:26 and astrosearch:27: 150 shared events (79% of the smaller set), 2:1 harmonic
- astrosearch:26 and astrosearch:28: 122 shared events (65% of the smaller set), 2:1 harmonic
- astrosearch:26 and astrosearch:32: 147 shared events (78% of the smaller set), 2:1 harmonic
- astrosearch:27 and astrosearch:28: 200 shared events (59% of the smaller set), same-period
- astrosearch:27 and astrosearch:32: 266 shared events (78% of the smaller set), same-period
- astrosearch:28 and astrosearch:29: 246 shared events (64% of the smaller set), same-period
- astrosearch:28 and astrosearch:31: 247 shared events (65% of the smaller set), same-period
- astrosearch:28 and astrosearch:33: 102 shared events (50% of the smaller set), 2:1 harmonic
- astrosearch:29 and astrosearch:31: 263 shared events (69% of the smaller set), same-period
- astrosearch:29 and astrosearch:33: 112 shared events (55% of the smaller set), 2:1 harmonic
- astrosearch:30 and astrosearch:33: 150 shared events (74% of the smaller set), 2:1 harmonic
- astrosearch:31 and astrosearch:33: 103 shared events (50% of the smaller set), 2:1 harmonic
- Depth check: **INCONSISTENT**
  - astrosearch:21: depth 2837.1 ppm but only 49.4 ppm at this ephemeris in the unmasked light curve (found in search pass 5 with 55.2% of samples masked); the dip depends on masking or detrending, not the data alone
  - astrosearch:22: depth 2694.9 ppm but only 177.3 ppm at this ephemeris in the unmasked light curve (found in search pass 5 with 55.2% of samples masked); the dip depends on masking or detrending, not the data alone
  - astrosearch:25: depth 2089.5 ppm but only 115.6 ppm at this ephemeris in the unmasked light curve (found in search pass 5 with 55.2% of samples masked); the dip depends on masking or detrending, not the data alone
  - astrosearch:26: depth 2947.1 ppm but only 438.6 ppm at this ephemeris in the unmasked light curve (found in search pass 5 with 55.2% of samples masked); the dip depends on masking or detrending, not the data alone
  - astrosearch:27: depth 2232.7 ppm but only 18.6 ppm at this ephemeris in the unmasked light curve (found in search pass 5 with 55.2% of samples masked); the dip depends on masking or detrending, not the data alone
  - astrosearch:28: depth 2090.5 ppm but only 141.0 ppm at this ephemeris in the unmasked light curve (found in search pass 5 with 55.2% of samples masked); the dip depends on masking or detrending, not the data alone
  - astrosearch:29: depth 2055.1 ppm but only 99.7 ppm at this ephemeris in the unmasked light curve (found in search pass 5 with 55.2% of samples masked); the dip depends on masking or detrending, not the data alone
  - astrosearch:30: depth 2277.3 ppm but only 315.7 ppm at this ephemeris in the unmasked light curve (found in search pass 4 with 36.5% of samples masked); the dip depends on masking or detrending, not the data alone
  - astrosearch:31: depth 2044.0 ppm but only 246.1 ppm at this ephemeris in the unmasked light curve (found in search pass 5 with 55.2% of samples masked); the dip depends on masking or detrending, not the data alone
  - astrosearch:32: depth 1947.8 ppm but only -18.6 ppm at this ephemeris in the unmasked light curve (found in search pass 5 with 55.2% of samples masked); the dip depends on masking or detrending, not the data alone
  - astrosearch:21: local depth 1.9 ppm is 0.1% of box depth 2837.1 ppm (dip no sharper than nearby variation)
  - astrosearch:22: local depth -115.7 ppm is -4.3% of box depth 2694.9 ppm (dip no sharper than nearby variation)
  - astrosearch:25: local depth 84.2 ppm is 4.0% of box depth 2089.5 ppm (dip no sharper than nearby variation)
  - astrosearch:26: local depth -587.5 ppm is -19.9% of box depth 2947.1 ppm (dip no sharper than nearby variation)
  - astrosearch:27: local depth -136.9 ppm is -6.1% of box depth 2232.7 ppm (dip no sharper than nearby variation)
  - astrosearch:28: local depth -1.1 ppm is -0.1% of box depth 2090.5 ppm (dip no sharper than nearby variation)
  - astrosearch:29: local depth 53.4 ppm is 2.6% of box depth 2055.1 ppm (dip no sharper than nearby variation)
  - astrosearch:30: local depth 365.1 ppm is 16.0% of box depth 2277.3 ppm (dip no sharper than nearby variation)
  - astrosearch:31: local depth 25.2 ppm is 1.2% of box depth 2044.0 ppm (dip no sharper than nearby variation)
  - astrosearch:32: local depth -186.5 ppm is -9.6% of box depth 1947.8 ppm (dip no sharper than nearby variation)
  - astrosearch:33: local depth 769.0 ppm is 23.4% of box depth 3286.6 ppm (dip no sharper than nearby variation)
  - reported depths in one event pattern span 1.9-2947.1 ppm (1557x); they come from different estimators: astrosearch:21 1.9 ppm [local (transit-masked detrending)]; astrosearch:22 2694.9 ppm [box (local not positive)]; astrosearch:25 84.2 ppm [local (transit-masked detrending)]; astrosearch:26 2947.1 ppm [box (local not positive)]; astrosearch:27 2232.7 ppm [box (local not positive)]; astrosearch:28 2090.5 ppm [box (local not positive)]; astrosearch:29 53.4 ppm [local (transit-masked detrending)]; astrosearch:30 365.1 ppm [local (transit-masked detrending)]; astrosearch:31 25.2 ppm [local (transit-masked detrending)]; astrosearch:32 1947.8 ppm [box (local not positive)]; astrosearch:33 769.0 ppm [local (transit-masked detrending)]

## Event pattern 2 (1 detections)

| Detection | Period (d) | Events | Reported depth (ppm) | Estimator | Box depth | Local depth |
|---|---|---|---|---|---|---|
| astrosearch:24 | 0.37544072 | 299 | 749.4 | local (transit-masked detrending) | 4145.4 | 749.4 |

- Depth check: **INCONSISTENT**
  - astrosearch:24: local depth 749.4 ppm is 18.1% of box depth 4145.4 ppm (dip no sharper than nearby variation)

## Event pattern 3 (3 detections)

| Detection | Period (d) | Events | Reported depth (ppm) | Estimator | Box depth | Local depth |
|---|---|---|---|---|---|---|
| tls:1 | 0.37544360 | 295 | 3864.6 | tls model depth |  |  |
| tls:2 | 0.37544065 | 291 | 3890.6 | tls model depth |  |  |
| tls:3 | 0.37544065 | 294 | 3854.0 | tls model depth |  |  |

- tls:1 and tls:2: 178 shared events (61% of the smaller set), same-period
- tls:1 and tls:3: 182 shared events (62% of the smaller set), same-period
- tls:2 and tls:3: 288 shared events (99% of the smaller set), same-period
- Depth check: **consistent**

## Variability conflicts

| Detection | Period (d) | Modulation period (d) | Depth before (ppm) | Depth after (ppm) | Retained | Verdict |
|---|---|---|---|---|---|---|
| astrosearch:21 | 0.18772020 | 0.37505136 | 49.4 | -91.1 | -184% | competing explanation: periodic modulation accounts for most of the dip (stellar variability, ellipsoidal or reflection modulation, or an alias of it) |
| astrosearch:22 | 0.18769245 | 0.37505136 | 177.3 | -53.7 | -30% | competing explanation: periodic modulation accounts for most of the dip (stellar variability, ellipsoidal or reflection modulation, or an alias of it) |
| astrosearch:24 | 0.37544072 | 0.37505136 | 2977.9 | 35.7 | 1% | competing explanation: periodic modulation accounts for most of the dip (stellar variability, ellipsoidal or reflection modulation, or an alias of it) |
| astrosearch:25 | 0.18775143 | 0.37505136 | 115.6 | -28.9 | -25% | competing explanation: periodic modulation accounts for most of the dip (stellar variability, ellipsoidal or reflection modulation, or an alias of it) |
| astrosearch:26 | 0.37524452 | 0.37505136 | 438.6 | -163.2 | -37% | competing explanation: periodic modulation accounts for most of the dip (stellar variability, ellipsoidal or reflection modulation, or an alias of it) |
| astrosearch:27 | 0.18762111 | 0.37505136 | 18.6 | -78.1 | -419% | competing explanation: periodic modulation accounts for most of the dip (stellar variability, ellipsoidal or reflection modulation, or an alias of it) |
| astrosearch:28 | 0.18767808 | 0.37505136 | 141.0 | -42.8 | -30% | competing explanation: periodic modulation accounts for most of the dip (stellar variability, ellipsoidal or reflection modulation, or an alias of it) |
| astrosearch:29 | 0.18774796 | 0.37505136 | 99.7 | -137.8 | -138% | competing explanation: periodic modulation accounts for most of the dip (stellar variability, ellipsoidal or reflection modulation, or an alias of it) |
| astrosearch:30 | 0.18771971 | 0.37505136 | 315.7 | 51.0 | 16% | competing explanation: periodic modulation accounts for most of the dip (stellar variability, ellipsoidal or reflection modulation, or an alias of it) |
| astrosearch:31 | 0.18776184 | 0.37505136 | 246.1 | -96.4 | -39% | competing explanation: periodic modulation accounts for most of the dip (stellar variability, ellipsoidal or reflection modulation, or an alias of it) |
| astrosearch:32 | 0.18762210 | 0.37505136 | -18.6 | -28.4 |  | no positive dip at this ephemeris before modulation removal |
| astrosearch:33 | 0.37564298 | 0.37505136 | 1782.2 | 168.5 | 9% | competing explanation: periodic modulation accounts for most of the dip (stellar variability, ellipsoidal or reflection modulation, or an alias of it) |
| tls:1 | 0.37544360 | 0.37505136 | 3544.2 | 40.5 | 1% | competing explanation: periodic modulation accounts for most of the dip (stellar variability, ellipsoidal or reflection modulation, or an alias of it) |
| tls:2 | 0.37544065 | 0.37505136 | 3843.1 | 335.6 | 9% | competing explanation: periodic modulation accounts for most of the dip (stellar variability, ellipsoidal or reflection modulation, or an alias of it) |
| tls:3 | 0.37544065 | 0.37505136 | 3674.4 | 164.5 | 4% | competing explanation: periodic modulation accounts for most of the dip (stellar variability, ellipsoidal or reflection modulation, or an alias of it) |

