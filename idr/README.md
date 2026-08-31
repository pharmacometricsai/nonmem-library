# Indirect response (IDR) model library

NONMEM control streams for the four basic indirect response models of
Dayneka, Garg & Jusko (1993). The canonical diagrams show only the
response turnover compartment and where the drug acts; they do not say
where the driving concentration comes from. These files supply it by
fitting **PK and PD simultaneously** in one control stream.

**4 IDR models x 6 PK backbones x 2 exposure links = 48 control streams.**

---

## The four models

All four share one turnover compartment - zero-order production `KIN`,
first-order loss `KOUT` - and differ only in which term the drug touches
and in which direction:

| Model | Drug acts on | ODE | Response |
|-------|--------------|-----|----------|
| I | production, inhibited | `dR/dt = KIN*(1 - IMAX*EFF) - KOUT*R` | decreases |
| II | elimination, inhibited | `dR/dt = KIN - KOUT*(1 - IMAX*EFF)*R` | increases |
| III | production, stimulated | `dR/dt = KIN*(1 + SMAX*EFF) - KOUT*R` | increases |
| IV | elimination, stimulated | `dR/dt = KIN - KOUT*(1 + SMAX*EFF)*R` | decreases |

with the fractional effect

```
EFF = (CD/XC50)**HILL / (1 + (CD/XC50)**HILL)
```

where `XC50` is `IC50` (models I, II) or `SC50` (models III, IV) and `CD`
is the driving concentration - plasma or effect-site, depending on the
link. `HILL` is estimable; **fix it to 1 to recover the plain Emax form**
drawn on the page.

### Telling them apart

Models I and IV both push the response **down**; models II and III both
push it **up**. Within each pair the distinguishing feature is not the
direction but the *return to baseline*:

* Models **I and III** act on `KIN`. `KOUT` is untouched, so the response
  returns to baseline at a dose-independent rate set by `KOUT`.
* Models **II and IV** act on `KOUT`. The apparent turnover slows (II) or
  accelerates (IV) while drug is present, so the return to baseline is
  **dose-dependent**.

That difference is only visible with more than one dose level and enough
sampling on the washout limb. Without both, model I and model IV (or II
and III) will fit a data set about equally well and the choice between
them is mechanistic, not statistical.

---

## Baseline - the one constraint that matters

`KIN` is **not** a THETA. It is computed as

```
KIN    = R0*KOUT
A_0(R) = R0
```

so the response starts at the baseline `R0` and, once drug is gone,
returns to it. This is what makes `R0` interpretable and keeps the model
from drifting off baseline in the pre-dose period. Do not add a separate
baseline THETA on top of it, and do not free `KIN` - if you do, `R0` stops
being the baseline and the initial condition becomes inconsistent with the
steady state.

---

## Contents

| # | File | IDR | PK | Route | Link |
|---|------|-----|----|-------|------|
| 1 | `idr1_1cmt_iv_direct.ctl` | I | 1-cmt | IV | direct |
| 2 | `idr1_1cmt_iv_effect.ctl` | I | 1-cmt | IV | effect cmt |
| 3 | `idr1_1cmt_oral_direct.ctl` | I | 1-cmt | oral | direct |
| 4 | `idr1_1cmt_oral_effect.ctl` | I | 1-cmt | oral | effect cmt |
| 5 | `idr1_2cmt_iv_direct.ctl` | I | 2-cmt | IV | direct |
| 6 | `idr1_2cmt_iv_effect.ctl` | I | 2-cmt | IV | effect cmt |
| 7 | `idr1_2cmt_oral_direct.ctl` | I | 2-cmt | oral | direct |
| 8 | `idr1_2cmt_oral_effect.ctl` | I | 2-cmt | oral | effect cmt |
| 9 | `idr1_3cmt_iv_direct.ctl` | I | 3-cmt | IV | direct |
| 10 | `idr1_3cmt_iv_effect.ctl` | I | 3-cmt | IV | effect cmt |
| 11 | `idr1_3cmt_oral_direct.ctl` | I | 3-cmt | oral | direct |
| 12 | `idr1_3cmt_oral_effect.ctl` | I | 3-cmt | oral | effect cmt |
| 13 | `idr2_1cmt_iv_direct.ctl` | II | 1-cmt | IV | direct |
| 14 | `idr2_1cmt_iv_effect.ctl` | II | 1-cmt | IV | effect cmt |
| 15 | `idr2_1cmt_oral_direct.ctl` | II | 1-cmt | oral | direct |
| 16 | `idr2_1cmt_oral_effect.ctl` | II | 1-cmt | oral | effect cmt |
| 17 | `idr2_2cmt_iv_direct.ctl` | II | 2-cmt | IV | direct |
| 18 | `idr2_2cmt_iv_effect.ctl` | II | 2-cmt | IV | effect cmt |
| 19 | `idr2_2cmt_oral_direct.ctl` | II | 2-cmt | oral | direct |
| 20 | `idr2_2cmt_oral_effect.ctl` | II | 2-cmt | oral | effect cmt |
| 21 | `idr2_3cmt_iv_direct.ctl` | II | 3-cmt | IV | direct |
| 22 | `idr2_3cmt_iv_effect.ctl` | II | 3-cmt | IV | effect cmt |
| 23 | `idr2_3cmt_oral_direct.ctl` | II | 3-cmt | oral | direct |
| 24 | `idr2_3cmt_oral_effect.ctl` | II | 3-cmt | oral | effect cmt |
| 25 | `idr3_1cmt_iv_direct.ctl` | III | 1-cmt | IV | direct |
| 26 | `idr3_1cmt_iv_effect.ctl` | III | 1-cmt | IV | effect cmt |
| 27 | `idr3_1cmt_oral_direct.ctl` | III | 1-cmt | oral | direct |
| 28 | `idr3_1cmt_oral_effect.ctl` | III | 1-cmt | oral | effect cmt |
| 29 | `idr3_2cmt_iv_direct.ctl` | III | 2-cmt | IV | direct |
| 30 | `idr3_2cmt_iv_effect.ctl` | III | 2-cmt | IV | effect cmt |
| 31 | `idr3_2cmt_oral_direct.ctl` | III | 2-cmt | oral | direct |
| 32 | `idr3_2cmt_oral_effect.ctl` | III | 2-cmt | oral | effect cmt |
| 33 | `idr3_3cmt_iv_direct.ctl` | III | 3-cmt | IV | direct |
| 34 | `idr3_3cmt_iv_effect.ctl` | III | 3-cmt | IV | effect cmt |
| 35 | `idr3_3cmt_oral_direct.ctl` | III | 3-cmt | oral | direct |
| 36 | `idr3_3cmt_oral_effect.ctl` | III | 3-cmt | oral | effect cmt |
| 37 | `idr4_1cmt_iv_direct.ctl` | IV | 1-cmt | IV | direct |
| 38 | `idr4_1cmt_iv_effect.ctl` | IV | 1-cmt | IV | effect cmt |
| 39 | `idr4_1cmt_oral_direct.ctl` | IV | 1-cmt | oral | direct |
| 40 | `idr4_1cmt_oral_effect.ctl` | IV | 1-cmt | oral | effect cmt |
| 41 | `idr4_2cmt_iv_direct.ctl` | IV | 2-cmt | IV | direct |
| 42 | `idr4_2cmt_iv_effect.ctl` | IV | 2-cmt | IV | effect cmt |
| 43 | `idr4_2cmt_oral_direct.ctl` | IV | 2-cmt | oral | direct |
| 44 | `idr4_2cmt_oral_effect.ctl` | IV | 2-cmt | oral | effect cmt |
| 45 | `idr4_3cmt_iv_direct.ctl` | IV | 3-cmt | IV | direct |
| 46 | `idr4_3cmt_iv_effect.ctl` | IV | 3-cmt | IV | effect cmt |
| 47 | `idr4_3cmt_oral_direct.ctl` | IV | 3-cmt | oral | direct |
| 48 | `idr4_3cmt_oral_effect.ctl` | IV | 3-cmt | oral | effect cmt |

Plus this `README.md`.

### Naming convention

```
idr<1-4>_<n>cmt_<route>_<link>.ctl
```

* `<route>`: `iv` (bolus or infusion) or `oral` (first-order depot)
* `<link>`: `direct` (plasma drives the effect) or `effect`
  (effect compartment with `KE0`)

---

## The two exposure links

**`direct`** - plasma concentration drives `EFF` immediately. Use this
first. An indirect response model already produces a delay between
concentration and effect, purely from the turnover of `R`.

**`effect`** - an effect compartment sits between plasma and the response:

```
dCe/dt = KE0*(C - Ce)
```

The effect compartment is **unscaled**: `A(effect)` holds the effect-site
concentration directly, not an amount, so it takes no drug mass out of the
central compartment. The response compartment is unscaled the same way -
`A(response)` is `R` itself.

A word of warning: turnover and effect-site equilibration are two ways of
producing the same delay, and `KE0` and `KOUT` are often poorly separated.
Reach for the `effect` variants only when a `direct` fit leaves obvious
hysteresis in the concentration-effect plot.

---

## PK backbones

| Structure | Central | Periph. 1 | Periph. 2 |
|-----------|---------|-----------|-----------|
| 1-compartment | `V` | - | - |
| 2-compartment | `V2` | `Q` / `V3` | - |
| 3-compartment | `V2` | `Q3` / `V3` | `Q4` / `V4` |

Because every file here is an ODE model, **these names do not shift**
between the IV and the oral versions - unlike the closed-form ADVANs in
`../1cmt`, `../2cmt` and `../3cmt`, where adding a depot renumbers
everything. The compartment map is written at the top of each file.

Oral versions add `KA` and a bioavailability `F1`, bounded on
(0.001, 0.999) with logit-scale IIV.

---

## Data set

Comma-delimited, at `../data/pkpddata.csv`, commented header row
(`IGNORE=@`), columns:

```
ID TIME AMT RATE DV EVID MDV CMT FLAG ADDL II WT
```

| Column | Meaning |
|--------|---------|
| `ID` | subject identifier |
| `TIME` | time of the event (h) |
| `AMT` | dose amount (mg) |
| `RATE` | 0 = bolus / instantaneous; > 0 = IV infusion rate |
| `DV` | observation - which analyte is set by `FLAG` |
| `EVID` | 0 = observation, 1 = dose |
| `MDV` | 0 = observation, 1 = missing dependent variable |
| `CMT` | see below |
| `FLAG` | 1 = plasma drug concentration, 2 = response |
| `ADDL`, `II` | additional-dose bookkeeping |
| `WT` | body weight (kg) - carried for covariate work, not yet used |

### `CMT`

* **Dose records**: `CMT=1` for the oral models (depot); the central
  compartment for the IV models (also 1, since there is no depot).
* **Observation records**: `CMT` = the central compartment for **both**
  analytes. `FLAG` chooses which prediction is compared to `DV`:

```
IF(FLAG.EQ.2) THEN
  IPRED = RESP                          ; response
  Y     = IPRED + IPRED*EPS(3) + EPS(4)
ELSE
  IPRED = CP                            ; plasma drug
  Y     = IPRED + IPRED*EPS(1) + EPS(2)
ENDIF
```

Drug and response get **separate residual error models** (`$SIGMA` 1-2 and
3-4). They are different assays in different units; sharing an error model
between them would be meaningless.

### Example records (2-compartment oral, IDR I)

```
C,ID,TIME,AMT,RATE,DV,EVID,MDV,CMT,FLAG,ADDL,II,WT
,1,0,0,0,98.4,0,0,2,2,0,0,70       <- pre-dose response, near R0
,1,0,100,0,.,1,1,1,0,0,0,70        <- 100 mg oral dose into the depot
,1,1,0,0,1.62,0,0,2,1,0,0,70       <- plasma drug, 1.62 mg/L
,1,6,0,0,71.3,0,0,2,2,0,0,70       <- response, suppressed
,1,24,0,0,93.8,0,0,2,2,0,0,70      <- response, recovering
```

A pre-dose response observation is worth including in every subject: it is
what actually informs `R0`.

---

## Output tables

Three tables per run:

* `<model>.tab` - per-record predictions and residuals.
* `<model>.prof` - the profiles that make an IDR fit interpretable:

  | Column | Meaning |
  |--------|---------|
  | `CP` | plasma drug concentration |
  | `CE` | effect-site concentration (`effect` link only) |
  | `EFF` | fractional effect, 0 to 1 |
  | `FACT` | the factor actually applied - `1 - IMAX*EFF` or `1 + SMAX*EFF` |
  | `RESP` | the response `R` |

  `FACT` is the useful one to plot: it shows directly how far `KIN` or
  `KOUT` was pushed, and when.
* `<model>.par` - one row per subject with the individual parameters,
  the derived `KIN`, and the ETAs.

---

## Parameterisation and variability

* Fixed effects are log-normally distributed: `P = THETA(n)*EXP(ETA(n))`.
* `IMAX` and `F1` are bounded on (0.001, 0.999) with **logit-scale** IIV,
  the same treatment `DF` gets in the absorption libraries. This matters
  for `IMAX`: an exponential ETA on a parameter capped at 1 will push
  individual values above the cap.
* `SMAX` is not capped - stimulation has no natural upper bound.
* IIV is **estimated** on `CL`, volumes, `Q`s, `KA`, `KE0`, `R0`, `KOUT`
  and the potency term `IC50`/`SC50` (initial variance 0.09, ~30% CV).
* IIV on `IMAX`, `SMAX`, `HILL` and `F1` is coded but set to
  `$OMEGA 0 FIX`. Maximum effect and Hill slope are usually properties of
  the system rather than the subject; potency is where between-subject
  variability normally lives.
* `$OMEGA` is diagonal. `R0` and `KOUT` often want a `$OMEGA BLOCK` once
  the diagonal model converges.

### Initial estimates

| Parameter | Value |
|-----------|-------|
| `CL` / central volume | 5 L/h, 50 L |
| `Q` / `Q3`, periph. 1 | 10 L/h, 100 L |
| `Q4`, periph. 2 | 2 L/h, 200 L |
| `KA` / `F1` | 1 /h, 0.7 |
| `KE0` | 0.5 /h |
| `R0` / `KOUT` | 100 response units, 0.1 /h (t-half ~7 h) |
| `IMAX` / `IC50` | 0.8, 1 mg/L |
| `SMAX` / `SC50` | 2, 1 mg/L |
| `HILL` | 1 |

The PK values match the other libraries in this repository. The PD values
are placeholders on an arbitrary response scale - **`R0`, `KOUT` and the
potency term must be replaced with values on your own response scale**,
and the additive response error (`$SIGMA` 4, initialised at 1.0) rescaled
with them.

---

## Practical notes

1. **Fit the PK first.** Even though these files estimate PK and PD
   together, start from PK parameters that already converged on their own;
   a joint fit launched from generic initial estimates will usually fail.
2. **Fix `HILL` to 1 to begin with.** Free it only if the
   concentration-effect relationship is visibly steeper or shallower than
   a plain Emax curve.
3. **`IMAX` at the bound.** If `IMAX` runs to 0.999, the data support
   complete inhibition - fix it to 1 and re-run rather than leaving it
   pinned.
4. `KOUT` is identified by how fast the response returns to baseline, not
   by how deep it goes. Without washout sampling it will be confounded
   with the potency term.
5. Each file carries a commented SAEM + IMP alternative to FOCE-I, which
   is often a more robust route for a joint PK/PD fit.
6. `$COVARIANCE` is enabled everywhere; drop it for exploratory runs.

## Provenance

The model forms follow Dayneka, Garg & Jusko (1993), Jusko & Ko (1994)
and Sharma & Jusko (1998); the effect compartment follows Sheiner et al.
(1979). Full citations are in the root `../README.md`.
