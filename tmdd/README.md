# Target-mediated drug disposition (TMDD) model library

NONMEM control streams for target-mediated drug disposition. The base
structures are one-, two- and three-compartment plasma disposition
wrapped around a single target-binding module. This library expands
that into the approximation ladder and the input routes you actually
need in practice:

**3 structures x 5 model forms x 3 input routes = 45 control streams.**

---

## !! Units - read this before anything else

TMDD models bind drug to target through terms like `KON*C*R`. That product
is only meaningful if **drug and target concentrations are in the same
molar unit**. These files assume:

| Quantity | Unit |
|----------|------|
| Concentrations (drug, target, complex) | nM |
| `AMT` | nmol |
| Volumes | L |
| Time | **days** |

Two consequences worth stating plainly:

1. Feeding a mg/L drug assay and an nM target assay into the same TMDD
   model is the classic silent failure - it will run, converge to
   something, and be wrong. Convert the dose and the drug assay to molar
   units first (nmol = mg / MW * 1e6 for MW in g/mol).
2. **Time is in days here, not hours.** The `1cmt/`, `2cmt/` and `3cmt/`
   libraries in this repository use hours. Do not mix a data set between
   them without rescaling.

---

## The five model forms

Written in order of decreasing complexity. Each row is a strictly stronger
set of assumptions than the one above it.

### 1. `full` - full TMDD (Mager & Jusko 2001)

States: free drug, free target, complex.

```
dA_drug/dt    = In - CL*C - distribution - (KON*C*R - KOFF*RC)*V
dA_target/dt  = KDEG*R0*V - KDEG*R*V - (KON*C*R - KOFF*RC)*V
dA_complex/dt = (KON*C*R - KOFF*RC)*V - KINT*RC*V
```

Parameters: `CL`, volumes, `Q`s, `KON`, `KOFF`, `KINT`, `KDEG`, `R0`.
`KSYN` is not a separate parameter - it is fixed at `KDEG*R0` so the
target sits at its baseline `R0` before the first dose, which is also the
initial condition (`A_0`).

### 2. `qe` - quasi-equilibrium / rapid binding (Mager & Jusko 2001)

Binding is instantaneous, so `C*R/RC = KD = KOFF/KON`. The states become
**total** drug and **total** target, and free drug comes from the
quadratic:

```
BB = Ctot - Rtot - KD
C  = 0.5*(BB + sqrt(BB^2 + 4*KD*Ctot))
RC = Ctot - C
R  = Rtot - RC
```

Parameters: `CL`, volumes, `Q`s, `KD`, `KINT`, `KDEG`, `R0`.

### 3. `qss` - quasi-steady-state (Gibiansky et al. 2008)

Identical equations to `qe`, with `KD` replaced by
`KSS = (KOFF + KINT)/KON`. QSS is the more general of the two and
degenerates to QE when `KINT << KOFF`. If you are unsure which to use,
start here.

Parameters: `CL`, volumes, `Q`s, `KSS`, `KINT`, `KDEG`, `R0`.

### 4. `wagner` - quasi-equilibrium with constant total target

Rapid binding **plus** the assumption that total target stays at `R0`.
There is no target-turnover state at all, so `KDEG` and `KSYN` drop out.

```
BB = Ctot - R0 - KD
C  = 0.5*(BB + sqrt(BB^2 + 4*KD*Ctot))
dA_drug/dt = In - CL*C - distribution - KINT*RC*V
```

Parameters: `CL`, volumes, `Q`s, `KD`, `KINT`, `R0`. Equivalent to the
constant-Rtot form when `KDEG = KINT`.

### 5. `mm` - Michaelis-Menten approximation

Target binding collapses into saturable elimination from plasma:

```
dA_drug/dt = In - CL*C - distribution - VM*C/(KM + C)
```

with `VM ~ KINT*R0*V` and `KM ~ KSS`. Adequate when drug concentration
substantially exceeds target concentration (occupancy near 100%).

**The target is not modelled in this form.** `Rtot` is held at `R0` and
the target readouts are derived (`RC = R0*C/(KM+C)`), not fitted. If you
are fitting an `mm` model to drug data only, add
`IGNORE(FLAG.EQ.2)` to the `$DATA` record.

---

## Contents

| # | File | Cmts | Form | Input |
|---|------|------|------|-------|
| 1 | `tmdd_1cmt_full_iv.ctl` | 1 | Full TMDD | IV input (bolus or infusion) |
| 2 | `tmdd_1cmt_full_sc.ctl` | 1 | Full TMDD | Subcutaneous, first-order absorption |
| 3 | `tmdd_1cmt_full_sclag.ctl` | 1 | Full TMDD | Subcutaneous, first-order absorption with lag time |
| 4 | `tmdd_1cmt_qe_iv.ctl` | 1 | Quasi-equilibrium (rapid binding) approximation | IV input (bolus or infusion) |
| 5 | `tmdd_1cmt_qe_sc.ctl` | 1 | Quasi-equilibrium (rapid binding) approximation | Subcutaneous, first-order absorption |
| 6 | `tmdd_1cmt_qe_sclag.ctl` | 1 | Quasi-equilibrium (rapid binding) approximation | Subcutaneous, first-order absorption with lag time |
| 7 | `tmdd_1cmt_qss_iv.ctl` | 1 | Quasi-steady-state approximation | IV input (bolus or infusion) |
| 8 | `tmdd_1cmt_qss_sc.ctl` | 1 | Quasi-steady-state approximation | Subcutaneous, first-order absorption |
| 9 | `tmdd_1cmt_qss_sclag.ctl` | 1 | Quasi-steady-state approximation | Subcutaneous, first-order absorption with lag time |
| 10 | `tmdd_1cmt_wagner_iv.ctl` | 1 | Wagner (quasi-equilibrium, constant total target) | IV input (bolus or infusion) |
| 11 | `tmdd_1cmt_wagner_sc.ctl` | 1 | Wagner (quasi-equilibrium, constant total target) | Subcutaneous, first-order absorption |
| 12 | `tmdd_1cmt_wagner_sclag.ctl` | 1 | Wagner (quasi-equilibrium, constant total target) | Subcutaneous, first-order absorption with lag time |
| 13 | `tmdd_1cmt_mm_iv.ctl` | 1 | Michaelis-Menten approximation | IV input (bolus or infusion) |
| 14 | `tmdd_1cmt_mm_sc.ctl` | 1 | Michaelis-Menten approximation | Subcutaneous, first-order absorption |
| 15 | `tmdd_1cmt_mm_sclag.ctl` | 1 | Michaelis-Menten approximation | Subcutaneous, first-order absorption with lag time |
| 16 | `tmdd_2cmt_full_iv.ctl` | 2 | Full TMDD | IV input (bolus or infusion) |
| 17 | `tmdd_2cmt_full_sc.ctl` | 2 | Full TMDD | Subcutaneous, first-order absorption |
| 18 | `tmdd_2cmt_full_sclag.ctl` | 2 | Full TMDD | Subcutaneous, first-order absorption with lag time |
| 19 | `tmdd_2cmt_qe_iv.ctl` | 2 | Quasi-equilibrium (rapid binding) approximation | IV input (bolus or infusion) |
| 20 | `tmdd_2cmt_qe_sc.ctl` | 2 | Quasi-equilibrium (rapid binding) approximation | Subcutaneous, first-order absorption |
| 21 | `tmdd_2cmt_qe_sclag.ctl` | 2 | Quasi-equilibrium (rapid binding) approximation | Subcutaneous, first-order absorption with lag time |
| 22 | `tmdd_2cmt_qss_iv.ctl` | 2 | Quasi-steady-state approximation | IV input (bolus or infusion) |
| 23 | `tmdd_2cmt_qss_sc.ctl` | 2 | Quasi-steady-state approximation | Subcutaneous, first-order absorption |
| 24 | `tmdd_2cmt_qss_sclag.ctl` | 2 | Quasi-steady-state approximation | Subcutaneous, first-order absorption with lag time |
| 25 | `tmdd_2cmt_wagner_iv.ctl` | 2 | Wagner (quasi-equilibrium, constant total target) | IV input (bolus or infusion) |
| 26 | `tmdd_2cmt_wagner_sc.ctl` | 2 | Wagner (quasi-equilibrium, constant total target) | Subcutaneous, first-order absorption |
| 27 | `tmdd_2cmt_wagner_sclag.ctl` | 2 | Wagner (quasi-equilibrium, constant total target) | Subcutaneous, first-order absorption with lag time |
| 28 | `tmdd_2cmt_mm_iv.ctl` | 2 | Michaelis-Menten approximation | IV input (bolus or infusion) |
| 29 | `tmdd_2cmt_mm_sc.ctl` | 2 | Michaelis-Menten approximation | Subcutaneous, first-order absorption |
| 30 | `tmdd_2cmt_mm_sclag.ctl` | 2 | Michaelis-Menten approximation | Subcutaneous, first-order absorption with lag time |
| 31 | `tmdd_3cmt_full_iv.ctl` | 3 | Full TMDD | IV input (bolus or infusion) |
| 32 | `tmdd_3cmt_full_sc.ctl` | 3 | Full TMDD | Subcutaneous, first-order absorption |
| 33 | `tmdd_3cmt_full_sclag.ctl` | 3 | Full TMDD | Subcutaneous, first-order absorption with lag time |
| 34 | `tmdd_3cmt_qe_iv.ctl` | 3 | Quasi-equilibrium (rapid binding) approximation | IV input (bolus or infusion) |
| 35 | `tmdd_3cmt_qe_sc.ctl` | 3 | Quasi-equilibrium (rapid binding) approximation | Subcutaneous, first-order absorption |
| 36 | `tmdd_3cmt_qe_sclag.ctl` | 3 | Quasi-equilibrium (rapid binding) approximation | Subcutaneous, first-order absorption with lag time |
| 37 | `tmdd_3cmt_qss_iv.ctl` | 3 | Quasi-steady-state approximation | IV input (bolus or infusion) |
| 38 | `tmdd_3cmt_qss_sc.ctl` | 3 | Quasi-steady-state approximation | Subcutaneous, first-order absorption |
| 39 | `tmdd_3cmt_qss_sclag.ctl` | 3 | Quasi-steady-state approximation | Subcutaneous, first-order absorption with lag time |
| 40 | `tmdd_3cmt_wagner_iv.ctl` | 3 | Wagner (quasi-equilibrium, constant total target) | IV input (bolus or infusion) |
| 41 | `tmdd_3cmt_wagner_sc.ctl` | 3 | Wagner (quasi-equilibrium, constant total target) | Subcutaneous, first-order absorption |
| 42 | `tmdd_3cmt_wagner_sclag.ctl` | 3 | Wagner (quasi-equilibrium, constant total target) | Subcutaneous, first-order absorption with lag time |
| 43 | `tmdd_3cmt_mm_iv.ctl` | 3 | Michaelis-Menten approximation | IV input (bolus or infusion) |
| 44 | `tmdd_3cmt_mm_sc.ctl` | 3 | Michaelis-Menten approximation | Subcutaneous, first-order absorption |
| 45 | `tmdd_3cmt_mm_sclag.ctl` | 3 | Michaelis-Menten approximation | Subcutaneous, first-order absorption with lag time |

Plus this `README.md`.

### Naming convention

```
tmdd_<n>cmt_<form>_<route>.ctl
```

* `<form>`: `full`, `qe`, `qss`, `wagner`, `mm`
* `<route>`: `iv` (bolus or infusion), `sc` (first-order depot),
  `sclag` (first-order depot with lag time)

---

## Volume and Q naming

Taken straight from the page-4 diagrams:

| Structure | Central | Periph. 1 | Periph. 2 |
|-----------|---------|-----------|-----------|
| 1-compartment | `V` | - | - |
| 2-compartment | `V2` | `Q` / `V3` | - |
| 3-compartment | `V2` | `Q3` / `V3` | `Q4` / `V4` |

Note the same `V2`/`V3` shift documented in `../2cmt/` and `../3cmt/`:
`V2` is the central volume in the 2- and 3-compartment TMDD models, but
the one-compartment model simply calls its central volume `V`.

**Only free drug distributes to the peripheral compartments.** In the
`qe`, `qss` and `wagner` forms the central state is *total* drug, so the
distribution term is built from the free concentration recovered from the
quadratic, not from the state itself.

---

## Data set

Comma-delimited, at `../data/tmdddata.csv`, commented header row
(`IGNORE=@`), columns:

```
ID TIME AMT RATE DV EVID MDV CMT FLAG ADDL II WT
```

`FLAG` is the addition relative to the other libraries in this repository.

| Column | Meaning |
|--------|---------|
| `ID` | subject identifier |
| `TIME` | time of the event (days) |
| `AMT` | dose amount (**nmol**) |
| `RATE` | 0 = bolus / instantaneous; > 0 = IV infusion rate |
| `DV` | observed concentration (nM) - which analyte is set by `FLAG` |
| `EVID` | 0 = observation, 1 = dose |
| `MDV` | 0 = observation, 1 = missing dependent variable |
| `CMT` | see below |
| `FLAG` | 1 = total drug, 2 = total target |
| `ADDL`, `II` | additional-dose bookkeeping |
| `WT` | body weight (kg) - carried for covariate work, not yet used |

### `CMT` on dose records

| Route | Dose `CMT` | `RATE` |
|-------|-----------|--------|
| `iv` | the plasma compartment (1) | 0 for a bolus, > 0 for an infusion |
| `sc`, `sclag` | 1 (SC depot) | 0 |

In the SC models every compartment index shifts by one, because the depot
takes slot 1. The compartment map is written at the top of each file.

### `CMT` on observation records

**All observations use `CMT` = the plasma compartment**, for both
analytes. `FLAG` chooses which prediction is compared to `DV`:

```
IF(FLAG.EQ.2) THEN
  IPRED = RTOT                          ; total target
  Y     = IPRED + IPRED*EPS(3) + EPS(4)
ELSE
  IPRED = CTOT                          ; total drug
  Y     = IPRED + IPRED*EPS(1) + EPS(2)
ENDIF
```

Drug and target therefore get **separate residual error models** -
`$SIGMA` 1-2 for drug, 3-4 for target. This is deliberate: the two assays
almost never have comparable precision.

### Example records (2-compartment QSS, IV)

```
C,ID,TIME,AMT,RATE,DV,EVID,MDV,CMT,FLAG,ADDL,II,WT
,1,0,700,0,.,1,1,1,0,0,0,70        <- 700 nmol IV bolus into plasma
,1,0.25,0,0,210,0,0,1,1,0,0,70     <- total drug, 210 nM
,1,0.25,0,0,1.9,0,0,1,2,0,0,70     <- total target, 1.9 nM
,1,7,0,0,44,0,0,1,1,0,0,70         <- total drug
```

---

## Output tables

Every `$TABLE` writes a `.tab` file. Each run produces three:

* `<model>.tab` - per-record predictions and residuals.
* `<model>.prof` - the analyte profiles, which is what makes a TMDD fit
  interpretable:

  | Column | Meaning |
  |--------|---------|
  | `CFR` | free drug concentration |
  | `CTOT` | total drug (free + complex) |
  | `RFR` | free target concentration |
  | `RCX` | drug-target complex concentration |
  | `RTOT` | total target (free + complex) |

  Every form produces all five, so profiles are directly comparable when
  you run the same data through `full`, `qss` and `mm`.
* `<model>.par` - one row per subject with individual parameters and ETAs.

---

## Parameterisation and variability

* Fixed effects are log-normally distributed: `P = THETA(n)*EXP(ETA(n))`.
* SC bioavailability `F1` is bounded on (0.001, 0.999) with logit-scale
  IIV, the same treatment the `DF` dose fraction gets elsewhere in this
  repository.
* IIV is **estimated** on `CL`, volumes, `Q`s, `KINT`, `KDEG`, `R0`, `VM`
  and `KA` (initial variance 0.09, ~30% CV).
* IIV on the binding constants (`KON`, `KOFF`, `KD`, `KSS`, `KM`) and on
  `F1` / `ALAG1` is coded but set to `$OMEGA 0 FIX`. Binding constants are
  properties of the molecule, not of the subject - free them only with a
  specific reason.
* `$OMEGA` is diagonal. `CL` and central volume usually want a
  `$OMEGA BLOCK` once the diagonal model converges.

### Initial estimates

Generic monoclonal-antibody-like placeholders:

| Parameter | Value |
|-----------|-------|
| `CL` | 0.2 L/day |
| central volume | 3 L |
| `Q` / `Q3`, periph. 1 | 0.5 L/day, 3 L |
| `Q4`, periph. 2 | 0.2 L/day, 5 L |
| `KON` | 1.0 1/(nM*day) |
| `KOFF` | 0.5 1/day |
| `KD` | 0.5 nM |
| `KSS` | 1.0 nM |
| `KINT` | 0.5 1/day |
| `KDEG` | 1.0 1/day |
| `R0` | 5 nM |
| `VM` / `KM` | 7.5 nmol/day, 1.0 nM |
| `KA` / `F1` / `ALAG1` | 0.25 1/day, 0.7, 0.1 day |

**Replace these with values for your compound before running.** TMDD
models are far more sensitive to initial estimates than the linear models
in the rest of this repository.

---

## Estimation

Every file ships with FOCE-I:

```
$ESTIMATION METHOD=1 INTER MAXEVAL=9999 NSIG=3 SIGL=9 PRINT=5 NOABORT POSTHOC
```

and carries a commented two-step SAEM + IMP alternative, which is often a
more robust route to the optimum for these models:

```
$ESTIMATION METHOD=SAEM INTERACTION NBURN=2000 NITER=1000 PRINT=50
$ESTIMATION METHOD=IMP INTERACTION EONLY=1 NITER=10 ISAMPLE=3000 PRINT=1
```

`TOL=9` is set on `$SUBROUTINES`. Full TMDD models are stiff around the
binding equilibrium; if a run crawls or fails, relax to `TOL=6` before
changing anything structural.

---

## Choosing a form

A practical order of attack:

1. **Start with `qss`.** It is the most generally applicable
   approximation and by far the most likely to converge.
2. Move **down** to `wagner` or `mm` if the target-turnover parameters
   (`KDEG`, and `R0` separately from `KINT`) are not identifiable.
3. Move **up** to `full` only when you have target and/or complex
   measurements over a range of doses. Without them `KON` and `KOFF` are
   individually unidentifiable - only their ratio is informed.
4. Compare `mm` against `qss` on the same data. If the fits are
   indistinguishable, the target is saturated throughout and `mm` is the
   honest choice.

## Caveats

1. **Identifiability is the whole game.** The full model has seven
   disposition/target parameters before covariates; drug-concentration
   data alone will not support it.
2. `KSYN` is never estimated - it is `KDEG*R0` by construction, which is
   what makes the pre-dose baseline consistent. Do not add it as a THETA
   without also removing that constraint.
3. The target is assumed to reside in the central volume. Models where the
   target sits in tissue need a different structure than anything here.
4. Negative-concentration guards (`IF(CTOT.LT.0) CTOT = 0`) are present
   before every square root. They protect against ODE solver noise near
   zero, not against a mis-specified model.
5. `$COVARIANCE` is enabled everywhere; drop it for exploratory runs.

## Provenance

Approximation forms follow Mager & Jusko (2001), Mager & Krzyzanski
(2005), Gibiansky et al. (2008) and the constant-total-target treatment
after Wagner (1973). Full citations are in the root `../README.md`.
