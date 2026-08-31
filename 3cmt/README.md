# Three-compartment NONMEM model library

NONMEM control streams (`.ctl`) for the **39 three-compartment models** on
the standard three-compartment PK model set. The absorption taxonomy is
identical to `../1cmt/` and `../2cmt/`; every model here has two peripheral
compartments hanging off central.

Each file is a complete, runnable template: structural model,
between-subject variability, combined residual error, estimation, covariance
step, output tables and a commented-out simulation block.

---

## Contents

| # | File | Absorption | Elimination | Solver |
|---|------|------------|-------------|--------|
| 1 | `3cmt_fo_nolag_lin.ctl` | First-order absorption without lag time | Linear elimination | ADVAN12 TRANS4 |
| 2 | `3cmt_fo_nolag_mm.ctl` | First-order absorption without lag time | Nonlinear elimination | ADVAN13 `$DES` |
| 3 | `3cmt_fo_nolag_mmlin.ctl` | First-order absorption without lag time | Nonlinear and linear elimination | ADVAN13 `$DES` |
| 4 | `3cmt_fo_lag_lin.ctl` | First-order absorption with lag time | Linear elimination | ADVAN12 TRANS4 |
| 5 | `3cmt_fo_lag_mm.ctl` | First-order absorption with lag time | Nonlinear elimination | ADVAN13 `$DES` |
| 6 | `3cmt_fo_lag_mmlin.ctl` | First-order absorption with lag time | Nonlinear and linear elimination | ADVAN13 `$DES` |
| 7 | `3cmt_sigmoid_dur_lin.ctl` | Sigmoid absorption with duration | Linear elimination | ADVAN12 TRANS4 |
| 8 | `3cmt_sigmoid_dur_mm.ctl` | Sigmoid absorption with duration | Nonlinear elimination | ADVAN13 `$DES` |
| 9 | `3cmt_sigmoid_dur_mmlin.ctl` | Sigmoid absorption with duration | Nonlinear and linear elimination | ADVAN13 `$DES` |
| 10 | `3cmt_sigmoid_rate_lin.ctl` | Sigmoid absorption with rate | Linear elimination | ADVAN12 TRANS4 |
| 11 | `3cmt_sigmoid_rate_mm.ctl` | Sigmoid absorption with rate | Nonlinear elimination | ADVAN13 `$DES` |
| 12 | `3cmt_sigmoid_rate_mmlin.ctl` | Sigmoid absorption with rate | Nonlinear and linear elimination | ADVAN13 `$DES` |
| 13 | `3cmt_zo_dur_lin.ctl` | Zero-order absorption with duration | Linear elimination | ADVAN11 TRANS4 |
| 14 | `3cmt_zo_dur_mm.ctl` | Zero-order absorption with duration | Nonlinear elimination | ADVAN13 `$DES` |
| 15 | `3cmt_zo_dur_mmlin.ctl` | Zero-order absorption with duration | Nonlinear and linear elimination | ADVAN13 `$DES` |
| 16 | `3cmt_zo_rate_lin.ctl` | Zero-order absorption with rate | Linear elimination | ADVAN11 TRANS4 |
| 17 | `3cmt_zo_rate_mm.ctl` | Zero-order absorption with rate | Nonlinear elimination | ADVAN13 `$DES` |
| 18 | `3cmt_zo_rate_mmlin.ctl` | Zero-order absorption with rate | Nonlinear and linear elimination | ADVAN13 `$DES` |
| 19 | `3cmt_transit_lin.ctl` | Transit compartment absorption | Linear elimination | ADVAN13 `$DES` |
| 20 | `3cmt_transit_mm.ctl` | Transit compartment absorption | Nonlinear elimination | ADVAN13 `$DES` |
| 21 | `3cmt_transit_mmlin.ctl` | Transit compartment absorption | Nonlinear and linear elimination | ADVAN13 `$DES` |
| 22 | `3cmt_parallel_fo_nolag_lin.ctl` | Parallel first-order absorption without lag time | Linear elimination | ADVAN13 `$DES` |
| 23 | `3cmt_parallel_fo_nolag_mm.ctl` | Parallel first-order absorption without lag time | Nonlinear elimination | ADVAN13 `$DES` |
| 24 | `3cmt_parallel_fo_nolag_mmlin.ctl` | Parallel first-order absorption without lag time | Nonlinear and linear elimination | ADVAN13 `$DES` |
| 25 | `3cmt_parallel_fo_lag_lin.ctl` | Parallel first-order absorption with lag time | Linear elimination | ADVAN13 `$DES` |
| 26 | `3cmt_parallel_fo_lag_mm.ctl` | Parallel first-order absorption with lag time | Nonlinear elimination | ADVAN13 `$DES` |
| 27 | `3cmt_parallel_fo_lag_mmlin.ctl` | Parallel first-order absorption with lag time | Nonlinear and linear elimination | ADVAN13 `$DES` |
| 28 | `3cmt_seq_fo_zo_dur_lin.ctl` | Sequential absorption with first-order process with lag time followed by zero-order with duration | Linear elimination | ADVAN12 TRANS4 |
| 29 | `3cmt_seq_fo_zo_dur_mm.ctl` | Sequential absorption with first-order process with lag time followed by zero-order with duration | Nonlinear elimination | ADVAN13 `$DES` |
| 30 | `3cmt_seq_fo_zo_dur_mmlin.ctl` | Sequential absorption with first-order process with lag time followed by zero-order with duration | Nonlinear and linear elimination | ADVAN13 `$DES` |
| 31 | `3cmt_seq_fo_zo_rate_lin.ctl` | Sequential absorption with first-order process with lag time followed by zero-order with rate | Linear elimination | ADVAN12 TRANS4 |
| 32 | `3cmt_seq_fo_zo_rate_mm.ctl` | Sequential absorption with first-order process with lag time followed by zero-order with rate | Nonlinear elimination | ADVAN13 `$DES` |
| 33 | `3cmt_seq_fo_zo_rate_mmlin.ctl` | Sequential absorption with first-order process with lag time followed by zero-order with rate | Nonlinear and linear elimination | ADVAN13 `$DES` |
| 34 | `3cmt_parallel_fo_zo_dur_lin.ctl` | Parallel absorption with first-order process and zero-order with duration | Linear elimination | ADVAN12 TRANS4 |
| 35 | `3cmt_parallel_fo_zo_dur_mm.ctl` | Parallel absorption with first-order process and zero-order with duration | Nonlinear elimination | ADVAN13 `$DES` |
| 36 | `3cmt_parallel_fo_zo_dur_mmlin.ctl` | Parallel absorption with first-order process and zero-order with duration | Nonlinear and linear elimination | ADVAN13 `$DES` |
| 37 | `3cmt_parallel_fo_zo_rate_lin.ctl` | Parallel absorption with first-order process and zero-order with rate | Linear elimination | ADVAN12 TRANS4 |
| 38 | `3cmt_parallel_fo_zo_rate_mm.ctl` | Parallel absorption with first-order process and zero-order with rate | Nonlinear elimination | ADVAN13 `$DES` |
| 39 | `3cmt_parallel_fo_zo_rate_mmlin.ctl` | Parallel absorption with first-order process and zero-order with rate | Nonlinear and linear elimination | ADVAN13 `$DES` |

Plus this `README.md`.

---

## Naming convention

```
3cmt_<absorption>_<elimination>.ctl
```

Model numbers 1-39 are listed in the Contents table above; they are
not part of the file name and do not appear inside the files.

* `<absorption>`:
  * `fo_nolag` - First-order absorption without lag time
  * `fo_lag` - First-order absorption with lag time
  * `sigmoid_dur` - Sigmoid absorption with duration
  * `sigmoid_rate` - Sigmoid absorption with rate
  * `zo_dur` - Zero-order absorption with duration
  * `zo_rate` - Zero-order absorption with rate
  * `transit` - Transit compartment absorption
  * `parallel_fo_nolag` - Parallel first-order absorption without lag time
  * `parallel_fo_lag` - Parallel first-order absorption with lag time
  * `seq_fo_zo_dur` - Sequential absorption with first-order process with lag time followed by zero-order with duration
  * `seq_fo_zo_rate` - Sequential absorption with first-order process with lag time followed by zero-order with rate
  * `parallel_fo_zo_dur` - Parallel absorption with first-order process and zero-order with duration
  * `parallel_fo_zo_rate` - Parallel absorption with first-order process and zero-order with rate
* `<elimination>`:
  * `lin` - linear elimination (CL)
  * `mm` - nonlinear Michaelis-Menten elimination (Vmax/Km)
  * `mmlin` - parallel nonlinear and linear elimination

---

## Volume and Q naming - read this before interpreting estimates

Names follow the NONMEM `TRANS4` convention, which is what the diagrams use,
so the **same symbol means different things in the two layouts**:

| Layout | Models | Central | Periph. 1 | Periph. 2 |
|--------|--------|---------|-----------|-----------|
| With a depot | 1-12, 19-39 | `V2` | `Q3` / `V3` | `Q4` / `V4` |
| No depot (dose into central) | 13-18 | `V1` | `Q2` / `V2` | `Q3` / `V3` |

This is the same shift that separates `ADVAN12` from `ADVAN11`, and it is
the single most common source of mix-ups when comparing runs across the
library. `V2` is the central volume in model 1 and a *peripheral* volume in
model 13.

### Distribution

Written in `$DES` with named intermediates so every line stays inside 80
columns:

```
CONC  = A(central)/Vc          ; central concentration
CP1   = A(periph1)/Vp1
CP2   = A(periph2)/Vp2
DIST1 = Qa*(CONC - CP1)        ; net flux central -> peripheral 1
DIST2 = Qb*(CONC - CP2)        ; net flux central -> peripheral 2
ELR   = <elimination>

dA_central/dt = input - DIST1 - DIST2 - ELR
dA_periph1/dt = DIST1
dA_periph2/dt = DIST2
```

### Elimination

| Variant | Elimination term `ELR` | Parameters |
|---------|------------------------|------------|
| `lin` | `CL*CONC` | `CL` |
| `mm` | `VM*CONC/(KM+CONC)` | `VM`, `KM` |
| `mmlin` | `CL*CONC + VM*CONC/(KM+CONC)` | `CL`, `VM`, `KM` |

`VM` is the NONMEM variable name used for Vmax.

---

## Absorption sub-models

**First-order absorption without lag time**  (`fo_nolag`)

```
Dose --(instantaneous)--> [Depot] --Ka--> [Central (V2)]
[Periph.1 (V3)] <--Q3--> [Central (V2)] <--Q4--> [Periph.2 (V4)]
                              |  elimination
```

* One dose record per administration: CMT=1 (depot), RATE=0.

**First-order absorption with lag time**  (`fo_lag`)

```
Dose at t=ALAG1 --> [Depot] --Ka--> [Central (V2)]
[Periph.1 (V3)] <--Q3--> [Central (V2)] <--Q4--> [Periph.2 (V4)]
                              |  elimination
```

* One dose record per administration: CMT=1 (depot), RATE=0.
* Absorption starts at TIME+ALAG1.

**Sigmoid absorption with duration**  (`sigmoid_dur`)

```
Dose --(zero-order, D1)--> [Depot] --Ka--> [Central (V2)]
[Periph.1 (V3)] <--Q3--> [Central (V2)] <--Q4--> [Periph.2 (V4)]
                              |  elimination
```

* One dose record per administration: CMT=1 (depot), RATE=-2
* (duration D1 modelled in $PK).  Zero-order input INTO the depot
* followed by first-order transfer to central gives the sigmoid
* (S-shaped) plasma profile.

**Sigmoid absorption with rate**  (`sigmoid_rate`)

```
Dose --(zero-order, R1)--> [Depot] --Ka--> [Central (V2)]
[Periph.1 (V3)] <--Q3--> [Central (V2)] <--Q4--> [Periph.2 (V4)]
                              |  elimination
```

* One dose record per administration: CMT=1 (depot), RATE=-1
* (rate R1 modelled in $PK).

**Zero-order absorption with duration**  (`zo_dur`)

```
Dose --(zero-order, D1)--> [Central (V1)]
[Periph.1 (V2)] <--Q2--> [Central (V1)] <--Q3--> [Periph.2 (V3)]
                              |  elimination
```

* One dose record per administration: CMT=1 (central), RATE=-2
* (duration D1 modelled in $PK).

**Zero-order absorption with rate**  (`zo_rate`)

```
Dose --(zero-order, R1)--> [Central (V1)]
[Periph.1 (V2)] <--Q2--> [Central (V1)] <--Q3--> [Periph.2 (V3)]
                              |  elimination
```

* One dose record per administration: CMT=1 (central), RATE=-1
* (rate R1 modelled in $PK).

**Transit compartment absorption**  (`transit`)

```
Dose --> [T1] --Ktr--> ... --Ktr--> [Tn] --Ktr--> [Depot] --Ka-->
     [Central (V2)]     N transit cmts,  Ktr = (N+1)/MTT
[Periph.1 (V3)] <--Q3--> [Central (V2)] <--Q4--> [Periph.2 (V4)]
                              |  elimination
```

* One dose record per administration: CMT=1 (depot), RATE=0.
* NOTE: the analytical transit solution tracks only the MOST RECENT
* dose.  It is exact for single-dose data; for multiple dosing either
* code the transit chain explicitly with a fixed N, or use the
* superposition form.

**Parallel first-order absorption without lag time**  (`parallel_fo_nolag`)

```
DF x Dose     --> [Depot1] --Ka1--\
                                   >--> [Central (V2)]
(1-DF) x Dose --> [Depot2] --Ka2--/
[Periph.1 (V3)] <--Q3--> [Central (V2)] <--Q4--> [Periph.2 (V4)]
                              |  elimination
```

* TWO dose records per administration, both with the full AMT:
  CMT=1 (depot 1), RATE=0   -> receives DF     x Dose via F1
  CMT=2 (depot 2), RATE=0   -> receives (1-DF) x Dose via F2

**Parallel first-order absorption with lag time**  (`parallel_fo_lag`)

```
DF x Dose     at t=ALAG1 --> [Depot1] --Ka1--\
                                              >--> [Central (V2)]
(1-DF) x Dose at t=ALAG2 --> [Depot2] --Ka2--/
[Periph.1 (V3)] <--Q3--> [Central (V2)] <--Q4--> [Periph.2 (V4)]
                              |  elimination
```

* TWO dose records per administration, both with the full AMT:
  CMT=1 (depot 1), RATE=0   -> DF     x Dose, delayed by ALAG1
  CMT=2 (depot 2), RATE=0   -> (1-DF) x Dose, delayed by ALAG2

**Sequential absorption with first-order process with lag time followed by zero-order with duration**  (`seq_fo_zo_dur`)

```
DF x Dose     --(zero-order, D2)--> [Central (V2)]
(1-DF) x Dose at t=D2 --> [Depot] --Ka--^
[Periph.1 (V3)] <--Q3--> [Central (V2)] <--Q4--> [Periph.2 (V4)]
                              |  elimination
```

* TWO dose records per administration, both with the full AMT:
  CMT=2 (central), RATE=-2  -> DF     x Dose, zero-order over D2
  CMT=1 (depot),   RATE=0   -> (1-DF) x Dose, released at t=D2
  (ALAG1 = D2) then absorbed with Ka

**Sequential absorption with first-order process with lag time followed by zero-order with rate**  (`seq_fo_zo_rate`)

```
DF x Dose     --(zero-order, R2)--> [Central (V2)]
(1-DF) x Dose at t=DF x Dose/R2 --> [Depot] --Ka--^
[Periph.1 (V3)] <--Q3--> [Central (V2)] <--Q4--> [Periph.2 (V4)]
                              |  elimination
```

* TWO dose records per administration, both with the full AMT:
  CMT=2 (central), RATE=-1  -> DF     x Dose, zero-order at rate R2
  CMT=1 (depot),   RATE=0   -> (1-DF) x Dose, released at
  t = DF x Dose / R2, absorbed with Ka

**Parallel absorption with first-order process and zero-order with duration**  (`parallel_fo_zo_dur`)

```
DF x Dose     --(zero-order, D2)--> [Central (V2)]
(1-DF) x Dose --> [Depot] --Ka--^
[Periph.1 (V3)] <--Q3--> [Central (V2)] <--Q4--> [Periph.2 (V4)]
                              |  elimination
```

* TWO dose records per administration, both with the full AMT:
  CMT=2 (central), RATE=-2  -> DF     x Dose, zero-order over D2
  CMT=1 (depot),   RATE=0   -> (1-DF) x Dose, absorbed with Ka
  (starts immediately - parallel)

**Parallel absorption with first-order process and zero-order with rate**  (`parallel_fo_zo_rate`)

```
DF x Dose     --(zero-order, R2)--> [Central (V2)]
(1-DF) x Dose --> [Depot] --Ka--^
[Periph.1 (V3)] <--Q3--> [Central (V2)] <--Q4--> [Periph.2 (V4)]
                              |  elimination
```

* TWO dose records per administration, both with the full AMT:
  CMT=2 (central), RATE=-1  -> DF     x Dose, zero-order at rate R2
  CMT=1 (depot),   RATE=0   -> (1-DF) x Dose, absorbed with Ka

---

## Solver choice

Closed-form (analytical) ADVANs are used wherever the model permits:

* `ADVAN11 TRANS4` (`CL`, `V1`, `Q2`, `V2`, `Q3`, `V3`) - three compartments,
  no depot (models 13, 16).
* `ADVAN12 TRANS4` (`CL`, `V2`, `Q3`, `V3`, `Q4`, `V4`, `KA`) - three
  compartments with a first-order depot
  (models 1, 4, 7, 10, 28, 31, 34, 37).

`ADVAN13` with an explicit `$DES` block is used whenever the closed form does
not apply:

* any model with **nonlinear (Michaelis-Menten) elimination**;
* **parallel first-order absorption** (two depot compartments);
* **transit compartment absorption** (gamma-density input term).

`TOL=9` is set on `$SUBROUTINES` for all ODE models. Three-compartment ODE
models with Michaelis-Menten elimination are the stiffest in this library and
also the slowest; if a run stalls or crawls, relax to `TOL=6` with
`$ESTIMATION ... NSIG=2 SIGL=6`, and consider `METHOD=SAEM` for the
`mm` / `mmlin` variants.

---

## Data set

All control streams expect a comma-delimited data set at
`../data/pkdata.csv` with a commented header row (`IGNORE=@`) and columns:

```
ID TIME AMT RATE DV EVID MDV CMT ADDL II WT
```

| Column | Meaning |
|--------|---------|
| `ID` | subject identifier |
| `TIME` | time of the event (h) |
| `AMT` | dose amount (mg); 0 or `.` on observation records |
| `RATE` | 0 = instantaneous; `-1` = rate modelled in `$PK`; `-2` = duration modelled in `$PK` |
| `DV` | observed concentration (mg/L) |
| `EVID` | 0 = observation, 1 = dose |
| `MDV` | 0 = observation, 1 = missing dependent variable |
| `CMT` | compartment the record refers to |
| `ADDL`, `II` | additional-dose bookkeeping |
| `WT` | body weight (kg) - carried for covariate work, not yet used |

One data set drives all three libraries: the `CMT` numbering on **dose**
records never changes, because peripheral compartments are always added
after the dosing compartments.

### `RATE` conventions

| Absorption model | `RATE` on the dose record |
|------------------|---------------------------|
| First-order, transit, parallel first-order | `0` |
| Sigmoid / zero-order **with duration** | `-2` (`D1` or `D2` estimated in `$PK`) |
| Sigmoid / zero-order **with rate** | `-1` (`R1` or `R2` estimated in `$PK`) |

### Models that need two dose records per administration

Models 22-39 split the dose between two input routes using bioavailability
fractions, so the data set must carry **two dose records at the same `TIME`,
each with the full `AMT`**, one per compartment; NONMEM applies `F1` and `F2`
to split them. Example for model 34 (parallel first-order + zero-order with
duration):

```
C,ID,TIME,AMT,RATE,DV,EVID,MDV,CMT,ADDL,II,WT
,1,0,100,0,.,1,1,1,0,0,70      <- (1-DF) x 100 mg into the depot
,1,0,100,-2,.,1,1,2,0,0,70     <- DF x 100 mg zero-order into central
,1,0.5,0,0,1.24,0,0,2,0,0,70   <- observation
```

---

## Parameterisation and variability

* All fixed effects are log-normally distributed:
  `P = THETA(n) * EXP(ETA(n))`.
* The dose fraction `DF` is bounded on (0.001, 0.999) with **logit-scale**
  IIV:

  ```
  LGTDF = LOG(THETA(n)/(1-THETA(n))) + ETA(n)
  DF    = EXP(LGTDF)/(1+EXP(LGTDF))
  ```

* IIV is **estimated** on `CL`, all three volumes, both `Q`s, `VM`, `KA`,
  `KA1`, `KA2` and `MTT` (initial variance 0.09, ~30% CV).
* IIV on `KM`, `D1`, `D2`, `R1`, `R2`, `ALAG1`, `ALAG2` and `DF` is coded but
  set to `$OMEGA 0 FIX`. Free them one at a time only if the data support it.
* A three-compartment model with IIV on all seven disposition parameters is
  over-parameterised for most data sets. Expect to fix several `$OMEGA`
  elements to zero - the second peripheral compartment's `Q` and `V` are
  usually the first to go.
* `$OMEGA` is diagonal; add a `$OMEGA BLOCK` on `CL`/central volume once the
  diagonal model converges.
* Residual error is combined proportional + additive:
  `Y = IPRED + IPRED*EPS(1) + EPS(2)`.

Initial estimates are generic placeholders (`CL` 5 L/h, central volume 50 L,
first `Q` 10 L/h and volume 100 L, second `Q` 2 L/h and volume 200 L, `KA`
1 /h, `Vmax` 100 mg/h, `Km` 5 mg/L). The second peripheral compartment is
deliberately initialised as the slower, larger one. **Replace all of these
with values appropriate to your compound and units before running.**

---

## Running a model

```bash
# plain NONMEM
nmfe75 3cmt_fo_nolag_lin.ctl 3cmt_fo_nolag_lin.lst

# PsN
execute 3cmt_fo_nolag_lin.ctl -dir=3cmt_fo_nolag_lin
```

Every `$TABLE` writes a `.tab` file. Each run produces two:

* `<model>.tab` - per-record predictions and residuals.
* `<model>.par` - one row per subject with individual parameters and ETAs.

## Simulating instead of estimating

Every file ends with a commented simulation block. Comment out
`$ESTIMATION` and `$COVARIANCE`, then un-comment it.

---

## Caveats

1. **Do not fit a three-compartment model without the data to support it.**
   Three phases need rich early sampling and a long terminal tail; otherwise
   the second peripheral compartment is unidentifiable and the run will
   either fail the covariance step or converge to a boundary.
2. **Transit models (19-21)** use the Savic analytical approximation, which
   tracks only the most recent dose. Exact for single-dose data; for
   multiple dosing, code the transit chain explicitly with a fixed `N`.
3. **Sequential models (28-33)** set `ALAG1` from the duration or from
   `DF*Dose/R2`, so the first-order process begins exactly when the
   zero-order input ends.
4. Nonlinear-elimination models are sensitive to initial estimates and to the
   dose range; without data spanning the saturable region, `VM` and `KM` will
   not be separately identifiable.
5. `$COVARIANCE` is enabled everywhere; drop it for exploratory runs.
6. Units are assumed consistent (mg, L, h). Check the `S1`/`S2`/`S3` scaling
   if your `DV` is in ng/mL while `AMT` is in mg.

## Provenance

The 39 structures are the standard three-compartment
absorption/elimination grid used throughout the pharmacometrics
literature; the root `../README.md` lists the background reading
behind each model family.
