# Two-compartment NONMEM model library

NONMEM control streams (`.ctl`) for the **39 two-compartment models** on
the standard two-compartment PK model set. The absorption taxonomy is
identical to the one-compartment set in `../1cmt/`; every model here adds a
peripheral compartment connected to central by inter-compartmental clearance
`Q`.

Each file is a complete, runnable template: structural model,
between-subject variability, combined residual error, estimation, covariance
step, output tables and a commented-out simulation block.

---

## Contents

| # | File | Absorption | Elimination | Solver |
|---|------|------------|-------------|--------|
| 1 | `2cmt_fo_nolag_lin.ctl` | First-order absorption without lag time | Linear elimination | ADVAN4 TRANS4 |
| 2 | `2cmt_fo_nolag_mm.ctl` | First-order absorption without lag time | Nonlinear elimination | ADVAN13 `$DES` |
| 3 | `2cmt_fo_nolag_mmlin.ctl` | First-order absorption without lag time | Nonlinear and linear elimination | ADVAN13 `$DES` |
| 4 | `2cmt_fo_lag_lin.ctl` | First-order absorption with lag time | Linear elimination | ADVAN4 TRANS4 |
| 5 | `2cmt_fo_lag_mm.ctl` | First-order absorption with lag time | Nonlinear elimination | ADVAN13 `$DES` |
| 6 | `2cmt_fo_lag_mmlin.ctl` | First-order absorption with lag time | Nonlinear and linear elimination | ADVAN13 `$DES` |
| 7 | `2cmt_sigmoid_dur_lin.ctl` | Sigmoid absorption with duration | Linear elimination | ADVAN4 TRANS4 |
| 8 | `2cmt_sigmoid_dur_mm.ctl` | Sigmoid absorption with duration | Nonlinear elimination | ADVAN13 `$DES` |
| 9 | `2cmt_sigmoid_dur_mmlin.ctl` | Sigmoid absorption with duration | Nonlinear and linear elimination | ADVAN13 `$DES` |
| 10 | `2cmt_sigmoid_rate_lin.ctl` | Sigmoid absorption with rate | Linear elimination | ADVAN4 TRANS4 |
| 11 | `2cmt_sigmoid_rate_mm.ctl` | Sigmoid absorption with rate | Nonlinear elimination | ADVAN13 `$DES` |
| 12 | `2cmt_sigmoid_rate_mmlin.ctl` | Sigmoid absorption with rate | Nonlinear and linear elimination | ADVAN13 `$DES` |
| 13 | `2cmt_zo_dur_lin.ctl` | Zero-order absorption with duration | Linear elimination | ADVAN3 TRANS4 |
| 14 | `2cmt_zo_dur_mm.ctl` | Zero-order absorption with duration | Nonlinear elimination | ADVAN13 `$DES` |
| 15 | `2cmt_zo_dur_mmlin.ctl` | Zero-order absorption with duration | Nonlinear and linear elimination | ADVAN13 `$DES` |
| 16 | `2cmt_zo_rate_lin.ctl` | Zero-order absorption with rate | Linear elimination | ADVAN3 TRANS4 |
| 17 | `2cmt_zo_rate_mm.ctl` | Zero-order absorption with rate | Nonlinear elimination | ADVAN13 `$DES` |
| 18 | `2cmt_zo_rate_mmlin.ctl` | Zero-order absorption with rate | Nonlinear and linear elimination | ADVAN13 `$DES` |
| 19 | `2cmt_transit_lin.ctl` | Transit compartment absorption | Linear elimination | ADVAN13 `$DES` |
| 20 | `2cmt_transit_mm.ctl` | Transit compartment absorption | Nonlinear elimination | ADVAN13 `$DES` |
| 21 | `2cmt_transit_mmlin.ctl` | Transit compartment absorption | Nonlinear and linear elimination | ADVAN13 `$DES` |
| 22 | `2cmt_parallel_fo_nolag_lin.ctl` | Parallel first-order absorption without lag time | Linear elimination | ADVAN13 `$DES` |
| 23 | `2cmt_parallel_fo_nolag_mm.ctl` | Parallel first-order absorption without lag time | Nonlinear elimination | ADVAN13 `$DES` |
| 24 | `2cmt_parallel_fo_nolag_mmlin.ctl` | Parallel first-order absorption without lag time | Nonlinear and linear elimination | ADVAN13 `$DES` |
| 25 | `2cmt_parallel_fo_lag_lin.ctl` | Parallel first-order absorption with lag time | Linear elimination | ADVAN13 `$DES` |
| 26 | `2cmt_parallel_fo_lag_mm.ctl` | Parallel first-order absorption with lag time | Nonlinear elimination | ADVAN13 `$DES` |
| 27 | `2cmt_parallel_fo_lag_mmlin.ctl` | Parallel first-order absorption with lag time | Nonlinear and linear elimination | ADVAN13 `$DES` |
| 28 | `2cmt_seq_fo_zo_dur_lin.ctl` | Sequential absorption with first-order process with lag time followed by zero-order with duration | Linear elimination | ADVAN4 TRANS4 |
| 29 | `2cmt_seq_fo_zo_dur_mm.ctl` | Sequential absorption with first-order process with lag time followed by zero-order with duration | Nonlinear elimination | ADVAN13 `$DES` |
| 30 | `2cmt_seq_fo_zo_dur_mmlin.ctl` | Sequential absorption with first-order process with lag time followed by zero-order with duration | Nonlinear and linear elimination | ADVAN13 `$DES` |
| 31 | `2cmt_seq_fo_zo_rate_lin.ctl` | Sequential absorption with first-order process with lag time followed by zero-order with rate | Linear elimination | ADVAN4 TRANS4 |
| 32 | `2cmt_seq_fo_zo_rate_mm.ctl` | Sequential absorption with first-order process with lag time followed by zero-order with rate | Nonlinear elimination | ADVAN13 `$DES` |
| 33 | `2cmt_seq_fo_zo_rate_mmlin.ctl` | Sequential absorption with first-order process with lag time followed by zero-order with rate | Nonlinear and linear elimination | ADVAN13 `$DES` |
| 34 | `2cmt_parallel_fo_zo_dur_lin.ctl` | Parallel absorption with first-order process and zero-order with duration | Linear elimination | ADVAN4 TRANS4 |
| 35 | `2cmt_parallel_fo_zo_dur_mm.ctl` | Parallel absorption with first-order process and zero-order with duration | Nonlinear elimination | ADVAN13 `$DES` |
| 36 | `2cmt_parallel_fo_zo_dur_mmlin.ctl` | Parallel absorption with first-order process and zero-order with duration | Nonlinear and linear elimination | ADVAN13 `$DES` |
| 37 | `2cmt_parallel_fo_zo_rate_lin.ctl` | Parallel absorption with first-order process and zero-order with rate | Linear elimination | ADVAN4 TRANS4 |
| 38 | `2cmt_parallel_fo_zo_rate_mm.ctl` | Parallel absorption with first-order process and zero-order with rate | Nonlinear elimination | ADVAN13 `$DES` |
| 39 | `2cmt_parallel_fo_zo_rate_mmlin.ctl` | Parallel absorption with first-order process and zero-order with rate | Nonlinear and linear elimination | ADVAN13 `$DES` |

Plus this `README.md`.

---

## Naming convention

```
2cmt_<absorption>_<elimination>.ctl
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

## Volume naming - important

Volume names follow the NONMEM `TRANS4` convention, which is also what the
diagrams use, so the **same symbol means different things in the two
layouts**:

| Layout | Models | Central | Peripheral |
|--------|--------|---------|------------|
| With a depot compartment | 1-12, 19-39 | `V2` | `V3` |
| No depot (dose straight into central) | 13-18 | `V1` | `V2` |

`Q` is the inter-compartmental clearance in every model.

### Distribution

```
dA_central/dt   = ... - Q*C_central + Q*C_periph - elimination
dA_periph/dt    =       Q*C_central - Q*C_periph
```

In the `$DES` models this is written with two intermediates, `CONC`
(central concentration) and `CP` (peripheral concentration), plus `ELR`
(the elimination rate), which keeps every line inside 80 columns.

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
[Central (V2)] <--Q--> [Periph. (V3)]  --elimination-->
```

* One dose record per administration: CMT=1 (depot), RATE=0.

**First-order absorption with lag time**  (`fo_lag`)

```
Dose at t=ALAG1 --> [Depot] --Ka--> [Central (V2)]
[Central (V2)] <--Q--> [Periph. (V3)]  --elimination-->
```

* One dose record per administration: CMT=1 (depot), RATE=0.
* Absorption starts at TIME+ALAG1.

**Sigmoid absorption with duration**  (`sigmoid_dur`)

```
Dose --(zero-order, D1)--> [Depot] --Ka--> [Central (V2)]
[Central (V2)] <--Q--> [Periph. (V3)]  --elimination-->
```

* One dose record per administration: CMT=1 (depot), RATE=-2
* (duration D1 modelled in $PK).  Zero-order input INTO the depot
* followed by first-order transfer to central gives the sigmoid
* (S-shaped) plasma profile.

**Sigmoid absorption with rate**  (`sigmoid_rate`)

```
Dose --(zero-order, R1)--> [Depot] --Ka--> [Central (V2)]
[Central (V2)] <--Q--> [Periph. (V3)]  --elimination-->
```

* One dose record per administration: CMT=1 (depot), RATE=-1
* (rate R1 modelled in $PK).

**Zero-order absorption with duration**  (`zo_dur`)

```
Dose --(zero-order, D1)--> [Central (V1)]
[Central (V1)] <--Q--> [Periph. (V2)]  --elimination-->
```

* One dose record per administration: CMT=1 (central), RATE=-2
* (duration D1 modelled in $PK).

**Zero-order absorption with rate**  (`zo_rate`)

```
Dose --(zero-order, R1)--> [Central (V1)]
[Central (V1)] <--Q--> [Periph. (V2)]  --elimination-->
```

* One dose record per administration: CMT=1 (central), RATE=-1
* (rate R1 modelled in $PK).

**Transit compartment absorption**  (`transit`)

```
Dose --> [T1] --Ktr--> ... --Ktr--> [Tn] --Ktr--> [Depot] --Ka-->
     [Central (V2)]     N transit cmts,  Ktr = (N+1)/MTT
[Central (V2)] <--Q--> [Periph. (V3)]  --elimination-->
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
[Central (V2)] <--Q--> [Periph. (V3)]  --elimination-->
```

* TWO dose records per administration, both with the full AMT:
  CMT=1 (depot 1), RATE=0   -> receives DF     x Dose via F1
  CMT=2 (depot 2), RATE=0   -> receives (1-DF) x Dose via F2

**Parallel first-order absorption with lag time**  (`parallel_fo_lag`)

```
DF x Dose     at t=ALAG1 --> [Depot1] --Ka1--\
                                              >--> [Central (V2)]
(1-DF) x Dose at t=ALAG2 --> [Depot2] --Ka2--/
[Central (V2)] <--Q--> [Periph. (V3)]  --elimination-->
```

* TWO dose records per administration, both with the full AMT:
  CMT=1 (depot 1), RATE=0   -> DF     x Dose, delayed by ALAG1
  CMT=2 (depot 2), RATE=0   -> (1-DF) x Dose, delayed by ALAG2

**Sequential absorption with first-order process with lag time followed by zero-order with duration**  (`seq_fo_zo_dur`)

```
DF x Dose     --(zero-order, D2)--> [Central (V2)]
(1-DF) x Dose at t=D2 --> [Depot] --Ka--^
[Central (V2)] <--Q--> [Periph. (V3)]  --elimination-->
```

* TWO dose records per administration, both with the full AMT:
  CMT=2 (central), RATE=-2  -> DF     x Dose, zero-order over D2
  CMT=1 (depot),   RATE=0   -> (1-DF) x Dose, released at t=D2
  (ALAG1 = D2) then absorbed with Ka

**Sequential absorption with first-order process with lag time followed by zero-order with rate**  (`seq_fo_zo_rate`)

```
DF x Dose     --(zero-order, R2)--> [Central (V2)]
(1-DF) x Dose at t=DF x Dose/R2 --> [Depot] --Ka--^
[Central (V2)] <--Q--> [Periph. (V3)]  --elimination-->
```

* TWO dose records per administration, both with the full AMT:
  CMT=2 (central), RATE=-1  -> DF     x Dose, zero-order at rate R2
  CMT=1 (depot),   RATE=0   -> (1-DF) x Dose, released at
  t = DF x Dose / R2, absorbed with Ka

**Parallel absorption with first-order process and zero-order with duration**  (`parallel_fo_zo_dur`)

```
DF x Dose     --(zero-order, D2)--> [Central (V2)]
(1-DF) x Dose --> [Depot] --Ka--^
[Central (V2)] <--Q--> [Periph. (V3)]  --elimination-->
```

* TWO dose records per administration, both with the full AMT:
  CMT=2 (central), RATE=-2  -> DF     x Dose, zero-order over D2
  CMT=1 (depot),   RATE=0   -> (1-DF) x Dose, absorbed with Ka
  (starts immediately - parallel)

**Parallel absorption with first-order process and zero-order with rate**  (`parallel_fo_zo_rate`)

```
DF x Dose     --(zero-order, R2)--> [Central (V2)]
(1-DF) x Dose --> [Depot] --Ka--^
[Central (V2)] <--Q--> [Periph. (V3)]  --elimination-->
```

* TWO dose records per administration, both with the full AMT:
  CMT=2 (central), RATE=-1  -> DF     x Dose, zero-order at rate R2
  CMT=1 (depot),   RATE=0   -> (1-DF) x Dose, absorbed with Ka

---

## Solver choice

Closed-form (analytical) ADVANs are used wherever the model permits:

* `ADVAN3 TRANS4` (`CL`, `V1`, `Q`, `V2`) - two compartments, no depot
  (models 13, 16).
* `ADVAN4 TRANS4` (`CL`, `V2`, `Q`, `V3`, `KA`) - two compartments with a
  first-order depot (models 1, 4, 7, 10, 28, 31, 34, 37).

`ADVAN13` with an explicit `$DES` block is used whenever the closed form does
not apply:

* any model with **nonlinear (Michaelis-Menten) elimination**;
* **parallel first-order absorption** (two depot compartments);
* **transit compartment absorption** (gamma-density input term).

`TOL=9` is set on `$SUBROUTINES` for all ODE models. Two-compartment ODE
models with Michaelis-Menten elimination are stiff at high concentrations; if
a run stalls, try `TOL=6` with `$ESTIMATION ... NSIG=2`, or switch to
`ADVAN13` with `SIGL=6`.

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

The same data set drives the one- and two-compartment libraries: the
`CMT` numbering for dose records is unchanged, because the peripheral
compartment is always added **after** the dosing compartments.

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

* IIV is **estimated** on `CL`, the central and peripheral volumes, `Q`,
  `VM`, `KA`, `KA1`, `KA2` and `MTT` (initial variance 0.09, ~30% CV).
* IIV on `KM`, `D1`, `D2`, `R1`, `R2`, `ALAG1`, `ALAG2` and `DF` is coded but
  set to `$OMEGA 0 FIX`. These random effects are frequently unidentifiable;
  free them one at a time only if the data support it.
* `$OMEGA` is diagonal. A two-compartment model usually needs a `$OMEGA BLOCK`
  on `CL`/central volume - add it once the diagonal model converges.
* The number of transit compartments `NN` is a THETA with no ETA.
* Residual error is combined proportional + additive:
  `Y = IPRED + IPRED*EPS(1) + EPS(2)`.

Initial estimates are generic placeholders (`CL` 5 L/h, central volume 50 L,
`Q` 10 L/h, peripheral volume 100 L, `KA` 1 /h, `Vmax` 100 mg/h, `Km` 5
mg/L). **Replace them with values appropriate to your compound and units
before running.**

---

## Running a model

```bash
# plain NONMEM
nmfe75 2cmt_fo_nolag_lin.ctl 2cmt_fo_nolag_lin.lst

# PsN
execute 2cmt_fo_nolag_lin.ctl -dir=2cmt_fo_nolag_lin
```

Every `$TABLE` writes a `.tab` file. Each run produces two:

* `<model>.tab` - per-record predictions and residuals.
* `<model>.par` - one row per subject with individual parameters and ETAs.

## Simulating instead of estimating

Every file ends with a commented simulation block. Comment out
`$ESTIMATION` and `$COVARIANCE`, then un-comment it.

---

## Caveats

1. **Transit models (19-21)** use the Savic analytical approximation, which
   tracks only the most recent dose. Exact for single-dose data; for
   multiple dosing, code the transit chain explicitly with a fixed `N`.
2. **Sequential models (28-33)** set `ALAG1` from the duration or from
   `DF*Dose/R2`, so the first-order process begins exactly when the
   zero-order input ends.
3. Nonlinear-elimination models are sensitive to initial estimates and to the
   dose range in the data; without data spanning the saturable region, `VM`
   and `KM` will not be separately identifiable.
4. Two-compartment models need a sampling schedule rich enough to resolve the
   distribution phase, otherwise `Q` and the peripheral volume collapse.
5. `$COVARIANCE` is enabled everywhere; drop it for exploratory runs.
6. Units are assumed consistent (mg, L, h). Check the `S1`/`S2`/`S3` scaling
   if your `DV` is in ng/mL while `AMT` is in mg.

## Provenance

The 39 structures are the standard two-compartment
absorption/elimination grid used throughout the pharmacometrics
literature; the root `../README.md` lists the background reading
behind each model family.
