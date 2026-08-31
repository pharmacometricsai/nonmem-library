# One-compartment NONMEM model library

NONMEM control streams (`.ctl`) for the **39 one-compartment models** shown on
the standard one-compartment PK model set. Each file is a
complete, runnable template: structural model, between-subject variability,
combined residual error, estimation, covariance step, output tables and a
commented-out simulation block.

---

## Contents

| # | File | Absorption | Elimination | Solver |
|---|------|------------|-------------|--------|
| 1 | `1cmt_fo_nolag_lin.ctl` | First-order absorption without lag time | Linear elimination | ADVAN2 TRANS2 |
| 2 | `1cmt_fo_nolag_mm.ctl` | First-order absorption without lag time | Nonlinear elimination | ADVAN13 `$DES` |
| 3 | `1cmt_fo_nolag_mmlin.ctl` | First-order absorption without lag time | Nonlinear and linear elimination | ADVAN13 `$DES` |
| 4 | `1cmt_fo_lag_lin.ctl` | First-order absorption with lag time | Linear elimination | ADVAN2 TRANS2 |
| 5 | `1cmt_fo_lag_mm.ctl` | First-order absorption with lag time | Nonlinear elimination | ADVAN13 `$DES` |
| 6 | `1cmt_fo_lag_mmlin.ctl` | First-order absorption with lag time | Nonlinear and linear elimination | ADVAN13 `$DES` |
| 7 | `1cmt_sigmoid_dur_lin.ctl` | Sigmoid absorption with duration | Linear elimination | ADVAN2 TRANS2 |
| 8 | `1cmt_sigmoid_dur_mm.ctl` | Sigmoid absorption with duration | Nonlinear elimination | ADVAN13 `$DES` |
| 9 | `1cmt_sigmoid_dur_mmlin.ctl` | Sigmoid absorption with duration | Nonlinear and linear elimination | ADVAN13 `$DES` |
| 10 | `1cmt_sigmoid_rate_lin.ctl` | Sigmoid absorption with rate | Linear elimination | ADVAN2 TRANS2 |
| 11 | `1cmt_sigmoid_rate_mm.ctl` | Sigmoid absorption with rate | Nonlinear elimination | ADVAN13 `$DES` |
| 12 | `1cmt_sigmoid_rate_mmlin.ctl` | Sigmoid absorption with rate | Nonlinear and linear elimination | ADVAN13 `$DES` |
| 13 | `1cmt_zo_dur_lin.ctl` | Zero-order absorption with duration | Linear elimination | ADVAN1 TRANS2 |
| 14 | `1cmt_zo_dur_mm.ctl` | Zero-order absorption with duration | Nonlinear elimination | ADVAN13 `$DES` |
| 15 | `1cmt_zo_dur_mmlin.ctl` | Zero-order absorption with duration | Nonlinear and linear elimination | ADVAN13 `$DES` |
| 16 | `1cmt_zo_rate_lin.ctl` | Zero-order absorption with rate | Linear elimination | ADVAN1 TRANS2 |
| 17 | `1cmt_zo_rate_mm.ctl` | Zero-order absorption with rate | Nonlinear elimination | ADVAN13 `$DES` |
| 18 | `1cmt_zo_rate_mmlin.ctl` | Zero-order absorption with rate | Nonlinear and linear elimination | ADVAN13 `$DES` |
| 19 | `1cmt_transit_lin.ctl` | Transit compartment absorption | Linear elimination | ADVAN13 `$DES` |
| 20 | `1cmt_transit_mm.ctl` | Transit compartment absorption | Nonlinear elimination | ADVAN13 `$DES` |
| 21 | `1cmt_transit_mmlin.ctl` | Transit compartment absorption | Nonlinear and linear elimination | ADVAN13 `$DES` |
| 22 | `1cmt_parallel_fo_nolag_lin.ctl` | Parallel first-order absorption without lag time | Linear elimination | ADVAN13 `$DES` |
| 23 | `1cmt_parallel_fo_nolag_mm.ctl` | Parallel first-order absorption without lag time | Nonlinear elimination | ADVAN13 `$DES` |
| 24 | `1cmt_parallel_fo_nolag_mmlin.ctl` | Parallel first-order absorption without lag time | Nonlinear and linear elimination | ADVAN13 `$DES` |
| 25 | `1cmt_parallel_fo_lag_lin.ctl` | Parallel first-order absorption with lag time | Linear elimination | ADVAN13 `$DES` |
| 26 | `1cmt_parallel_fo_lag_mm.ctl` | Parallel first-order absorption with lag time | Nonlinear elimination | ADVAN13 `$DES` |
| 27 | `1cmt_parallel_fo_lag_mmlin.ctl` | Parallel first-order absorption with lag time | Nonlinear and linear elimination | ADVAN13 `$DES` |
| 28 | `1cmt_seq_fo_zo_dur_lin.ctl` | Sequential absorption with first-order process with lag time followed by zero-order with duration | Linear elimination | ADVAN2 TRANS2 |
| 29 | `1cmt_seq_fo_zo_dur_mm.ctl` | Sequential absorption with first-order process with lag time followed by zero-order with duration | Nonlinear elimination | ADVAN13 `$DES` |
| 30 | `1cmt_seq_fo_zo_dur_mmlin.ctl` | Sequential absorption with first-order process with lag time followed by zero-order with duration | Nonlinear and linear elimination | ADVAN13 `$DES` |
| 31 | `1cmt_seq_fo_zo_rate_lin.ctl` | Sequential absorption with first-order process with lag time followed by zero-order with rate | Linear elimination | ADVAN2 TRANS2 |
| 32 | `1cmt_seq_fo_zo_rate_mm.ctl` | Sequential absorption with first-order process with lag time followed by zero-order with rate | Nonlinear elimination | ADVAN13 `$DES` |
| 33 | `1cmt_seq_fo_zo_rate_mmlin.ctl` | Sequential absorption with first-order process with lag time followed by zero-order with rate | Nonlinear and linear elimination | ADVAN13 `$DES` |
| 34 | `1cmt_parallel_fo_zo_dur_lin.ctl` | Parallel absorption with first-order process and zero-order with duration | Linear elimination | ADVAN2 TRANS2 |
| 35 | `1cmt_parallel_fo_zo_dur_mm.ctl` | Parallel absorption with first-order process and zero-order with duration | Nonlinear elimination | ADVAN13 `$DES` |
| 36 | `1cmt_parallel_fo_zo_dur_mmlin.ctl` | Parallel absorption with first-order process and zero-order with duration | Nonlinear and linear elimination | ADVAN13 `$DES` |
| 37 | `1cmt_parallel_fo_zo_rate_lin.ctl` | Parallel absorption with first-order process and zero-order with rate | Linear elimination | ADVAN2 TRANS2 |
| 38 | `1cmt_parallel_fo_zo_rate_mm.ctl` | Parallel absorption with first-order process and zero-order with rate | Nonlinear elimination | ADVAN13 `$DES` |
| 39 | `1cmt_parallel_fo_zo_rate_mmlin.ctl` | Parallel absorption with first-order process and zero-order with rate | Nonlinear and linear elimination | ADVAN13 `$DES` |

Plus this `README.md`.

---

## Naming convention

```
1cmt_<absorption>_<elimination>.ctl
```

Model numbers 1-39 are listed in the Contents table above; they are
not part of the file name and do not appear inside the files.

* `<absorption>` – absorption sub-model:
  * `fo_nolag` – First-order absorption without lag time
  * `fo_lag` – First-order absorption with lag time
  * `sigmoid_dur` – Sigmoid absorption with duration
  * `sigmoid_rate` – Sigmoid absorption with rate
  * `zo_dur` – Zero-order absorption with duration
  * `zo_rate` – Zero-order absorption with rate
  * `transit` – Transit compartment absorption
  * `parallel_fo_nolag` – Parallel first-order absorption without lag time
  * `parallel_fo_lag` – Parallel first-order absorption with lag time
  * `seq_fo_zo_dur` – Sequential absorption with first-order process with lag time followed by zero-order with duration
  * `seq_fo_zo_rate` – Sequential absorption with first-order process with lag time followed by zero-order with rate
  * `parallel_fo_zo_dur` – Parallel absorption with first-order process and zero-order with duration
  * `parallel_fo_zo_rate` – Parallel absorption with first-order process and zero-order with rate
* `<elimination>`:
  * `lin` – linear elimination (CL/V)
  * `mm` – nonlinear Michaelis-Menten elimination (Vmax/Km)
  * `mmlin` – parallel nonlinear and linear elimination (CL/V + Vmax/Km)

---

## Model structure

### Elimination

| Variant | Differential / analytical form | Parameters |
|---------|-------------------------------|------------|
| `lin` | `-CL*C` | `CL`, `V` |
| `mm` | `-Vmax*C/(Km+C)` | `VM`, `KM`, `V` |
| `mmlin` | `-CL*C - Vmax*C/(Km+C)` | `CL`, `V`, `VM`, `KM` |

`VM` is used as the NONMEM variable name for Vmax (`VMAX` exceeds the
safe 4-character abbreviation habit and `VM` avoids any clash with `V`).

### Absorption

**First-order absorption without lag time**  (`fo_nolag`)

```
Dose --(instantaneous)--> [Depot] --Ka--> [Central (V)] --elimination-->
```

* One dose record per administration: CMT=1 (depot), RATE=0.

**First-order absorption with lag time**  (`fo_lag`)

```
Dose at t=ALAG1 --> [Depot] --Ka--> [Central (V)] --elimination-->
```

* One dose record per administration: CMT=1 (depot), RATE=0.
* Absorption starts at TIME+ALAG1.

**Sigmoid absorption with duration**  (`sigmoid_dur`)

```
Dose --(zero-order, D1)--> [Depot] --Ka--> [Central (V)] --elimination-->
```

* One dose record per administration: CMT=1 (depot), RATE=-2
* (duration D1 modelled in $PK).  Zero-order input INTO the depot
* followed by first-order transfer to central gives the sigmoid
* (S-shaped) plasma profile.

**Sigmoid absorption with rate**  (`sigmoid_rate`)

```
Dose --(zero-order, R1)--> [Depot] --Ka--> [Central (V)] --elimination-->
```

* One dose record per administration: CMT=1 (depot), RATE=-1
* (rate R1 modelled in $PK).

**Zero-order absorption with duration**  (`zo_dur`)

```
Dose --(zero-order, D1)--> [Central (V)] --elimination-->
```

* One dose record per administration: CMT=1 (central), RATE=-2
* (duration D1 modelled in $PK).

**Zero-order absorption with rate**  (`zo_rate`)

```
Dose --(zero-order, R1)--> [Central (V)] --elimination-->
```

* One dose record per administration: CMT=1 (central), RATE=-1
* (rate R1 modelled in $PK).

**Transit compartment absorption**  (`transit`)

```
Dose --> [T1] --Ktr--> ... --Ktr--> [Tn] --Ktr--> [Depot] --Ka--> [Central (V)] --elimination-->
        N transit compartments,  Ktr = (N+1)/MTT
```

* One dose record per administration: CMT=1 (depot), RATE=0.
* NOTE: the analytical transit solution tracks only the MOST RECENT
* dose.  It is exact for single-dose data; for multiple dosing either
* code the transit chain explicitly with a fixed N, or use the
* superposition form.

**Parallel first-order absorption without lag time**  (`parallel_fo_nolag`)

```
DF x Dose     --> [Depot1] --Ka1--\
                                   >--> [Central (V)] --elimination-->
(1-DF) x Dose --> [Depot2] --Ka2--/
```

* TWO dose records per administration, both with the full AMT:
  CMT=1 (depot 1), RATE=0   -> receives DF     x Dose via F1
  CMT=2 (depot 2), RATE=0   -> receives (1-DF) x Dose via F2

**Parallel first-order absorption with lag time**  (`parallel_fo_lag`)

```
DF x Dose     at t=ALAG1 --> [Depot1] --Ka1--\
                                              >--> [Central (V)] --elimination-->
(1-DF) x Dose at t=ALAG2 --> [Depot2] --Ka2--/
```

* TWO dose records per administration, both with the full AMT:
  CMT=1 (depot 1), RATE=0   -> DF     x Dose, delayed by ALAG1
  CMT=2 (depot 2), RATE=0   -> (1-DF) x Dose, delayed by ALAG2

**Sequential absorption with first-order process with lag time followed by zero-order with duration**  (`seq_fo_zo_dur`)

```
DF x Dose     --(zero-order, D2)--> [Central (V)] --elimination-->
(1-DF) x Dose at t=D2 --> [Depot] --Ka--^
```

* TWO dose records per administration, both with the full AMT:
  CMT=2 (central), RATE=-2  -> DF     x Dose, zero-order over D2
  CMT=1 (depot),   RATE=0   -> (1-DF) x Dose, released at t=D2
  (ALAG1 = D2) then absorbed with Ka

**Sequential absorption with first-order process with lag time followed by zero-order with rate**  (`seq_fo_zo_rate`)

```
DF x Dose     --(zero-order, R2)--> [Central (V)] --elimination-->
(1-DF) x Dose at t=DF x Dose/R2 --> [Depot] --Ka--^
```

* TWO dose records per administration, both with the full AMT:
  CMT=2 (central), RATE=-1  -> DF     x Dose, zero-order at rate R2
  CMT=1 (depot),   RATE=0   -> (1-DF) x Dose, released at
  t = DF x Dose / R2, absorbed with Ka

**Parallel absorption with first-order process and zero-order with duration**  (`parallel_fo_zo_dur`)

```
DF x Dose     --(zero-order, D2)--> [Central (V)] --elimination-->
(1-DF) x Dose --> [Depot] --Ka--^
```

* TWO dose records per administration, both with the full AMT:
  CMT=2 (central), RATE=-2  -> DF     x Dose, zero-order over D2
  CMT=1 (depot),   RATE=0   -> (1-DF) x Dose, absorbed with Ka
  (starts immediately - parallel)

**Parallel absorption with first-order process and zero-order with rate**  (`parallel_fo_zo_rate`)

```
DF x Dose     --(zero-order, R2)--> [Central (V)] --elimination-->
(1-DF) x Dose --> [Depot] --Ka--^
```

* TWO dose records per administration, both with the full AMT:
  CMT=2 (central), RATE=-1  -> DF     x Dose, zero-order at rate R2
  CMT=1 (depot),   RATE=0   -> (1-DF) x Dose, absorbed with Ka

---

## Solver choice

Closed-form (analytical) ADVANs are used wherever the model permits, because
they are an order of magnitude faster and numerically exact:

* `ADVAN1 TRANS2` – one compartment, no depot (models 13, 16).
* `ADVAN2 TRANS2` – one compartment with a first-order depot
  (models 1, 4, 7, 10, 28, 31, 34, 37).

`ADVAN13` with an explicit `$DES` block is used whenever the closed form does
not apply:

* any model with **nonlinear (Michaelis-Menten) elimination**;
* **parallel first-order absorption** (two depot compartments);
* **transit compartment absorption** (gamma-density input term).

`TOL=9` is set on `$SUBROUTINES` for all ODE models.

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
| `WT` | body weight (kg) – carried for covariate work, not yet used |

### `RATE` conventions used here

| Absorption model | `RATE` on the dose record |
|------------------|---------------------------|
| First-order, transit, parallel first-order | `0` |
| Sigmoid / zero-order **with duration** | `-2` (`D1` or `D2` estimated in `$PK`) |
| Sigmoid / zero-order **with rate** | `-1` (`R1` or `R2` estimated in `$PK`) |

### Models that need two dose records per administration

Models 22-39 split the dose between two input routes using bioavailability
fractions. The data set must therefore carry **two dose records at the same
`TIME`, each with the full `AMT`**, one per compartment; NONMEM applies `F1`
and `F2` to split them. Example for model 34
(parallel first-order + zero-order with duration):

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
* The dose fraction `DF` is bounded on (0, 1) and given **logit-scale** IIV:

  ```
  LGTDF = LOG(THETA(n)/(1-THETA(n))) + ETA(n)
  DF    = EXP(LGTDF)/(1+EXP(LGTDF))
  ```

* IIV is **estimated** on `CL`, `V`, `VM`, `KA`, `KA1`, `KA2` and `MTT`
  (initial variance 0.09, i.e. ~30% CV).
* IIV on `KM`, `D1`, `D2`, `R1`, `R2`, `ALAG1`, `ALAG2` and `DF` is present in
  the code but **`$OMEGA 0 FIX`**. These random effects are frequently
  unidentifiable; free them one at a time only if the data support it.
* The number of transit compartments `NN` is a THETA with no ETA.
* Residual error is combined proportional + additive:
  `Y = IPRED + IPRED*EPS(1) + EPS(2)`.

Initial estimates are generic placeholders (`CL` 5 L/h, `V` 50 L, `KA` 1 /h,
`Vmax` 100 mg/h, `Km` 5 mg/L). **Replace them with values appropriate to your
compound and units before running.**

---

## Running a model

```bash
# plain NONMEM
nmfe75 1cmt_fo_nolag_lin.ctl 1cmt_fo_nolag_lin.lst

# PsN
execute 1cmt_fo_nolag_lin.ctl -dir=1cmt_fo_nolag_lin
```

Every `$TABLE` writes a `.tab` file. Each run produces two:

* `<model>.tab` – per-record predictions and residuals
  (`IPRED IRES IWRES CWRES PRED RES WRES`).
* `<model>.par` – one row per subject with the individual parameters and ETAs.

## Simulating instead of estimating

Every file ends with a commented simulation block. Comment out
`$ESTIMATION` and `$COVARIANCE`, then un-comment:

```
$SIMULATION (20260830) (20260831 NORMAL) ONLYSIM SUBPROBLEMS=200
$TABLE ID TIME AMT RATE EVID MDV CMT DV IPRED NOAPPEND ONEHEADER NOPRINT FILE=<model>_sim.tab
```

---

## Caveats

1. **Transit models (19-21)** use the Savic analytical approximation, which
   tracks only the most recent dose. It is exact for single-dose data; for
   multiple dosing, code the transit chain explicitly with a fixed `N`.
2. **Sequential models (28-33)** set `ALAG1` from the duration or from
   `DF*Dose/R2`, so the first-order process begins exactly when the
   zero-order input ends. `DOSE` is captured in `$PK` via `NEWIND`/`AMT`
   bookkeeping.
3. Nonlinear-elimination models are sensitive to initial estimates and to the
   dose range in the data; without data spanning the saturable region, `VM`
   and `KM` will not be separately identifiable.
4. `$COVARIANCE` is enabled everywhere; drop it for exploratory runs.
5. Units are assumed consistent (mg, L, h). Check `S1`/`S2`/`S3` scaling if
   your `DV` is in ng/mL while `AMT` is in mg.

## Provenance

The 39 structures are the standard one-compartment absorption/elimination
grid used throughout the pharmacometrics literature; the root `../README.md`
lists the background reading behind each model family.
