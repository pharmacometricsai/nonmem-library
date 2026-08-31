# NONMEM model library

## Experimental, please check before use.

A library of **210 runnable NONMEM control streams** covering the standard
structural models used in population PK and PK/PD analysis: compartmental
disposition with thirteen absorption sub-models, target-mediated drug
disposition, and indirect response.

Every file is a complete template - structural model, between-subject
variability, residual error, estimation, covariance step, output tables and a
commented-out simulation block. Nothing needs to be assembled before a run;
point `$DATA` at your data set, replace the initial estimates, and go.

---

## Layout

| Folder | Files | Contents |
|--------|-------|----------|
| [`1cmt/`](1cmt/) | 39 | One-compartment disposition x 13 absorption models x 3 elimination forms |
| [`2cmt/`](2cmt/) | 39 | The same grid with one peripheral compartment |
| [`3cmt/`](3cmt/) | 39 | The same grid with two peripheral compartments |
| [`tmdd/`](tmdd/) | 45 | Target-mediated drug disposition: 3 structures x 5 approximations x 3 input routes |
| [`idr/`](idr/) | 48 | Indirect response models I-IV x 6 PK backbones x 2 exposure links |

Each folder has its own `README.md` with a full contents table, the equations,
the data-set requirements and the caveats specific to that family. **Read the
folder README before running anything from it** - the parameter naming and the
unit conventions differ between families in ways that matter.

---

## Quick start

```bash
cd 1cmt
nmfe75 1cmt_fo_nolag_lin.ctl 1cmt_fo_nolag_lin.lst

# or with PsN
execute 1cmt_fo_nolag_lin.ctl -dir=1cmt_fo_nolag_lin
```

Before the first run of any file:

1. Point `$DATA` at your data set (the paths here are placeholders under
   `../data/`).
2. Replace every initial estimate. The values shipped are generic
   placeholders, not priors.
3. Check the unit table below against your data.

---

## File naming

```
<family>_<structure>_<variant>.ctl
```

| Folder | Pattern | Example |
|--------|---------|---------|
| `1cmt/`, `2cmt/`, `3cmt/` | `<n>cmt_<absorption>_<elimination>` | `2cmt_transit_mmlin.ctl` |
| `tmdd/` | `tmdd_<n>cmt_<form>_<route>` | `tmdd_3cmt_qss_sc.ctl` |
| `idr/` | `idr<1-4>_<n>cmt_<route>_<link>` | `idr2_2cmt_oral_effect.ctl` |

Model numbers appear only in each folder's Contents table - not in the file
name, and not inside the files. Output tables are named after the control
stream that wrote them, so `2cmt_transit_mmlin.ctl` produces
`2cmt_transit_mmlin.tab`, `2cmt_transit_mmlin.par`, and - where the family has
one - `2cmt_transit_mmlin.prof`.

---

## Shared conventions

These hold across all five folders unless a folder README says otherwise.

### Units - check this first

| Folder | Concentration | Amount | Time |
|--------|---------------|--------|------|
| `1cmt/`, `2cmt/`, `3cmt/`, `idr/` | mg/L | mg | **hours** |
| `tmdd/` | **nM** | **nmol** | **days** |

The TMDD library is deliberately different. Its binding terms (`KON*C*R`) are
only meaningful if drug and target are in the same molar unit, and antibody
kinetics are conventionally expressed per day. **Do not share a data set
between `tmdd/` and the other folders without rescaling.**

### Data set

```
ID TIME AMT RATE DV EVID MDV CMT ADDL II WT          (1cmt, 2cmt, 3cmt)
ID TIME AMT RATE DV EVID MDV CMT FLAG ADDL II WT     (tmdd, idr)
```

Comma-delimited with a commented header row (`IGNORE=@`).

`RATE` carries the input mechanism:

| Value | Meaning |
|-------|---------|
| `0` | instantaneous - IV bolus, or a dose into a depot |
| `> 0` | IV infusion at that rate |
| `-1` | rate modelled in `$PK` (`R1` or `R2`) |
| `-2` | duration modelled in `$PK` (`D1` or `D2`) |

`FLAG` appears in the two families that observe more than one analyte. All
observation records use `CMT` = the central drug compartment, and `FLAG`
selects which prediction is compared to `DV`:

| Folder | `FLAG=1` | `FLAG=2` |
|--------|----------|----------|
| `tmdd/` | total drug in plasma | total target |
| `idr/` | plasma drug concentration | response |

Each analyte gets its own residual error model (`$SIGMA` 1-2 and 3-4). They
are different assays in different units; sharing an error model between them
would be meaningless.

### Two dose records per administration

Models that split a dose between two input routes - the `parallel_*` and
`seq_*` absorption models in `1cmt/`, `2cmt/` and `3cmt/` - need **two dose
records at the same `TIME`, each carrying the full `AMT`**, one per
compartment. NONMEM applies `F1` and `F2` to split them. The folder READMEs
carry worked examples.

### Parameterisation

* Fixed effects are log-normally distributed: `P = THETA(n)*EXP(ETA(n))`.
* Parameters bounded on (0, 1) - bioavailability `F1`, dose fraction `DF`,
  maximum inhibition `IMAX` - use **logit-scale** IIV. An exponential ETA on a
  parameter capped at 1 pushes individual values above the cap.
* IIV is estimated on the parameters that usually carry it (clearances,
  volumes, absorption and turnover rates, potency) and coded but fixed at
  `$OMEGA 0 FIX` on the ones that usually do not (binding constants, Hill
  slopes, lag times, maximum effects). Free them one at a time, with a reason.
* `$OMEGA` is diagonal throughout. A `$OMEGA BLOCK` on CL and central volume
  is normally the first thing to add once a diagonal model converges.
* Residual error is combined proportional + additive.
* Structural constraints that keep a model interpretable are **built in, not
  estimated**: `KSYN = KDEG*R0` in `tmdd/`, `KIN = R0*KOUT` in `idr/`. Both
  are paired with an `A_0` initial condition so the system starts at its
  pre-dose steady state. Breaking either constraint makes the baseline
  parameter stop meaning what its name says.

### Solvers

Closed-form ADVANs are used wherever the model admits them, because they are
an order of magnitude faster and numerically exact:

| ADVAN | Used for |
|-------|----------|
| `ADVAN1`/`ADVAN2 TRANS2` | 1-compartment, linear elimination |
| `ADVAN3`/`ADVAN4 TRANS4` | 2-compartment, linear elimination |
| `ADVAN11`/`ADVAN12 TRANS4` | 3-compartment, linear elimination |
| `ADVAN13` + `$DES`, `TOL=9` | everything else |

`ADVAN13` is used whenever the closed form does not apply: any Michaelis-Menten
elimination, parallel first-order absorption (two depots), transit chains, and
every model in `tmdd/` and `idr/`.

**A naming consequence worth knowing.** The closed-form ADVANs impose the
`TRANS4` parameter names, and those shift when a depot is added: `V2` is the
central volume in `2cmt_fo_nolag_lin.ctl` but a *peripheral* volume in
`2cmt_zo_dur_lin.ctl`, which has no depot. The `tmdd/` and `idr/` files are all
ODE models and therefore use one consistent naming regardless of route. This
is the single most common source of mix-ups when comparing estimates across
folders; each folder README states its own table explicitly.

### Estimation

Every file ships with FOCE-I:

```
$ESTIMATION METHOD=1 INTER MAXEVAL=9999 NSIG=3 SIGL=9 PRINT=5 NOABORT POSTHOC
```

The `tmdd/` and `idr/` files also carry a commented two-step SAEM + IMP
alternative, which is often a more robust route to the optimum for a
nonlinear or joint PK/PD fit. `$COVARIANCE` is enabled everywhere; drop it for
exploratory runs.

### Output tables

Every `$TABLE` in this library writes a `.tab`-family file named after its
control stream:

| File | Scope | Contents |
|------|-------|----------|
| `<model>.tab` | one row per record | `IPRED IRES IWRES CWRES PRED RES WRES`, the `ET` and `SG` columns below, and the input items |
| `<model>.par` | one row per subject (`FIRSTONLY`) | the individual structural parameters and the `ET` columns |
| `<model>.prof` | one row per record | analyte/response profiles - `tmdd/` and `idr/` only |
| `<model>_sim.tab` | one row per record | written by the commented-out simulation block |

`CWRES` is requested in every main table. It is a built-in NONMEM 7 item, so
no `$ABBREVIATED COMRES` machinery is needed.

Two sets of derived columns are exported so that diagnostics can be built
without going back to the `.ext` file:

* **`ET1`, `ET2`, ...** - the individual random effects, defined explicitly in
  `$PK` as `ETn = ETA(n)`, one per `$OMEGA`. Writing them out as named
  variables rather than relying on the `ETAn` table alias keeps them portable
  across NONMEM versions.
* **`SG1`, `SG2`, ...** - the residual error variances, defined in `$ERROR` as
  `SGn = SIGMA(n,n)`, one per `$SIGMA`. `SIGMA(i,j)` is a reserved variable in
  NM-TRAN abbreviated code from NONMEM 7 onward; on an older version, delete
  the two or four lines under the "residual error variances" comment and the
  matching `SG` entries in `$TABLE`.

In `1cmt/`, `2cmt/` and `3cmt/` there are two `SG` columns (proportional,
additive). In `tmdd/` and `idr/` there are four, because each of the two
analytes carries its own error model.

---

## The five sub-libraries

### `1cmt/`, `2cmt/`, `3cmt/` - compartmental PK

The same 39-model grid at three levels of disposition complexity: **13
absorption sub-models x 3 elimination forms**.

Absorption: first-order (with and without lag), sigmoid (zero-order input into
a depot followed by first-order transfer, with duration or rate), zero-order
directly into central (duration or rate), transit compartments, parallel
first-order from two depots (with and without lag), sequential first-order
after zero-order (duration or rate), and parallel first-order plus zero-order
(duration or rate).

Elimination: linear (`CL/V`), Michaelis-Menten (`Vmax/Km`), or both in
parallel.

The three folders share their absorption taxonomy exactly, so
`1cmt_transit_mmlin.ctl`, `2cmt_transit_mmlin.ctl` and
`3cmt_transit_mmlin.ctl` differ only in disposition. That is what makes them
usable as a structural search: fix the absorption model, walk the disposition
ladder, and compare.

Two things to know before running them:

* The **transit** models use the Savic analytical approximation, which tracks
  only the most recent dose. It is exact for single-dose data; for multiple
  dosing, code the chain explicitly with a fixed `N`.
* The **sequential** models set `ALAG1` from the zero-order duration or from
  `DF*Dose/R2`, so the first-order process begins exactly when the zero-order
  input ends.

### `tmdd/` - target-mediated drug disposition

Three plasma structures (1-, 2-, 3-compartment) wrapped around one
target-binding module, expanded across **5 model forms x 3 input routes**.

The five forms are an approximation ladder, each row a strictly stronger set of
assumptions than the one above:

| Form | States | Assumption |
|------|--------|------------|
| `full` | free drug, free target, complex | none beyond the mechanism |
| `qe` | total drug, total target | binding at equilibrium, `KD = KOFF/KON` |
| `qss` | total drug, total target | complex at quasi-steady state, `KSS = (KOFF+KINT)/KON` |
| `wagner` | total drug | quasi-equilibrium **and** total target held at `R0` |
| `mm` | drug only | binding collapses to `VM*C/(KM+C)` |

Input routes: IV (bolus or infusion), subcutaneous first-order, subcutaneous
with lag.

**Start at `qss`.** It is the most generally applicable approximation and by
far the most likely to converge. Move down to `wagner` or `mm` when the
target-turnover parameters are not identifiable; move up to `full` only when
you have target or complex measurements across a dose range - without them
`KON` and `KOFF` are individually unidentifiable and only their ratio is
informed.

### `idr/` - indirect response

The four basic indirect response models, each paired with **6 PK backbones**
(1/2/3-compartment x IV/oral) and **2 exposure links** (plasma-driven, or an
effect compartment with `KE0`), fitted as simultaneous PK/PD.

| Model | Drug acts on | Response |
|-------|--------------|----------|
| I | production, inhibited | decreases |
| II | elimination, inhibited | increases |
| III | production, stimulated | increases |
| IV | elimination, stimulated | decreases |

Models I and IV both push the response down; II and III both push it up.
Within each pair, the discriminator is not the direction but the **return to
baseline**: I and III act on `KIN`, so recovery happens at a dose-independent
rate set by `KOUT`; II and IV act on `KOUT` itself, so recovery is
dose-dependent. Distinguishing them needs more than one dose level and enough
sampling on the washout limb. Without both, the choice is mechanistic rather
than statistical.

The drug effect is a sigmoid `EFF = (CD/XC50)**HILL / (1 + (CD/XC50)**HILL)`.
Fix `HILL` to 1 to recover the plain Emax form and free it only if the
concentration-effect relationship is visibly steeper or shallower.

---

## Background reading

Citations are given as author, title, journal and year - enough to find each
paper, without page numbers that are easy to get wrong.

### Compartmental disposition and absorption

* Gibaldi M, Perrier D. *Pharmacokinetics*, 2nd ed. Marcel Dekker, 1982. The
  standard treatment of linear compartmental disposition and the closed-form
  solutions behind `ADVAN1`-`ADVAN4`, `ADVAN11` and `ADVAN12`.
* Beal SL, Sheiner LB, Boeckmann AJ, Bauer RJ (eds). *NONMEM Users Guides*,
  ICON plc. The authority on the `ADVAN`/`TRANS` definitions, `$DES`, `A_0`
  initial conditions and the `RATE = -1 / -2` conventions used throughout this
  library.
* Savic RM, Jonker DM, Kerbusch T, Karlsson MO. Implementation of a transit
  compartment model for describing drug absorption in pharmacokinetic studies.
  *Journal of Pharmacokinetics and Pharmacodynamics*, 2007. The transit chain,
  `KTR = (N+1)/MTT`, and the analytic gamma-density input used in the
  `transit` models here.
* Zhou H. Pharmacokinetic strategies in deciphering atypical drug absorption
  profiles. *Journal of Clinical Pharmacology*, 2003. Parallel and sequential
  zero-/first-order absorption and double-peak profiles - the background to
  the `parallel_*` and `seq_*` models.

### Nonlinear elimination

* Michaelis L, Menten ML. Die Kinetik der Invertinwirkung. *Biochemische
  Zeitschrift*, 1913. The origin of the saturable rate law that the `mm` and
  `mmlin` variants apply to drug elimination.

### Target-mediated drug disposition

* Levy G. Pharmacologic target-mediated drug disposition. *Clinical
  Pharmacology and Therapeutics*, 1994. The paper that named the phenomenon.
* Mager DE, Jusko WJ. General pharmacokinetic model for drugs exhibiting
  target-mediated drug disposition. *Journal of Pharmacokinetics and
  Pharmacodynamics*, 2001. The full model implemented in the `full` variants.
* Mager DE, Krzyzanski W. Quasi-equilibrium pharmacokinetic model for drugs
  exhibiting target-mediated drug disposition. *Pharmaceutical Research*,
  2005. The rapid-binding (QE) approximation.
* Gibiansky L, Gibiansky E, Kakkar T, Ma P. Approximations of the
  target-mediated drug disposition model and identifiability of model
  parameters. *Journal of Pharmacokinetics and Pharmacodynamics*, 2008. The
  QSS and Michaelis-Menten approximations, and the identifiability analysis
  that should govern which form you pick.
* Wagner JG, 1973. The constant-total-target quasi-equilibrium treatment of
  saturable binding, conventionally referred to as the Wagner model and
  implemented here as the `wagner` variant.
* Peletier LA, Gabrielsson J. Dynamics of target-mediated drug disposition:
  characteristic profiles and parameter identification. *Journal of
  Pharmacokinetics and Pharmacodynamics*, 2012. The mathematical anatomy of
  the characteristic four-phase TMDD profile.
* Dua P, Hawkins E, van der Graaf PH. A tutorial on target-mediated drug
  disposition (TMDD) models. *CPT: Pharmacometrics and Systems Pharmacology*,
  2015. The best entry point if this family is new to you.

### Indirect response and PK/PD linking

* Dayneka NL, Garg V, Jusko WJ. Comparison of four basic models of indirect
  pharmacodynamic responses. *Journal of Pharmacokinetics and
  Biopharmaceutics*, 1993. The source of models I-IV.
* Jusko WJ, Ko HC. Physiologic indirect response models characterize diverse
  types of pharmacodynamic effects. *Clinical Pharmacology and Therapeutics*,
  1994.
* Sharma A, Jusko WJ. Characteristics of indirect pharmacodynamic models and
  applications to clinical drug responses. *British Journal of Clinical
  Pharmacology*, 1998. How the four models actually behave, and what data are
  needed to tell them apart.
* Sheiner LB, Stanski DR, Vozeh S, Miller RD, Ham J. Simultaneous modeling of
  pharmacokinetics and pharmacodynamics: application to d-tubocurarine.
  *Clinical Pharmacology and Therapeutics*, 1979. The effect compartment and
  `KE0`, used in the `_effect` variants.
* Hill AV. The possible effects of the aggregation of the molecules of
  haemoglobin on its dissociation curves. *Journal of Physiology*, 1910. The
  origin of the `HILL` coefficient.

---

## Caveats

1. **The initial estimates are placeholders.** Every file ships with generic
   values chosen to be dimensionally sensible, not to be priors. Nonlinear
   models - Michaelis-Menten elimination, all of `tmdd/`, all of `idr/` - are
   sensitive to them, and a joint PK/PD or TMDD fit launched from generic
   starting values will usually fail. Fit the PK first and carry those
   estimates in.
2. **Structural richness is not identifiability.** A three-compartment model
   needs rich early sampling and a long terminal tail; the full TMDD model
   needs target measurements across a dose range; models II and IV in `idr/`
   need washout data at more than one dose. Fitting a structure the data
   cannot support produces a converged run and a meaningless parameter.
3. **Check the unit table** above before mixing data between folders.
4. **Nothing here has been run against real data.** These are structurally
   verified templates - parameter/ETA index consistency, compartment counts
   against `$DES`, correct `ADVAN` parameter names, records within 80 columns -
   not fitted models.
5. `$COVARIANCE` is enabled everywhere. Drop it for exploratory runs; a failed
   covariance step on a first pass tells you very little.
