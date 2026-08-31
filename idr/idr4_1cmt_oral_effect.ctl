;; ==========================================================================
;; IDR4_1CMT_ORAL_EFFECT
;; Indirect response model, simultaneous PK/PD
;;   PD   : IDR model IV - Stimulation of elimination
;;   PK   : 1-compartment, Oral, first-order absorption
;;   Link : Effect compartment (KE0), for counter-clockwise hysteresis
;; --------------------------------------------------------------------------
;; Structure
;;   Dose --> [Depot] --KA--> [Central (V)]
;;   [Central (V)] --CL--> eliminated
;;   plasma --KE0--> [Effect site (Ce)] --> drives the response
;;   KIN --> [Response (R)] --KOUT-->
;;   dR/dt = KIN - KOUT*(1 + SMAX*EFF)*R
;;   EFF = (CD/SC50)**HILL / (1 + (CD/SC50)**HILL)
;;   where CD = effect-site concentration Ce
;;
;; Behaviour
;;   The drug stimulates the first-order loss of the response.
;;   Response falls below baseline.  Because the drug acts on KOUT,
;;   the apparent turnover accelerates while drug is present - the
;;   return to baseline is dose-dependent, which is what separates
;;   model IV from model I.
;;   Net direction of the response: decrease from baseline.
;;
;; Baseline
;;   KIN is NOT a separate parameter: KIN = R0*KOUT, and the response
;;   compartment is initialised at A_0 = R0.  The response therefore
;;   starts at, and returns to, R0 by construction.  Do not add a
;;   separate baseline THETA on top of this.
;;
;; Data set / dosing requirements
;;   Dose records: CMT=1 (depot), RATE=0.
;;   Observation records: CMT=2 for BOTH analytes; the FLAG column
;;   selects which one is being observed:
;;     FLAG=1  plasma drug concentration (mg/L)
;;     FLAG=2  response R (response units)
;;
;; Units: mg, L, h - the same as ../1cmt, ../2cmt and ../3cmt.
;;   (../tmdd uses nM and days; do not share a data set with it.)
;;
;; Solver: general non-linear ODE (ADVAN13)
;; ==========================================================================

$PROBLEM IDR4 1cmt oral effect

$INPUT   ID TIME AMT RATE DV EVID MDV CMT FLAG ADDL II WT
$DATA    ../data/pkpddata.csv IGNORE=@

$SUBROUTINES ADVAN13 TOL=9
$MODEL  NCOMPARTMENTS=4
        COMP=(DEPOT, DEFDOSE)
        COMP=(CENTRAL, DEFOBS)
        COMP=(EFFECT)
        COMP=(RESPONSE)

; compartment 1 = oral depot
; compartment 2 = central drug
; compartment 3 = effect site (holds Ce directly)
; compartment 4 = response R (holds R directly)
; The effect and response compartments are unscaled: A() holds the
; effect-site concentration and the response itself, not an amount.

$PK
; ---- typical values and between-subject variability ------------------
  CL    = THETA(1)*EXP(ETA(1))
  V     = THETA(2)*EXP(ETA(2))
  KA    = THETA(3)*EXP(ETA(3))
  LGTF1 = LOG(THETA(4)/(1-THETA(4))) + ETA(4)
  F1    = EXP(LGTF1)/(1+EXP(LGTF1))
  KE0   = THETA(5)*EXP(ETA(5))
  R0    = THETA(6)*EXP(ETA(6))
  KOUT  = THETA(7)*EXP(ETA(7))
  SMAX  = THETA(8)*EXP(ETA(8))
  SC50  = THETA(9)*EXP(ETA(9))
  HILL  = THETA(10)*EXP(ETA(10))

; ---- random effects exported to $TABLE --------------------------------
  ET1   = ETA(1)
  ET2   = ETA(2)
  ET3   = ETA(3)
  ET4   = ETA(4)
  ET5   = ETA(5)
  ET6   = ETA(6)
  ET7   = ETA(7)
  ET8   = ETA(8)
  ET9   = ETA(9)
  ET10  = ETA(10)

; ---- response baseline -----------------------------------------------
;     KIN = R0*KOUT, so the response sits at R0 before the first dose.
  KIN   = R0*KOUT
  A_0(4) = R0

; ---- scaling (concentration = amount / volume) ------------------------
  S2 = V

$DES
  CP   = A(2)/V
  IF(CP.LT.0) CP = 0
  DADT(1) = -KA*A(1)
  DADT(2) = KA*A(1) - CL*CP
  DADT(3) = KE0*(CP - A(3))
  CD   = A(3)
  EFF  = 0
  IF(CD.GT.0) THEN
    XX  = (CD/SC50)**HILL
    EFF = XX/(1+XX)
  ENDIF
  DADT(4) = KIN - KOUT*(1 + SMAX*EFF)*A(4)

$ERROR
  CP   = A(2)/V
  IF(CP.LT.0) CP = 0
  CE   = A(3)
  CD   = A(3)
  EFF  = 0
  IF(CD.GT.0) THEN
    XX  = (CD/SC50)**HILL
    EFF = XX/(1+XX)
  ENDIF
  FACT = 1 + SMAX*EFF        ; factor applied to KOUT
  RESP = A(4)
  IF(FLAG.EQ.2) THEN
    IPRED = RESP                        ; response
    IF(IPRED.LE.0) IPRED = 1.0E-10
    Y     = IPRED + IPRED*EPS(3) + EPS(4)
  ELSE
    IPRED = CP                          ; plasma drug
    IF(IPRED.LE.0) IPRED = 1.0E-10
    Y     = IPRED + IPRED*EPS(1) + EPS(2)
  ENDIF
  IRES  = DV - IPRED
  IWRES = IRES/IPRED

; ---- residual error variances exported to $TABLE ----------------------
  SG1   = SIGMA(1,1)   ; proportional, analyte 1
  SG2   = SIGMA(2,2)   ; additive, analyte 1
  SG3   = SIGMA(3,3)   ; proportional, analyte 2
  SG4   = SIGMA(4,4)   ; additive, analyte 2

; ---- initial estimates -------------------------------------------------
$THETA  (0, 5.0)             ; 1 CL    clearance (L/h)
$THETA  (0, 50.0)            ; 2 V     central volume (L)
$THETA  (0, 1.0)             ; 3 KA    absorption rate constant (1/h)
$THETA  (0.001, 0.7, 0.999)  ; 4 F1    oral bioavailability (logit-scale IIV)
$THETA  (0, 0.5)             ; 5 KE0   effect-site equilibration rate (1/h)
$THETA  (0, 100.0)           ; 6 R0    baseline response (response units)
$THETA  (0, 0.1)             ; 7 KOUT  response turnover rate (1/h)
$THETA  (0, 2.0)             ; 8 SMAX  max. fractional stimulation
$THETA  (0, 1.0)             ; 9 SC50  conc. for 50% of SMAX (mg/L)
$THETA  (0.01, 1.0)          ; 10 HILL  Hill coefficient (fix to 1 for Emax)

$OMEGA  0.09             ; 1 IIV CL
$OMEGA  0.09             ; 2 IIV V
$OMEGA  0.09             ; 3 IIV KA
$OMEGA  0 FIX            ; 4 IIV F1    (free it if the data support it)
$OMEGA  0.09             ; 5 IIV KE0
$OMEGA  0.09             ; 6 IIV R0
$OMEGA  0.09             ; 7 IIV KOUT
$OMEGA  0 FIX            ; 8 IIV SMAX  (free it if the data support it)
$OMEGA  0.09             ; 9 IIV SC50
$OMEGA  0 FIX            ; 10 IIV HILL  (free it if the data support it)

$SIGMA  0.04             ; 1 proportional RUV, drug     (FLAG=1)
$SIGMA  0.01             ; 2 additive RUV, drug
$SIGMA  0.04             ; 3 proportional RUV, response (FLAG=2)
$SIGMA  1.0              ; 4 additive RUV, response

$ESTIMATION METHOD=1 INTER MAXEVAL=9999 NSIG=3 SIGL=9 PRINT=5 NOABORT POSTHOC
; A joint PK/PD fit is often easier to reach in two steps; the usual
; alternative to the FOCE-I line above is:
; $ESTIMATION METHOD=SAEM INTERACTION NBURN=2000 NITER=1000 PRINT=50
; $ESTIMATION METHOD=IMP INTERACTION EONLY=1 NITER=10 ISAMPLE=3000 PRINT=1
$COVARIANCE PRINT=E UNCONDITIONAL

$TABLE   ID TIME AMT RATE EVID MDV CMT FLAG DV IPRED IRES IWRES CWRES
         PRED RES WRES ET1 ET2 ET3 ET4 ET5 ET6 ET7 ET8 ET9 ET10 SG1 SG2
         SG3 SG4
         ONEHEADER NOPRINT FILE=idr4_1cmt_oral_effect.tab
$TABLE   ID TIME CP CE EFF FACT RESP
         ONEHEADER NOPRINT FILE=idr4_1cmt_oral_effect.prof
$TABLE   ID CL V KA F1 KE0 R0 KOUT SMAX SC50 HILL KIN ET1 ET2 ET3 ET4
         ET5 ET6 ET7 ET8 ET9 ET10
         FIRSTONLY ONEHEADER NOPRINT FILE=idr4_1cmt_oral_effect.par

;; ----------------------------------------------------------------------
;; SIMULATION BLOCK
;; To simulate instead of estimate: comment out $ESTIMATION and
;; $COVARIANCE above, then un-comment the two lines below.
;; ----------------------------------------------------------------------
; $SIMULATION (20260830) (20260831 NORMAL) ONLYSIM SUBPROBLEMS=200
; $TABLE ID TIME AMT EVID MDV CMT FLAG DV IPRED CP EFF RESP
;        NOAPPEND ONEHEADER NOPRINT FILE=idr4_1cmt_oral_effect_sim.tab

