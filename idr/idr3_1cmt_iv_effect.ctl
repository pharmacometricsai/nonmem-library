;; ==========================================================================
;; IDR3_1CMT_IV_EFFECT
;; Indirect response model, simultaneous PK/PD
;;   PD   : IDR model III - Stimulation of production
;;   PK   : 1-compartment, IV input (bolus or infusion)
;;   Link : Effect compartment (KE0), for counter-clockwise hysteresis
;; --------------------------------------------------------------------------
;; Structure
;;   Dose --(IV)--> [Central (V)]
;;   [Central (V)] --CL--> eliminated
;;   plasma --KE0--> [Effect site (Ce)] --> drives the response
;;   KIN --> [Response (R)] --KOUT-->
;;   dR/dt = KIN*(1 + SMAX*EFF) - KOUT*R
;;   EFF = (CD/SC50)**HILL / (1 + (CD/SC50)**HILL)
;;   where CD = effect-site concentration Ce
;;
;; Behaviour
;;   The drug stimulates the zero-order production of the response.
;;   Response rises above baseline; the return to baseline is governed
;;   by KOUT alone, so its rate is dose-independent.
;;   Net direction of the response: increase from baseline.
;;
;; Baseline
;;   KIN is NOT a separate parameter: KIN = R0*KOUT, and the response
;;   compartment is initialised at A_0 = R0.  The response therefore
;;   starts at, and returns to, R0 by construction.  Do not add a
;;   separate baseline THETA on top of this.
;;
;; Data set / dosing requirements
;;   Dose records: CMT=1 (central).
;;   RATE=0 gives an IV bolus; RATE>0 an IV infusion.
;;   Observation records: CMT=1 for BOTH analytes; the FLAG column
;;   selects which one is being observed:
;;     FLAG=1  plasma drug concentration (mg/L)
;;     FLAG=2  response R (response units)
;;
;; Units: mg, L, h - the same as ../1cmt, ../2cmt and ../3cmt.
;;   (../tmdd uses nM and days; do not share a data set with it.)
;;
;; Solver: general non-linear ODE (ADVAN13)
;; ==========================================================================

$PROBLEM IDR3 1cmt iv effect

$INPUT   ID TIME AMT RATE DV EVID MDV CMT FLAG ADDL II WT
$DATA    ../data/pkpddata.csv IGNORE=@

$SUBROUTINES ADVAN13 TOL=9
$MODEL  NCOMPARTMENTS=3
        COMP=(CENTRAL, DEFDOSE, DEFOBS)
        COMP=(EFFECT)
        COMP=(RESPONSE)

; compartment 1 = central drug
; compartment 2 = effect site (holds Ce directly)
; compartment 3 = response R (holds R directly)
; The effect and response compartments are unscaled: A() holds the
; effect-site concentration and the response itself, not an amount.

$PK
; ---- typical values and between-subject variability ------------------
  CL    = THETA(1)*EXP(ETA(1))
  V     = THETA(2)*EXP(ETA(2))
  KE0   = THETA(3)*EXP(ETA(3))
  R0    = THETA(4)*EXP(ETA(4))
  KOUT  = THETA(5)*EXP(ETA(5))
  SMAX  = THETA(6)*EXP(ETA(6))
  SC50  = THETA(7)*EXP(ETA(7))
  HILL  = THETA(8)*EXP(ETA(8))

; ---- random effects exported to $TABLE --------------------------------
  ET1   = ETA(1)
  ET2   = ETA(2)
  ET3   = ETA(3)
  ET4   = ETA(4)
  ET5   = ETA(5)
  ET6   = ETA(6)
  ET7   = ETA(7)
  ET8   = ETA(8)

; ---- response baseline -----------------------------------------------
;     KIN = R0*KOUT, so the response sits at R0 before the first dose.
  KIN   = R0*KOUT
  A_0(3) = R0

; ---- scaling (concentration = amount / volume) ------------------------
  S1 = V

$DES
  CP   = A(1)/V
  IF(CP.LT.0) CP = 0
  DADT(1) = -CL*CP
  DADT(2) = KE0*(CP - A(2))
  CD   = A(2)
  EFF  = 0
  IF(CD.GT.0) THEN
    XX  = (CD/SC50)**HILL
    EFF = XX/(1+XX)
  ENDIF
  DADT(3) = KIN*(1 + SMAX*EFF) - KOUT*A(3)

$ERROR
  CP   = A(1)/V
  IF(CP.LT.0) CP = 0
  CE   = A(2)
  CD   = A(2)
  EFF  = 0
  IF(CD.GT.0) THEN
    XX  = (CD/SC50)**HILL
    EFF = XX/(1+XX)
  ENDIF
  FACT = 1 + SMAX*EFF        ; factor applied to KIN
  RESP = A(3)
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
$THETA  (0, 0.5)             ; 3 KE0   effect-site equilibration rate (1/h)
$THETA  (0, 100.0)           ; 4 R0    baseline response (response units)
$THETA  (0, 0.1)             ; 5 KOUT  response turnover rate (1/h)
$THETA  (0, 2.0)             ; 6 SMAX  max. fractional stimulation
$THETA  (0, 1.0)             ; 7 SC50  conc. for 50% of SMAX (mg/L)
$THETA  (0.01, 1.0)          ; 8 HILL  Hill coefficient (fix to 1 for Emax)

$OMEGA  0.09             ; 1 IIV CL
$OMEGA  0.09             ; 2 IIV V
$OMEGA  0.09             ; 3 IIV KE0
$OMEGA  0.09             ; 4 IIV R0
$OMEGA  0.09             ; 5 IIV KOUT
$OMEGA  0 FIX            ; 6 IIV SMAX  (free it if the data support it)
$OMEGA  0.09             ; 7 IIV SC50
$OMEGA  0 FIX            ; 8 IIV HILL  (free it if the data support it)

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
         PRED RES WRES ET1 ET2 ET3 ET4 ET5 ET6 ET7 ET8 SG1 SG2 SG3 SG4
         ONEHEADER NOPRINT FILE=idr3_1cmt_iv_effect.tab
$TABLE   ID TIME CP CE EFF FACT RESP
         ONEHEADER NOPRINT FILE=idr3_1cmt_iv_effect.prof
$TABLE   ID CL V KE0 R0 KOUT SMAX SC50 HILL KIN ET1 ET2 ET3 ET4 ET5 ET6
         ET7 ET8
         FIRSTONLY ONEHEADER NOPRINT FILE=idr3_1cmt_iv_effect.par

;; ----------------------------------------------------------------------
;; SIMULATION BLOCK
;; To simulate instead of estimate: comment out $ESTIMATION and
;; $COVARIANCE above, then un-comment the two lines below.
;; ----------------------------------------------------------------------
; $SIMULATION (20260830) (20260831 NORMAL) ONLYSIM SUBPROBLEMS=200
; $TABLE ID TIME AMT EVID MDV CMT FLAG DV IPRED CP EFF RESP
;        NOAPPEND ONEHEADER NOPRINT FILE=idr3_1cmt_iv_effect_sim.tab

