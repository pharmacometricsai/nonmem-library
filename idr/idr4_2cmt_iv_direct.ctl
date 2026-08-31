;; ==========================================================================
;; IDR4_2CMT_IV_DIRECT
;; Indirect response model, simultaneous PK/PD
;;   PD   : IDR model IV - Stimulation of elimination
;;   PK   : 2-compartment, IV input (bolus or infusion)
;;   Link : Direct link (plasma concentration drives the effect)
;; --------------------------------------------------------------------------
;; Structure
;;   Dose --(IV)--> [Central (V2)]
;;   [Central (V2)] <--Q--> [Periph.1 (V3)]
;;   [Central (V2)] --CL--> eliminated
;;   plasma concentration drives the response directly
;;   KIN --> [Response (R)] --KOUT-->
;;   dR/dt = KIN - KOUT*(1 + SMAX*EFF)*R
;;   EFF = (CD/SC50)**HILL / (1 + (CD/SC50)**HILL)
;;   where CD = plasma concentration
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

$PROBLEM IDR4 2cmt iv direct

$INPUT   ID TIME AMT RATE DV EVID MDV CMT FLAG ADDL II WT
$DATA    ../data/pkpddata.csv IGNORE=@

$SUBROUTINES ADVAN13 TOL=9
$MODEL  NCOMPARTMENTS=3
        COMP=(CENTRAL, DEFDOSE, DEFOBS)
        COMP=(PERIPH1)
        COMP=(RESPONSE)

; compartment 1 = central drug
; compartment 2 = peripheral 1
; compartment 3 = response R (holds R directly)
; The effect and response compartments are unscaled: A() holds the
; effect-site concentration and the response itself, not an amount.

$PK
; ---- typical values and between-subject variability ------------------
  CL    = THETA(1)*EXP(ETA(1))
  V2    = THETA(2)*EXP(ETA(2))
  Q     = THETA(3)*EXP(ETA(3))
  V3    = THETA(4)*EXP(ETA(4))
  R0    = THETA(5)*EXP(ETA(5))
  KOUT  = THETA(6)*EXP(ETA(6))
  SMAX  = THETA(7)*EXP(ETA(7))
  SC50  = THETA(8)*EXP(ETA(8))
  HILL  = THETA(9)*EXP(ETA(9))

; ---- response baseline -----------------------------------------------
;     KIN = R0*KOUT, so the response sits at R0 before the first dose.
  KIN   = R0*KOUT
  A_0(3) = R0

; ---- scaling (concentration = amount / volume) ------------------------
  S1 = V2

$DES
  CP   = A(1)/V2
  IF(CP.LT.0) CP = 0
  CP1  = A(2)/V3
  DIST1= Q*(CP - CP1)
  DADT(1) = -CL*CP - DIST1
  DADT(2) =  DIST1
  CD   = CP
  EFF  = 0
  IF(CD.GT.0) THEN
    XX  = (CD/SC50)**HILL
    EFF = XX/(1+XX)
  ENDIF
  DADT(3) = KIN - KOUT*(1 + SMAX*EFF)*A(3)

$ERROR
  CP   = A(1)/V2
  IF(CP.LT.0) CP = 0
  CD   = CP
  EFF  = 0
  IF(CD.GT.0) THEN
    XX  = (CD/SC50)**HILL
    EFF = XX/(1+XX)
  ENDIF
  FACT = 1 + SMAX*EFF        ; factor applied to KOUT
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

; ---- initial estimates -------------------------------------------------
$THETA  (0, 5.0)             ; 1 CL    clearance (L/h)
$THETA  (0, 50.0)            ; 2 V2    central volume (L)
$THETA  (0, 10.0)            ; 3 Q     inter-cmt clearance, periph.1 (L/h)
$THETA  (0, 100.0)           ; 4 V3    peripheral volume 1 (L)
$THETA  (0, 100.0)           ; 5 R0    baseline response (response units)
$THETA  (0, 0.1)             ; 6 KOUT  response turnover rate (1/h)
$THETA  (0, 2.0)             ; 7 SMAX  max. fractional stimulation
$THETA  (0, 1.0)             ; 8 SC50  conc. for 50% of SMAX (mg/L)
$THETA  (0.01, 1.0)          ; 9 HILL  Hill coefficient (fix to 1 for Emax)

$OMEGA  0.09             ; 1 IIV CL
$OMEGA  0.09             ; 2 IIV V2
$OMEGA  0.09             ; 3 IIV Q
$OMEGA  0.09             ; 4 IIV V3
$OMEGA  0.09             ; 5 IIV R0
$OMEGA  0.09             ; 6 IIV KOUT
$OMEGA  0 FIX            ; 7 IIV SMAX  (free it if the data support it)
$OMEGA  0.09             ; 8 IIV SC50
$OMEGA  0 FIX            ; 9 IIV HILL  (free it if the data support it)

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

$TABLE   ID TIME AMT RATE EVID MDV CMT FLAG DV IPRED IRES IWRES
         CWRES PRED RES WRES
         ONEHEADER NOPRINT FILE=idr4_2cmt_iv_direct.tab
$TABLE   ID TIME CP EFF FACT RESP
         ONEHEADER NOPRINT FILE=idr4_2cmt_iv_direct.prof
$TABLE   ID CL V2 Q V3 R0 KOUT SMAX SC50 HILL KIN ETA1 ETA2 ETA3 ETA4
         ETA5 ETA6 ETA7 ETA8 ETA9
         FIRSTONLY ONEHEADER NOPRINT FILE=idr4_2cmt_iv_direct.par

;; ----------------------------------------------------------------------
;; SIMULATION BLOCK
;; To simulate instead of estimate: comment out $ESTIMATION and
;; $COVARIANCE above, then un-comment the two lines below.
;; ----------------------------------------------------------------------
; $SIMULATION (20260830) (20260831 NORMAL) ONLYSIM SUBPROBLEMS=200
; $TABLE ID TIME AMT EVID MDV CMT FLAG DV IPRED CP EFF RESP
;        NOAPPEND ONEHEADER NOPRINT FILE=idr4_2cmt_iv_direct.tab

