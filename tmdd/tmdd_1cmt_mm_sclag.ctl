;; ==========================================================================
;; TMDD_1CMT_MM_SCLAG
;; Target-mediated drug disposition - 1-compartment model
;;   Form : Michaelis-Menten approximation
;;   Input: Subcutaneous, first-order absorption with lag time
;; --------------------------------------------------------------------------
;; Structure
;;   Dose --> [SC depot] --KA--> [Plasma (V)]
;;   [Plasma (V)] --CL--> eliminated
;;   target binding lumped into VM*C/(KM+C) from plasma
;;
;; Assumptions
;;   Target binding collapses to saturable elimination:
;;   VM*C/(KM + C), with VM ~ KINT*R0*V and KM ~ KSS.  Valid when the
;;   drug concentration substantially exceeds the target concentration
;;   (target occupancy near 100%).  The target is NOT modelled: Rtot is
;;   held at R0 and the target readouts are derived, not fitted.
;;
;; Data set / dosing requirements
;;   Dose records: CMT=1 (SC depot), RATE=0.
;;   Absorption starts at TIME+ALAG1.
;;   Observation records: CMT=2 for BOTH analytes; the FLAG column
;;   selects which one is being observed:
;;     FLAG=1  total drug concentration in plasma (nM)
;;     FLAG=2  total target concentration (nM)
;;
;; UNITS - read this first
;;   Drug and target concentrations MUST be in the SAME molar unit,
;;   because KON*C*R only makes sense that way.  These files assume
;;   nM for concentrations, nmol for AMT, L for volumes and days for
;;   time.  Feeding mg/L drug and nM target into a TMDD model is the
;;   classic silent failure - convert the dose to nmol first.
;;
;; Solver: general non-linear ODE (ADVAN13)
;; ==========================================================================

$PROBLEM TMDD 1cmt mm | sclag

$INPUT   ID TIME AMT RATE DV EVID MDV CMT FLAG ADDL II WT
$DATA    ../data/tmdddata.csv IGNORE=@

$SUBROUTINES ADVAN13 TOL=9
$MODEL  NCOMPARTMENTS=2
        COMP=(DEPOT, DEFDOSE)
        COMP=(CENTRAL, DEFOBS)

; compartment 1 = SC depot
; compartment 2 = central drug

$PK
; ---- typical values and between-subject variability ------------------
  CL    = THETA(1)*EXP(ETA(1))
  V     = THETA(2)*EXP(ETA(2))
  VM    = THETA(3)*EXP(ETA(3))
  KM    = THETA(4)*EXP(ETA(4))
  R0    = THETA(5)*EXP(ETA(5))
  KA    = THETA(6)*EXP(ETA(6))
  LGTF1 = LOG(THETA(7)/(1-THETA(7))) + ETA(7)
  F1    = EXP(LGTF1)/(1+EXP(LGTF1))
  ALAG1 = THETA(8)*EXP(ETA(8))

; ---- random effects exported to $TABLE --------------------------------
  ET1   = ETA(1)
  ET2   = ETA(2)
  ET3   = ETA(3)
  ET4   = ETA(4)
  ET5   = ETA(5)
  ET6   = ETA(6)
  ET7   = ETA(7)
  ET8   = ETA(8)

; ---- scaling (concentration = amount / volume) ------------------------
  S2 = V

$DES
  C    = A(2)/V
  IF(C.LT.0) C = 0
  DADT(1) = -KA*A(1)
  DADT(2) = KA*A(1) - CL*C - VM*C/(KM+C)

$ERROR
  CFR  = A(2)/V
  IF(CFR.LT.0) CFR = 0
  CTOT = CFR
  RCX  = R0*CFR/(KM+CFR)   ; derived, constant-target
  RTOT = R0
  RFR  = R0 - RCX
  IF(FLAG.EQ.2) THEN
    IPRED = RTOT                        ; total target
    IF(IPRED.LE.0) IPRED = 1.0E-10
    Y     = IPRED + IPRED*EPS(3) + EPS(4)
  ELSE
    IPRED = CTOT                        ; total drug
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
$THETA  (0, 0.2)             ; 1 CL    clearance (L/day)
$THETA  (0, 3.0)             ; 2 V     central volume (L)
$THETA  (0, 7.5)             ; 3 VM    max. target-mediated rate (nmol/day)
$THETA  (0, 1.0)             ; 4 KM    Michaelis constant (nM)
$THETA  (0, 5.0)             ; 5 R0    baseline target conc., readout only (nM)
$THETA  (0, 0.25)            ; 6 KA    absorption rate constant (1/day)
$THETA  (0.001, 0.7, 0.999)  ; 7 F1    SC bioavailability (logit-scale IIV)
$THETA  (0, 0.1)             ; 8 ALAG1 absorption lag time (day)

$OMEGA  0.09             ; 1 IIV CL
$OMEGA  0.09             ; 2 IIV V
$OMEGA  0.09             ; 3 IIV VM
$OMEGA  0 FIX            ; 4 IIV KM    (free it if the data support it)
$OMEGA  0.09             ; 5 IIV R0
$OMEGA  0.09             ; 6 IIV KA
$OMEGA  0 FIX            ; 7 IIV F1    (free it if the data support it)
$OMEGA  0 FIX            ; 8 IIV ALAG1 (free it if the data support it)

$SIGMA  0.04             ; 1 proportional RUV, drug   (FLAG=1)
$SIGMA  0.01             ; 2 additive RUV, drug
$SIGMA  0.04             ; 3 proportional RUV, target (FLAG=2)
$SIGMA  0.01             ; 4 additive RUV, target

$ESTIMATION METHOD=1 INTER MAXEVAL=9999 NSIG=3 SIGL=9 PRINT=5 NOABORT POSTHOC
; TMDD models often need a more robust route to the optimum; a common
; two-step alternative to the FOCE-I line above is:
; $ESTIMATION METHOD=SAEM INTERACTION NBURN=2000 NITER=1000 PRINT=50
; $ESTIMATION METHOD=IMP INTERACTION EONLY=1 NITER=10 ISAMPLE=3000 PRINT=1
$COVARIANCE PRINT=E UNCONDITIONAL

$TABLE   ID TIME AMT RATE EVID MDV CMT FLAG DV IPRED IRES IWRES CWRES
         PRED RES WRES ET1 ET2 ET3 ET4 ET5 ET6 ET7 ET8 SG1 SG2 SG3 SG4
         ONEHEADER NOPRINT FILE=tmdd_1cmt_mm_sclag.tab
$TABLE   ID TIME CFR CTOT RFR RCX RTOT
         ONEHEADER NOPRINT FILE=tmdd_1cmt_mm_sclag.prof
$TABLE   ID CL V VM KM R0 KA F1 ALAG1 ET1 ET2 ET3 ET4 ET5 ET6 ET7 ET8
         FIRSTONLY ONEHEADER NOPRINT FILE=tmdd_1cmt_mm_sclag.par

;; ----------------------------------------------------------------------
;; SIMULATION BLOCK
;; To simulate instead of estimate: comment out $ESTIMATION and
;; $COVARIANCE above, then un-comment the two lines below.
;; ----------------------------------------------------------------------
; $SIMULATION (20260830) (20260831 NORMAL) ONLYSIM SUBPROBLEMS=200
; $TABLE ID TIME AMT RATE EVID MDV CMT FLAG DV IPRED CFR CTOT RFR RCX
;        NOAPPEND ONEHEADER NOPRINT FILE=tmdd_1cmt_mm_sclag_sim.tab

