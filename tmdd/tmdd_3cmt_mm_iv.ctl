;; ==========================================================================
;; TMDD_3CMT_MM_IV
;; Target-mediated drug disposition - 3-compartment model
;;   Form : Michaelis-Menten approximation
;;   Input: IV input (bolus or infusion)
;; --------------------------------------------------------------------------
;; Structure
;;   Dose --(IV)--> [Plasma (V2)]
;;   [Plasma (V2)] <--Q3--> [Periph.1 (V3)]
;;   [Plasma (V2)] <--Q4--> [Periph.2 (V4)]
;;   [Plasma (V2)] --CL--> eliminated
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
;;   Dose records: CMT=1 (plasma).
;;   RATE=0 gives an IV bolus; RATE>0 an IV infusion.
;;   Observation records: CMT=1 for BOTH analytes; the FLAG column
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

$PROBLEM TMDD 3cmt mm | iv

$INPUT   ID TIME AMT RATE DV EVID MDV CMT FLAG ADDL II WT
$DATA    ../data/tmdddata.csv IGNORE=@

$SUBROUTINES ADVAN13 TOL=9
$MODEL  NCOMPARTMENTS=3
        COMP=(CENTRAL, DEFDOSE, DEFOBS)
        COMP=(PERIPH1)
        COMP=(PERIPH2)

; compartment 1 = central drug
; compartment 2 = peripheral 1
; compartment 3 = peripheral 2

$PK
; ---- typical values and between-subject variability ------------------
  CL    = THETA(1)*EXP(ETA(1))
  V2    = THETA(2)*EXP(ETA(2))
  Q3    = THETA(3)*EXP(ETA(3))
  V3    = THETA(4)*EXP(ETA(4))
  Q4    = THETA(5)*EXP(ETA(5))
  V4    = THETA(6)*EXP(ETA(6))
  VM    = THETA(7)*EXP(ETA(7))
  KM    = THETA(8)*EXP(ETA(8))
  R0    = THETA(9)*EXP(ETA(9))

; ---- scaling (concentration = amount / volume) ------------------------
  S1 = V2

$DES
  C    = A(1)/V2
  IF(C.LT.0) C = 0
  CP1  = A(2)/V3
  CP2  = A(3)/V4
  DIST1= Q3*(C - CP1)
  DIST2= Q4*(C - CP2)
  DADT(1) = -CL*C - DIST1 - DIST2 - VM*C/(KM+C)
  DADT(2) =  DIST1
  DADT(3) =  DIST2

$ERROR
  CFR  = A(1)/V2
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

; ---- initial estimates -------------------------------------------------
$THETA  (0, 0.2)             ; 1 CL    clearance (L/day)
$THETA  (0, 3.0)             ; 2 V2    central volume (L)
$THETA  (0, 0.5)             ; 3 Q3    inter-cmt clearance, periph.1 (L/day)
$THETA  (0, 3.0)             ; 4 V3    peripheral volume 1 (L)
$THETA  (0, 0.2)             ; 5 Q4    inter-cmt clearance, periph.2 (L/day)
$THETA  (0, 5.0)             ; 6 V4    peripheral volume 2 (L)
$THETA  (0, 7.5)             ; 7 VM    max. target-mediated rate (nmol/day)
$THETA  (0, 1.0)             ; 8 KM    Michaelis constant (nM)
$THETA  (0, 5.0)             ; 9 R0    baseline target conc., readout only (nM)

$OMEGA  0.09             ; 1 IIV CL
$OMEGA  0.09             ; 2 IIV V2
$OMEGA  0.09             ; 3 IIV Q3
$OMEGA  0.09             ; 4 IIV V3
$OMEGA  0.09             ; 5 IIV Q4
$OMEGA  0.09             ; 6 IIV V4
$OMEGA  0.09             ; 7 IIV VM
$OMEGA  0 FIX            ; 8 IIV KM    (free it if the data support it)
$OMEGA  0.09             ; 9 IIV R0

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

$TABLE   ID TIME AMT RATE EVID MDV CMT FLAG DV IPRED IRES IWRES
         CWRES PRED RES WRES
         ONEHEADER NOPRINT FILE=tmdd_3cmt_mm_iv.tab
$TABLE   ID TIME CFR CTOT RFR RCX RTOT
         ONEHEADER NOPRINT FILE=tmdd_3cmt_mm_iv.prof
$TABLE   ID CL V2 Q3 V3 Q4 V4 VM KM R0 ETA1 ETA2 ETA3 ETA4 ETA5 ETA6
         ETA7 ETA8 ETA9
         FIRSTONLY ONEHEADER NOPRINT FILE=tmdd_3cmt_mm_iv.par

;; ----------------------------------------------------------------------
;; SIMULATION BLOCK
;; To simulate instead of estimate: comment out $ESTIMATION and
;; $COVARIANCE above, then un-comment the two lines below.
;; ----------------------------------------------------------------------
; $SIMULATION (20260830) (20260831 NORMAL) ONLYSIM SUBPROBLEMS=200
; $TABLE ID TIME AMT RATE EVID MDV CMT FLAG DV IPRED CFR CTOT RFR RCX
;        NOAPPEND ONEHEADER NOPRINT FILE=tmdd_3cmt_mm_iv.tab

