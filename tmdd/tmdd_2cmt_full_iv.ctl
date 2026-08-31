;; ==========================================================================
;; TMDD_2CMT_FULL_IV
;; Target-mediated drug disposition - 2-compartment model
;;   Form : Full TMDD
;;   Input: IV input (bolus or infusion)
;; --------------------------------------------------------------------------
;; Structure
;;   Dose --(IV)--> [Plasma (V2)]
;;   [Plasma (V2)] <--Q--> [Periph.1 (V3)]
;;   [Plasma (V2)] --CL--> eliminated
;;   drug + [Free target (R)] <--KON/KOFF--> [Complex (RC)]
;;   target: KSYN = KDEG*R0 in, KDEG out; complex out via KINT
;;
;; Assumptions
;;   Free drug, free target and complex are all explicit states.
;;   Binding is described by KON / KOFF; the complex is internalised
;;   with KINT.  Target turnover is zero-order synthesis and
;;   first-order degradation, with KSYN = KDEG*R0 so that the target
;;   sits at R0 before the first dose (Mager & Jusko 2001).
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

$PROBLEM TMDD 2cmt full | iv

$INPUT   ID TIME AMT RATE DV EVID MDV CMT FLAG ADDL II WT
$DATA    ../data/tmdddata.csv IGNORE=@

$SUBROUTINES ADVAN13 TOL=9
$MODEL  NCOMPARTMENTS=4
        COMP=(CENTRAL, DEFDOSE, DEFOBS)
        COMP=(PERIPH1)
        COMP=(RECEPT)
        COMP=(COMPLEX)

; compartment 1 = central drug
; compartment 2 = peripheral 1
; compartment 3 = free target
; compartment 4 = drug-target complex

$PK
; ---- typical values and between-subject variability ------------------
  CL    = THETA(1)*EXP(ETA(1))
  V2    = THETA(2)*EXP(ETA(2))
  Q     = THETA(3)*EXP(ETA(3))
  V3    = THETA(4)*EXP(ETA(4))
  KON   = THETA(5)*EXP(ETA(5))
  KOFF  = THETA(6)*EXP(ETA(6))
  KINT  = THETA(7)*EXP(ETA(7))
  KDEG  = THETA(8)*EXP(ETA(8))
  R0    = THETA(9)*EXP(ETA(9))

; ---- pre-dose target baseline ----------------------------------------
;     KSYN = KDEG*R0, so the target starts at steady state R0.
  A_0(3) = R0*V2

; ---- scaling (concentration = amount / volume) ------------------------
  S1 = V2

$DES
  C    = A(1)/V2
  IF(C.LT.0) C = 0
  CP1  = A(2)/V3
  RF   = A(3)/V2
  RC   = A(4)/V2
  BIND = KON*C*RF*V2 - KOFF*RC*V2
  DIST1= Q*(C - CP1)
  DADT(1) = -CL*C - DIST1 - BIND
  DADT(2) =  DIST1
  DADT(3) =  KDEG*R0*V2 - KDEG*RF*V2 - BIND
  DADT(4) =  BIND - KINT*RC*V2

$ERROR
  CFR  = A(1)/V2
  RFR  = A(3)/V2
  RCX  = A(4)/V2
  CTOT = CFR + RCX
  RTOT = RFR + RCX
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
$THETA  (0, 0.5)             ; 3 Q     inter-cmt clearance, periph.1 (L/day)
$THETA  (0, 3.0)             ; 4 V3    peripheral volume 1 (L)
$THETA  (0, 1.0)             ; 5 KON   association rate (1/(nM*day))
$THETA  (0, 0.5)             ; 6 KOFF  dissociation rate (1/day)
$THETA  (0, 0.5)             ; 7 KINT  complex internalisation rate (1/day)
$THETA  (0, 1.0)             ; 8 KDEG  free target degradation rate (1/day)
$THETA  (0, 5.0)             ; 9 R0    baseline total target conc. (nM)

$OMEGA  0.09             ; 1 IIV CL
$OMEGA  0.09             ; 2 IIV V2
$OMEGA  0.09             ; 3 IIV Q
$OMEGA  0.09             ; 4 IIV V3
$OMEGA  0 FIX            ; 5 IIV KON   (free it if the data support it)
$OMEGA  0 FIX            ; 6 IIV KOFF  (free it if the data support it)
$OMEGA  0.09             ; 7 IIV KINT
$OMEGA  0.09             ; 8 IIV KDEG
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
         ONEHEADER NOPRINT FILE=tmdd_2cmt_full_iv.tab
$TABLE   ID TIME CFR CTOT RFR RCX RTOT
         ONEHEADER NOPRINT FILE=tmdd_2cmt_full_iv.prof
$TABLE   ID CL V2 Q V3 KON KOFF KINT KDEG R0 ETA1 ETA2 ETA3 ETA4 ETA5
         ETA6 ETA7 ETA8 ETA9
         FIRSTONLY ONEHEADER NOPRINT FILE=tmdd_2cmt_full_iv.par

;; ----------------------------------------------------------------------
;; SIMULATION BLOCK
;; To simulate instead of estimate: comment out $ESTIMATION and
;; $COVARIANCE above, then un-comment the two lines below.
;; ----------------------------------------------------------------------
; $SIMULATION (20260830) (20260831 NORMAL) ONLYSIM SUBPROBLEMS=200
; $TABLE ID TIME AMT RATE EVID MDV CMT FLAG DV IPRED CFR CTOT RFR RCX
;        NOAPPEND ONEHEADER NOPRINT FILE=tmdd_2cmt_full_iv.tab

