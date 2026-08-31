;; ==========================================================================
;; TMDD_3CMT_FULL_SCLAG
;; Target-mediated drug disposition - 3-compartment model
;;   Form : Full TMDD
;;   Input: Subcutaneous, first-order absorption with lag time
;; --------------------------------------------------------------------------
;; Structure
;;   Dose --> [SC depot] --KA--> [Plasma (V2)]
;;   [Plasma (V2)] <--Q3--> [Periph.1 (V3)]
;;   [Plasma (V2)] <--Q4--> [Periph.2 (V4)]
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

$PROBLEM TMDD 3cmt full | sclag

$INPUT   ID TIME AMT RATE DV EVID MDV CMT FLAG ADDL II WT
$DATA    ../data/tmdddata.csv IGNORE=@

$SUBROUTINES ADVAN13 TOL=9
$MODEL  NCOMPARTMENTS=6
        COMP=(DEPOT, DEFDOSE)
        COMP=(CENTRAL, DEFOBS)
        COMP=(PERIPH1)
        COMP=(PERIPH2)
        COMP=(RECEPT)
        COMP=(COMPLEX)

; compartment 1 = SC depot
; compartment 2 = central drug
; compartment 3 = peripheral 1
; compartment 4 = peripheral 2
; compartment 5 = free target
; compartment 6 = drug-target complex

$PK
; ---- typical values and between-subject variability ------------------
  CL    = THETA(1)*EXP(ETA(1))
  V2    = THETA(2)*EXP(ETA(2))
  Q3    = THETA(3)*EXP(ETA(3))
  V3    = THETA(4)*EXP(ETA(4))
  Q4    = THETA(5)*EXP(ETA(5))
  V4    = THETA(6)*EXP(ETA(6))
  KON   = THETA(7)*EXP(ETA(7))
  KOFF  = THETA(8)*EXP(ETA(8))
  KINT  = THETA(9)*EXP(ETA(9))
  KDEG  = THETA(10)*EXP(ETA(10))
  R0    = THETA(11)*EXP(ETA(11))
  KA    = THETA(12)*EXP(ETA(12))
  LGTF1 = LOG(THETA(13)/(1-THETA(13))) + ETA(13)
  F1    = EXP(LGTF1)/(1+EXP(LGTF1))
  ALAG1 = THETA(14)*EXP(ETA(14))

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
  ET11  = ETA(11)
  ET12  = ETA(12)
  ET13  = ETA(13)
  ET14  = ETA(14)

; ---- pre-dose target baseline ----------------------------------------
;     KSYN = KDEG*R0, so the target starts at steady state R0.
  A_0(5) = R0*V2

; ---- scaling (concentration = amount / volume) ------------------------
  S2 = V2

$DES
  C    = A(2)/V2
  IF(C.LT.0) C = 0
  CP1  = A(3)/V3
  CP2  = A(4)/V4
  RF   = A(5)/V2
  RC   = A(6)/V2
  BIND = KON*C*RF*V2 - KOFF*RC*V2
  DIST1= Q3*(C - CP1)
  DIST2= Q4*(C - CP2)
  DADT(1) = -KA*A(1)
  DADT(2) = KA*A(1) - CL*C - DIST1 - DIST2 - BIND
  DADT(3) =  DIST1
  DADT(4) =  DIST2
  DADT(5) =  KDEG*R0*V2 - KDEG*RF*V2 - BIND
  DADT(6) =  BIND - KINT*RC*V2

$ERROR
  CFR  = A(2)/V2
  RFR  = A(5)/V2
  RCX  = A(6)/V2
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

; ---- residual error variances exported to $TABLE ----------------------
  SG1   = SIGMA(1,1)   ; proportional, analyte 1
  SG2   = SIGMA(2,2)   ; additive, analyte 1
  SG3   = SIGMA(3,3)   ; proportional, analyte 2
  SG4   = SIGMA(4,4)   ; additive, analyte 2

; ---- initial estimates -------------------------------------------------
$THETA  (0, 0.2)             ; 1 CL    clearance (L/day)
$THETA  (0, 3.0)             ; 2 V2    central volume (L)
$THETA  (0, 0.5)             ; 3 Q3    inter-cmt clearance, periph.1 (L/day)
$THETA  (0, 3.0)             ; 4 V3    peripheral volume 1 (L)
$THETA  (0, 0.2)             ; 5 Q4    inter-cmt clearance, periph.2 (L/day)
$THETA  (0, 5.0)             ; 6 V4    peripheral volume 2 (L)
$THETA  (0, 1.0)             ; 7 KON   association rate (1/(nM*day))
$THETA  (0, 0.5)             ; 8 KOFF  dissociation rate (1/day)
$THETA  (0, 0.5)             ; 9 KINT  complex internalisation rate (1/day)
$THETA  (0, 1.0)             ; 10 KDEG  free target degradation rate (1/day)
$THETA  (0, 5.0)             ; 11 R0    baseline total target conc. (nM)
$THETA  (0, 0.25)            ; 12 KA    absorption rate constant (1/day)
$THETA  (0.001, 0.7, 0.999)  ; 13 F1    SC bioavailability (logit-scale IIV)
$THETA  (0, 0.1)             ; 14 ALAG1 absorption lag time (day)

$OMEGA  0.09             ; 1 IIV CL
$OMEGA  0.09             ; 2 IIV V2
$OMEGA  0.09             ; 3 IIV Q3
$OMEGA  0.09             ; 4 IIV V3
$OMEGA  0.09             ; 5 IIV Q4
$OMEGA  0.09             ; 6 IIV V4
$OMEGA  0 FIX            ; 7 IIV KON   (free it if the data support it)
$OMEGA  0 FIX            ; 8 IIV KOFF  (free it if the data support it)
$OMEGA  0.09             ; 9 IIV KINT
$OMEGA  0.09             ; 10 IIV KDEG
$OMEGA  0.09             ; 11 IIV R0
$OMEGA  0.09             ; 12 IIV KA
$OMEGA  0 FIX            ; 13 IIV F1    (free it if the data support it)
$OMEGA  0 FIX            ; 14 IIV ALAG1 (free it if the data support it)

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
         PRED RES WRES ET1 ET2 ET3 ET4 ET5 ET6 ET7 ET8 ET9 ET10 ET11
         ET12 ET13 ET14 SG1 SG2 SG3 SG4
         ONEHEADER NOPRINT FILE=tmdd_3cmt_full_sclag.tab
$TABLE   ID TIME CFR CTOT RFR RCX RTOT
         ONEHEADER NOPRINT FILE=tmdd_3cmt_full_sclag.prof
$TABLE   ID CL V2 Q3 V3 Q4 V4 KON KOFF KINT KDEG R0 KA F1 ALAG1 ET1 ET2
         ET3 ET4 ET5 ET6 ET7 ET8 ET9 ET10 ET11 ET12 ET13 ET14
         FIRSTONLY ONEHEADER NOPRINT FILE=tmdd_3cmt_full_sclag.par

;; ----------------------------------------------------------------------
;; SIMULATION BLOCK
;; To simulate instead of estimate: comment out $ESTIMATION and
;; $COVARIANCE above, then un-comment the two lines below.
;; ----------------------------------------------------------------------
; $SIMULATION (20260830) (20260831 NORMAL) ONLYSIM SUBPROBLEMS=200
; $TABLE ID TIME AMT RATE EVID MDV CMT FLAG DV IPRED CFR CTOT RFR RCX
;        NOAPPEND ONEHEADER NOPRINT FILE=tmdd_3cmt_full_sclag_sim.tab

