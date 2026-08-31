;; ==========================================================================
;; TMDD_2CMT_WAGNER_SC
;; Target-mediated drug disposition - 2-compartment model
;;   Form : Wagner (quasi-equilibrium, constant total target)
;;   Input: Subcutaneous, first-order absorption
;; --------------------------------------------------------------------------
;; Structure
;;   Dose --> [SC depot] --KA--> [Plasma (V2)]
;;   [Plasma (V2)] <--Q--> [Periph.1 (V3)]
;;   [Plasma (V2)] --CL--> eliminated
;;   drug + target in rapid equilibrium (KD); total target held
;;   at R0; complex cleared by KINT
;;
;; Assumptions
;;   Quasi-equilibrium binding PLUS the assumption that the TOTAL
;;   target concentration stays at its baseline R0, so there is no
;;   target-turnover state at all.  KSYN and KDEG drop out; the only
;;   target parameters are R0 and KD.  Equivalent to the constant-Rtot
;;   form when KDEG = KINT (Wagner 1973).
;;
;; Data set / dosing requirements
;;   Dose records: CMT=1 (SC depot), RATE=0.
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

$PROBLEM TMDD 2cmt wagner | sc

$INPUT   ID TIME AMT RATE DV EVID MDV CMT FLAG ADDL II WT
$DATA    ../data/tmdddata.csv IGNORE=@

$SUBROUTINES ADVAN13 TOL=9
$MODEL  NCOMPARTMENTS=3
        COMP=(DEPOT, DEFDOSE)
        COMP=(CENTRAL, DEFOBS)
        COMP=(PERIPH1)

; compartment 1 = SC depot
; compartment 2 = central drug
; compartment 3 = peripheral 1

$PK
; ---- typical values and between-subject variability ------------------
  CL    = THETA(1)*EXP(ETA(1))
  V2    = THETA(2)*EXP(ETA(2))
  Q     = THETA(3)*EXP(ETA(3))
  V3    = THETA(4)*EXP(ETA(4))
  KD    = THETA(5)*EXP(ETA(5))
  KINT  = THETA(6)*EXP(ETA(6))
  R0    = THETA(7)*EXP(ETA(7))
  KA    = THETA(8)*EXP(ETA(8))
  LGTF1 = LOG(THETA(9)/(1-THETA(9))) + ETA(9)
  F1    = EXP(LGTF1)/(1+EXP(LGTF1))

; ---- scaling (concentration = amount / volume) ------------------------
  S2 = V2

$DES
  CTOT = A(2)/V2
  IF(CTOT.LT.0) CTOT = 0
  CP1  = A(3)/V3
  BB   = CTOT - R0 - KD
  CF   = 0.5*(BB + SQRT(BB*BB + 4*KD*CTOT))
  RC   = CTOT - CF
  DIST1= Q*(CF - CP1)
  DADT(1) = -KA*A(1)
  DADT(2) = KA*A(1) - CL*CF - DIST1 - KINT*RC*V2
  DADT(3) =  DIST1

$ERROR
  CTOT = A(2)/V2
  IF(CTOT.LT.0) CTOT = 0
  BB   = CTOT - R0 - KD
  CFR  = 0.5*(BB + SQRT(BB*BB + 4*KD*CTOT))
  RCX  = CTOT - CFR
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
$THETA  (0, 0.5)             ; 3 Q     inter-cmt clearance, periph.1 (L/day)
$THETA  (0, 3.0)             ; 4 V3    peripheral volume 1 (L)
$THETA  (0, 0.5)             ; 5 KD    equilibrium dissoc. constant (nM)
$THETA  (0, 0.5)             ; 6 KINT  complex internalisation rate (1/day)
$THETA  (0, 5.0)             ; 7 R0    baseline total target conc. (nM)
$THETA  (0, 0.25)            ; 8 KA    absorption rate constant (1/day)
$THETA  (0.001, 0.7, 0.999)  ; 9 F1    SC bioavailability (logit-scale IIV)

$OMEGA  0.09             ; 1 IIV CL
$OMEGA  0.09             ; 2 IIV V2
$OMEGA  0.09             ; 3 IIV Q
$OMEGA  0.09             ; 4 IIV V3
$OMEGA  0 FIX            ; 5 IIV KD    (free it if the data support it)
$OMEGA  0.09             ; 6 IIV KINT
$OMEGA  0.09             ; 7 IIV R0
$OMEGA  0.09             ; 8 IIV KA
$OMEGA  0 FIX            ; 9 IIV F1    (free it if the data support it)

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
         ONEHEADER NOPRINT FILE=tmdd_2cmt_wagner_sc.tab
$TABLE   ID TIME CFR CTOT RFR RCX RTOT
         ONEHEADER NOPRINT FILE=tmdd_2cmt_wagner_sc.prof
$TABLE   ID CL V2 Q V3 KD KINT R0 KA F1 ETA1 ETA2 ETA3 ETA4 ETA5 ETA6
         ETA7 ETA8 ETA9
         FIRSTONLY ONEHEADER NOPRINT FILE=tmdd_2cmt_wagner_sc.par

;; ----------------------------------------------------------------------
;; SIMULATION BLOCK
;; To simulate instead of estimate: comment out $ESTIMATION and
;; $COVARIANCE above, then un-comment the two lines below.
;; ----------------------------------------------------------------------
; $SIMULATION (20260830) (20260831 NORMAL) ONLYSIM SUBPROBLEMS=200
; $TABLE ID TIME AMT RATE EVID MDV CMT FLAG DV IPRED CFR CTOT RFR RCX
;        NOAPPEND ONEHEADER NOPRINT FILE=tmdd_2cmt_wagner_sc.tab

