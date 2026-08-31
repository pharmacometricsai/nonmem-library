;; ==========================================================================
;; TMDD_1CMT_WAGNER_IV
;; Target-mediated drug disposition - 1-compartment model
;;   Form : Wagner (quasi-equilibrium, constant total target)
;;   Input: IV input (bolus or infusion)
;; --------------------------------------------------------------------------
;; Structure
;;   Dose --(IV)--> [Plasma (V)]
;;   [Plasma (V)] --CL--> eliminated
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

$PROBLEM TMDD 1cmt wagner | iv

$INPUT   ID TIME AMT RATE DV EVID MDV CMT FLAG ADDL II WT
$DATA    ../data/tmdddata.csv IGNORE=@

$SUBROUTINES ADVAN13 TOL=9
$MODEL  NCOMPARTMENTS=1
        COMP=(CENTRAL, DEFDOSE, DEFOBS)

; compartment 1 = central drug

$PK
; ---- typical values and between-subject variability ------------------
  CL    = THETA(1)*EXP(ETA(1))
  V     = THETA(2)*EXP(ETA(2))
  KD    = THETA(3)*EXP(ETA(3))
  KINT  = THETA(4)*EXP(ETA(4))
  R0    = THETA(5)*EXP(ETA(5))

; ---- random effects exported to $TABLE --------------------------------
  ET1   = ETA(1)
  ET2   = ETA(2)
  ET3   = ETA(3)
  ET4   = ETA(4)
  ET5   = ETA(5)

; ---- scaling (concentration = amount / volume) ------------------------
  S1 = V

$DES
  CTOT = A(1)/V
  IF(CTOT.LT.0) CTOT = 0
  BB   = CTOT - R0 - KD
  CF   = 0.5*(BB + SQRT(BB*BB + 4*KD*CTOT))
  RC   = CTOT - CF
  DADT(1) = -CL*CF - KINT*RC*V

$ERROR
  CTOT = A(1)/V
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

; ---- residual error variances exported to $TABLE ----------------------
  SG1   = SIGMA(1,1)   ; proportional, analyte 1
  SG2   = SIGMA(2,2)   ; additive, analyte 1
  SG3   = SIGMA(3,3)   ; proportional, analyte 2
  SG4   = SIGMA(4,4)   ; additive, analyte 2

; ---- initial estimates -------------------------------------------------
$THETA  (0, 0.2)             ; 1 CL    clearance (L/day)
$THETA  (0, 3.0)             ; 2 V     central volume (L)
$THETA  (0, 0.5)             ; 3 KD    equilibrium dissoc. constant (nM)
$THETA  (0, 0.5)             ; 4 KINT  complex internalisation rate (1/day)
$THETA  (0, 5.0)             ; 5 R0    baseline total target conc. (nM)

$OMEGA  0.09             ; 1 IIV CL
$OMEGA  0.09             ; 2 IIV V
$OMEGA  0 FIX            ; 3 IIV KD    (free it if the data support it)
$OMEGA  0.09             ; 4 IIV KINT
$OMEGA  0.09             ; 5 IIV R0

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
         PRED RES WRES ET1 ET2 ET3 ET4 ET5 SG1 SG2 SG3 SG4
         ONEHEADER NOPRINT FILE=tmdd_1cmt_wagner_iv.tab
$TABLE   ID TIME CFR CTOT RFR RCX RTOT
         ONEHEADER NOPRINT FILE=tmdd_1cmt_wagner_iv.prof
$TABLE   ID CL V KD KINT R0 ET1 ET2 ET3 ET4 ET5
         FIRSTONLY ONEHEADER NOPRINT FILE=tmdd_1cmt_wagner_iv.par

;; ----------------------------------------------------------------------
;; SIMULATION BLOCK
;; To simulate instead of estimate: comment out $ESTIMATION and
;; $COVARIANCE above, then un-comment the two lines below.
;; ----------------------------------------------------------------------
; $SIMULATION (20260830) (20260831 NORMAL) ONLYSIM SUBPROBLEMS=200
; $TABLE ID TIME AMT RATE EVID MDV CMT FLAG DV IPRED CFR CTOT RFR RCX
;        NOAPPEND ONEHEADER NOPRINT FILE=tmdd_1cmt_wagner_iv_sim.tab

