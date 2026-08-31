;; ==========================================================================
;; TMDD_1CMT_QSS_SCLAG
;; Target-mediated drug disposition - 1-compartment model
;;   Form : Quasi-steady-state approximation
;;   Input: Subcutaneous, first-order absorption with lag time
;; --------------------------------------------------------------------------
;; Structure
;;   Dose --> [SC depot] --KA--> [Plasma (V)]
;;   [Plasma (V)] --CL--> eliminated
;;   drug + target in quasi-steady state (KSS) binding; complex cleared by
;;     KINT
;;   target: KSYN = KDEG*R0 in, KDEG out
;;
;; Assumptions
;;   As the QE form, but the complex is assumed to be at quasi-steady
;;   state rather than at equilibrium, so KD is replaced by
;;   KSS = (KOFF + KINT)/KON.  QSS degenerates to QE when KINT << KOFF
;;   (Gibiansky et al. 2008).
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

$PROBLEM TMDD 1cmt qss | sclag

$INPUT   ID TIME AMT RATE DV EVID MDV CMT FLAG ADDL II WT
$DATA    ../data/tmdddata.csv IGNORE=@

$SUBROUTINES ADVAN13 TOL=9
$MODEL  NCOMPARTMENTS=3
        COMP=(DEPOT, DEFDOSE)
        COMP=(CENTRAL, DEFOBS)
        COMP=(TARGET)

; compartment 1 = SC depot
; compartment 2 = central drug
; compartment 3 = total target

$PK
; ---- typical values and between-subject variability ------------------
  CL    = THETA(1)*EXP(ETA(1))
  V     = THETA(2)*EXP(ETA(2))
  KSS   = THETA(3)*EXP(ETA(3))
  KINT  = THETA(4)*EXP(ETA(4))
  KDEG  = THETA(5)*EXP(ETA(5))
  R0    = THETA(6)*EXP(ETA(6))
  KA    = THETA(7)*EXP(ETA(7))
  LGTF1 = LOG(THETA(8)/(1-THETA(8))) + ETA(8)
  F1    = EXP(LGTF1)/(1+EXP(LGTF1))
  ALAG1 = THETA(9)*EXP(ETA(9))

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

; ---- pre-dose target baseline ----------------------------------------
;     KSYN = KDEG*R0, so the target starts at steady state R0.
  A_0(3) = R0*V

; ---- scaling (concentration = amount / volume) ------------------------
  S2 = V

$DES
  CTOT = A(2)/V
  IF(CTOT.LT.0) CTOT = 0
  RTOT = A(3)/V
  BB   = CTOT - RTOT - KSS
  CF   = 0.5*(BB + SQRT(BB*BB + 4*KSS*CTOT))
  RC   = CTOT - CF
  RF   = RTOT - RC
  DADT(1) = -KA*A(1)
  DADT(2) = KA*A(1) - CL*CF - KINT*RC*V
  DADT(3) =  KDEG*R0*V - KDEG*RF*V - KINT*RC*V

$ERROR
  CTOT = A(2)/V
  IF(CTOT.LT.0) CTOT = 0
  RTOT = A(3)/V
  BB   = CTOT - RTOT - KSS
  CFR  = 0.5*(BB + SQRT(BB*BB + 4*KSS*CTOT))
  RCX  = CTOT - CFR
  RFR  = RTOT - RCX
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
$THETA  (0, 1.0)             ; 3 KSS   quasi-steady-state constant (nM)
$THETA  (0, 0.5)             ; 4 KINT  complex internalisation rate (1/day)
$THETA  (0, 1.0)             ; 5 KDEG  free target degradation rate (1/day)
$THETA  (0, 5.0)             ; 6 R0    baseline total target conc. (nM)
$THETA  (0, 0.25)            ; 7 KA    absorption rate constant (1/day)
$THETA  (0.001, 0.7, 0.999)  ; 8 F1    SC bioavailability (logit-scale IIV)
$THETA  (0, 0.1)             ; 9 ALAG1 absorption lag time (day)

$OMEGA  0.09             ; 1 IIV CL
$OMEGA  0.09             ; 2 IIV V
$OMEGA  0 FIX            ; 3 IIV KSS   (free it if the data support it)
$OMEGA  0.09             ; 4 IIV KINT
$OMEGA  0.09             ; 5 IIV KDEG
$OMEGA  0.09             ; 6 IIV R0
$OMEGA  0.09             ; 7 IIV KA
$OMEGA  0 FIX            ; 8 IIV F1    (free it if the data support it)
$OMEGA  0 FIX            ; 9 IIV ALAG1 (free it if the data support it)

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
         PRED RES WRES ET1 ET2 ET3 ET4 ET5 ET6 ET7 ET8 ET9 SG1 SG2 SG3
         SG4
         ONEHEADER NOPRINT FILE=tmdd_1cmt_qss_sclag.tab
$TABLE   ID TIME CFR CTOT RFR RCX RTOT
         ONEHEADER NOPRINT FILE=tmdd_1cmt_qss_sclag.prof
$TABLE   ID CL V KSS KINT KDEG R0 KA F1 ALAG1 ET1 ET2 ET3 ET4 ET5 ET6
         ET7 ET8 ET9
         FIRSTONLY ONEHEADER NOPRINT FILE=tmdd_1cmt_qss_sclag.par

;; ----------------------------------------------------------------------
;; SIMULATION BLOCK
;; To simulate instead of estimate: comment out $ESTIMATION and
;; $COVARIANCE above, then un-comment the two lines below.
;; ----------------------------------------------------------------------
; $SIMULATION (20260830) (20260831 NORMAL) ONLYSIM SUBPROBLEMS=200
; $TABLE ID TIME AMT RATE EVID MDV CMT FLAG DV IPRED CFR CTOT RFR RCX
;        NOAPPEND ONEHEADER NOPRINT FILE=tmdd_1cmt_qss_sclag_sim.tab

