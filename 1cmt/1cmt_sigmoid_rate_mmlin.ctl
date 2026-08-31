;; ==========================================================================
;; 1CMT_SIGMOID_RATE_MMLIN
;;   Sigmoid absorption with rate
;;   Nonlinear and linear elimination
;; --------------------------------------------------------------------------
;; Structure
;;   Dose --(zero-order, R1)--> [Depot] --Ka--> [Central (V)] --CL +
;;     Vmax/Km-->
;;
;; Data set / dosing requirements
;;   One dose record per administration: CMT=1 (depot), RATE=-1
;;   (rate R1 modelled in $PK).
;;
;; Solver: general non-linear ODE (ADVAN13)
;; ==========================================================================

$PROBLEM 1CMT sigmoid_rate | Nonlinear and linear elimination

$INPUT   ID TIME AMT RATE DV EVID MDV CMT ADDL II WT
$DATA    ../data/pkdata.csv IGNORE=@

$SUBROUTINES ADVAN13 TOL=9
$MODEL  NCOMPARTMENTS=2
        COMP=(DEPOT, DEFDOSE)
        COMP=(CENTRAL, DEFOBS)

; compartment 1 = depot
; compartment 2 = central (observation)

$PK
; ---- typical values and between-subject variability ------------------
  CL    = THETA(1)*EXP(ETA(1))
  V     = THETA(2)*EXP(ETA(2))
  VM    = THETA(3)*EXP(ETA(3))
  KM    = THETA(4)*EXP(ETA(4))
  KA    = THETA(5)*EXP(ETA(5))
  R1    = THETA(6)*EXP(ETA(6))

; ---- random effects exported to $TABLE --------------------------------
  ET1   = ETA(1)
  ET2   = ETA(2)
  ET3   = ETA(3)
  ET4   = ETA(4)
  ET5   = ETA(5)
  ET6   = ETA(6)

; ---- scaling (concentration = amount / volume) ------------------------
  S2 = V

$DES
  CONC    = A(2)/V
  DADT(1) = -KA*A(1)
  DADT(2) =  KA*A(1) - (CL*CONC + VM*CONC/(KM+CONC))

$ERROR
  IPRED = F
  IF(IPRED.LE.0) IPRED = 1.0E-10
  Y     = IPRED + IPRED*EPS(1) + EPS(2)   ; combined prop. + add. error
  IRES  = DV - IPRED
  IWRES = IRES/IPRED

; ---- residual error variances exported to $TABLE ----------------------
  SG1   = SIGMA(1,1)          ; proportional
  SG2   = SIGMA(2,2)          ; additive

; ---- initial estimates -------------------------------------------------
$THETA  (0, 5.0)             ; 1 CL     linear clearance (L/h)
$THETA  (0, 50.0)            ; 2 V      central volume (L)
$THETA  (0, 100.0)           ; 3 VM     Vmax (mg/h)
$THETA  (0, 5.0)             ; 4 KM     Km (mg/L)
$THETA  (0, 1.0)             ; 5 KA     absorption rate constant (1/h)
$THETA  (0, 50.0)            ; 6 R1     zero-order rate into cmt 1 (mg/h)

$OMEGA  0.09             ; 1 IIV CL
$OMEGA  0.09             ; 2 IIV V
$OMEGA  0.09             ; 3 IIV VM
$OMEGA  0 FIX            ; 4 IIV KM    (free to estimate if supported)
$OMEGA  0.09             ; 5 IIV KA
$OMEGA  0 FIX            ; 6 IIV R1    (free to estimate if supported)

$SIGMA  0.04             ; 1 proportional residual error (CV ~20%)
$SIGMA  0.01             ; 2 additive residual error (variance)

$ESTIMATION METHOD=1 INTER MAXEVAL=9999 NSIG=3 SIGL=9 PRINT=5 NOABORT POSTHOC
$COVARIANCE PRINT=E UNCONDITIONAL

$TABLE   ID TIME AMT RATE EVID MDV CMT DV IPRED IRES IWRES CWRES PRED
         RES WRES ET1 ET2 ET3 ET4 ET5 ET6 SG1 SG2
         ONEHEADER NOPRINT FILE=1cmt_sigmoid_rate_mmlin.tab
$TABLE   ID CL V VM KM KA R1 ET1 ET2 ET3 ET4 ET5 ET6
         FIRSTONLY ONEHEADER NOPRINT FILE=1cmt_sigmoid_rate_mmlin.par

;; ----------------------------------------------------------------------
;; SIMULATION BLOCK
;; To simulate instead of estimate: comment out $ESTIMATION and
;; $COVARIANCE above, then un-comment the two lines below.
;; ----------------------------------------------------------------------
; $SIMULATION (20260830) (20260831 NORMAL) ONLYSIM SUBPROBLEMS=200
; $TABLE ID TIME AMT RATE EVID MDV CMT DV IPRED
;        NOAPPEND ONEHEADER NOPRINT FILE=1cmt_sigmoid_rate_mmlin_sim.tab

