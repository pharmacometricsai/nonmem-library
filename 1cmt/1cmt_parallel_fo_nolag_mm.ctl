;; ==========================================================================
;; 1CMT_PARALLEL_FO_NOLAG_MM
;;   Parallel first-order absorption without lag time
;;   Nonlinear elimination
;; --------------------------------------------------------------------------
;; Structure
;;   DF x Dose     --> [Depot1] --Ka1--\
;;                                      >--> [Central (V)] --Vmax/Km-->
;;   (1-DF) x Dose --> [Depot2] --Ka2--/
;;
;; Data set / dosing requirements
;;   TWO dose records per administration, both with the full AMT:
;;     CMT=1 (depot 1), RATE=0   -> receives DF     x Dose via F1
;;     CMT=2 (depot 2), RATE=0   -> receives (1-DF) x Dose via F2
;;
;; Solver: general non-linear ODE (ADVAN13)
;; ==========================================================================

$PROBLEM 1CMT parallel_fo_nolag | Nonlinear elimination

$INPUT   ID TIME AMT RATE DV EVID MDV CMT ADDL II WT
$DATA    ../data/pkdata.csv IGNORE=@

$SUBROUTINES ADVAN13 TOL=9
$MODEL  NCOMPARTMENTS=3
        COMP=(DEPOT1, DEFDOSE)
        COMP=(DEPOT2)
        COMP=(CENTRAL, DEFOBS)

; compartment 1 = depot 1
; compartment 2 = depot 2
; compartment 3 = central (observation)

$PK
; ---- typical values and between-subject variability ------------------
  VM    = THETA(1)*EXP(ETA(1))
  KM    = THETA(2)*EXP(ETA(2))
  V     = THETA(3)*EXP(ETA(3))
  KA1   = THETA(4)*EXP(ETA(4))
  KA2   = THETA(5)*EXP(ETA(5))
  LGTDF = LOG(THETA(6)/(1-THETA(6))) + ETA(6)
  DF    = EXP(LGTDF)/(1+EXP(LGTDF))

; ---- random effects exported to $TABLE --------------------------------
  ET1   = ETA(1)
  ET2   = ETA(2)
  ET3   = ETA(3)
  ET4   = ETA(4)
  ET5   = ETA(5)
  ET6   = ETA(6)

; ---- structural / input specification ---------------------------------
  F1 = DF
  F2 = 1 - DF

; ---- scaling (concentration = amount / volume) ------------------------
  S3 = V

$DES
  CONC    = A(3)/V
  DADT(1) = -KA1*A(1)
  DADT(2) = -KA2*A(2)
  DADT(3) =  KA1*A(1) + KA2*A(2) - (VM*CONC/(KM+CONC))

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
$THETA  (0, 100.0)           ; 1 VM     Vmax (mg/h)
$THETA  (0, 5.0)             ; 2 KM     Km (mg/L)
$THETA  (0, 50.0)            ; 3 V      central volume (L)
$THETA  (0, 1.5)             ; 4 KA1    absorption rate constant, depot 1 (1/h)
$THETA  (0, 0.3)             ; 5 KA2    absorption rate constant, depot 2 (1/h)
$THETA  (0.001, 0.5, 0.999)  ; 6 DF     dose fraction (logit-scale IIV)

$OMEGA  0.09             ; 1 IIV VM
$OMEGA  0 FIX            ; 2 IIV KM    (free to estimate if supported)
$OMEGA  0.09             ; 3 IIV V
$OMEGA  0.09             ; 4 IIV KA1
$OMEGA  0.09             ; 5 IIV KA2
$OMEGA  0 FIX            ; 6 IIV DF    (free to estimate if supported)

$SIGMA  0.04             ; 1 proportional residual error (CV ~20%)
$SIGMA  0.01             ; 2 additive residual error (variance)

$ESTIMATION METHOD=1 INTER MAXEVAL=9999 NSIG=3 SIGL=9 PRINT=5 NOABORT POSTHOC
$COVARIANCE PRINT=E UNCONDITIONAL

$TABLE   ID TIME AMT RATE EVID MDV CMT DV IPRED IRES IWRES CWRES PRED
         RES WRES ET1 ET2 ET3 ET4 ET5 ET6 SG1 SG2
         ONEHEADER NOPRINT FILE=1cmt_parallel_fo_nolag_mm.tab
$TABLE   ID VM KM V KA1 KA2 DF ET1 ET2 ET3 ET4 ET5 ET6
         FIRSTONLY ONEHEADER NOPRINT FILE=1cmt_parallel_fo_nolag_mm.par

;; ----------------------------------------------------------------------
;; SIMULATION BLOCK
;; To simulate instead of estimate: comment out $ESTIMATION and
;; $COVARIANCE above, then un-comment the two lines below.
;; ----------------------------------------------------------------------
; $SIMULATION (20260830) (20260831 NORMAL) ONLYSIM SUBPROBLEMS=200
; $TABLE ID TIME AMT RATE EVID MDV CMT DV IPRED
;        NOAPPEND ONEHEADER NOPRINT FILE=1cmt_parallel_fo_nolag_mm_sim.tab

