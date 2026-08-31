;; ==========================================================================
;; 1CMT_PARALLEL_FO_LAG_LIN
;; #25 of the one-compartment set
;;   Parallel first-order absorption with lag time
;;   Linear elimination
;; --------------------------------------------------------------------------
;; Structure
;;   DF x Dose     at t=ALAG1 --> [Depot1] --Ka1--\
;;                                                 >--> [Central (V)] --CL-->
;;   (1-DF) x Dose at t=ALAG2 --> [Depot2] --Ka2--/
;;
;; Data set / dosing requirements
;;   TWO dose records per administration, both with the full AMT:
;;     CMT=1 (depot 1), RATE=0   -> DF     x Dose, delayed by ALAG1
;;     CMT=2 (depot 2), RATE=0   -> (1-DF) x Dose, delayed by ALAG2
;;
;; Solver: general non-linear ODE (ADVAN13)
;; ==========================================================================

$PROBLEM 1CMT parallel_fo_lag | Linear elimination

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
  CL    = THETA(1)*EXP(ETA(1))
  V     = THETA(2)*EXP(ETA(2))
  KA1   = THETA(3)*EXP(ETA(3))
  KA2   = THETA(4)*EXP(ETA(4))
  LGTDF = LOG(THETA(5)/(1-THETA(5))) + ETA(5)
  DF    = EXP(LGTDF)/(1+EXP(LGTDF))
  ALAG1 = THETA(6)*EXP(ETA(6))
  ALAG2 = THETA(7)*EXP(ETA(7))

; ---- structural / input specification ---------------------------------
  F1 = DF
  F2 = 1 - DF

; ---- scaling (concentration = amount / volume) ------------------------
  S3 = V

$DES
  CONC    = A(3)/V
  DADT(1) = -KA1*A(1)
  DADT(2) = -KA2*A(2)
  DADT(3) =  KA1*A(1) + KA2*A(2) - (CL*CONC)

$ERROR
  IPRED = F
  IF(IPRED.LE.0) IPRED = 1.0E-10
  Y     = IPRED + IPRED*EPS(1) + EPS(2)   ; combined prop. + add. error
  IRES  = DV - IPRED
  IWRES = IRES/IPRED

; ---- initial estimates -------------------------------------------------
$THETA  (0, 5.0)             ; 1 CL     clearance (L/h)
$THETA  (0, 50.0)            ; 2 V      central volume (L)
$THETA  (0, 1.5)             ; 3 KA1    absorption rate constant, depot 1 (1/h)
$THETA  (0, 0.3)             ; 4 KA2    absorption rate constant, depot 2 (1/h)
$THETA  (0.001, 0.5, 0.999)  ; 5 DF     dose fraction (logit-scale IIV)
$THETA  (0, 0.25)            ; 6 ALAG1  lag time, compartment 1 (h)
$THETA  (0, 1.0)             ; 7 ALAG2  lag time, compartment 2 (h)

$OMEGA  0.09             ; 1 IIV CL
$OMEGA  0.09             ; 2 IIV V
$OMEGA  0.09             ; 3 IIV KA1
$OMEGA  0.09             ; 4 IIV KA2
$OMEGA  0 FIX            ; 5 IIV DF    (free to estimate if supported)
$OMEGA  0 FIX            ; 6 IIV ALAG1 (free to estimate if supported)
$OMEGA  0 FIX            ; 7 IIV ALAG2 (free to estimate if supported)

$SIGMA  0.04             ; 1 proportional residual error (CV ~20%)
$SIGMA  0.01             ; 2 additive residual error (variance)

$ESTIMATION METHOD=1 INTER MAXEVAL=9999 NSIG=3 SIGL=9 PRINT=5 NOABORT POSTHOC
$COVARIANCE PRINT=E UNCONDITIONAL

$TABLE   ID TIME AMT RATE EVID MDV CMT DV IPRED IRES IWRES CWRES
         PRED RES WRES
         ONEHEADER NOPRINT FILE=1cmt_parallel_fo_lag_lin.tab
$TABLE   ID CL V KA1 KA2 DF ALAG1 ALAG2 ETA1 ETA2 ETA3 ETA4 ETA5 ETA6
         ETA7
         FIRSTONLY ONEHEADER NOPRINT FILE=1cmt_parallel_fo_lag_lin.par

;; ----------------------------------------------------------------------
;; SIMULATION BLOCK
;; To simulate instead of estimate: comment out $ESTIMATION and
;; $COVARIANCE above, then un-comment the two lines below.
;; ----------------------------------------------------------------------
; $SIMULATION (20260830) (20260831 NORMAL) ONLYSIM SUBPROBLEMS=200
; $TABLE ID TIME AMT RATE EVID MDV CMT DV IPRED
;        NOAPPEND ONEHEADER NOPRINT FILE=1cmt_parallel_fo_lag_lin.tab

